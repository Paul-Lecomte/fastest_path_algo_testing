// Comparison of Dijkstra and A* algorithms on a public-transport graph (GTFS-like data)
// Support for either simulated GTFS data or loading GTFS files from a folder
// Prototype originally considered in JS; now implemented in C++ for performance testing

#include <algorithm>
#include <iostream>
#include <stdexcept>
#include <vector>
#include <string>
#include <unordered_map>
#include <queue>
#include <climits>
#include <random>
#include <chrono>
#include <cmath>
#include <fstream>
#include <sstream>
#include <functional>
using namespace std;

// safeStoi: helper that validates a numeric string before calling stoi
// This prevents exceptions for empty or non-numeric strings and produces clearer error messages.
int safeStoi(const string& s) {
    if (s.empty()) throw std::invalid_argument("Empty string passed to stoi");
    size_t i = 0;
    if (s[0] == '-' || s[0] == '+') i = 1;
    for (; i < s.size(); ++i) {
        if (!isdigit(s[i]))
            throw std::invalid_argument("Non-numeric string for stoi: " + s);
    }
    return std::stoi(s);
}

// GTFS-like data structures (simplified versions of common GTFS entities)
// These structs hold strings for ids and basic attributes; numeric times are handled in StopTime.
struct Agency {
    string id, name, url, timezone;
};
struct Calendar {
    string service_id;
    bool monday, tuesday, wednesday, thursday, friday, saturday, sunday;
    string start_date, end_date;
};
struct CalendarDate {
    string service_id, date;
    int exception_type;
};
struct FeedInfo {
    string publisher_name, publisher_url, lang, start_date, end_date;
};
struct Route {
    string id, agency_id, short_name, long_name, type;
};
struct Stop {
    string id, name, desc, lat, lon;
};
struct Trip {
    string id, route_id, service_id, headsign;
};

// StopTime stores the scheduled arrival and departure times (converted to minutes).
struct StopTime {
    string trip_id, stop_id;
    int arrival_time, departure_time, stop_sequence;
};
struct Transfer {
    string from_stop_id, to_stop_id;
    int min_transfer_time; // minimum required transfer time in minutes
};

// Node represents an event vertex in a time-expanded graph: a stop at a particular time.
// Two nodes are equal if they reference the same stop and the same time.
struct Node {
    string stop_id;
    int time;
    bool operator==(const Node& o) const { return stop_id == o.stop_id && time == o.time; }
};
struct NodeHash {
    size_t operator()(const Node& n) const { return hash<string>()(n.stop_id) ^ hash<int>()(n.time); }
};

// Simple heuristic for A*.
// Currently uses numeric difference of stop ids when they are formatted like 'S123' -> extracts 123.
// If parsing fails, it returns 0 which reduces A* to Dijkstra behavior.
int heuristic(const Node& n, const string& goal_id) {
    int id1 = 0, id2 = 0;
    try {
        id1 = safeStoi(n.stop_id.substr(1));
        id2 = safeStoi(goal_id.substr(1));
    } catch (...) {
        return 0; // fallback to admissible heuristic 0
    }
    return abs(id1 - id2);
}

// Generic CSV reader for GTFS files. It reads a file line-by-line, skips the header
// and calls a caller-provided parser lambda to convert CSV fields to a struct T.
// Returns true if the file could be opened and processed (parsing errors are reported but do not abort).
template<typename T>
bool readCSV(const string& filename, vector<T>& container, function<T(const vector<string>&)> parser) {
    ifstream file(filename);
    if (!file.is_open()) return false;
    string line;
    getline(file, line); // skip header line
    while (getline(file, line)) {
        vector<string> fields;
        stringstream ss(line);
        string field;
        while (getline(ss, field, ',')) {
            // Trim surrounding quotes if present
            if (!field.empty() && field.front() == '"' && field.back() == '"')
                field = field.substr(1, field.length() - 2);
            fields.push_back(field);
        }
        if (!fields.empty()) {
            try {
                container.push_back(parser(fields));
            } catch (const std::exception& e) {
                cerr << "Error parsing line: " << line << " (" << e.what() << ")\n";
            }
        }
    }
    return true;
}

// Parsers for each GTFS-like file. These expect fields in a specific order; they are minimal and don't validate field counts.
Agency parseAgency(const vector<string>& f) { return {f[0], f[1], f[2], f[3]}; }
Route parseRoute(const vector<string>& f) { return {f[0], f[1], f[2], f[3], f[4]}; }
Stop parseStop(const vector<string>& f) { return {f[0], f[1], f[2], f[3], f[4]}; }
Trip parseTrip(const vector<string>& f) { return {f[0], f[1], f[2], f[3]}; }

// Improved StopTime parser: converts HH:MM:SS (or HH:MM) to integer minutes.
StopTime parseStopTime(const vector<string>& f) {
    auto timeToMin = [](const string& t) {
        int h=0, m=0, s=0;
        // sscanf will parse H:M or H:M:S; if it reads at least H and M, use minutes
        if (sscanf(t.c_str(), "%d:%d:%d", &h, &m, &s) < 2) return 0;
        // convert to minutes, rounding seconds >= 30 up to the next minute
        return h*60 + m + (s>=30?1:0);
    };
    int stop_sequence = 0;
    if (f.size() > 4) {
        try {
            stop_sequence = safeStoi(f[4]);
        } catch (...) {
            stop_sequence = 0;
        }
    }
    return {f[0], f[3], timeToMin(f[1]), timeToMin(f[2]), stop_sequence};
}

// Transfer parser: some GTFS transfer files include a minimum transfer time field; parse it if present.
Transfer parseTransfer(const vector<string>& f) {
    int min_transfer_time = 0;
    if (f.size() > 3) {
        try {
            min_transfer_time = safeStoi(f[3]);
        } catch (...) {
            min_transfer_time = 0;
        }
    }
    return {f[0], f[1], min_transfer_time};
}

// Load GTFS-like data from a folder by reading standard GTFS filenames. Returns true if all mandatory files could be opened.
bool loadGTFSData(
    const string& gtfsFolder,
    vector<Agency>& agencies, vector<Route>& routes, vector<Trip>& trips,
    vector<Stop>& stops, vector<StopTime>& stopTimes, vector<Transfer>& transfers
) {
    bool ok = true;
    ok &= readCSV(gtfsFolder + "/agency.txt", agencies, std::function<Agency(const vector<string>&)>(parseAgency));
    ok &= readCSV(gtfsFolder + "/routes.txt", routes, std::function<Route(const vector<string>&)>(parseRoute));
    ok &= readCSV(gtfsFolder + "/stops.txt", stops, std::function<Stop(const vector<string>&)>(parseStop));
    ok &= readCSV(gtfsFolder + "/trips.txt", trips, std::function<Trip(const vector<string>&)>(parseTrip));
    ok &= readCSV(gtfsFolder + "/stop_times.txt", stopTimes, std::function<StopTime(const vector<string>&)>(parseStopTime));
    // transfers are optional: try to open and parse if present
    ifstream f(gtfsFolder + "/transfers.txt");
    if (f.is_open()) { f.close(); ok &= readCSV(gtfsFolder + "/transfers.txt", transfers, std::function<Transfer(const vector<string>&)>(parseTransfer)); }
    return ok;
}

// Generate synthetic GTFS-like data for testing and benchmarking. Default sizes are large to approximate a realistic dataset.
// The function randomly creates sequences of stop_times for trips and random transfers between stops.
void generateData(
    vector<Agency>& agencies, vector<Route>& routes, vector<Trip>& trips,
    vector<Stop>& stops, vector<StopTime>& stopTimes, vector<Transfer>& transfers,
    int N_AGENCIES = 467, int N_ROUTES = 5000, int N_TRIPS = 1800000, int N_STOPS = 93000, int N_TRANSFERS = 112000
) {
    // Create agencies with simple identifiers and URLs
    for (int i = 0; i < N_AGENCIES; ++i)
        agencies.push_back({"A" + to_string(i), "Agency " + to_string(i), "https://agency" + to_string(i) + ".fr", "Europe/Paris"});
    // Create stops with synthetic lat/lon strings (not real coordinates) and a short description
    for (int i = 0; i < N_STOPS; ++i)
        stops.push_back({"S" + to_string(i), "Stop " + to_string(i), "Desc", "48." + to_string(800 + i), "2." + to_string(200 + i)});
    // Create routes and assign them to agencies in round-robin fashion
    for (int i = 0; i < N_ROUTES; ++i)
        routes.push_back({"R" + to_string(i), agencies[i % N_AGENCIES].id, "L" + to_string(i), "Ligne " + to_string(i), "3"});
    // Create trips associated with routes; service_id is set to 'WD' (weekday) for simplicity
    for (int i = 0; i < N_TRIPS; ++i)
        trips.push_back({"T" + to_string(i), routes[i % N_ROUTES].id, "WD", "Towards centre"});

    // Create stop times for each trip using random lengths, starting hours and spacing between stops
    random_device rd; mt19937 gen(rd());
    uniform_int_distribution<> stopDist(0, N_STOPS - 1), tripLenDist(5, 10), timeDist(2, 8), hourDist(6, 21);
    for (auto& trip : trips) {
        int tripLen = tripLenDist(gen);
        int startHour = hourDist(gen), t = startHour * 60 + uniform_int_distribution<>(0, 59)(gen);
        int stopIdx = stopDist(gen);
        for (int j = 0; j < tripLen; ++j) {
            string stop_id = stops[(stopIdx + j) % N_STOPS].id;
            stopTimes.push_back({trip.id, stop_id, t, t, j});
            // Add a random transit time to the next stop
            t += timeDist(gen);
        }
    }
    // Create random transfers between stops with random minimum transfer times
    uniform_int_distribution<> transferTimeDist(1, 10);
    for (int i = 0; i < N_TRANSFERS; ++i) {
        int fromIdx = stopDist(gen), toIdx = stopDist(gen);
        if (fromIdx != toIdx)
            transfers.push_back({stops[fromIdx].id, stops[toIdx].id, transferTimeDist(gen)});
    }
}

// Build a time-expanded directed graph from trips, stop times and transfers. The graph maps Node -> list of (Node, cost) edges.
// Nodes represent events (stop + time). Edges include: trip edges (ride a vehicle), transfer edges (walk between stops), and waiting edges at the same stop.
void buildGraph(
    const vector<Trip>& trips, const vector<StopTime>& stopTimes, const vector<Transfer>& transfers,
    unordered_map<Node, vector<pair<Node, int>>, NodeHash>& graph
) {
    // Index stop_times by trip to create edges along each vehicle trip (from one stop event to the next)
    unordered_map<string, vector<const StopTime*>> tripStops;
    for (const auto& st : stopTimes)
        tripStops[st.trip_id].push_back(&st);
    // Create trip edges: for each consecutive pair in a trip, create an edge from departure to arrival with travel time cost
    for (const auto& [trip_id, stopsVec] : tripStops) {
        auto stopsCopy = stopsVec;
        sort(stopsCopy.begin(), stopsCopy.end(), [](const StopTime* a, const StopTime* b) { return a->stop_sequence < b->stop_sequence; });
        for (size_t i = 1; i < stopsCopy.size(); ++i) {
            Node from{stopsCopy[i-1]->stop_id, stopsCopy[i-1]->departure_time};
            Node to{stopsCopy[i]->stop_id, stopsCopy[i]->arrival_time};
            int cost = to.time - from.time; // travel time between the events
            graph[from].push_back({to, cost});
        }
    }

    // Index events by stop to build waiting and transfer edges
    unordered_map<string, vector<const StopTime*>> stopEvents;
    for (const auto& st : stopTimes)
        stopEvents[st.stop_id].push_back(&st);

    // Create transfer edges: for each transfer record, link appropriate departure events at the origin stop
    // to the earliest reachable arrival event at the destination stop that occurs after the minimum transfer time.
    for (const auto& tr : transfers) {
        auto& fromEvents = stopEvents[tr.from_stop_id];
        auto& toEvents = stopEvents[tr.to_stop_id];
        for (const auto* st : fromEvents) {
            int t = st->departure_time + tr.min_transfer_time;
            auto it = lower_bound(toEvents.begin(), toEvents.end(), t,
                [](const StopTime* st, int t) { return st->arrival_time < t; });
            if (it != toEvents.end()) {
                int t2 = (*it)->arrival_time;
                Node from{tr.from_stop_id, st->departure_time};
                Node to{tr.to_stop_id, t2};
                graph[from].push_back({to, t2 - st->departure_time});
            }
        }
    }

    // Create waiting edges at the same stop between consecutive departure events so an agent can wait for the next vehicle.
    for (auto& [stop_id, evs] : stopEvents) {
        sort(evs.begin(), evs.end(), [](const StopTime* a, const StopTime* b) { return a->departure_time < b->departure_time; });
        for (size_t i = 1; i < evs.size(); ++i) {
            int t1 = evs[i-1]->departure_time, t2 = evs[i]->departure_time;
            if (t2 > t1)
                graph[{stop_id, t1}].push_back({{stop_id, t2}, t2 - t1});
        }
    }
}

// Classic Dijkstra search on the time-expanded graph. startNodes are event nodes representing the possible starting times at the origin stop.
// The function uses a priority queue keyed by f-score (cost + heuristic) to allow reuse with the same queue structure for A* as well.
// Returns the minimum travel time (in minutes) and fills outPath with the sequence of event nodes (stop+time). Returns -1 if no path.
int dijkstra(
    const unordered_map<Node, vector<pair<Node, int>>, NodeHash>& graph,
    const vector<Node>& startNodes, const string& endStopId, vector<Node>& outPath
) {
    unordered_map<Node, int, NodeHash> dist;
    unordered_map<Node, Node, NodeHash> prev;
    auto cmp = [&](const pair<int, Node>& a, const pair<int, Node>& b) { return a.first > b.first; };
    priority_queue<pair<int, Node>, vector<pair<int, Node>>, decltype(cmp)> pq(cmp);
    for (auto& node : startNodes) { dist[node] = 0; pq.push({0, node}); }
    Node bestEnd; int minCost = INT_MAX;
    while (!pq.empty()) {
        auto [fscore, u] = pq.top(); pq.pop();
        int cost = dist[u];
        // If we reached the destination stop (any event at that stop), record the best one
        if (u.stop_id == endStopId && cost < minCost) { minCost = cost; bestEnd = u; }
        if (dist[u] < cost) continue; // ignore outdated queue entries
        auto it = graph.find(u);
        if (it != graph.end()) {
            for (auto& [v, edgeCost] : it->second) {
                int newCost = cost + edgeCost;
                if (!dist.count(v) || newCost < dist[v]) {
                    dist[v] = newCost; prev[v] = u;
                    // push using f = g + h where heuristic is allowed (here heuristic may be zero)
                    int f = newCost + heuristic(v, endStopId);
                    pq.push({f, v});
                }
            }
        }
    }
    if (minCost == INT_MAX) return -1;
    // Reconstruct path from bestEnd backwards using the prev map
    for (Node at = bestEnd; dist.count(at); at = prev[at]) {
        outPath.push_back(at); if (!prev.count(at)) break;
    }
    reverse(outPath.begin(), outPath.end());
    return minCost;
}

// A* search on the same graph. Implementation is very similar to dijkstra above but starts the priority queue with heuristic values.
// The code uses the same dist/prev containers and reconstructs the path the same way.
int astar(
    const unordered_map<Node, vector<pair<Node, int>>, NodeHash>& graph,
    const vector<Node>& startNodes, const string& endStopId, vector<Node>& outPath
) {
    unordered_map<Node, int, NodeHash> dist;
    unordered_map<Node, Node, NodeHash> prev;
    auto cmp = [&](const pair<int, Node>& a, const pair<int, Node>& b) { return a.first > b.first; };
    priority_queue<pair<int, Node>, vector<pair<int, Node>>, decltype(cmp)> pq(cmp);
    for (auto& node : startNodes) { dist[node] = 0; pq.push({heuristic(node, endStopId), node}); }
    Node bestEnd; int minCost = INT_MAX;
    while (!pq.empty()) {
        auto [fscore, u] = pq.top(); pq.pop();
        int cost = dist[u];
        if (u.stop_id == endStopId && cost < minCost) { minCost = cost; bestEnd = u; }
        if (dist[u] < cost) continue;
        auto it = graph.find(u);
        if (it != graph.end()) {
            for (auto& [v, edgeCost] : it->second) {
                int newCost = cost + edgeCost;
                if (!dist.count(v) || newCost < dist[v]) {
                    dist[v] = newCost; prev[v] = u;
                    int f = newCost + heuristic(v, endStopId);
                    pq.push({f, v});
                }
            }
        }
    }
    if (minCost == INT_MAX) return -1;
    for (Node at = bestEnd; dist.count(at); at = prev[at]) {
        outPath.push_back(at); if (!prev.count(at)) break;
    }
    reverse(outPath.begin(), outPath.end());
    return minCost;
}

// Dijkstra on a time-expanded graph built entirely inside this function. This provides an alternative implementation
// that constructs its own EventNode type and graph and demonstrates the same technique used in buildGraph.
int dijkstra_time_expanded(
    const vector<Stop>& stops,
    const vector<Trip>& trips,
    const vector<StopTime>& stopTimes,
    const vector<Transfer>& transfers,
    const string& startStopId, int startTime,
    const string& endStopId,
    vector<pair<string, int>>& outPath
) {
    // 1) Index events by stop and sort by departure time so we can build waiting edges quickly
    unordered_map<string, vector<const StopTime*>> events;
    for (const auto& st : stopTimes)
        events[st.stop_id].push_back(&st);
    for (auto& [stop, evs] : events)
        sort(evs.begin(), evs.end(), [](const StopTime* a, const StopTime* b) { return a->departure_time < b->departure_time; });

    // EventNode and hash for the local graph representation
    struct EventNode {
        string stop_id;
        int time;
        bool operator==(const EventNode& o) const { return stop_id == o.stop_id && time == o.time; }
    };
    struct EventNodeHash {
        size_t operator()(const EventNode& n) const { return hash<string>()(n.stop_id) ^ hash<int>()(n.time); }
    };

    // 2) Build the local graph: trip edges, waiting edges and transfers (same logic as buildGraph)
    unordered_map<EventNode, vector<pair<EventNode, int>>, EventNodeHash> graph;
    // a) Trip edges by trip_id
    unordered_map<string, vector<const StopTime*>> tripStops;
    for (const auto& st : stopTimes)
        tripStops[st.trip_id].push_back(&st);
    for (const auto& [trip_id, stopsVec] : tripStops) {
        auto stopsCopy = stopsVec;
        sort(stopsCopy.begin(), stopsCopy.end(), [](const StopTime* a, const StopTime* b) { return a->stop_sequence < b->stop_sequence; });
        for (size_t i = 1; i < stopsCopy.size(); ++i) {
            EventNode from{stopsCopy[i-1]->stop_id, stopsCopy[i-1]->departure_time};
            EventNode to{stopsCopy[i]->stop_id, stopsCopy[i]->arrival_time};
            int cost = to.time - from.time;
            graph[from].push_back({to, cost});
        }
    }
    // b) Waiting edges at same stop
    for (auto& [stop, evs] : events) {
        for (size_t i = 1; i < evs.size(); ++i) {
            int t1 = evs[i-1]->departure_time, t2 = evs[i]->departure_time;
            if (t2 > t1)
                graph[{stop, t1}].push_back({{stop, t2}, t2 - t1});
        }
    }
    // c) Transfer edges between stops using transfers dataset
    for (const auto& tr : transfers) {
        auto& fromEvents = events[tr.from_stop_id];
        auto& toEvents = events[tr.to_stop_id];
        for (const auto* st : fromEvents) {
            int t = st->departure_time + tr.min_transfer_time;
            auto it = lower_bound(toEvents.begin(), toEvents.end(), t,
                [](const StopTime* st, int t) { return st->arrival_time < t; });
            if (it != toEvents.end()) {
                int t2 = (*it)->arrival_time;
                graph[{tr.from_stop_id, st->departure_time}].push_back({{tr.to_stop_id, t2}, t2 - st->departure_time});
            }
        }
    }

    // 3) Run standard Dijkstra on the constructed EventNode graph
    unordered_map<EventNode, int, EventNodeHash> dist;
    unordered_map<EventNode, EventNode, EventNodeHash> prev;
    auto cmp = [&](const pair<int, EventNode>& a, const pair<int, EventNode>& b) { return a.first > b.first; };
    priority_queue<pair<int, EventNode>, vector<pair<int, EventNode>>, decltype(cmp)> pq(cmp);

    // Source: all events at startStop with departure_time >= startTime
    vector<const StopTime*> startEvents;
    for (const auto* st : events[startStopId])
        if (st->departure_time >= startTime)
            startEvents.push_back(st);
    if (startEvents.empty()) return -1; // no feasible starting event after startTime
    for (const auto* st : startEvents) {
        EventNode n{startStopId, st->departure_time};
        dist[n] = 0;
        pq.push({0, n});
    }

    // Dijkstra main loop: find the earliest arrival event at endStop
    EventNode bestEnd; int minCost = INT_MAX;
    while (!pq.empty()) {
        auto [cost, u] = pq.top(); pq.pop();
        if (u.stop_id == endStopId && cost < minCost) { minCost = cost; bestEnd = u; }
        if (dist[u] < cost) continue;
        for (auto& [v, edgeCost] : graph[u]) {
            int newCost = cost + edgeCost;
            if (!dist.count(v) || newCost < dist[v]) {
                dist[v] = newCost; prev[v] = u;
                pq.push({newCost, v});
            }
        }
    }
    if (minCost == INT_MAX) return -1;

    // Reconstruct the found path (stop_id, time) pairs
    for (EventNode at = bestEnd; dist.count(at); at = prev[at]) {
        outPath.push_back({at.stop_id, at.time});
        if (!prev.count(at)) break;
    }
    reverse(outPath.begin(), outPath.end());
    return minCost;
}

int main() {
    // Containers to hold loaded or generated GTFS-like data
    vector<Agency> agencies;
    vector<Route> routes;
    vector<Trip> trips;
    vector<Stop> stops;
    vector<StopTime> stopTimes;
    vector<Transfer> transfers;

    cout << "Data mode:\n1. Generate data\n2. Load GTFS\nYour choice: ";
    int choix; cin >> choix;
    if (choix == 2) {
        string gtfsFolder;
        cout << "GTFS folder path (e.g.: ./gtfs): ";
        cin >> gtfsFolder;
        if (!loadGTFSData(gtfsFolder, agencies, routes, trips, stops, stopTimes, transfers)) {
            cout << "Error loading GTFS, falling back to generated data.\n";
            // For fallback in interactive use, use much smaller sizes for quick runs
            generateData(agencies, routes, trips, stops, stopTimes, transfers, 2, 10, 50, 30, 40);
        }
    } else {
        // Default: generate small synthetic dataset suitable for quick tests
        generateData(agencies, routes, trips, stops, stopTimes, transfers, 2, 10, 50, 30, 40);
    }

    // Build an index of stops by id for readable output
    unordered_map<string, Stop> stopById;
    for (const auto& s : stops) stopById[s.id] = s;

    // Build the time-expanded graph used by Dijkstra and A*
    unordered_map<Node, vector<pair<Node, int>>, NodeHash> graph;
    buildGraph(trips, stopTimes, transfers, graph);

    // Option to force long-distance origin/destination (disabled by default)
    bool longDistance = false; // set to true if you want start = first stop and end = last stop

    string startStopId, endStopId;
    if (longDistance) {
        startStopId = stops.front().id;
        endStopId = stops.back().id;
    } else {
        // Choose two distinct random stops as origin and destination
        random_device rd; mt19937 gen(rd());
        uniform_int_distribution<> stopDist(0, stops.size() - 1);
        int startStopIdx = stopDist(gen), endStopIdx = stopDist(gen);
        while (endStopIdx == startStopIdx) endStopIdx = stopDist(gen);
        startStopId = stops[startStopIdx].id;
        endStopId = stops[endStopIdx].id;
    }

    // Collect all event nodes corresponding to departure events at the chosen start stop
    vector<Node> startNodes;
    for (const auto& st : stopTimes)
        if (st.stop_id == startStopId)
            startNodes.push_back({st.stop_id, st.departure_time});

    // Run Dijkstra and A* and measure execution time (microseconds)
    vector<Node> pathD, pathA;
    auto t1 = chrono::high_resolution_clock::now();
    int costD = dijkstra(graph, startNodes, endStopId, pathD);
    auto t2 = chrono::high_resolution_clock::now();
    int costA = astar(graph, startNodes, endStopId, pathA);
    auto t3 = chrono::high_resolution_clock::now();

    auto dijkstraTime = chrono::duration_cast<chrono::microseconds>(t2 - t1).count();
    auto astarTime = chrono::duration_cast<chrono::microseconds>(t3 - t2).count();

    cout << "Origin: " << stopById[startStopId].name << ", Destination: " << stopById[endStopId].name << endl;
    cout << "-----------------------------" << endl;
    cout << "Dijkstra: ";
    if (costD == -1) cout << "No path found.\n";
    else {
        for (auto& n : pathD) {
            int h = n.time / 60, m = n.time % 60;
            cout << stopById[n.stop_id].name << " (" << (h < 10 ? "0" : "") << h << ":" << (m < 10 ? "0" : "") << m << ") ";
        }
        cout << "\nTotal duration: " << costD << " min, Stops: " << pathD.size()
             << ", Algorithm time: " << dijkstraTime << " us\n";
    }
    cout << "-----------------------------" << endl;
    cout << "A*: ";
    if (costA == -1) cout << "No path found.\n";
    else {
        for (auto& n : pathA) {
            int h = n.time / 60, m = n.time % 60;
            cout << stopById[n.stop_id].name << " (" << (h < 10 ? "0" : "") << h << ":" << (m < 10 ? "0" : "") << m << ") ";
        }
        cout << "\nTotal duration: " << costA << " min, Stops: " << pathA.size()
             << ", Algorithm time: " << astarTime << " us\n";
    }

    // Dijkstra on an independently built time-expanded graph (demonstration / alternative)
    vector<pair<string, int>> pathTE;
    auto t4 = chrono::high_resolution_clock::now();
    int costTE = dijkstra_time_expanded(stops, trips, stopTimes, transfers, startStopId, 6*60, endStopId, pathTE);
    auto t5 = chrono::high_resolution_clock::now();
    auto teTime = chrono::duration_cast<chrono::microseconds>(t5 - t4).count();

    cout << "-----------------------------" << endl;
    cout << "Dijkstra Time-Expanded: ";
    if (costTE == -1) cout << "No path found.\n";
    else {
        for (auto& [sid, t] : pathTE) {
            int h = t / 60, m = t % 60;
            cout << stopById[sid].name << " (" << (h < 10 ? "0" : "") << h << ":" << (m < 10 ? "0" : "") << m << ") ";
        }
        cout << "\nTotal duration: " << costTE << " min, Stops : " << pathTE.size()
             << ", Algorithm time: " << teTime << " us\n";
    }
    return 0;
}