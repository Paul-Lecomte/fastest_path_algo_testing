# `algo_compare.cpp` - Dijkstra / A* comparison

Summary
-------
Small C++ program that compares three path-finding approaches on GTFS-like data (real or synthetic): Dijkstra on a time-event graph, A\* (simple heuristic), and a time-expanded Dijkstra.

Features
--------
- Basic GTFS file reading (`agency.txt`, `routes.txt`, `stops.txt`, `trips.txt`, `stop_times.txt`, optional `transfers.txt`) or synthetic data generation.
- Builds a time-expanded graph (stop/time events) with trip, transfer and waiting edges.
- Runs and compares:
    - Dijkstra (on event graph),
    - A\* (same graph, simple heuristic),
    - Dijkstra time-expanded (explicit event graph).
- Prints found path (stops + times), total duration and algorithm runtime in microseconds.

Build and run
-------------
- Requirement: modern C++ compiler (C++17).
- On Windows with MinGW/g++ (example):
    - `g++ -std=c++17 -O2 algo_compare.cpp -o algo_compare.exe`
    - `.\algo_compare.exe`
- In CLion: add `algo_compare.cpp` to a target, build and run.

Program behavior
----------------
- On start the program asks:
    - Mode: generated data (default) or load GTFS folder (`Path to GTFS folder`, e.g. `./gtfs`).
- Synthetic data is produced by `generateData(...)`. Modify sizes in source to test different loads.
- Option `longDistance` in `main()` forces start = first stop and end = last stop; otherwise stops are chosen randomly.
- A\* heuristic: numeric difference of stop IDs (e.g. `S123`). Only meaningful if IDs reflect spatial distance.

CSV parsing and formats
-----------------------
- CSV parser is minimal (simple quote handling); may fail for complex CSVs.
- `stop_times.txt` times are parsed to minutes; stop sequences and transfer times are handled with basic fallbacks.

Limitations and notes
---------------------
- Heuristic for A\* is very basic; replace by geographic distance for real GTFS.
- Memory usage can be large for big time-expanded graphs.
- Times are stored in minutes; algorithm runtimes are shown in microseconds.
- Parser and robustness are intentionally simple for demonstration.

Expected GTFS files
-------------------
- `agency.txt`, `routes.txt`, `stops.txt`, `trips.txt`, `stop_times.txt`
- Optional: `transfers.txt`

License
-------
Example code; reuse freely (MIT recommended).
