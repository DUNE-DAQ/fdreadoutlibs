# DAPHNE scan checks

`rate_producer.c` is the board-local paced MMIO producer used in the scans.
Arguments are frequency in Hz and maximum writes; signals stop it cleanly and
an alarm bounds execution to 120 seconds. It requires gateware `14f56c3`,
software-only mode, continuation disabled, and `/dev/mem` access.

`board_counters.py` saves coherent snapshots of the 32-channel counters.
`tp_ta_replay.cxx` reads native TP-v2 records, sorts them, and runs the existing
six-channel/10-us coincidence and prompt-light builders. Arguments are input
binary and output TA JSONL paths; stdout reports timing and counts.

The scan configurations and evidence are recorded in
[the throughput report](../../docs/pds-throughput-check.md).
