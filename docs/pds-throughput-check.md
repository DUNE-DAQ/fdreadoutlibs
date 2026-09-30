# DAPHNE throughput checks

Area: `/nfs/home/marroyav/workareas/daq/daphne/pds-trigger-dev-20260930`.
Evidence: `checks/daq-bandwidth/`; gateware `14f56c3`.
Continuation disabled; software-only triggers; channels 0–15 connected.
DAQ uses independent 10 ms random windows at mean rate 10 Hz.
Each rate has a one-second settling period and a ten-second measurement.

## Raw acquisition

Run 45330 scanned 40–130 kHz in 10 kHz steps.
Connected-channel board counters showed no rejections through 120 kHz.
At 130 kHz, buffer-full rejections appeared; accepted frames plateaued near
124 kHz per channel. This counter test does not establish NIC losslessness.
[Raw counters](benchmarks/daphne/raw-scan.csv).

## Live descriptor TPs

Runs 45331 and 45332 enabled the descriptor TP hook with integral cut zero,
`SimplePDSChannelMap` (single-board test channels 100–115), and TA processing off.
Threshold 8 generated noise descriptors and descriptor overflows, which the TP
hook rejects. The original threshold 64 was restored after each test.
With noise descriptors, sampled NIC misses were 0% at 40/50 kHz, 13% at
60 kHz, 30% at 70 kHz, 16% at 80 kHz, and 62% at 120 kHz.
Quiet input at 120 kHz still showed 37% NIC misses. These short tests have
variable descriptor populations; they do not define a universal frequency limit.
TP-vector send failures were zero. Run 45331 persisted 3,854,763 TPs;
producer monitoring counted 3,871,986, leaving 17,223 unaccounted for in storage.
[TP and NIC counters](benchmarks/daphne/tp-scan.csv).

## TA replay

Sorted replay of the captured TPs produced 93,126 TAs using six distinct channels
within 10 us. Builder time was 0.764 s; sorting plus building took 0.914 s.
This is an offline replay, not a commissioned live TA chain or worst-case load test.
[Replay counters](benchmarks/daphne/ta-replay.json).

## Receive path

Captured UDP payloads contain 16 frames in 8192 bytes. The grouped 10 Gb/s link
ceiling is about 151 kHz for one frame per connected channel per trigger.
The receiver uses MTU 9000, ring 4096, and burst 2048. One board source IP maps
to one RX queue/core. Its frame callbacks invoke preprocessing synchronously;
the enabled TP extractor therefore runs on the receive thread.
[Grouping evidence](benchmarks/daphne/packet-grouping.json).
