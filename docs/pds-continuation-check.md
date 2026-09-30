# DAPHNE continuation check, 2026-09-30

Board: `np04-daphne-015`, gateware `14f56c3`, ABI `0x20000`.
Software tested: `fdreadoutlibs` commit `75e5311` in the development area
specified in [the scaffold guide](pds-descriptor-trigger-scaffold.md).
Only physical fiber 0, carrying hardware channels 0–15, was tested.

## Working capture configuration

- ADC register 3: `0x2000` (16-bit serialization).
- ADC register 4: `0x0008` (14-bit resolution, offset binary, LSB first).
- Negative-pulse polarity; invert-enable mask zero.
- Software trigger selector `0`; correlation threshold remains 16.
- Connected-channel continuation register: `0x80200040` (enabled, 32 quiet
  samples, residual/activity threshold 64 ADC counts).
- Software triggers at 100 Hz; receiver running before link enable.

The initial ADC register 4 was `0x0018`, selecting MSB first. The deployed
receiver implementation (`ip_repo/daphne_ip/rtl/frontend/febit3.vhd`) requires
LSB first. Clearing bit 4 produced ordinary baseline values near 8000 with
small sample fluctuations, replacing the apparent large, repeated excursions.
The bit definition is documented in the
[AFE5808A register map](https://www.ti.com/lit/ds/symlink/afe5808a.pdf).
Changing descriptor thresholds alone would have hidden this configuration error.

The first 20-trigger capture contained 960 frames: 320 complete channel-event
chains of three 256-sample fragments, with consecutive timestamps separated by
256 ticks. All descriptors and overflow bits matched independent reconstruction
from the packed ADC samples. There were no overflows, busy/full counter increases
on connected channels, or peaks above threshold 64. Strict native replay passed
and emitted zero TPs/TAs, as expected from those descriptors.

The final 100-trigger confirmation produced 4,800 frames: 1,600 complete
three-fragment channel-event chains, including 3,200 continuation frames.
All timestamps were contiguous, all waveform/descriptor checks passed, and
there were zero overflows or connected-channel busy/full counter increases.
The receiver reported zero missed packets and input errors. Strict replay
passed with zero TPs/TAs. PCAP SHA-256:
`322bab1f9a145bedcb6cdae86993ebc2681aa4db063d9525ad0843456cfee6a3`.

An exploratory threshold of 8 produced 966 frames, including chains of four and
five fragments. Every chain remained contiguous and all descriptor fields
matched the samples. However, two frames overflowed; strict replay rejected them.
The explicit incomplete-input replay produced 36 TPs and no six-channel TAs.
Threshold 8 is therefore **not** a qualified setting from this test.

These are bounded software-trigger capture tests and offline native TP/TA replay.
They do not validate a light-pulse response, six-channel coincidence efficiency,
calibrated prompt light, autonomous self-trigger selection, or a live DAQ TA
module. A controlled optical/electrical stimulus is needed for those checks.

## Evidence and restoration

Evidence is in the development area's `checks/continuation-lsb-64` and
`checks/continuation-lsb-8`, with the final confirmation in
`checks/continuation-lsb-64-confirm`: raw PCAP, concatenated native frames, register
snapshots, trigger reports, channel map, replay output, and `result.json`.
The single-board diagnostic map is `slot*100 + hardware_channel`.

The first capture harness incorrectly assumed the records counter counted
initial triggers; it includes continuation fragments. The raw capture was
recovered and independently audited. Subsequent captures check frame counters
and validate the number of initial events and complete chains from frame data.

Each test restores the original ADC and continuation settings and disables both
links. The tested LSB-first correction must be made in the acquisition
configuration before reuse; restoring the original state also restores its
MSB-first mismatch. No detector configuration database was changed.
