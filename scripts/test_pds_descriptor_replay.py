#!/usr/bin/env python3
"""Independent synthetic wire-format acceptance tests for pds_descriptor_replay."""
import json
from pathlib import Path
import struct
import subprocess
import sys
import tempfile

executable = str(Path(sys.argv[1]).resolve())
checks = 0

with tempfile.TemporaryDirectory(prefix='pds-replay-') as temporary:
    root = Path(temporary)
    mapping = root / 'channels.csv'
    mapping.write_text(''.join(f'2,1,0,0,{ch},{1000 + ch}\n' for ch in range(6)))
    data = root / 'frames.bin'

    def frame(channel, timestamp=100000, start=0, version=4, overflow=False):
        daq_header = 1 | (2 << 6) | (1 << 12) | (63 << 52)
        metadata = (channel << 56) | (version << 52) | (1 << 50) | (int(overflow) << 49)
        peak = 250 | (100 << 22) | (3 << 36) | (1 << 44) | (start << 52) | (1 << 60)
        return struct.pack('<8Q', daq_header, timestamp, metadata, peak, 0, 0, 0, 0) + bytes(448)

    def run(blob, *options, success=True):
        global checks
        data.write_bytes(blob)
        result = subprocess.run([executable, '--input', str(data), '--channel-map', str(mapping), *options],
                                text=True, capture_output=True)
        if (result.returncode == 0) != success:
            raise AssertionError(f'Unexpected status {result.returncode}: {result.stderr}')
        checks += 1
        return [json.loads(line) for line in result.stdout.splitlines()] if success else []

    frames = [frame(ch, start=2 * ch) for ch in range(6)]
    records = run(b''.join(reversed(frames)))  # Ordering across input frames must be restored.
    tps = [r for r in records if r['kind'] == 'tp']
    activities = [r for r in records if r['kind'] == 'pds_ta']
    light = [r for r in records if r['kind'] == 'descriptor_light_proxy']
    assert [r['channel'] for r in tps] == list(range(1000, 1006))
    assert [r['start'] for r in tps] == list(range(100000, 100012, 2))
    assert len(activities) == 1 and activities[0]['distinct_channels'] == 6
    assert activities[0]['integral'] == 1500
    assert len(light) == 1 and light[0]['prompt_integral'] == 1000 and light[0]['total_integral'] == 1500
    assert abs(light[0]['fraction'] - 2 / 3) < 1.e-6

    assert not any(r['kind'] == 'pds_ta' for r in run(b''.join(frames[:5])))
    outside = b''.join(frames[:5]) + frame(5, timestamp=100626)
    assert not any(r['kind'] == 'pds_ta' for r in run(outside))
    # At the inclusive boundary all original five channels are still retained.
    boundary = b''.join(frames[:5]) + frame(5, timestamp=100625)
    assert sum(r['kind'] == 'pds_ta' for r in run(boundary)) == 1
    assert run(b'') == []
    maximum = (1 << 64) - 1
    near_wrap = run(frame(0, timestamp=maximum - 1500) + frame(1, timestamp=maximum - 700),
                    '--total-us', '16')
    final_light = [r for r in near_wrap if r['kind'] == 'descriptor_light_proxy']
    assert len(final_light) == 1 and final_light[0]['total_integral'] == 500
    assert final_light[0]['end_exclusive'] == maximum - 500
    run(frames[0][:-1], success=False)
    run(frames[0][8:], success=False)  # Firmware-only 504-byte payload is not a DAQ frame.
    run(frames[0] * 2, success=False)
    run(frame(0, version=3), success=False)
    run(frame(0, overflow=True), success=False)
    incomplete = run(frame(0, overflow=True), '--allow-descriptor-overflow')
    assert incomplete and all(r['incomplete_input'] for r in incomplete)
    run(b''.join(frames), '--max-tps', '5', success=False)
    run(b''.join(frames), '--min-channels', '0', success=False)
    mapping.write_text('2,1,0,0,0,1000\n')
    run(b''.join(frames), success=False)
    mapping.write_text('2,1,0,0,0,1000\n2,1,0,0,1,1000\n')
    run(frames[0], success=False)

print(f'{checks} binary replay acceptance checks passed')
