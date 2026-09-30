# PDS descriptor trigger scaffold

Development area on `np04-srv-017`:
`/nfs/home/marroyav/workareas/daq/daphne/pds-trigger-dev-20260930`.

The area was created with DBT v8.14.0 against `fddaq-v5.6.2-a9-1` and copied from
the adjacent `fddaq-v5.6.2-a9-1` area. Eight repositories were compared by SHA-256
before editing. `provenance/source-manifest.json` records exact source commits and
file hashes; `provenance/source-recipe.yaml` preserves the original recipe.
Build and install directories are fresh. Python uses the release's shared venv
(`dbt-create -q`); no local Python development packages were installed.

Changes are confined to `fdreadoutlibs`, on branch
`marroyav/pds-descriptor-trigger-scaffold`. No DAQ session or hardware was started.

## Components

* `include/fdreadoutlibs/pds/DescriptorToTP.hpp`: format-4 descriptor conversion
  into this release's native `trgdataformats::TriggerPrimitive` (version 2).
* `src/daphneeth/DAPHNEEthFrameProcessor.cpp`: restores the previously commented
  TP extraction path. Existing post-processing enablement gates registration.
  A configured PDS channel map and `TriggerPrimitiveVector` output are required.
  Existing offline channel masks and minimum integral cut apply. Invalid frames
  produce an ERS warning, limited to one per monitoring interval; the rejection
  count is logged at bookkeeping debug level. Queue send failures use existing
  failed-TP monitoring.
* `include/fdreadoutlibs/pds/ActivityBuilder.hpp`: independently usable, streaming
  coincidence and light-window builders. These are not yet a registered
  `TriggerActivityMaker` plugin or DAQ activity module.
* `test/apps/pds_descriptor_replay.cxx`: bounded finite-file replay producing
  JSON lines for TPs, PDS TAs, and descriptor-light summaries.
* `unittest/PDSDescriptorToTP_test.cxx`: conversion and malformed-frame tests.
* `unittest/PDSDescriptorPipeline_test.cxx`: window boundaries, invalid input,
  configuration, time ordering, capacity, and end-to-end synthetic tests.

## Data and timing contract

The pinned `DAPHNEEthFrame.hpp` represents a **512-byte DAQ frame**: 16-byte
`DAQEthHeader`, 8-byte common DAPHNE header, five 8-byte descriptor slots, and
448 waveform bytes. Firmware's 504-byte payload replaces the 16-byte DAQ header
with an 8-byte timestamp. Replay accepts only concatenated 512-byte DAQ frames,
not firmware-only payloads, PCAP packets, or HDF5 containers. Use the DAQ format
definition from the pinned source; do not assume the two layouts are interchangeable.

Descriptor `sample_start` is relative to frame sample zero; `time_peak` is
relative to the descriptor start. Duration is `duration_minus_one + 1`, including
256-sample excursions. Integral and peak are already positive, baseline-subtracted
quantities. Conversion does not subtract baseline a second time.

The current hardware/readout convention uses a 62.5 MHz timestamp clock and one
timestamp tick per ADC sample. The standalone converter/builders expose
`ticks_per_sample` for other integral ratios; the DAQ hook and replay deliberately
use the verified current convention. Microsecond windows round **up** to whole
ticks: 10 us = 625 ticks, 0.1 us = 7 ticks = 112 ns. Defaults are illustrative,
not a physics tuning recommendation. Clock changes require coordinated changes
to conversion, window configuration, and downstream consumers.

Offline channel IDs must be globally unique within an activity region. Mapping
uses detector, crate, slot, stream, and hardware channel; board-local channel
numbers alone are insufficient. Replay requires an explicit CSV map. Its demo
uses synthetic single-board channels only.

Absent descriptors must be zero; found descriptors must have a positive peak.
Wrong versions, non-fragment-local format,
invalid bounds, impossible integrals, absent-slot gaps, overlapping descriptors,
invalid channel maps, and timestamp wrap are rejected. A frame is validated
atomically before any of its TPs can be emitted. Descriptor-overflow frames are
rejected by default because they contain incomplete light information. The
standalone converter can explicitly accept them, but then the caller owns the
incompleteness accounting. Continuation fragments are accepted: a boundary-crossing
excursion produces fragment-local TPs, with no cross-fragment pulse stitching.

Replay provides `--allow-descriptor-overflow` for diagnostic inspection of retained
peaks. If any such frames occur, **every output record** is conservatively marked
`incomplete_input: true`; stderr reports their count. This does not recover missing
peaks. The DAQ TP hook retains the rejection policy.

A read-only audit of the existing `test-256-20260929/board-rate-fiber.pcap` found
6,400 format-4 frames and 32,000 structurally valid descriptor slots; all frames
had descriptor overflow. Evidence and the capture hash are recorded in
`provenance/capture-descriptor-audit.json`. This is structural validation, not a
physics light-yield or live-TA validation.

## Coincidence policy

At TP time `t`, retain inputs in the inclusive sliding interval `[t - W, t]`.
Emit a native PDS TA when at least `minimum_channels` **distinct offline channels**
are present. The default is six, implementing “more than five channels”. Multiple
peaks on a channel contribute charge but count once toward channel multiplicity.
After emitting, consume all contributing inputs; another TA requires new inputs.
The triggering sixth TP closes the TA immediately. Later TPs with the same
timestamp can therefore belong to the next candidate; this is a threshold-crossing
policy, not a complete event clustering algorithm or a rising-edge holdoff policy.

TA inputs preserve native TPs. `time_start` is the first contributing start,
`time_end` is the exclusive end of the last-ending excursion, `time_activity` is
the threshold-crossing TP start, and `time_peak` belongs to the largest individual
peak (first on ties). `adc_integral` sums all contributing descriptor integrals;
`adc_peak` is the largest individual peak, not a simultaneous channel sum.
The release has no PDS coincidence algorithm enum, so native algorithm remains
`kUnknown`; replay labels it explicitly as a prototype. A production algorithm
ID must be agreed before downstream selection uses it.

## Prompt-light proxy

The first TP starts a non-overlapping total-light gate `[t0, t0 + total_ticks)`.
Sum whole descriptor integrals whose **start times** are in that gate. The prompt
sum uses `[t0, t0 + prompt_ticks)`. Report `prompt_integral / total_integral`, or
no fraction for zero total charge. A TP exactly on the total boundary starts a
new gate. This gate is seeded independently of the coincidence trigger.

This is a thresholded **descriptor-integral proxy**, in ADC-sample units. It is
not calibrated photoelectrons or exact prompt charge. A descriptor crossing a
gate boundary cannot be split from its integral alone; waveform samples are
needed for exact integration. Below-threshold light and dropped/overflowed frames
are absent. Continuation pieces are counted separately by their starts. No
uniform-charge interpolation or fabricated PE conversion is applied.

## Ordering and lifecycle

Both builders require an at-most-once stream ordered by TP start time across all
participating channels in one detector/region. Instantiate separate builders for
separate regions. They reject late TPs and mixed detector IDs. Network arrival
order and independently ordered per-link batches do not satisfy this contract.
Replay sorts a bounded finite file; it rejects repeated channel/start pairs.
It is not an online reorder buffer.

`advance(watermark)` asserts that all future TPs start at or after the watermark.
Advance using a common watermark from every participating source, including
quiet-source progress. This expires stale coincidence inputs and completes
quiet prompt windows. Reset by constructing new builders on a new run/timing
epoch. At end of a finite replay, `finish()` closes the last light gate at its
actual end, including gates near the timestamp limit; repeated calls produce
no duplicate output. Capacity and numeric overflow raise errors instead of truncating.
The reference coincidence implementation scans the retained window per TP;
it has not been throughput qualified.

## Build and run

```bash
ssh np04-srv-017
cd /nfs/home/marroyav/workareas/daq/daphne/pds-trigger-dev-20260930
source env.sh
dbt-build -j 4
ctest --test-dir build/fdreadoutlibs --output-on-failure --no-tests=error
./build/fdreadoutlibs/test/apps/pds_descriptor_replay --demo
python3 sourcecode/fdreadoutlibs/scripts/test_pds_descriptor_replay.py \
  build/fdreadoutlibs/test/apps/pds_descriptor_replay
```

For DAQ-frame replay, supply an explicit map with one row per hardware channel:

```text
# det,crate,slot,stream,hardware_channel,offline_channel
2,1,0,0,0,1000
2,1,0,0,1,1001
```

```bash
./build/fdreadoutlibs/test/apps/pds_descriptor_replay \
  --input frames.bin --channel-map channels.csv \
  --min-channels 6 --window-us 10 --prompt-us 0.1 --total-us 10 \
  --minimum-integral 0 --max-tps 1000000 > replay.jsonl
```

## Next integration steps

1. Validate with captured format-4 DAQ frames and the actual PDS offline map.
   Compare firmware descriptors against waveform-based integrals and peak timing.
2. Add a DAQ activity module/plugin and OKS configuration for region membership,
   windows, capacity, and the agreed trigger algorithm identifier.
3. Add the cross-link TP merger with bounded lateness, idle-source watermarks,
   duplicate suppression, loss accounting, and run-transition handling.
4. Publish descriptor rejection, overflow, late-input, occupancy, TA-rate, and
   light-summary monitoring; avoid per-frame warning floods under persistent faults.
5. Choose event/TA seeding, overlap/holdoff policy, calibration, masks, threshold
   treatment, and waveform boundary handling for a physics prompt-light observable.
6. Measure throughput and latency before enabling TPG and TA production in a run.

The TP output hook is implemented and buildable, but the complete live
descriptor-to-TA DAQ chain is not configured or commissioned by this scaffold.
