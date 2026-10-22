# Compact PDS descriptors after latency-buffer insertion

This branch targets fddaq-v5.7.1-a9-1. PDS follows the charge processor's
TaskRawDataProcessorModel architecture: timestamp checks run before insertion;
extract_tps is registered as a post-processing task before base configuration
creates its pointer queue and dedicated worker. DPDK forwards complete frames
and has no descriptor-only queue, worker pool, or new configuration schema.

PDS retains its SkipList buffer because channels can arrive out of timestamp
order. The immediate-insertion pointer fix from datahandlinglibs PR 106 is
required (version 3.2.4 onward; 5.7.1 supplies 3.2.7). This workarea additionally
patches stop ordering in datahandlinglibs: queued processing must drain before
SkipList flush releases frame storage. The regression test holds a worker while
three out-of-order frames arrive, verifies one TP per original frame on another
thread, and checks storage survives until the worker drains.

The companion fddetdataformats branch carries the existing firmware's compact
format-4 contract: 512 bytes, 256 samples, and up to five peak descriptors.
It retains the release's FrameConcept and modern Python bindings. Rebuild
fdreadoutmodules and dpdklibs with this adapter; their source has no TPG worker
changes. This format remains a development wire-format extension, not an
upstream 5.7.1 default. Geometry/channel mapping, masks, integrated-ADC cut and
native TriggerPrimitiveVector output remain in fdreadoutlibs. Malformed or
overflowed descriptors are rejected, with a rate-limited warning.

## Verification

The workarea's checks/post-insertion directory contains build logs and package
CTest logs: 9 fddetdataformats, 5 fdreadoutlibs, 4 datahandlinglibs and 5 dpdklibs
tests passed. DAPHNE-14 configured and started run 45433, but the endpoint was
in state 6 with invalid timestamps; no hardware data qualification is claimed.
The alternative M2-004 configuration deliberately uses a plain DaqApplication
containing only HermesModule. It omits DaphneV3ControllerModule, whose normal
configuration writes analog bias settings. Timing/endpoint operations and
network configuration can therefore be tested without sending a bias request.

On 2026-10-22, the operator requested a time update through the coldbox TLU.
BOREAS_TLU_NP04 was synced from host UTC with kSoftware, and GIB was reloaded
from the upstream TLU with kUpstream. Numeric timestamp comparison verified
both within about 1 ms of host UTC. M2-004 still reported endpoint state 6;
that remains a blocker for valid fresh acquisition.
