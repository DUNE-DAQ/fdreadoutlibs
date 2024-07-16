# fdreadoutlibs - Far Detector readout libraries
Collection of Far Detector FrontEnd specific readout specializations. This includes type definitions to be used with the implementations in `fdreadoutmodules` and frontend specific specializations (i.e. frame processors or software hit finding). 

## Building and setting up the workarea

How to clone and build DUNE DAQ packages, including this one, is covered in [the daq-buildtools instructions](https://dune-daq-sw.readthedocs.io/en/latest/packages/daq-buildtools/).

## Frontends and features provided by `fdreadoutlibs`
The following frontends and features are provided by this package:
* *Daphne*: Frame processor and request handler to be used with the SkipList latency buffer.
* *SSP*: Only frame processor
* *WIB*: Frame processor for `WIB`, software and hardware `WIB` TPs. Implementation of avx based software hit finding (software tpg) for `WIB`
* *WIB2*: Frame processor for `WIB2`, implementation of avx based software hit finding (software tpg) for `WIB2`
* *WIBETH*: Frame processor for `WIBETH`, implementation of avx based software hit finding (software tpg) for `WIBETH`
