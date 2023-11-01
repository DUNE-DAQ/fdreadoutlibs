# TPG Applications
Here is a short summary of the scripts available in this directory. Refer to the code for further details. 

* `plot_trigprim_output_data.py.py` plots the Trigger Primitive output file obtained through `wibeth_tpg_algorithms_emulator` (when `save_trigprim` flag is enabled) and produces a plot called `output_trigger_primitives.png` . This script requires the use of `matplotlib`. To use the script run the following command: 
```sh
python3 plot_trigprim_output_data.py  -f TP_OUTPUT.TXT
```

#### Setup matplotlib on NP04 machines (e.g. `np04-srv-019`)
To use the `matplotlib` python module run the following command on a console where the DUNE-DAQ software area has not been sourced:
```sh
pip install --prefix=$PREFIX_PATH matplotlib
```



