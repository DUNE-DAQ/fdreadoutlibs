`hdf5_converter.py` converts the HDF5 files from the DUNE DAQ to a more readable format.
You can change the output format by changing the value of the variable `FORMAT`. The possible values are:
  - `txt`: the output is a text file
  - `npy`: the output is a numpy array
  - `img_groups`: TPs are grouped, and an image is produced for each group
  - `img_all`: an image with all the tps is produced. 

You can simply run the script by using:
```sh
python hdf5_converter.py --input_file <input-file> --output_folder <output-folder> --drift_direction <0 or 1>  --make_fixed_size
```
More parameters are tunable, if needed. The main ones are:
* --num_records: you can set the number of records to convert
* --time_start/--time_end: you can set an interval save only tps with ```time_start``` within the interval.
* --drift_direction: you must set the correct drift direction, in order to have a correct division between planes. '0 for horizontal drift, 1 for vertical drift'.
* --format : you can choose the format. By default it will save the tps in all possible formats.

This script requires the use of `matplotlib` to produce the images.