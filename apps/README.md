# TPG Applications
Here is a short summary of the applications available in this directory. Refer to the code for further details. 

`wibeth_tpg_algorithms_emulator` is a testing tool for validating different TPG algorithms, either in a naive or in AVX implementation. To use the tool use the following. Refer to the code for further details.

```sh
$ wibeth_tpg_algorithms_emualator --help 
Test TPG algorithms
Usage: wibeth_tpg_algorithms_emulator [OPTIONS]

Options:
  -h,--help                   Print this help message and exit
  -f,--frame_file_path TEXT   Path to the input frame file
  -a,--algorithm TEXT         TPG Algorithm (SimpleThreshold / AbsRS)
  -i,--implementation TEXT    TPG implementation (AVX / NAIVE)
  -d,--duration_test INT      Duration (in seconds) to run the test
  -n,--num_frames_to_read INT Number of frames to read. Default: select all frames.
  -t,--swtpg_threshold INT    Value of the TPG threshold
  --save_adc_data             Save ADC data
  --save_trigprim             Save trigger primitive data
```

The command line option `save_adc_data` allows to save the raw ADC values in a txt file after the 14-bit to 16-bit expansion. The command line option `save_trigprim`  allows to save the in a file the Trigger Primitive object information in a txt file. 

Example of usage: 
```sh
$ ./TestTPGAlgorithmsWIBEth --frame_file_path FRAMES_FILE --algorithm SimpleThreshold --implementation NAIVE --save_adc_data 1
$ ./TestTPGAlgorithmsWIBEth.cxx --frame_file_path FRAMES_FILE --algorithm SimpleThreshold --implementation AVX  --save_trigprim 1
```


`wibeth_binary_frame_reader`: reads a WIBEth frame file (`.bin` file) and prints all the ADC values. Usage `WIBEthBinwibeth_binary_frame_readeraryFrameReader <input_file_name>`.  

`wibeth_binary_frame_modifier` is used to create a custom WIBEth frame file suitable for testing different patterns. The application will produce an output file `wibeth_output.bin`. There are no command line options, please refer to the code for further details (e.g. what ADC value to set, which time frame to use, etc.). 


## Note
- The repository also contains tools for WIB2 frames but they are not kept up to date. Please refer to the code or ask mainteners of the repository for help. 
