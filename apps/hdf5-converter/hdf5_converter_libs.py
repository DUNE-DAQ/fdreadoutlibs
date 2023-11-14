'''
Description:
This file contains the function that converts a HDF5 file containing TPStream data from DUNE DAQ to different formats.

The function is called by main.py.

The function tpstream_hdf5_converter returns the tps ordered by time_start.

The function tpstream_hdf5_converter_as_cpp returns the tps ordered by fragment, as in the cpp implementation.
'''

import daqdataformats
import detdataformats
import fddetdataformats
import trgdataformats
import detchannelmaps
from hdf5libs import HDF5RawDataFile

import argparse
import sys
import numpy as np

def tp_to_numpy(tp):
    return np.array([tp.time_start, tp.time_over_threshold, tp.time_peak, tp.channel, tp.adc_integral, tp.adc_peak, tp.detid, tp.flag, tp.version])


def tpstream_hdf5_converter(filename, num_records=-1, out_format=["txt"]):
    h5_file = HDF5RawDataFile(filename)

    # Get all records (TimeSlice)
    records = h5_file.get_all_record_ids()
    # Set the number of records to process
    print(f'Input file: {filename}')
    print(f'Total number of records: {len(records)}')
    if num_records == -1:
        num_records = len(records)
        print(f'Number of records to process: {num_records}')
    elif num_records > len(records):
        num_records = len(records)
        print(f'Number of records to process is greater than the number of records in the file. Setting number of records to process to: {num_records}')
    else:
        print(f'Number of records to process: {num_records}')

    all_tps = []

    # Get the size of a single TP in memory, useful to iterate over the TPs in the fragments
    one_tp_size_in_memory = trgdataformats.TriggerPrimitive.sizeof()
    for record in records[:num_records]:
        print(f'Processing record {record}')
        # Get all data sources
        datasets = h5_file.get_fragment_dataset_paths(record)
        n_tps = []
        n_tps_left = []
        next_tp_in_datasets = []
        fragments = []
        for dataset in datasets:
            # Get the fragment
            frag = h5_file.get_frag(dataset)
            fragments.append(frag)
            # Get the number of TPs
            n_tps_left.append(int(frag.get_data_size()/trgdataformats.TriggerPrimitive.sizeof()))
            n_tps.append(int(frag.get_data_size()/trgdataformats.TriggerPrimitive.sizeof()))
            # Get the data
            next_tp_in_datasets.append(trgdataformats.TriggerPrimitive(frag.get_data(0)).time_start)
        while len(n_tps_left) > 0:
            index = np.argmin(next_tp_in_datasets)
            tp = trgdataformats.TriggerPrimitive(fragments[index].get_data(one_tp_size_in_memory*(n_tps[index]-n_tps_left[index])))
            all_tps.append(tp_to_numpy(tp))
            n_tps_left[index] -= 1
            if n_tps_left[index] == 0:
                del n_tps_left[index]
                del next_tp_in_datasets[index]
                del n_tps[index]
                del fragments[index]
            else:
                next_tp_in_datasets[index] = trgdataformats.TriggerPrimitive(fragments[index].get_data(one_tp_size_in_memory*(n_tps[index]-n_tps_left[index]))).time_start
                


    all_tps = np.array(all_tps)
    print(f"Final shape: {all_tps.shape}")

    return all_tps

def tpstream_hdf5_converter_as_cpp(filename, num_records=-1, out_format=["txt"]):

    h5_file = HDF5RawDataFile(filename)

    # Get all records (TimeSlice)
    records = h5_file.get_all_record_ids()
    # Set the number of records to process
    print(f'Input file: {filename}')
    print(f'Total number of records: {len(records)}')
    if num_records == -1:
        num_records = len(records)
        print(f'Number of records to process: {num_records}')
    elif num_records > len(records):
        num_records = len(records)
        print(f'Number of records to process is greater than the number of records in the file. Setting number of records to process to: {num_records}')
    else:
        print(f'Number of records to process: {num_records}')

    all_tps = []

    # Get the size of a single TP in memory, useful to iterate over the TPs in the fragments
    one_tp_size_in_memory = trgdataformats.TriggerPrimitive.sizeof()
    for record in records[:num_records]:
        print(f'Processing record {record}')
        # Get all data sources
        datasets = h5_file.get_fragment_dataset_paths(record)

        for dataset in datasets:
            # Get the fragment
            frag = h5_file.get_frag(dataset)
            # Get the number of TPs
            n_tps=int(frag.get_data_size()/trgdataformats.TriggerPrimitive.sizeof())

            for i in range(n_tps):
                tp = trgdataformats.TriggerPrimitive(frag.get_data(one_tp_size_in_memory*i))
                # print_tp(tp)
                tp_np = tp_to_numpy(tp)
                all_tps.append(tp_np)



    all_tps = np.array(all_tps)
    print(f"Final shape: {all_tps.shape}")

    return all_tps


