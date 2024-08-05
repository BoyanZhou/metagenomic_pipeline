"""
This script trims the 3' end of reads, given the directory of sequencing files (XXX_R1.gz, XXX_R2.gz)
module add trimgalore/0.6.6 (include python/cpu/3.6.5)

command line example:
python step1_trim_galore_20240715.py input_directory_path output_directory_path
"""


import os
import pipeline_metagenomic as pipe_meta
import sys


def build_sample_dict_under_dir(data_dir):
    # all paired gz files, data_dir = "XXX/raw_data"
    all_gz_list = [i for i in os.listdir(data_dir) if i.endswith(".gz")]
    all_gz_list = [i for i in all_gz_list if "combined" in i]   # only used combined
    all_gz_list = sorted(all_gz_list)
    sample_name_fq_dict = {}
    # must be sorted file
    for gz_file in sorted(all_gz_list):
        sample_name = gz_file.split("_")[0]
        if sample_name not in sample_name_fq_dict:
            sample_name_fq_dict.update({sample_name: [os.path.join(data_dir, gz_file)]})
        else:
            sample_name_fq_dict[sample_name].append(os.path.join(data_dir, gz_file))
    for sample_name, fqs in sample_name_fq_dict.items():
        print(f"The number of fqs for {sample_name} is {len(fqs)}")
    return sample_name_fq_dict


if __name__ == "__main__":
    gz_file_dir, output_dir = sys.argv[1:3]
    sample_fq_dict = build_sample_dict_under_dir(gz_file_dir)  # {"sample1": ["path/XXX_R1.fq.gz", "path/XXX_R2.fq.gz"]}
    for sample, fq_path_list in sample_fq_dict.items():
        pipe_meta.trim_galore_process_step1(output_dir, fq_path_list, quality_threshold=25, thread=8)
