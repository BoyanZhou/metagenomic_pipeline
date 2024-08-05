"""
This script generates the functional profiling by humann3

module add python/cpu/3.7.2
"""


import os
import pipeline_metagenomic as pipe_meta
import sys


def build_gz_sample_dict_under_kneaddata(data_kneaddata_dir):
    """
    Currently kneaddata/fastq has been gzipped
    :param data_kneaddata_dir:
    :return:
    """
    all_fq_list = [i for i in os.listdir(data_kneaddata_dir) if i.endswith("kneaddata_paired_1.fastq.gz") or i.endswith("kneaddata_paired_2.fastq.gz")]
    sample_name_fq_dict = {}
    # must be sorted file
    for gz_file in sorted(all_fq_list):
        sample_name = gz_file.split("_")[0]
        if sample_name not in sample_name_fq_dict:
            sample_name_fq_dict.update({sample_name: [os.path.join(data_kneaddata_dir, gz_file)]})
        else:
            sample_name_fq_dict[sample_name].append(os.path.join(data_kneaddata_dir, gz_file))
    for sample_name, fqs in sample_name_fq_dict.items():
        print(f"The number of fqs for {sample_name} is {len(fqs)}")
    return sample_name_fq_dict


# main function: humann3
def humann3_process(sample_fq_gz_dict, total_output_dir, metaphlan_dir, metaphlan_index, humann3_nucleotide_db, humann3_protein_db):
    # create folder under the total output folder for each sample
    for sample_id in list(sample_fq_gz_dict.keys()):
        sample_dir = os.path.join(total_output_dir, sample_id)
        if not os.path.exists(sample_dir):
            os.system(f"mkdir -p {sample_dir}")
        else:
            # sample result directory exist, do not over write
            print(f"Sample {sample_id} has been processed. Skip this one!")
            continue
        # combine previous pair reads
        os.chdir(sample_dir)
        fq_path_list = sample_fq_gz_dict[sample_id]     # [1.fastq.gz, 2.fastq.gz]
        combined_fq_path = os.path.join(sample_dir, f"{sample_id}_kneaddata_combined.fastq.gz")
        os.system(f"cat {fq_path_list[0]} {fq_path_list[1]} > {combined_fq_path}")
        pipe_meta.humann_analysis_step5(combined_fq_path, sample_dir, sample_id, metaphlan_dir, metaphlan_index, humann3_nucleotide_db, humann3_protein_db)
        os.system(f"rm {combined_fq_path}")


if __name__ == "__main__":
    kneaddata_dir_as_input, humann3_dir_as_output, metaphlan_dir, metaphlan_index, humann3_nucleotide_db, humann3_protein_db = sys.argv[1:7]
    """
    kneaddata_dir = "user_path/processed_data/kneaddata"
    humann3_output_dir = "user_path/processed_data/humann3_result"
    humann3_nucleotide_db_dir = "user_path/software/HUMAnN3_database/chocophlan"
    humann3_protein_db_dir = "user_path/software/HUMAnN3_database/uniref"
    """
    # build dict of sample ID and fastq files
    sample_fq_gz_dict = build_gz_sample_dict_under_kneaddata(kneaddata_dir_as_input)  # {"sample1": [fastq.gz, fastq.gz]}
    humann3_process(sample_fq_gz_dict, humann3_dir_as_output, metaphlan_dir, metaphlan_index, humann3_nucleotide_db, humann3_protein_db)

