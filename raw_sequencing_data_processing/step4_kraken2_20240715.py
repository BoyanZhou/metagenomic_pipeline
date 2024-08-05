"""
This script generates the taxonomic profiling by Kraken2: can use either NCBI or GTDB database

module add kraken/2.0.8
module add python/cpu/3.6.5
"""


import os
import pipeline_metagenomic as pipe_meta
import sys


def build_sample_dict_under_dir(data_after_knead_dir):
    """
    The input directory is the output directory of Kneaddata in Step3
    :param data_after_knead_dir:
    :return:
    """
    # all fq file, data_dir = "/gpfs/data/lilab/home/zhoub03/Lama_oxalate/Lama_three_studies/processed_data/kneaddata"
    all_fq_list = [i for i in os.listdir(data_after_knead_dir) if i.endswith("kneaddata_paired_1.fastq.gz") or i.endswith("kneaddata_paired_2.fastq.gz")]
    sample_name_fq_dict = {}
    # must be sorted file
    for gz_file in sorted(all_fq_list):
        sample_name = gz_file.split("_")[0]
        if sample_name not in sample_name_fq_dict:
            sample_name_fq_dict.update({sample_name: [os.path.join(data_after_knead_dir, gz_file)]})
        else:
            sample_name_fq_dict[sample_name].append(os.path.join(data_after_knead_dir, gz_file))
    for sample_name, fqs in sample_name_fq_dict.items():
        print(f"The number of fqs for {sample_name} is {len(fqs)}")
    return sample_name_fq_dict


# main function: Kraken2 process
def kraken2_process(sample_fq_dict, kraken2_output_dir, kraken2_db_dir):
    sample_id_list = list(sample_fq_dict.keys())
    all_sample_fq_path_list = [sample_fq_dict[i] for i in sample_id_list]
    sample_output_folder_list = [kraken2_output_dir for i in sample_id_list]    # all sample to the same folder
    print(f"sample_id_list {sample_id_list}")
    print(f"all_sample_fq_path_list {all_sample_fq_path_list}")
    print(f"sample_output_folder_list {sample_output_folder_list}")
    print(f"kraken2_db_dir {kraken2_db_dir}")
    pipe_meta.kraken_process_step4(sample_id_list, all_sample_fq_path_list, sample_output_folder_list,
                                   kraken2_db_dir, whether_use_mpa_style=True)


if __name__ == "__main__":
    """
    # fq_data_dir = "user/raw_data"
    # trim_galore_dir = "user/trim_galore"
    # kneaddata_dir = "user/kneaddata"
    # Example of host genome
    mouse_bowtie_ref = "user/software/Kneaddata/contamination_database/mouse_C57BL_6NJ"
    human_bowtie_ref = "user/software/human_database/ftp.broadinstitute.org/bundle/hg38/hg38"
    # Example of Kraken2 database
    kraken2_gtdb_dir = "user/software/kraken2_GTDB_download"    # 150G database
    kraken2_ncbi_dir = "user/software/kraken2/NCBI_standard"
    """
    kneaddata_dir_as_input, kraken2_dir_as_output, kraken2_database_dir = sys.argv[1:4]
    # build dict of sample ID and fastq files
    sample_fq_dict = build_sample_dict_under_dir(kneaddata_dir_as_input)  # {"sample1": ["path/XXX_R1.fastq.gz", "path/XXX_R2.fastq.gz"]}

    kraken2_process(sample_fq_dict, kraken2_dir_as_output, kraken2_database_dir)
