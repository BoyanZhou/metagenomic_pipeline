"""
This script removes potential host reads, given the directory of sequencing files after Trim Galore

module add python/cpu/3.7.2
"""


import os
import pipeline_metagenomic as pipe_meta
import sys


def build_sample_dict_under_dir(data_after_trim_galore_dir):
    """
    The input directory is the output directory of trim_galore in Step1
    :param data_after_trim_galore_dir:
    :return:
    """
    all_gz_list = [i for i in os.listdir(data_after_trim_galore_dir) if i.endswith(".gz")]
    sample_name_fq_dict = {}
    # must be sorted file
    for gz_file in sorted(all_gz_list):
        sample_name = gz_file.split("_")[0]
        if sample_name not in sample_name_fq_dict:
            sample_name_fq_dict.update({sample_name: [os.path.join(data_after_trim_galore_dir, gz_file)]})
        else:
            sample_name_fq_dict[sample_name].append(os.path.join(data_after_trim_galore_dir, gz_file))
    for sample_name, fqs in sample_name_fq_dict.items():
        print(f"The number of fqs for {sample_name} is {len(fqs)}")
    return sample_name_fq_dict


# main function: Kneaddata
def knead_process(sample_fq_dict, output_knead_dir, bowtie_ref):
    """
    :param sample_fq_dict: result of build_sample_dict_under_dir()
    :param output_knead_dir: output directory
    :param bowtie_ref: the bowtie ref (built from genome .fas) of host (e.g. human, mouse)
    :return:
    """
    for sample_id in sample_fq_dict.keys():
        fq_path_list = sample_fq_dict[sample_id]
        pipe_meta.kneaddata_clean_step3(output_knead_dir, fq_path_list, bowtie_ref, thread=4)


if __name__ == "__main__":
    """
    # fq_data_dir = "user/raw_data"
    # trim_galore_dir = "user/trim_galore"
    # kneaddata_dir = "user/kneaddata"
    # Example of host genome
    mouse_bowtie_ref = "user/software/Kneaddata/contamination_database/mouse_C57BL_6NJ"
    human_bowtie_ref = "user/software/human_database/ftp.broadinstitute.org/bundle/hg38/hg38"
    """
    trim_galore_dir, knead_dir, host_bowtie_ref = sys.argv[1:4]
    sample_name_fq_dict1 = build_sample_dict_under_dir(trim_galore_dir)
    knead_process(sample_name_fq_dict1, knead_dir, host_bowtie_ref)
