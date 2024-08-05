import os
import subprocess
import numpy as np
import pandas as pd
import sys
import time
import logging
import logging.handlers
import datetime


"""
# metagenomic requires
module add kraken/2.0.8
module add python/cpu/3.6.5
module add bowtie2/2.3.5.1

update log:
date: 05/26/2023
change: add metaphlan database path and index to the function humann_analysis_step5
"""


#######################################################
#######################################################
# #                                                 # #
# #                  Common Functions               # #
# #                                                 # #
#######################################################
#######################################################
def sample_list_divide(sample_list, number_of_parts):
    """
    :param sample_list: [1,2,3,4,5]
    :param number_of_parts: 3
    :return:
    """
    part_length = int(len(sample_list)/number_of_parts)
    sample_divided_list = []
    part_start = 0
    for i in range(number_of_parts-1):
        sample_divided_list.append(sample_list[part_start: (part_start + part_length)])
        part_start += part_length
    sample_divided_list.append(sample_list[part_start:])
    return sample_divided_list


def generate_whole_pipeline_sbatch(path_sbatch):
    with open(path_sbatch, "w") as sbatch_f:
        sbatch_f.write(f"#!/bin/bash\n"
                       f"#SBATCH --partition=cpu_medium\n"
                       f"#SBATCH --job-name=metagenomic\n"
                       f"#SBATCH --mem-per-cpu=20G\n"
                       f"#SBATCH --time=72:00:00\n"
                       f"#SBATCH --tasks=1\n"
                       f"#SBATCH --cpus-per-task=4\n"
                       f"#SBATCH --nodes=1\n\n")


##########################################################################
##########################################################################
# #                                                                    # #
# #                  STEP 1 : QC and host genome removal               # #
# #                                                                    # #
##########################################################################
##########################################################################


def trim_galore_process_step1(output_dir, fq_path_list, quality_threshold, thread):
    """
    trim the low quality bases at 3' reads end of input fq, just paired-end reads of one sample
    :param output_dir:
    :param fq_path_list: [XXX_R1.fq, XXX_R2.fq]
    :param quality_threshold:
    :return:
    """
    if not os.path.exists(output_dir):
        os.system(f"mkdir -p {output_dir}")
    step1a_cmd = f"trim_galore --paired --gzip --output_dir {output_dir} --cores {thread} --clip_R1 1 --clip_R2 1 " \
                 f"--three_prime_clip_R1 1 --three_prime_clip_R2 1 --quality {quality_threshold} " \
                 f"{fq_path_list[0]} {fq_path_list[1]} --no_report_file"
    process = subprocess.Popen(step1a_cmd.split(), stdout=subprocess.PIPE)
    output, error = process.communicate()
    print(output, error)


def fastqc_report_step2(qc_output_dir, fq_path_list, thread):
    """
    optional step, generate quality report for fq file, not really process the fq
    :param qc_output_dir:
    :param fq_path_list:
    :param thread:
    :return:
    """
    if not os.path.exists(qc_output_dir):
        os.system(f"mkdir -p {qc_output_dir}")
    for fq_path in fq_path_list:
        step1b_cmd = f"fastqc {fq_path} -t {thread} -o {qc_output_dir}"
        os.system(step1b_cmd)


def kneaddata_clean_step3(clean_output_dir, fq_path_list, db_path, thread, remove_temp_files=True):
    """
    :param clean_output_dir: the output dir of knead data
    :param fq_path_list: fq_path_list = ["PATH/XXX_R1.fq.gz", "PATH/XXX_R2.fq.gz"]
    :param db_path: path of contamination bowtie2 database, path/ftp.broadinstitute.org/bundle/hg38/hg38
    :param thread:
    :param remove_temp_files: whether remove temp files (only reserve useful final output)
    :return: For pair end reads, the useful file names are XXX_R1_kneaddata_paired_1.fastq,
    XXX_R1_kneaddata_paired_2.fastq
    """
    if not os.path.exists(clean_output_dir):
        os.system(f"mkdir -p {clean_output_dir}")
    step1c_cmd = f"kneaddata --input {fq_path_list[0]} --input {fq_path_list[1]} " \
                 f"--reference-db {db_path} --output {clean_output_dir} --threads {thread}"
    process = subprocess.Popen(step1c_cmd.split(), stdout=subprocess.PIPE)
    output, error = process.communicate()
    print(output, error)
    # ------------------------------------------------------------------------------------------------------------------
    if remove_temp_files:
        # get the prefix of all files
        fq_name = os.path.split(fq_path_list[0])[-1]
        if fq_name.endswith("gz"):
            prefix_of_all_files = ".".join(fq_name.split(".")[:-2])
        else:
            # XXX.fq or XXX.fastq
            prefix_of_all_files = ".".join(fq_name.split(".")[:-1])
        print(prefix_of_all_files)
        os.chdir(clean_output_dir)
        all_output_files = [i for i in os.listdir(clean_output_dir) if i.startswith(prefix_of_all_files)]
        file_to_remove = []
        for file_name in all_output_files:
            if file_name.endswith("removed.unmatched.1.fastq") or file_name.endswith("removed.unmatched.2.fastq"):
                file_to_remove.append(file_name)
            if file_name.endswith("contam_1.fastq") or file_name.endswith("contam_2.fastq"):
                file_to_remove.append(file_name)
            if file_name.endswith("kneaddata.trimmed.1.fastq") or file_name.endswith("kneaddata.trimmed.2.fastq"):
                file_to_remove.append(file_name)
            if file_name.endswith("unmatched_1_contam.fastq") or file_name.endswith("unmatched_2_contam.fastq"):
                file_to_remove.append(file_name)
            if file_name.endswith("trimmed.single.1.fastq") or file_name.endswith("trimmed.single.2.fastq"):
                file_to_remove.append(file_name)
            if file_name.endswith("_unmatched_1.fastq") or file_name.endswith("_unmatched_2.fastq"):
                file_to_remove.append(file_name)
            if file_name.endswith("repeats.removed.1.fastq") or file_name.endswith("repeats.removed.2.fastq"):
                file_to_remove.append(file_name)
        print(prefix_of_all_files)
        for file_name in file_to_remove:
            os.system(f"rm {file_name}")


##########################################################################
##########################################################################
# #                                                                    # #
# #                     STEP 2 : Kraken2 profiling                     # #
# #                                                                    # #
##########################################################################
##########################################################################
def kraken_process_step4(sample_name_list, all_sample_fq_path_list, sample_output_folder_list, kraken_database_path,
                         whether_use_mpa_style=True):
    """
    Get processed fq and relative abundance report by Kraken
    :param sample_name_list: ["sample1", "sample2"]
    :param all_sample_fq_path_list: [["path/x.fq"], ["path/x_R1.fq", "path/x_R2.fq"]]
    :param sample_output_folder_list: ["output_path/subject1/sample1", "output_path/subject1/sample2"]
    :param kraken_database_path: "user_path/software/kraken2/NCBI_standard"
    :param whether_use_mpa_style: whether output metaphlan style output
    :return: sample_output.txt, sample_report.txt, for single,
    """

    use_mpa_style = ""
    if whether_use_mpa_style:
        use_mpa_style = "--use-mpa-style "
    for sample_index in range(len(sample_name_list)):
        sample_name = sample_name_list[sample_index]
        sample_fq_path_list = all_sample_fq_path_list[sample_index]
        sample_output_folder = sample_output_folder_list[sample_index]
        if not os.path.exists(sample_output_folder):
            os.system(f"mkdir -p {sample_output_folder}")

        # paired-end or single-end
        if len(sample_fq_path_list) == 1:
            # single
            command_kraken2 = f"kraken2 --db {kraken_database_path} --report {sample_output_folder}/{sample_name}_" \
                              f"report.txt {use_mpa_style}--use-names --report-zero-counts --classified-out " \
                              f"{sample_output_folder}/{sample_name}_#.fq {sample_fq_path_list[0]} " \
                              f"--output {sample_output_folder}/{sample_name}_output.txt"
        elif len(sample_fq_path_list) == 2:
            # paied
            command_kraken2 = f"kraken2 --db {kraken_database_path} --report {sample_output_folder}/{sample_name}_" \
                              f"report.txt {use_mpa_style}--use-names --report-zero-counts --paired --classified-out " \
                              f"{sample_output_folder}/{sample_name}_#.fq {sample_fq_path_list[0]} " \
                              f"{sample_fq_path_list[1]} --output {sample_output_folder}/{sample_name}_output.txt"
        else:
            print(f"Error! {sample_name} has incorrect number of fqs: {sample_fq_path_list} Terminated.")
            sys.exit()

        print(f"Processing {sample_name} by Kraken: {command_kraken2}")
        print(command_kraken2)
        os.system(command_kraken2)


##########################################################################
##########################################################################
# #                                                                    # #
# #                    STEP 3 : metaphlan3 and humann3                 # #
# #                                                                    # #
##########################################################################
##########################################################################
def humann_analysis_step5(input_fq_path, humann_output_dir, output_basename, metaphlan_dir, metaphlan_index, chocophlan_dir, uniref_dir, thread=16):
    """
    HUMAnN3 analysis command
    Module required: module add python/cpu/3.7.2
    :param input_fq_path: only single file, pair-end reads should be concatenated
    :param humann_output_dir:
    :param output_basename: sample name prefix
    :param metaphlan_dir: user_path/software/metaphlan4_database_vJan21
    :param metaphlan_index: mpa_vJan21_CHOCOPhlAnSGB_202103
    :param chocophlan_dir: metaphlan database, nucleotide-database
    user_path/software/HUMAnN3_database/chocophlan
    :param uniref_dir: protein-database, user_path/software/HUMAnN3_database/uniref
    :param thread: Ze used 16
    :return:
    """
    command_humann3 = f"humann --input {input_fq_path} --output {humann_output_dir} -" \
                      f"-output-basename {output_basename} --metaphlan-options '--bowtie2db {metaphlan_dir} --index {metaphlan_index} --read_min_len 10 -t marker_ab_table --add_viruses' " \
                      f"--nucleotide-database {chocophlan_dir} --protein-database {uniref_dir} " \
                      f"--prescreen-threshold 0.001 --threads {thread} --verbose"
    print(f"The command of HUMAnN3 is {command_humann3}")
    os.system(command_humann3)


def humann_pathway_report_combine(humann_pathway_report_path_list, combined_report_path):
    df0 = pd.read_table(humann_pathway_report_path_list[0])
    for humann_pathway_report_path in humann_pathway_report_path_list[1:]:
        df_add = pd.read_table(humann_pathway_report_path)
        df0 = pd.merge(df0, df_add, on='# Pathway', how="outer")
    df0.fillna(0)
    df0.to_csv(combined_report_path, sep="\t", index=False)
    # return df0


if __name__ == "__main__":
    """
    logging.basicConfig(filename="longstrain_test.log", format='%(asctime)s\t%(levelname)s\t%(name)s: %(message)s',
                        level=logging.DEBUG)
    my_logger = logging.getLogger('mylogger')
    my_logger.setLevel(logging.DEBUG)
    # rf_handler = logging.handlers.TimedRotatingFileHandler(sys.argv[-1], when='midnight', interval=1, backupCount=7, atTime=datetime.time(0, 0, 0, 0))
    rf_handler = logging.FileHandler(sys.argv[-1])
    rf_handler.setFormatter(logging.Formatter("%(asctime)s - %(levelname)s - %(message)s"))
    my_logger.addHandler(rf_handler)
    my_logger.info("Processing command line is: " + " ".join(sys.argv))
    """
    pass
