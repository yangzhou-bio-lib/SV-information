#!/usr/bin/env python3
# -*-coding:utf-8 -*-

"""
# File       : target_vertion1.py
# Time       : 2021/9/6 16:11
# Author     : "yanglv_yoyo"
# version    : python 3.7
# Description: Based on the deletion list, multiple bam files are provided for joint genotyping.
Provide a list of bam file paths and a bed file containing deletion breakpoints. Each line in the bam file list is a bam file path, separated by the newline character '/n'; 
the bed file containing deletion breakpoints only needs to provide the coloring number, start position, and end position.
This method detects the typing category of each deletion locus for each bam file, and finally merges the detection results of multiple bam files to generate VCF files.
"""

import pysam
import pandas as pd
import numpy as np
import time
from multiprocessing import Pool
import sys
import getopt
import os
import re   


def get_cnv_region(file_path=str):
    cnv_candidate_regions = pd.read_table(file_path, sep='\t', header=None, names=['chr', 'start', 'end'])
    return cnv_candidate_regions


def get_pysam_AlignmentFile(bam_file_path=str):
    bam_file = pysam.AlignmentFile(bam_file_path, 'rb')
    return bam_file


def main(bam_file, cnv_region_file, outputfile):
    # bam_file_path = '/sdb1/usrdata/yanglv/cattledata/African_111cattleAnd_2buffalo/African_111cattleAnd_2buffalo_bam/SRR12452190/SRR12452190_mkdup.bam'
    bam_file_path = bam_file
    bam_file = get_pysam_AlignmentFile(bam_file_path)
    # cnv_candidate_regions = get_cnv_region(
    #     '/home/huyan/data/yanglv_project/genotyping/SRR12452190_chr1_CNV_candidate.txt')
    cnv_candidate_regions = get_cnv_region(cnv_region_file)
    SM_ID = bam_file.header.to_dict()['RG'][0]['SM']
    outputfile_path = outputfile + SM_ID + '.bed'

    with open(outputfile_path, 'w') as result:
        # vcf ：#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SRR8588260
        # result.write('#CHROM\tSTART\tEND\tINFO_L\tINFO_R\n')
        result.write('#CHROM\tSTART\tEND\tL_READS_COUNT\tL_READS_INFO\tR_READS_COUNT\tR_READS_INFO\tCOMMON_READS\n')
        result.close()

    for i in range(0, len(cnv_candidate_regions)):
        contig = cnv_candidate_regions.iloc[i, 0]
        # contig = 'CM008168.2'
        region_start = cnv_candidate_regions.iloc[i, 1] - 1
        region_end = cnv_candidate_regions.iloc[i, 2] - 1

        left_start = int(region_start - 100)
        left_end = int(region_start + 100)

        right_start = int(region_end - 100)
        right_end = int(region_end + 100)

        region_left_reads = bam_file.fetch(contig=contig, start=left_start, end=left_end)
        region_right_reads = bam_file.fetch(contig=contig, start=right_start, stop=right_end)

        cnv_segment_dir = {'chr': 0, 'start': 0, 'end': 0}
        left_reads_dir = {'count': 0, 'support_count': 0, 'over_count': 0, 'clip_reads_count': 0,
                          'absolute_read_count': 0, 'read_info': [], 'clip_reads_id': []}
        # right_reads_dir = {'count': 0, 'read_breakpoint': dict(), 'read_id': [], 'MQ': []}
        right_reads_dir = {'count': 0, 'support_count': 0, 'over_count': 0, 'clip_reads_count': 0,
                           'absolute_read_count': 0, 'read_info': [], 'clip_reads_id': []}
        cnv_segment_dir['chr'] = contig
        cnv_segment_dir['start'] = region_start
        cnv_segment_dir['end'] = region_end

        read_info_left = []
        for read in region_left_reads:
            left_read_cigar_stat = {'M': 0, 'S': 0, 'H': 0}
            if read.mapq >= 30 and read.get_tag('NM') <= 1:
                left_reads_dir['count'] = left_reads_dir.get('count') + 1
                read_start = read.reference_start  
                read_end = read.reference_end  

                cigartuples = read.cigartuples  
                for iterm in cigartuples:
                    if iterm[0] == 0:
                        left_read_cigar_stat['M'] = iterm[1]
                    if iterm[0] == 4:
                        left_read_cigar_stat['S'] = iterm[1]
                    if iterm[0] == 5:
                        left_read_cigar_stat['H'] = iterm[1]

                read_breakpoint = int(read_start) + left_read_cigar_stat.get('M')
                if read_breakpoint == region_start:
                    # left_reads_dir['support_count'] = left_reads_dir.get('support_count') + 1
                    if (left_read_cigar_stat.get('S') + left_read_cigar_stat.get('H')) >= 20:
                        left_reads_dir['clip_reads_count'] = left_reads_dir.get('clip_reads_count') + 1
                        left_reads_dir['clip_reads_id'].append(read.query_name)
                        left_reads_dir['support_count'] = left_reads_dir.get('support_count') + 1

                    # elif (left_read_cigar_stat.get('S') + left_read_cigar_stat.get('H')) == 0:
                    #     left_reads_dir['absolute_read_count'] = left_reads_dir.get('absolute_read_count') + 1
                    #     left_reads_dir['support_count'] = left_reads_dir.get('support_count') + 1

                elif read_breakpoint > region_start:
                    if read_start <= region_start <= read_end:
                        left_reads_dir['over_count'] = left_reads_dir.get('over_count') + 1
                        read_info_left.append(min(abs(region_start - read_end), abs(region_start - read_start)))

                left_reads_dir['read_info'].append((read.query_name, read.cigartuples, read_breakpoint, read.mapq))

        read_info_right = []
        for read_r in region_right_reads:
            right_read_cigar_stat = {'M': 0, 'S': 0, 'H': 0}
            if read_r.mapq >= 30 and read_r.get_tag('NM') <= 1:
                # right_reads_dir['count'] = right_reads_dir.get('count') + 1
                read_start = read_r.reference_start 
                read_end = read_r.reference_end  
                cigartuples = read_r.cigartuples  

                for iterm in cigartuples:
                    if iterm[0] == 0:
                        right_read_cigar_stat['M'] = iterm[1]
                        # if len(cigartuples) == 1:
                        #     right_read_cigar_stat['ONLYM'] = 1
                    elif iterm[0] == 4:
                        right_read_cigar_stat['S'] = iterm[1]
                    elif iterm[0] == 5:
                        right_read_cigar_stat['H'] = iterm[1]

                read_breakpoint = int(read_end) - right_read_cigar_stat.get('M') - 1
                if read_breakpoint == region_end:
                    # right_reads_dir['support_count'] = right_reads_dir.get('support_count') + 1
                    if (right_read_cigar_stat.get('S') + right_read_cigar_stat.get('H')) >= 20:
                        # right_reads_dir['clip_reads_count'] = right_reads_dir.get('clip_reads_count') + 1
                        right_reads_dir['clip_reads_id'].append(read_r.query_name)
                        right_reads_dir['support_count'] = right_reads_dir.get('support_count') + 1

                elif read_breakpoint < region_end:
                    if read_start <= region_end <= read_end:
                        right_reads_dir['over_count'] = right_reads_dir.get('over_count') + 1
                        read_info_right.append(min(abs(region_end - read_end), abs(region_end - read_start)))

                right_reads_dir['read_info'].append(
                    (read_r.query_name, read_r.cigartuples, read_breakpoint, read_r.mapq))

        # read_info_left = ','.join(str(left_reads_dir['read_info'][i][2] - cnv_segment_dir['start']) for i in
        #                           range(len(left_reads_dir['read_info'])))

        # read_info_right = ','.join(str(right_reads_dir['read_info'][i][2] - cnv_segment_dir['end']) for i in
        #                            range(len(right_reads_dir['read_info'])))

        # read_info_right = ','.join(str(right_reads_dir['read_info'][i][2] - cnv_segment_dir['end']) for i in
        #                            range(len(right_reads_dir['read_info'])))
        #
        # read_info_right = right_reads_dir['read_info']

        set_left = set(left_reads_dir['clip_reads_id'][i] for i in range(len(left_reads_dir['clip_reads_id'])))
        # print(set_left)
        set_right = set(right_reads_dir['clip_reads_id'][i] for i in range(len(right_reads_dir['clip_reads_id'])))
        inset = set_left.intersection(set_right)

        with open(outputfile_path, 'a+') as result:
            result.write("%s\t%s\t%s\t%s\t" % (cnv_segment_dir['chr'], cnv_segment_dir['start'], cnv_segment_dir['end'],
                                               "L_READS_COUNT=%s,%s,%s\tRead_info=%s" % (
                                                   left_reads_dir['support_count'] + left_reads_dir['over_count'],
                                                   left_reads_dir['support_count'],
                                                   left_reads_dir['over_count'],
                                                   ",".join(str(i) for i in read_info_left))))
            result.write("R_READS_COUNT=%s,%s,%s\tRead_info=%s\tCOMMON_READS:%s\n" % (
                right_reads_dir['support_count'] + right_reads_dir['over_count'], right_reads_dir['support_count'],
                right_reads_dir['over_count'], ",".join(str(i) for i in read_info_right),
                str(len(inset))))
            result.close()
    bam_file.close()

def get_bam_list(bam_list_path=None):
    with open(bam_list_path, 'r') as bam_file:
        bam_list = bam_file.readlines()

    return bam_list


class Vcf:
    """
    Merge bed files, output vcf file.
    """
    def __init__(self, bam_file_list=None, cnv_region_file=None, output_path=None):
        self.file_format = 'VCFv4.2'
        self.file_data = time.strftime("%Y%m%d")
        self.bam_file_list = bam_file_list
        contig_name_list, contig_len_list = self.get_contig_INFO(bam_file_list)
        self.header = self.get_header(contig_name_list, contig_len_list)
        self.frist_nine_col = self.get_frist_nine_col(cnv_region_file)
        self.genotype_result = self.get_genotype()
        self.vcf_content = pd.merge(self.frist_nine_col, self.genotype_result, left_index=True, right_index=True, how="left")
        self.vcf_output(output_path)


    def get_header(self, contig_name_list = list, contig_len_list = list):
        """
        header information
        """
        # ##contig=<ID=1,length=157905821>
        contig_ID_LEN = zip(contig_name_list, contig_len_list)
        header_txt = "##fileformat=" + self.file_format + "\n"
        header_txt = header_txt + "##fileDate=" + self.file_data + "\n"
        for contig_ID, contig_LEN in contig_ID_LEN:
            header_txt = header_txt + "##contig=<ID=" + str(contig_ID) + ",length=" + str(contig_LEN) + ">" + "\n"

        header_txt = header_txt + \
        """##INFO=<ID=PR,Number=0,Type=Flag,Description="Provisional reference allele, may not be based on real reference genome">\n""" + \
        """##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">\n""" + \
        """##INFO=<ID=SVTYPE,Number=2,Type=String,Description="Type of structural variant">\n""" + \
        """##ALT=<ID=2,Description="Deletion">\n""" + \
        """##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n"""
        return header_txt

    def get_contig_INFO(self, bam_file_list):
        """
        get contig infomation from bam.
        """
        bam_file = pysam.AlignmentFile(bam_file_list[0].strip(), 'rb')
        contig_name_list = bam_file.header.references
        contig_len_list = bam_file.header.lengths
        del bam_file
        return contig_name_list, contig_len_list

    def get_frist_nine_col(self, cnv_region_file):
        """
        #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT
        """
        deletion_data = pd.read_table(cnv_region_file,names=["#CHROM", "POS", "END"])
        first_nine_col = deletion_data[["#CHROM", "POS"]]
        first_nine_col["ID"] = deletion_data["#CHROM"].astype(str) + "_" + deletion_data["POS"].astype(str) + "_" + deletion_data["END"].astype(str)
        n_raw = len(deletion_data)
        first_nine_col["REF"] = [1] * n_raw
        first_nine_col["ALT"] = [2] * n_raw
        first_nine_col["QUAL"] = ["."] * n_raw
        first_nine_col["FILTER"] = ["."] * n_raw
        INFO_pd = pd.DataFrame(index=first_nine_col.index)
        INFO_pd["SVTYPE"] = ["SVTYPE=DEL;"] * n_raw
        INFO_pd["END"] = ["END="] * n_raw
        INFO_pd["END_num"] = deletion_data["END"].values
        INFO_pd["PR"] = [";PR"] * n_raw
        first_nine_col["INFO"] = INFO_pd["SVTYPE"].astype(str) + INFO_pd["END"].astype(str) + INFO_pd["END_num"].astype(str) + INFO_pd["PR"].astype(str)
        first_nine_col["FORMAT"] = ["GT"] * n_raw
        first_nine_col.index = first_nine_col["ID"].values
        return first_nine_col

    def get_genotype(self):
        """
        genotyping
        """
        file_list = os.listdir(temp_dir)
        result = pd.DataFrame(index=self.frist_nine_col["ID"].values)
        for sample in file_list:
            # file = 'SRR934437.bed'
            self.get_array(sample,result)
        # print(result)
        genotype_result = result.replace(0,'0/0').replace(1,'0/1').replace(2,'1/1').replace(np.nan, './.')
        return genotype_result

    def get_array(self,sample, result):
        file = temp_dir + r'/' + sample
        file_df = pd.read_table(file)
        new_df = pd.DataFrame(index=file_df['#CHROM'])
        L_READS_COUNT = file_df['L_READS_COUNT'].str.split('=|,', expand=True)
        new_df['L_all'] = L_READS_COUNT[1].values
        new_df['L_clip'] = L_READS_COUNT[2].values
        new_df['L_absolute'] = L_READS_COUNT[3].values
        L_READS_INFO = file_df['L_READS_INFO'].str.split('=', expand=True)[1]
        new_df['L_info'] = L_READS_INFO.values
        R_READS_COUNT = file_df['R_READS_COUNT'].str.split('=|,', expand=True)
        new_df['R_all'] = R_READS_COUNT[1].values
        new_df['R_clip'] = R_READS_COUNT[2].values
        new_df['R_absolute'] = R_READS_COUNT[3].values
        R_READS_INFO = file_df['R_READS_INFO'].str.split('=', expand=True)[1]
        new_df['R_info'] = R_READS_INFO.values
        new_df['COMMON_READS'] = file_df['COMMON_READS'].str.split(':', expand=True)[1].values
        my_result = new_df.apply(
            lambda x: get_genotype(x.L_all, x.L_clip, x.L_absolute, x.L_info, x.R_all, x.R_clip, x.R_absolute, x.R_info,
                                   x.COMMON_READS), axis=1)
        print("**************")
        print(my_result)
        print("***************")
        sample_ID = re.sub('\.bed', '', sample)
        result[sample_ID] = my_result.values

    def vcf_output(self,output_path='_deletions.vcf'):
        """
        输出vcf文件
        """
        output_path = output_path + '_deletions.vcf'
        with open(output_path,'w+') as f:
            f.write(self.header)
        self.vcf_content.to_csv(output_path,mode='a', index=None, sep='\t')


def get_genotype(L_all, L_clip, L_absolute, L_info, R_all, R_clip, R_absolute, R_info, COMMON_READS):
    L_info_list = re.findall('\d\d|\d', L_info)
    L_info_list1 = []
    for iterm in L_info_list:
        if int(iterm) < 10:
            L_info_list1.append(iterm)
    n1 = len(L_info_list1)  #  

    R_info_list = re.findall('\d\d|\d', R_info)
    R_info_list1 = []
    for iterm in R_info_list:
        if int(iterm) < 10:
            R_info_list1.append(iterm)

    n2 = len(R_info_list1)
    RD = 5
    if ((int(L_all) - n1) >= RD or (int(R_all) - n2) >= RD) and int(L_clip) == 0 and int(R_clip) == 0:
        genotype = int(0)
    elif ((int(L_all) - n1) >= RD or (int(R_all) - n2) >= RD) and (int(L_absolute) - n1) == 0 and (int(R_absolute) - n2) == 0:
        genotype = int(2)
    elif ((int(L_all) - n1) >= RD or (int(R_all) - n2) >= RD) and (int(L_absolute) > 0 and int(L_clip) > 0) or (int(R_absolute) > 0 and int(R_clip) > 0):
        genotype = int(1)
    else:
        genotype = np.NaN
    return genotype

if __name__ == '__main__':
    start_time = time.time()
    t = 1
    pwd_path = os.getcwd()
    global temp_dir
    temp_dir = pwd_path + r'/temp_target_bed_files/'
    try:
        os.mkdir(temp_dir)
    except:
        pass
    
    try:
        opts, args = getopt.getopt(sys.argv[1:], "hb:l:o:t:v", ["help", "bamfile_list=", "deletion_list=", "output=", "thread=", "version="])
        for opt, arg in opts:
            if opt in ('-h', "--help"):
                print(r"""
usage: GGDTRS.py [-h] -b input_bam_list -o output_vcf_file
                    -l deletion_list [-t num of thread]

The script software joint genotype for the provided multiple BAM files according to the deletion list.
The BAM file path list of each line is a BAM file path and a bed file containing deletion breakpoint was provided, the bed file that contains deletion breakpoints only needs to provide the chromosome numbers, START positions, and END positions.  
In this script, the classification type of each DELETION site was detected for each BAM file. Finally, the detection results of multiple BAM files were merged to generate a VCF file.  


optional arguments:
-h, --help            show this help message and exit
-b, --bamfile_list
                    List file of input BAM files. Must be indexed.
-l, --deletion_list
                    Bed file of DELETION SV. 
-o, --outfile
                    Prefix for output filenames (same as the input BAM
                    filename without the extension by default)
-t, --thread
                    The number of thread(default=1).
-v, --version         show program's version number and exit    
            """)
                sys.exit(2)
            elif opt in ("-v", "--version"):
                print("GDTRS.py 1.0")
                sys.exit(2)

            if not opts:
                print('ERROR: -b -l -o are necessary!\n'
                    'usage: /[PATH]/python3 target.py -b <input_bam_list> -l <deletion_list> -o <output_vcf_file> -t <num of thread>')
                sys.exit(2)
            

    except getopt.GetoptError:
        print('usage: /[PATH]/python3 target.py -b <input_bam_list> -l <deletion_list> -o <output_vcf_file> -t <num of thread>')
        sys.exit(2)

    for opt, arg in opts:
        if opt in ("-b", "--bamfile_list"):
            bam_file_list = arg
        elif opt in ("-o", "--outfile"):
            outputfile = arg
        elif opt in ("-l", "--deletion_list"):
            cnv_region_file = arg
        elif opt in ("-t", "--thread"):
            t = arg

    if 'bam_file_list' in dir() and 'outputfile' in dir() and 'cnv_region_file' in dir():
        pass
    else:
        print('ERROR: -b -l -o are necessary!\n'
        'usage: /[PATH]/python3 target.py -b <input_bam_list> -l <deletion_list> -o <output_vcf_file> -t <num of thread>')
        sys.exit(2)

    print("Deletion list is " + str(cnv_region_file))
    bam_lists = get_bam_list(bam_file_list)
    p = Pool(int(t))
    for sample in bam_lists:
        print("Genotyping for " + sample)
        p.apply_async(main, args=(sample.strip(),cnv_region_file, temp_dir,))
        
    p.close()
    p.join()
    #main(sample.strip(), cnv_region_file, temp_dir)
    vcf = Vcf(bam_file_list = bam_lists, cnv_region_file = cnv_region_file, output_path=outputfile)
    end_time = time.time()
    print('run time %s s' % (end_time - start_time))
