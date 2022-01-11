# script-SV-genotyping
This is a pipeline to identify the exact deletion boundaries based on split reads within the target region to overcome differences among the deletion boundaries identified by different software.
##  Workflow description 
First, the target regions were defined as the largest intervals of the boundaries reported by different SV tools. The reads located in each target region were evaluated and re-filtered using the following thresholds: mapped reads quality ≥ 30, only allow ≤ 1 mutation by considering the diversity of different cattle individuals and breeds, split reads were separated into two parts that mapped to two candidate regions with at least 30 bp. Deletions were separated into three categories according to the mapped reads: breakpoints with microhomology repeats, breakpoints with short insertions, and perfectly support the deletion breakpoints by split reads. Thus, the exact boundaries of the deletion were defined. The number of reads adjacent or across the deletion boundaries was first evaluated to genotype the deletion for each animal. Only the deletion boundaries supported by at least five reads were considered for genotyping. The deletion with both normal reads and split reads were genotyped as heterozygous (0/1). Otherwise, the deletion with only normal or split reads was genotyped as normal (0/0) or total deletion (1/1), respectively. 
![image](https://github.com/yangzhou-bio-lib/SV-information/raw/main/SV%20flow%20chart.png)
  
# 1. Deletion breakpoints and its marker sequnce identification.
## 1.1 step1.sh
### Input files
Input files need to be prepared as shown below:  
`Deletion tag, candidate regions for 5` of the deletion breakpoint, candidate regions for 3` of the deletion breakpoint`  
Like this:
```
deletion_1,chr1:383081-383294,chr1:391608-391818
deletion_2,chr1:483102-483337,chr1:491647-491859
```
### Output
We look for evidence of slipt reads in the original bam file and output the eligible reads to the `${sample}.${a}.splitters.reads.comm` file

## 1.2 step2.R
Based on the result file obtained from the `step1.sh` script, further determine the location of the SV breakpoint and get the SV map.

# 2. Genotyping
## GGDTRS.py
The script joint genotype for the provided multiple BAM files according to the deletion list.
The BAM file path list of each line is a BAM file path and a bed file containing deletion breakpoint was provided, the bed file that contains deletion breakpoints only needs to provide the chromosome numbers, START positions, and END positions.  
In this script, the classification type of each DELETION site was detected for each BAM file. Finally, the detection results of multiple BAM files were merged to generate a VCF file.
### Rely on third-party libraries
pysam (0.16.0.1); pandas (1.1.4); numpy (1.19.4)
### optional arguments:
-h, --help            show this help message and exit  
-b, --bamfile_list    List file of input BAM files. Must be indexed.  
-l, --deletion_list   Bed file of DELETION SV.   
-o, --outfile         Prefix for output filenames (same as the input BAM filename without the extension by default)  
-t, --thread          The number of thread(default=1).  
-v, --version         Show program's version number and exit  
