# script-SV-genotyping
These are the scripts used for the project of cattle SV catalogues.

#  The description of workflow
The Deletion list needs to be constructed. We used the results of four software, retaining 50bp to 50 Mbp, at least two software-supported SV (overlapping 80%).  
First,  The result files generated by the four software programs was normalized to extract the breakpoint information of SV.  
Then, integrate the results of multiple software, keeping SVs that are supported by at least two software. Extract all deletion breakpoints to get the deletion list.  
Finally, using the bam files and the deletion list file as input files, the obtained deletion joint genotyping with the GGDTRS.py script.  

# Integrated results from multiple software to obtain SV maps.
## 1. smooth_break_point.sh
### Input files
Input files need to be prepared as shown below:
`Common SV regions, Upstream breakpoint regions detected by different software, Downstream breakpoint regions detected by different software`  
Like this:
```
chr1_383181_391717,CM008168.2:383081-383294,CM008168.2:391608-391818
chr1_483202_491758,CM008168.2:483102-483337,CM008168.2:491647-491859
```
### Output
We look for evidence of slipt reads in the original bam file and output the eligible reads to the `${sample}.${a}.splitters.reads.comm` file

## 2. Identify_breakpoint.R
Based on the result file obtained from the `smooth_break_point.sh` script, further determine the location of the SV breakpoint and get the SV map.

# GGDTRS.py
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
