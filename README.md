# script-SV-genotyping
These are the scripts used for the project of ‘Assembly of a pan-genome for global cattle reveals missing sequence and novel structural variation, providing new insights into their diversity and evolution history’

# GGDTRS.py
The script joint genotype for the provided multiple BAM files according to the deletion list.
The BAM file path list of each line is a BAM file path and a bed file containing deletion breakpoint was provided, the bed file that contains deletion breakpoints only needs to provide the chromosome numbers, START positions, and END positions.  
In this script, the classification type of each DELETION site was detected for each BAM file. Finally, the detection results of multiple BAM files were merged to generate a VCF file.
#### optional arguments:
-h, --help            show this help message and exit  
-b, --bamfile_list    List file of input BAM files. Must be indexed.  
-l, --deletion_list   Bed file of DELETION SV.   
-o, --outfile         Prefix for output filenames (same as the input BAM filename without the extension by default)  
-t, --thread          The number of thread(default=1).  
-v, --version         Show program's version number and exit  
