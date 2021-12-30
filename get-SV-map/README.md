# Integrated results from multiple software to obtain SV maps.
## 1. step1.sh
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

## 2. step2.R
Based on the result file obtained from the `step1.sh` script, further determine the location of the SV breakpoint and get the SV map.
