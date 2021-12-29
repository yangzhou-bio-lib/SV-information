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
Based on the result file obtained from the `smooth_break_point.sh` script, further determine the location of the SV breakpoint.
