#!/bin/bash

samplelist=""
helpmode=0


while getopts s:h option
do
  case "${option}"
  in
    s) samplelist=${OPTARG};;
    h) helpmode=1;;
  esac
done

for sample_path in $(less $samplelist);
    do
      index=$(basename $sample_path)
      sample=$(echo $index | awk -F "_m|_d|_c|.mkd" '{print $1}')
      mkdir sample
      cd sample
      cat $PATH/${sample}.DEL.support_arrange_tag |while IFS="," read a b c
      do

        samtools view ${sample_path} $b -q 30 |awk '$12=="NM:i:0"|$12=="NM:i:1"{print $0}'| $PATH/extractSplitReads_BwaMem -i stdin > ${sample}.${a}.splitters.reads1.temp
        samtools view ${sample_path} $c -q 30 |awk '$12=="NM:i:0"|$12=="NM:i:1"{print $0}'| $PATH/extractSplitReads_BwaMem -i stdin > ${sample}.${a}.splitters.reads2.temp
        awk '{print $1}' ${sample}.${a}.splitters.reads1.temp|sort > ${sample}.${a}.splitters.tag1.temp
        awk '{print $1}' ${sample}.${a}.splitters.reads2.temp|sort > ${sample}.${a}.splitters.tag2.temp
        comm -12 ${sample}.${a}.splitters.tag1.temp ${sample}.${a}.splitters.tag2.temp|uniq > ${sample}.${a}.tag

        if [ -s ${sample}.${a}.tag ]
          then
          grep -f ${sample}.${a}.tag ${sample}.${a}.splitters.reads1.temp|sort > ${sample}.${a}.splitters.reads1.comm
          grep -f ${sample}.${a}.tag ${sample}.${a}.splitters.reads2.temp|sort > ${sample}.${a}.splitters.reads2.comm
          ln1=$(wc -l ${sample}.${a}.splitters.reads1.comm)
          ln2=$(wc -l ${sample}.${a}.splitters.reads2.comm)
          awk -v OFS='\t' '{print $4,$6}' ${sample}.${a}.splitters.reads1.comm > ${sample}.${a}.splitters.reads1.comm.temp
          awk -v OFS='\t' '{print $4,$6}' ${sample}.${a}.splitters.reads2.comm > ${sample}.${a}.splitters.reads2.comm.temp
          paste -d "\t" ${sample}.${a}.splitters.reads1.comm.temp ${sample}.${a}.splitters.reads2.comm.temp > ${sample}.${a}.splitters.reads.comm
          rm *.temp
          rm *.tag
          fi
       done
done






