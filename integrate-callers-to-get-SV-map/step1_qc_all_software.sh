#!/bin/bash

dellypath=""
pindelpath=""
breakdancerpath=""
lumpypath=""
samplelist=""
outpath=""
svtype=""
helpmode=0


while getopts d:b:l:p:s:o:t:h option
do
  case "${option}"
  in
    d) dellypath=${OPTARG};;
    b) breakdancerpath=${OPTARG};;
    l) lumpypath=${OPTARG};;
    p) pindelpath=${OPTARG};;
    s) samplelist=${OPTARG};;
    o) outpath="${OPTARG}";;
    t) svtype=${OPTARG};;
    h) helpmode=1;;
  esac
done

for samname in $(less $samplelist);
    do
        file1=${samname}.pindel_${svtype}.txt
	      file2=${samname}.delly_${svtype}.txt
        file3=${samname}.breakdancer_${svtype}.txt
        file4=${samname}.lumpy_${svtype}.txt

awk -v OFS='\t' '$3-$2<=5000000&& $3-$2>50 {print $0}' ${pindelpath}/${file1} > ${outpath}/${samname}.pindel_${svtype}.temp
awk -v OFS='\t' '$3-$2<=5000000&& $3-$2>50 {print $0}' ${dellypath}/${file2} > ${outpath}/${samname}.delly_${svtype}.temp
awk -v OFS='\t' '$3-$2<=5000000&& $3-$2>50 {print $0}' ${breakdancerpath}/${file3} > ${outpath}/${samname}.breakdancer_${svtype}.temp
awk -v OFS='\t' '$3-$2<=5000000&& $3-$2>50 {print $0}' ${lumpypath}/${file4} > ${outpath}/${samname}.lumpy_${svtype}.temp


sort -k1,1 -k2,2n ${outpath}/${samname}.pindel_${svtype}.temp|uniq > ${outpath}/${samname}.pindel_${svtype}_sort
bedtools merge -i ${outpath}/${samname}.pindel_${svtype}_sort -c 1 -o count > ${outpath}/${samname}.pindel_${svtype}_merge
awk -v OFS='\t' '$4==1 {print $0}' ${outpath}/${samname}.pindel_${svtype}_merge > ${outpath}/${samname}.pindel_${svtype}

sort -k1,1 -k2,2n ${outpath}/${samname}.breakdancer_${svtype}.temp|uniq > ${outpath}/${samname}.breakdancer_${svtype}_sort
bedtools merge -i ${outpath}/${samname}.breakdancer_${svtype}_sort -c 1 -o count > ${outpath}/${samname}.breakdancer_${svtype}_merge
awk -v OFS='\t' '$4==1 {print $0}' ${outpath}/${samname}.breakdancer_${svtype}_merge > ${outpath}/${samname}.breakdancer_${svtype}

sort -k1,1 -k2,2n ${outpath}/${samname}.delly_${svtype}.temp|uniq > ${outpath}/${samname}.delly_${svtype}_sort
bedtools merge -i ${outpath}/${samname}.delly_${svtype}_sort -c 1 -o count > ${outpath}/${samname}.delly_${svtype}_merge
awk -v OFS='\t' '$4==1 {print $0}' ${outpath}/${samname}.delly_${svtype}_merge > ${outpath}/${samname}.delly_${svtype}

sort -k1,1 -k2,2n ${outpath}/${samname}.lumpy_${svtype}.temp|uniq > ${outpath}/${samname}.lumpy_${svtype}_sort
bedtools merge -i ${outpath}/${samname}.lumpy_${svtype}_sort -c 1 -o count > ${outpath}/${samname}.lumpy_${svtype}_merge
awk -v OFS='\t' '$4==1 {print $0}' ${outpath}/${samname}.lumpy_${svtype}_merge > ${outpath}/${samname}.lumpy_${svtype}

rm ${outpath}/*.temp
rm ${outpath}/*_sort
rm ${outpath}/*_merge
done
