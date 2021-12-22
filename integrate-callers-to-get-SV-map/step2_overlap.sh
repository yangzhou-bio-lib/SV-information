#!/bin/bash

inputpath=""
samplelist=""
outpath=""
svtype=""
helpmode=0


while getopts i:s:o:t:h option
do
  case "${option}"
  in
    i) inputpath=${OPTARG};;
    s) samplelist=${OPTARG};;
    o) outpath="${OPTARG}";;
    t) svtype=${OPTARG};;
    h) helpmode=1;;
  esac
done

for samname in $(less $samplelist);
    do
        file1=${samname}.pindel_${svtype}
        file2=${samname}.delly_${svtype}
        file3=${samname}.breakdancer_${svtype}
        file4=${samname}.lumpy_${svtype}
        
bedtools intersect -f 0.8 -r  -a ${inputpath}/${file1} -b ${inputpath}/${file2} > ${outpath}/${samname}.${svtype}.temp0
bedtools intersect -f 0.8 -r  -a ${inputpath}/${file1} -b ${inputpath}/${file3} > ${outpath}/${samname}.${svtype}.temp1
bedtools intersect -f 0.8 -r  -a ${inputpath}/${file1} -b ${inputpath}/${file4} > ${outpath}/${samname}.${svtype}.temp2
bedtools intersect -f 0.8 -r  -a ${inputpath}/${file2} -b ${inputpath}/${file3} > ${outpath}/${samname}.${svtype}.temp3
bedtools intersect -f 0.8 -r  -a ${inputpath}/${file2} -b ${inputpath}/${file4} > ${outpath}/${samname}.${svtype}.temp4
bedtools intersect -f 0.8 -r  -a ${inputpath}/${file3} -b ${inputpath}/${file4} > ${outpath}/${samname}.${svtype}.temp5
cat ${outpath}/${samname}.${svtype}.temp0 ${outpath}/${samname}.${svtype}.temp1 ${outpath}/${samname}.${svtype}.temp2 ${outpath}/${samname}.${svtype}.temp3 ${outpath}/${samname}.${svtype}.temp4 ${outpath}/${samname}.${svtype}.temp5 > ${outpath}/${samname}.${svtype}.temp
sort -k1,1 -k2,2n ${outpath}/${samname}.${svtype}.temp|uniq > ${outpath}/${samname}.${svtype}_sort
bedtools merge -i ${outpath}/${samname}.${svtype}_sort -c 1 -o count > ${outpath}/${samname}.${svtype}

rm ${outpath}/${samname}.${svtype}.temp*
rm ${outpath}/${samname}.${svtype}_sort
done
