function my_cmd()
{
    run_sample=$1
    sample_name=`echo $run_sample|cut -d "." -f1`
	cat "/Pindel_vcf/$run_sample" | \
        awk '$1!~/^#/{
            if($7=="PASS"){
                if($10!~/^0\/0/ && $10!~/\.\/\./){
                    split($8,x,";");
                    for(i in x){
                            n=split(x[i],y,"=");
                            if(n==2)INFO[y[1]]= y[2];
                    }
                    if(INFO["SVTYPE"]=="INS") printf("%s\t%d\t%d\t%s;%s\n", $1,$2,$2+1,INFO["SVTYPE"],INFO["SVLEN"]+1);
                    else if(INFO["SVTYPE"]=="DEL") printf("%s\t%d\t%d\t%s;%s\n", $1,$2,INFO["END"],INFO["SVTYPE"],-(INFO["SVLEN"])+1);
                    else if(INFO["SVTYPE"]=="RPL") printf("%s\t%d\t%d\t%s;%s:%s\n", $1,$2,INFO["END"],INFO["SVTYPE"],-(INFO["SVLEN"])+1,INFO["NTLEN"]);
                    else  printf("%s\t%d\t%d\t%s;%s\n", $1,$2,INFO["END"],INFO["SVTYPE"],INFO["SVLEN"]+1);
                    delete INFO;
                }
            }
        }' - > /Pindel_vcf/${sample_name}.txt
}
export -f my_cmd
parallel -j 10 my_cmd :::: /Pindel_vcf/pindel_path.list
