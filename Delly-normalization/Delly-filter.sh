function my_cmd()
{
    run_sample=$1
    sample_name=`echo $run_sample|cut -d "." -f1`
    cat "./Delly/$run_sample" | \
    awk '$1!~/^#/{
        if($7=="PASS" && $10!~/^0\/0/ && $10!~/^\.\/\./){
            split($8,x,";"); 
            for(i in x){ 
                n=split(x[i],y,"="); 
                if(n==2)INFO[y[1]]= y[2];
            }
        if(INFO["SVTYPE"]=="INS") printf("%s\t%d\t%d\t%s;%s\n", $1,$2,INFO["END"],INFO["SVTYPE"],INFO["INSLEN"]);
        else if(INFO["SVTYPE"]=="TRA") printf("%s\t%d\t%s\t%d\t%s\n", $1,$2,INFO["CHR2"],INFO["END"],INFO["SVTYPE"]);
        else  printf("%s\t%d\t%d\t%s;%d\n", $1,$2,INFO["END"],INFO["SVTYPE"],INFO["END"]-$2+1);
        delete INFO;
        }
    }' - > ./Delly/${sample_name}.txt
}
export -f my_cmd
parallel -j 20 -N 1 my_cmd :::: ./Delly/delly_path.list
