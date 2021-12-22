function my_cmd()
{
    run_sample=$1
    sample_name=`echo $run_sample|cut -d "." -f1`
    cat "/lumpy_path/$run_sample" | \
	awk '$1!~/^#/{
        split($8,x,";");
        for(i in x){
                n=split(x[i],y,"=");
                if(n==2)INFO[y[1]]= y[2];
        }
        if(INFO["SVTYPE"]=="DEL") printf("%s\t%d\t%d\t%s;%s\n", $1,$2,INFO["END"],INFO["SVTYPE"],-(INFO["SVLEN"])+1);
        else  printf("%s\t%d\t%d\t%s;%d\n", $1,$2,INFO["END"],INFO["SVTYPE"],INFO["SVLEN"]+1);
        delete INFO;
   	}' - > /lumpy_path/${sample_name}.txt
}
export -f my_cmd
parallel -j 80 -N 1 my_cmd :::: /lumpy_path/lumpy_path.list
