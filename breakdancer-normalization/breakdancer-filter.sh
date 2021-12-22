function my_cmd()
{
    run_sample=$1
    sample_name=`echo $run_sample|cut -d "." -f1`
    cat "./Breakdancer/$run_sample" | \
	awk '$1!~/^#/{
        printf("%s\t%d\t%d\t%s;%d\n", $1,$2,$5,$7,$8);
   	}' - > ./Breakdancer/${sample_name}.txt
}

export -f my_cmd
parallel -j 10 -N 1 my_cmd ::::./Breakdancer/Breakdancer_path.list
