# -*- coding: utf-8 -*-
from multiprocessing.pool import Pool
from time import sleep, time
import pandas as pd
import os
import re
import io
import sys


sys.stdout = io.TextIOWrapper(sys.stdout.buffer,encoding='utf-8')
os.chdir("/lumpy_path/")


def my_cmd(i):
    print("%s子进程开始，进程ID：%d" % (i, os.getpid()))
    start = time()
    data=pd.read_csv(i,sep='\t',header=None)
    data.columns=['a','b','c','d']
    data=data.sort_values(by='a')
    lis=["%s%s"%("chr",i) for i in range(1,31)]
    times=pd.value_counts(data.iloc[:,0])[list(data.iloc[:,0].drop_duplicates().sort_values())]
    data.iloc[:,0]=sum(([x]*y for x,y in zip(lis, times)),[])
    data.to_csv("/lumpy_path/"+i.split("_L.")[0]+".lumpy.txt",sep='\t',header=0,index=None)
    end = time()
    print("%s子进程结束，进程ID：%d。耗时0.2%f" % (i, os.getpid(), end-start))


if __name__ == "__main__":
    print("父进程开始")
    p = Pool(32)
    for i in [i for i in os.listdir() if i[-3:].endswith("txt")]:
        p.apply_async(my_cmd, args=(i,))
    p.close()
    p.join()
    print("父进程结束。")
