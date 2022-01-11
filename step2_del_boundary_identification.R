setwd("$PATH")
path<-c("$PATH")

list1<-list.files("$PATH")
i=1
j=1
for(i in 1:length(list1)){
    all_result<-NULL
    list2<-list.files(paste("$PATH",list1[i],sep=""),pattern="*.reads.comm")
    for(j in 1:length(list2)){
        data<-read.table(paste(path,list1[i],list2[j],sep="/"),fill=TRUE)
        data.temp<-na.omit(data)
        if(dim(data)[1]!=dim(data.temp)[1]){
            result<-cbind(list2[j],"N","N","N","N")
            }
        else{
            #filter the reads in same location
            data<-data[which(data[,1]!=data[,3]),]
            #filter the reads located in differnet region
            if(dim(data)[1]==0){
                result<-cbind(list2[j],"N","N","N","N")
                }
            else{
                data<-data[which(nchar(data[,2])<=7&nchar(data[,4])<=7),]
                if(dim(data)[1]==0){
                    result<-cbind(list2[j],"N","N","N","N")
                    }
                else{
                    data1<-data.frame(gsub("S","H",data[,2]))
                    data2<-data.frame(gsub("S","H",data[,4]))
                    data1<-do.call(rbind,strsplit(as.character(data1[,1]),"H"))
                    data1<-do.call(rbind,strsplit(as.character(data1[,1]),"M"))
                    data2<-do.call(rbind,strsplit(as.character(data2[,1]),"M"))
                    data2<-do.call(rbind,strsplit(as.character(data2[,1]),"H"))
                    data3<-cbind(data[,c(1,3)],data1,data2)
                    data3<-data3[which(data3[,3]>=30&data3[,4]>=30),]
                    if(dim(data3)[2]!=6|dim(data3)[1]==0){
                        result<-cbind(list2[j],"N","N","N","N")
                        }

                    else{
                        start<-unique(as.numeric(data3[,1])+as.numeric(data3[,3]))
                        end<-unique(as.numeric(data3[,2]))-1
                        repeats<-as.numeric(data3[,3])+as.numeric(data3[,6])
                        reads<-as.numeric(data3[,3])+as.numeric(data3[,4])

                        if(length(start)!=1|length(end)!=1){
                            start<-data.frame(table(as.numeric(data3[,1])+as.numeric(data3[,3])))
                            end<-data.frame(table(as.numeric(data3[,2])-1))
                            start<-as.character(start[which(start[,2]==max(start[,2])),][1,1])
                            end<-as.character(end[which(end[,2]==max(end[,2])),][1,1])
                            repeats<-repeats-reads
                            repeats<-data.frame(table(repeats))
                            repeats_real<-as.character(repeats[which(repeats[,2]==max(repeats[,2])),][1,1])
                            result<-cbind(list2[j],start,end,repeats_real,max(repeats[,2]))
                            }
                        else{
                            start<-unique(as.numeric(data3[,1])+as.numeric(data3[,3]))
                            end<-unique(as.numeric(data3[,2]))-1
                            repeats<-unique(repeats-reads)
                            result<-cbind(list2[j],start,end,repeats,dim(data3)[1])
                            }
                        }
                    }
                }
            }
            all_result<-rbind(all_result,result)
        }
    all_result<-cbind(do.call(rbind,strsplit(all_result[,1],"-"))[,2],all_result[,c(2:5)])
    write.table(all_result,paste("$PATH",list1[i],"_breakpoint",sep=""),sep="\t",col.names=FALSE,row.names=FALSE,quote=FALSE)
}
