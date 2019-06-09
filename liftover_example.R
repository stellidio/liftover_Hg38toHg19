#liftover example

#installation of 'rtracklayer' package and attachment
# source("http://www.bioconductor.org/biocLite.R")
# biocLite("rtracklayer")
library("rtracklayer")


#download and unzip the file "hg38ToHg19.over.chain" from http://hgdownload.cse.ucsc.edu/gbdb/hg38/liftOver/
#import chain file
chain_file=import.chain("hg38ToHg19.over.chain")

liftover_function<-function(csv_position_GRCh38) {
  
  #read csv file with position and create a dataframe
  readcsvfun<-function(csvfile) {
    readcsv_df<-as.data.frame(read.csv(csvfile, header=T,fill=T, row.names=NULL, 
                                       check.names=FALSE, colClasses = c("factor")))
    return(readcsv_df)
  }
  position_file<-readcsvfun(csv_position_GRCh38) 
  
  #create an empty dataframe to be filled with the positions
  #first column will be the clear hg38 position
  #second column will be the clear hg38 position -1 
  #the third column will be the clear hg19 position
  new_hg38_pos<-data.frame(ClearPosition_hg38=character(0), 
                             ClearPosition_hg38_Minus1=character(0),
                             ClearPosition_hg19=character(0),
                             stringsAsFactors=FALSE)
  
  #function that substract 1 from hg38 position
  number<-function(charact) {
    a<-as.numeric(charact)-1
  }
  

  for(i in 1:nrow(position_file)) {
    
    #remove the letter to keep a clear position number
    #first column of new_hg38_pos is clear hg38 position
    new_hg38_pos[i,1]<-lapply(position_file$Position_GRCh38[[i]], 
                               function(x) as.numeric(gsub("[^[:digit:]]", "", x)))
    
    #slit the clear positions and apply number function to create thwo columns in the df
    #second column of new_hg38_pos is hg38 position -1
    new_hg38_pos[i,2]<-lapply(new_hg38_pos[i,1], number)
    
    #2nd and 1st columns combined make a genomic range
    #make a df with ranges 
    rangesdf<-data.frame(chr=as.character(position_file$Chromosome[i]),
                         start=as.numeric(new_hg38_pos$ClearPosition_hg38_Minus1[i]),
                         end=as.numeric(new_hg38_pos$ClearPosition_hg38[i]))
    
    #change df to GRanges and apply liftover
    make_ranges<-makeGRangesFromDataFrame(rangesdf, TRUE)
    
    change_assembly<-liftOver(make_ranges, chain_file)
    hg19_df<-as.data.frame(change_assembly)
    
    #3rd column of new_hg38_pos is clear hg19 position
    new_hg38_pos[i,3]<-paste(hg19_df$end)

    
  }
  #combine clear hg38 and hg19 positions with the initial initial(input) csv file 
  new_position_file<-cbind(position_file,
                           Position_hg38=new_hg38_pos$ClearPosition_hg38, 
                           Position_hg19=new_hg38_pos$ClearPosition_hg19)
  
  #write a new csv with all the positions
  #initial csv including 2 more columns(clear hg38 position & hg19 position)
  write.csv(new_position_file, file = paste("lifted_position_GRCh38",".csv", sep = ""), row.names=FALSE)
  
}

input_file<-liftover_function("position_GRCh38.txt")




