#setwd('./RobustClone/example')
rm(list=ls())
library(R.matlab)
library(reticulate)
library(rlist)
#library(blockcluster)
library(igraph)
library(leiden)
path_to_python <- "/home/akanksha/anaconda3/bin/python"
#path_to_python <- "/usr/bin/python"
use_python(path_to_python, required = T)
py_config()
source('/mnt/storage/Akanksha/SingleCell/June_2021_analysis/RobustClone-master/matlab_and_R_scripts/Clustering_EvolutionaryTree_function.R')
files <- list.files("/mnt/storage/Akanksha/SingleCell/June_2021_analysis/sim_data_June2021/vary_sites_n/GT_500x2000_10/DT_NOISY_RGL/")
## example for SNV data
for (i in 44:length(files)){  
  setwd("/mnt/storage/Akanksha/SingleCell/June_2021_analysis/sim_data_June2021/vary_sites_n/GT_500x2000_10/DT_NOISY_RGL/")
  temp <- files[i]
  datamat <- as.matrix(readMat(temp)[[1]])
  str_folder <- strsplit(temp,".mat")[[1]][1]
  dir.create(str_folder)
  pathstr <- paste("/mnt/storage/Akanksha/SingleCell/June_2021_analysis/sim_data_June2021/vary_sites_n/GT_500x2000_10/DT_NOISY_RGL/",str_folder,"/", sep="")
  setwd(pathstr)
  AA <- datamat
  AA <- t(AA)
  S_str <- strsplit(temp,"_denoised")[[1]][1]
  S_str_combine <- paste("/mnt/storage/Akanksha/SingleCell/June_2021_analysis/sim_data_June2021/vary_sites_n/GT_500x2000_10/DT_NOISY_RGL_Smatrices/", S_str,"_S.mat",sep="")
  Smat <- as.matrix(readMat(S_str_combine)[[1]])
  start_time <- Sys.time()
  #out<-leiden(Smat)
  out<-leiden(Smat, resolution_parameter=0.02)
  print(length(unique(out)))
  uni_class_cell <-unique(out)
  robust_clone_cell <-{}
  for (un in 1:length(uni_class_cell)){
    idx_cell <- which(out==uni_class_cell[un])
    robust_clone_cell <- c(robust_clone_cell,list(idx_cell))
  }
  robust_clone <- robust_clone_cell
  clone_gety <- subclone_GTM(AA, robust_clone, 'SNV')
  # obtain subclonal GTM
  MST <- plot_MST(clone_gety, robust_clone, 'SNV', 'exampledata') # calculate and plot clonal MST
  clones_mt <- new_variant_site(clone_gety, MST, 'SNV', 'parent') # obtain the variant SNV loci each subclone compared with its parent subclone
  
  end_time <- Sys.time()
  diffTimecal <- end_time - start_time
  print("Time Elapsed")
  print(diffTimecal)
  #Saving outputfiles
  save(robust_clone, file=paste("LVcluster_",str_folder,".Rdata",sep=""));
  write.table(clone_gety, file="sub_clonal_GTM", row.names = F, col.names = F, sep=",");
  write.table(MST, file="MST", row.names = F, col.names = F, sep=",");
  save(clones_mt, file=paste("SNVloci_",str_folder,".Rdata",sep=""));
}




