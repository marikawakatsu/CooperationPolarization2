################################################################################
#
# Script for merging datasets
# Last updated: 26 Jul 2021
#
################################################################################

rm(list = ls())
source("rsrc/utils/functions.R")

# data parameters -- fixed for all plots
epsilon <- 1
N       <- 40
u       <- 0.001 # fixed, for now
gens    <- 20000000

# set directories
setwd("~/.julia/dev/CooperationPolarization2/") # to be updated later
file_dir  <- sprintf( "data/gens_%s/", format(gens, scientific = FALSE) )

###############################################
# Functions for merging and exporting data
###############################################
# function to load and merge data from multiple simulations
merge_data  <- function(file_dir, file_list){
  
  simdata   <- data.frame() # initialize data frame
  
  for (i in 1:length(file_list)){
    temp_data    <- read.csv( paste0(file_dir, file_list[i]), header = TRUE)
    temp_data$id <- paste0("run",i) # add id to identify data source
    
    if ( dim(temp_data[temp_data$N == 0,])[1] == 0 ){
      if (i == 1){
        simdata    <- rbind(simdata, temp_data) # bind new data to simdata
      }else if (i > 1){
        commoncols <- intersect(colnames(temp_data), colnames(simdata))   # select common columns
        simdata    <- rbind(simdata[,commoncols], temp_data[,commoncols]) # bind new data to simdata
      }
    }
  }
  print( sprintf("Done! Read %s files", format(length(file_list), scientific = FALSE) ) )
  return(simdata)
}
# function to check number of runs per case and export
export_data <- function(simdata, file_dir, file_name){

  casecount <- simdata %>% group_by(M, K, v, p1, β) %>% summarize(COUNT = n())
  
  if ( length(unique(casecount$COUNT)) == 1 ){
    write.csv( simdata, paste0(file_dir, file_name) )
    print("data successfully exported, no thresholding needed")
  }else if ( min(casecount$COUNT) == 150 ){
    simdata <- simdata %>% group_by(M, K, v, p1) %>% slice_head( n = min(casecount$COUNT) ) 
    write.csv( simdata, paste0(file_dir, file_name) )
    print( sprintf("data successfully exported, with thresholding to %d runs", min(casecount$COUNT)) )
  }else{
    warning("!!! check casecount !!!")
  }
  return(casecount)
}

##############################
# For Fig. 2 and Fig. S4
##############################
pattern   <- sprintf( "run_multi|vsweep|other" ) 
file_list <- list.files(path = file_dir, pattern = pattern)
simdata   <- merge_data(file_dir, file_list)
casecount <- export_data(simdata, file_dir, "della_all_merged.csv")

##############################
# For Fig. 3, Fig. 4, Fig. S2, Fig. S3
##############################
pattern   <- sprintf( "vsweep" ) 
file_list <- list.files(path = file_dir, pattern = pattern)
simdata   <- merge_data(file_dir, file_list)
casecount <- export_data(simdata, file_dir, "della_vsweep_merged.csv")

####################################
# For Fig. S5
####################################
pattern   <- sprintf( "run_multi" ) 
file_list <- list.files(path = file_dir, pattern = pattern)
simdata   <- merge_data(file_dir, file_list)
casecount <- export_data(simdata, file_dir, "della_run_multi_merged.csv")

####################################################
# For Fig. S6 (neutral simulations with M = 1) 
####################################################
# select file type
pattern   <- sprintf( "neutral_M1" )
file_list <- list.files(path = file_dir, pattern = pattern)
simdata   <- merge_data(file_dir, file_list)
casecount <- export_data(simdata, file_dir, "della_neutral_merged_M1.csv")

######################################################################
# For Fig. S7 (neutral simulations with M = 3, K = 2, plus more) 
######################################################################
# select file type
pattern   <- sprintf( "neutral_M2|neutral_M3" )
file_list <- list.files(path = file_dir, pattern = pattern)
simdata   <- merge_data(file_dir, file_list)
casecount <- export_data(simdata, file_dir, "della_neutral_merged_M2M3.csv")
