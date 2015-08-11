rm(list=ls())
setwd('/data/michaloa/S3')
#                 ~ remove above from the final script ~                

#====> Q: How to best code the paths and file names

# PATHS:

#. to input 
    
    dir_in <- 'data/raw/drug_response'
    dir_meta_in <- 'data/raw/metadata/'
    
    file_in_samples_ids <- file.path(dir_meta_in, 'samples_ids_drug_response.csv')

#. to output
    
    dir_out <- 'data/processed/drug_response'
        if (!dir.exists(dir_out))  dir.create(dir_out, recursive = T)
    
    dir_out_drc <- 'data/processed/drug_response/drc'
        if (!dir.exists(dir_out_drc))  dir.create(dir_out_drc, recursive = T)
    
    dir_out_response <- 'data/processed/drug_response/response'
        if (!dir.exists(dir_out_response))  dir.create(dir_out_response, recursive = T)
    
    dir_out_dose <- 'data/processed/drug_response/dose'
        if (!dir.exists(dir_out_dose))  dir.create(dir_out_dose, recursive = T)
    
    dir_meta_out <- 'data/processed/metadata'
        if (!dir.exists(dir_meta_out))  dir.create(dir_meta_out, recursive = T)
    
    file_out_compounds_ids <- file.path(dir_meta_out, 'compounds_ids_drug_response.csv')
    

# HARD CODED:

#. to get compound unique id
row_var <- 'SID'

#. to split into data types
compounds_var <- c('SID', 'readout', 'name', 'target', 'smi')
drc_var <- c('NPT', 'CCLASS', 'CCLASS2', 'HILL', 'INF', 'ZERO', 'MAXR', 'TAUC', 'FAUC', 'LAC50')
doses_var = paste('C', 0:10, sep='')

response_var = paste('DATA', 0:10, sep='')

#. vars to assign AC50
bad_curve_class2 = 4
curve_class2 = 'CCLASS2'
replace_ac50 = 'C10'
log_ac50 = 'LAC50'


# FUNCTIONS:

#. calculate AC50: antilog LAC50
calc_ac50 <- function(lac50){
    return(10^6*10^lac50)
}

#. calculate LAC50: log10(AC50/10^6)
calc_lac50 <- function(ac50, negative=TRUE){
    if(negative){
        return(-log10(ac50/10^6))
        }
        else {
            return(log10(ac50/10^6))
        }
}

# RUN:

# read sample ids and files
samples_ids <- read.csv(file_in_samples_ids)
samples_list <- samples_ids$SAMPLE_FILE
sample_id <- samples_ids$SAMPLE_ID

# read sample data to list
dat_in_list <- list()
for (i in 1:length(samples_list)){
    dat_in_list[[i]] <- read.csv(file.path(dir_in, samples_list[i]))
}

# filter out incomplete and align data rows (compounds)
rows_any <- sapply(dat_in_list, function(x) paste(x[,row_var]))
temp <- table(rows_any)
rows_in_all <- names(temp)[temp == length(dat_in_list)]
for (i in 1:length(dat_in_list)){
    temp <- dat_in_list[[i]]
    temp <- temp[ match(rows_in_all, temp[,row_var]), ]
    dat_in_list[[i]] <- temp
}

# write out compound ids
temp = dat_in_list[[1]]
temp = temp[, match(compounds_var, colnames(temp))]
write.csv(temp, file_out_compounds_ids, row.names=F, quote=T)

# write out drug doses
for (i in 1:length(dat_in_list)){
    temp = dat_in_list[[i]]
    temp = temp[ ,match(doses_var, colnames(temp))]
    temp = data.frame(SID = rows_in_all, temp)
    write.csv(temp, file.path(dir_out_dose, paste(sample_id[i], 'dose.csv', sep='_')),
        row.names=F, quote=T)
}
# write out response (normalized percentage of viable cells) - multiple files or ?RData
for (i in 1:length(dat_in_list)){
    temp = dat_in_list[[i]]
    temp = temp[ ,match(response_var, colnames(temp))]
    temp = data.frame(SID = rows_in_all, temp)
    write.csv(temp, file.path(dir_out_response, paste(sample_id[i], 'response.csv', sep='_')),
        row.names=F, quote=T)
}

# write out drc (dose response curve fit and sensitivity parameters)
for (i in 1:length(dat_in_list)){
    temp <- dat_in_list[[i]]
    max_dose <- temp[ ,replace_ac50]
    temp = temp[ ,match(drc_var, colnames(temp))]
    AC50 <- calc_ac50(temp[,log_ac50])
    iAC50 <- ifelse(temp[ ,curve_class2]==is_bad & is.na(temp[ ,log_ac50]), max_dose, AC50)
    iLAC50 = calc_lac50(iAC50, negative=TRUE)
    temp = data.frame(SID = rows_in_all, temp, AC50, iAC50, iLAC50)
    write.csv(temp, file.path(dir_out_drc, paste(sample_id[i], 'drc.csv', sep='_')),
        row.names=F, quote=T)
}
