library(dplyr)
###################
## Downsample
##################

## data: An input data frame. Must contain the following columns: CDR3.amino.acid.sequence, Read.count. Each row is a unique CDR3 amino acid sequence.
## B: Number of downsampled experiments
## target_library_size: Downsampled to this target_library_size
## seed: A random seed

downsample <- function(data, B = 100, target_library_size, seed = 12345){
    downsample.ls <- NULL
    for (b in 1:B) {
        
        set.seed(seed+b-1)
        downsampled_data <- table(sample(data$CDR3.amino.acid.sequence,
                                         target_library_size,
                                         prob = data$Read.count/sum(data$Read.count),replace = T))
        downsample.ls[[b]] <- as.data.frame(downsampled_data) %>% 
            dplyr::rename(CDR3.amino.acid.sequence = Var1,Read.count = Freq) %>%
            mutate(CDR3.amino.acid.sequence = as.character(CDR3.amino.acid.sequence))
        
    }
    downsample.ls
}
