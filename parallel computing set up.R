install.packages("doParalle")
library(doParallel)
number_of_cores=detectCores()
myCluster <- makeCluster(number_of_cores-1, # number of cores to use
                         type = "PSOCK") # type of cluster
#notice that we need to leave one core for computer system
registerDoParallel(myCluster)
