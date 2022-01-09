#  Test: Parallelization of evaluations
chromosome <- c(paste0("chr", 1:19), "chrX", "chrY", "chrM")
reorderRnameLevels <- function(variable, chromosome) {
    command <- paste0("<- forcats::fct_relevel(", variable, "$rname, chromosome)")
    operation <- makeOperation(paste0(variable, "$rname"), command)
    return(eval(parse(text = operation), envir = globalenv()))
}

#  Test: Does this work?
# reorderRnameLevels(variable = variable, chromosome = chromosome)
# GB.alt$rname %>% head(., n = 1)  # It works!

#  Test: Does this work?
# lapply(list(variable), reorderRnameLevels, chromosome) %>% invisible()
# GB.alt$rname %>% head(., n = 1)  # It works!

#  Test: Does this work?
# cores <- parallel::detectCores(logical = FALSE) - (parallel::detectCores(logical = FALSE) / 2)
# parallel::mclapply(
#     list(variable), reorderRnameLevels, chromosome, mc.cores = cores
# ) %>% invisible()
# GB.alt$rname %>% head(., n = 1)  # It works!
# GB.ref$rname %>% head(., n = 1)  # It works!


# # -----------------------------------------------------------------------------
# #  Load tibbles from .rds files -----------------------------------------------
# # command <- paste0(variable, " <- loadRDS(", variable, ", file = \"", variable, ".chr1.rds\")")
# # eval(parse(text = command))
# 
# loadTibbleFromRDS <- function(variable) {
#     command <- paste0(variable, " <- loadRDS(", variable, ", file = \"", variable, ".chr1.rds\")") %>% as.list()
#     return(eval(parse(text = command), envir = globalenv()))
# }
# 
# # loadTibbleFromRDS(variable = variable)  # 26.43 seconds
# 
# cores <- parallel::detectCores(logical = FALSE) - (parallel::detectCores(logical = FALSE) / 2)
# # parallel::mclapply(list(variable), loadTibbleFromRDS, mc.cores = 4)  # 24.78 seconds
# parallel::mclapply(variable, loadTibbleFromRDS, mc.cores = 4)  # 29.58 seconds
# -----------------------------------------------------------------------------


# -----------------------------------------------------------------------------
cores <- 4
cluster <- parallel::makeCluster(
    cores, 
    type = "FORK"
)
doParallel::registerDoParallel(cl = cluster)

# print(cluster)

# foreach::getDoParRegistered()
# foreach::getDoParWorkers()

foreach::foreach(i = 1:length(file)) %dopar% {
    loadTibbleFromRDS(variable = variable[i], file = file[i])
}

parallel::stopCluster(cluster)
# -----------------------------------------------------------------------------
