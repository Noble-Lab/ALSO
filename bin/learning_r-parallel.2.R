#!/usr/bin/env Rscript

#  Following the tutorial at this URL:
#+ blasbenito.com/post/02_parallelizing_loops_with_r/
list.of.packages <- c(
    "foreach",
    "doParallel",
    "ranger",
    "palmerpenguins",
    "tidyverse",
    "kableExtra"
)

# Loading packages
for(package.i in list.of.packages){
    suppressPackageStartupMessages(
        library(
            package.i, 
            character.only = TRUE
        )
    )
}

#  Beyond for: building loops with foreach ------------------------------------
# --
x <- vector()
for(i in 1:10){
    x[i] <- sqrt(i)
}
x

# --
x <- foreach(i = 1:10) %do% {
    sqrt(i)
}
x

# --
x <- foreach(i = 1:10, .combine = 'c') %do% {
    sqrt(i)
}
x

x <- foreach(
    i = 1:3, 
    j = 1:3, 
    k = 1:3, 
    .combine = 'c'
) %do% {
    i + j + k
}
x


#  Running foreach loops in parallel ------------------------------------------
# --
x <- foreach(i = 1:10, .combine = 'c') %dopar% {
    sqrt(i)
}  ## Warning: executing %dopar% sequentially: no parallel backend registered


#  Setup of a parallel backend ------------------------------------------------
# Setup for a single computer
parallel::detectCores()
n.cores <- parallel::detectCores() - 1

# Create the cluster
my.cluster <- parallel::makeCluster(
    n.cores, 
    type = "PSOCK"
)

# Check cluster definition (optional)
print(my.cluster)

# Register the cluster to be used by %dopar%
doParallel::registerDoParallel(cl = my.cluster)

# Check that it is indeed registered (optional)
foreach::getDoParRegistered()

# Check how many workers are available (optional)
foreach::getDoParWorkers()

# Run a set of tasks in parallel
x <- foreach(i = 1:10, .combine = 'c') %dopar% {
    sqrt(i)
}
x  # It works, but no gain in speed b/c the operation is very simple

# Stop the cluster when we are done working with it
parallel::stopCluster(cl = my.cluster)


#  Setup for a Beowulf cluster ------------------------------------------------
# Skip


#  Practical examples ---------------------------------------------------------
# Load example data
data("penguins")

penguins <- as.data.frame(
    na.omit(
        penguins[, c(
            "species",
            "bill_length_mm",
            "bill_depth_mm",
            "flipper_length_mm",
            "body_mass_g"
        )]
    )
)

# Fit classification model
m <- ranger::ranger(
    data = penguins,
    dependent.variable.name = "species",
    importance = "permutation"
)

# Summary
m

# Variable importance
m$variable.importance


#  Tuning random forest hyperparameters ---------------------------------------
sensitivity.df <- expand.grid(
    num.trees = c(500, 1000, 1500),
    mtry = 2:4,
    min.node.size = c(1, 10, 20)
)

# Create and register cluster
my.cluster <- parallel::makeCluster(n.cores)
doParallel::registerDoParallel(cl = my.cluster)

# Fit each rf model with different hyperparameters
prediction.error <- foreach(
    num.trees = sensitivity.df$num.trees,
    mtry = sensitivity.df$mtry,
    min.node.size = sensitivity.df$min.node.size,
    .combine = 'c', 
    .packages = "ranger"
) %dopar% {
    # Fit model
    m.i <- ranger::ranger(
        data = penguins,
        dependent.variable.name = "species",
        num.trees = num.trees,
        mtry = mtry,
        min.node.size = min.node.size
    )
    
    # Return prediction error as percentage
    return(m.i$prediction.error * 100)
}

# Add a prediction error column
sensitivity.df$prediction.error <- prediction.error

# Plot the results
ggplot2::ggplot(data = sensitivity.df) + 
    ggplot2::aes(
        x = mtry,
        y = as.factor(min.node.size),
        fill = prediction.error
    ) + 
    ggplot2::facet_wrap(as.factor(sensitivity.df$num.trees)) +
    ggplot2::geom_tile() + 
    ggplot2::scale_y_discrete(breaks = c(1, 10, 20)) +
    ggplot2::scale_fill_viridis_c() + 
    ggplot2::ylab("min.node.size")

best.hyperparameters <- sensitivity.df %>% 
    dplyr::arrange(prediction.error) %>% 
    dplyr::slice(1)


#  Confidence intervals of variable importance scores -------------------------
importance_to_df <- function(model){
    #  Transform the vector of importance scores returned by ranger into a data
    #+ frame (of one row)
    x <- as.data.frame(model$variable.importance)
    x$variable <- rownames(x)
    colnames(x)[1] <- "importance"
    rownames(x) <- NULL
    return(x)
}

# No need to create the cluster; it is still up
print(my.cluster)  ## socket cluster with 7 nodes on host 'localhost'

# Asses execution time
system.time(
    # Perform 1000 iterations in parallel
    importance.scores <- foreach(
        i = 1:1000, 
        .combine = 'rbind', 
        .packages = "ranger"
    ) %dopar% {
        # Fit model
        m.i <- ranger::ranger(
            data = penguins,
            dependent.variable.name = "species",
            importance = "permutation",
            mtry = best.hyperparameters$mtry,
            num.trees = best.hyperparameters$num.trees,
            min.node.size = best.hyperparameters$min.node.size
        )
        
        #format importance
        m.importance.i <- importance_to_df(model = m.i)
        
        #returning output
        return(m.importance.i)
    }
)
#  user  system elapsed 
# 0.396   0.065   8.905

ggplot2::ggplot(data = importance.scores) + 
    ggplot2::aes(
        y = reorder(variable, importance), 
        x = importance
    ) +
    ggplot2::geom_boxplot() + 
    ggplot2::ylab("")

# --
# List to save results
importance.scores.list <- list()

# Performing 1000 iterations sequentially
system.time(
    for(i in 1:1000){
        # Fit model
        m.i <- ranger::ranger(
            data = penguins,
            dependent.variable.name = "species",
            importance = "permutation",
            seed = i,
            num.threads = parallel::detectCores() - 1
        )
        
        # Format importance
        importance.scores.list[[i]] <- importance_to_df(model = m.i)
    }
)

parallel::stopCluster(cl = my.cluster)
