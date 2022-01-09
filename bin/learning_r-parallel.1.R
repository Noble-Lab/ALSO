
lapply(  # Input: 1,2,3,4,5; output: square of input
  1:5,
  function(x) x^2
)

# We can also generate multiple outputs, say x^2 and x^3
lapply(  # Output: square and cube of input
  1:5,
  function(x) c(x^2, x^3)
)
# This output is very large but sometimes a list is useful. 

#  An alternative is to use the sapply() function which generates a vector,
#+ matrix or array output.
sapply(  # This output is a vector
  1:5,
  function(x) x^2
)

sapply(  # This outputs a matrix
  1:5,
  function(x) c(x^2, x^3)
)

#  sapply() also provides two additional parameters: simplify and USE.NAMES.
#+ If both are kept FALSE, the output generated is the same as lappy()
sapply(  #  Output is same as for lapply()
  1:5,
  function(x) x^2,
  simplify = FALSE,
  USE.NAMES = FALSE
) 


# -------------------------------------
#  Include the parallel library
library(magrittr)
library(parallel)
#  The “parallel” packageR can perform tasks in parallel by providing the
#+ ability to allocate cores to R.
#+ 
#+ The workflow involves finding the number of cores in the system and
#+ allocating all of them or a subset to make a cluster.
#+ 
#+ We can then use the parallel version of various functions and run them by
#+ passing the cluster as an additional argument.
#+ 
#+ A word of caution: It is important to close the cluster at the end of
#+ execution step so that core memory is released.


#  Use the detectCores() function to find the number of cores in system
no_cores <- detectCores()
clust <- makeCluster(no_cores - 4)

#  The parallel version of lapply() is parLapply(), and it needs an additional
#+ cluster argument
parLapply(
  clust,
  1:5,
  function(x) c(x^2,x^3)
)

stopCluster(clust)

#  If we want a similar output but with sapply(), we use parSapply(). First,
#+ set a base variable 
base <- 4
#  This line is required so that all cores in the cluster have this variable
#+ available

no_cores <- detectCores()
clust <- makeCluster(no_cores - 4)

#  This line is required so that all cores in cluster have this variable available
clusterExport(clust, "base")
#  clusterExport() is a special command needed to change the variable scope in
#+ parallel execution. Normally, variables such as "base" have a scope which
#+ does not allow them to be accessible at all cores. We need to use the
#+ clusterExport() function and send the variable to all the assigned cores in
#+ the cluster. This is why we pass both the cluster variable as well as the
#+ variable we need to export. Changing the base variable after export will
#+ have no effect as all the cores will not see that change.

#  Using the parSapply() function
parSapply(
  clust,
  1:5,
  function(exponent) base^exponent
)

stopCluster(clust)
