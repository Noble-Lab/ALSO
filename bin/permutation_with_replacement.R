#install if necessary
install.packages('gtools')
#load library
library(gtools)
#urn with 3 balls
x <- c('red', 'blue', 'black')
#pick 2 balls from the urn with replacement
#get all permutations
permutations(n=3,r=2,v=x,repeats.allowed=T)
#      [,1]    [,2]   
# [1,] "black" "black"
# [2,] "black" "blue" 
# [3,] "black" "red"  
# [4,] "blue"  "black"
# [5,] "blue"  "blue" 
# [6,] "blue"  "red"  
# [7,] "red"   "black"
# [8,] "red"   "blue" 
# [9,] "red"   "red"
#number of permutations
nrow(permutations(n=3,r=2,v=x,repeats.allowed=T))
#[1] 9