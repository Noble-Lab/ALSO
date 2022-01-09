#install if necessary
install.packages('gtools')
#load library
library(gtools)
#urn with 3 balls
x <- c('red', 'blue', 'black')
#pick 2 balls from the urn with replacement
#get all permutations
permutations(n=3,r=2,v=x)
#     [,1]    [,2]   
#[1,] "black" "blue" 
#[2,] "black" "red"  
#[3,] "blue"  "black"
#[4,] "blue"  "red"  
#[5,] "red"   "black"
#[6,] "red"   "blue"
#number of permutations
nrow(permutations(n=3,r=2,v=x))
#[1] 6