comb_with_replacement <- function(n, r){
  return( factorial(n + r - 1) / (factorial(r) * factorial(n - 1)) )
}

#have 3 elements, choosing 3
comb_with_replacement(3,3)
#[1] 10