perm_without_replacement <- function(n, r){
  return(factorial(n)/factorial(n - r))
}

#16 choices, choose 16
perm_without_replacement(16,16)
#[1] 2.092279e+13

#16 choices, choose 3
perm_without_replacement(16,3)
#[1] 3360