#use chr22 as an example
#how many entries on the negative strand of chr22?
table(bam_df$rname == 'chr22' & bam_df$flag == 16)
# FALSE    TRUE 
#3875997   24413

#function for checking negative strand
check_neg <- function(x){
  if (intToBits(x)[5] == 1){
    return(T)
  } else {
    return(F)
  }
}

#test neg function with subset of chr22
test <- subset(bam_df, rname == 'chr22')
dim(test)
#[1] 56426    13
table(apply(as.data.frame(test$flag), 1, check_neg))
#number same as above
#FALSE  TRUE 
#32013 24413

#function for checking positive strand
check_pos <- function(x){
  if (intToBits(x)[3] == 1){
    return(F)
  } else if (intToBits(x)[5] != 1){
    return(T)
  } else {
    return(F)
  }
}

#check pos function
table(apply(as.data.frame(test$flag), 1, check_pos))
#looks OK
#FALSE  TRUE 
#24413 32013

#store the mapped positions on the plus and minus strands
chr22_neg <- bam_df[bam_df$rname == 'chr22' &
                    apply(as.data.frame(bam_df$flag), 1, check_neg),
                    'pos'
                   ]
length(chr22_neg)
#[1] 24413
chr22_pos <- bam_df[bam_df$rname == 'chr22' &
                    apply(as.data.frame(bam_df$flag), 1, check_pos),
                    'pos'
                   ]
length(chr22_pos)
#[1] 32013
 
#calculate the densities
chr22_neg_density <- density(chr22_neg)
chr22_pos_density <- density(chr22_pos)
 
#display the negative strand with negative values
chr22_neg_density$y <- chr22_neg_density$y * -1
 
plot(chr22_pos_density,
     ylim = range(c(chr22_neg_density$y, chr22_pos_density$y)),
     main = "Coverage plot of mapped CAGE reads",
     xlab = "Chromosome 22",
     col = 'blue',
     lwd=2.5)
lines(chr22_neg_density, lwd=2.5, col = 'red')