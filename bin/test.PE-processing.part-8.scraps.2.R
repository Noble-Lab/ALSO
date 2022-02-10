
cigar <- bam$cigar.KA.odd
cigar2 <- c(bam$cigar.KA.odd, bam$cigar.KA.even)

headTail <- function(x) {
    y <- x %>% head(10) %>% tail()
    return(y)
}

explodeCigarOps(cigar2) %>% headTail()
explodeCigarOps(cigar2, ops = c("I", "S")) %>% headTail()
explodeCigarOps(cigar2, ops = c("=")) %>% headTail()
explodeCigarOps(cigar2, ops = c("X")) %>% headTail()
explodeCigarOps(cigar2, ops = c("D")) %>% headTail()

explodeCigarOpLengths(cigar2) %>% headTail()
explodeCigarOpLengths(cigar2, ops = c("I", "S")) %>% headTail()
explodeCigarOpLengths(cigar2, ops = c("=")) %>% headTail()
explodeCigarOpLengths(cigar2, ops = c("X")) %>% headTail()
explodeCigarOpLengths(cigar2, ops = c("D")) %>% headTail()

explodeCigarOps(cigar2, ops = CIGAR_OPS) %>% headTail()
explodeCigarOpLengths(cigar2, ops = CIGAR_OPS) %>% headTail()
tmp <- cigarToRleList(cigar2) %>% headTail()


bam$cigar.KA.odd %>% headTail()
bam$read_sequence.KA.odd %>% headTail()
bam$matches.KA.odd %>% headTail()
bam$reference_sequence.KA.odd %>% headTail()
