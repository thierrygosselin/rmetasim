#These are functions that manipulate the landscape to produce summary statistics
#implemented in the 'ape' or 'pegas' packages on CRAN.
#pegas and ape must be installed and loaded for these to work.
#interface to theta.h
landscape.theta.h <- function(rland)
  {
    retval <- matrix(0,ncol=length(rland$loci),nrow=rland$intparam$habitats)
    for (i in 1:rland$intparam$habitats)
      {
        rland.tmp <- rland
        rland.tmp$individuals <- rland.tmp$individuals[landscape.populations(rland.tmp)==i,]
        for (j in 1:length(rland$loci))
          {
            alleledist <- as.factor(landscape.locus(rland.tmp,lnum=j)[,c(-1:-(landscape.democol()))])
            if (length(unique(alleledist))>1)
              retval[i,j] <- pegas::theta.h(alleledist)
            else
              retval[i,j] <- NA
          }
      }
    retval
  }

#interface to theta.k
landscape.theta.k <- function(rland)
  {
    retval <- matrix(0,ncol=length(rland$loci),nrow=rland$intparam$habitats)
    for (i in 1:rland$intparam$habitats)
      {
        rland.tmp <- rland
        rland.tmp$individuals <- rland.tmp$individuals[landscape.populations(rland.tmp)==i,]
        for (j in 1:length(rland$loci))
          {
            alleledist <- as.factor(landscape.locus(rland.tmp,lnum=j)[,c(-1:-(landscape.democol()))])
            if (length(unique(alleledist))>1)
              retval[i,j] <- pegas::theta.k(alleledist)
            else
              retval[i,j] <- NA
          }
      }
    retval
  }

# Extract every sampled gene copy, not just unique allele states.
.landscape.sequence.summary <- function(rland, tajima = FALSE) {
    retval <- matrix(NA_real_, nrow = rland$intparam$habitats,
                     ncol = length(rland$loci))
    populations <- rland$individuals[, 1] %/% rland$intparam$stages + 1L
    locusvec <- landscape.locusvec(rland)
    for (j in seq_along(rland$loci)) {
        if (rland$loci[[j]]$type != 253) next
        alleles <- rland$loci[[j]]$alleles
        indices <- vapply(alleles, function(a) as.character(a$aindex), character(1))
        states <- vapply(alleles, function(a) as.character(a$state), character(1))
        cols <- which(locusvec == j) + landscape.democol()
        for (i in seq_len(rland$intparam$habitats)) {
            sampled <- as.vector(rland$individuals[populations == i, cols, drop = FALSE])
            n <- length(sampled)
            if (n < 2L) next
            matched <- match(as.character(sampled), indices)
            if (anyNA(matched)) stop("Sampled sequence allele index is not in the allele table")
            sequences <- toupper(states[matched])
            widths <- nchar(sequences)
            if (anyNA(sequences) || length(unique(widths)) != 1L ||
                any(widths == 0L) || any(grepl("[^ACGT]", sequences))) {
                stop("Sequence summaries require equal-length, non-empty A/C/G/T sequences")
            }
            dna <- do.call(rbind, strsplit(sequences, "", fixed = TRUE))
            S <- sum(apply(dna, 2L, function(site) length(unique(site)) > 1L))
            if (!tajima) {
                # Watterson's theta per locus: n is the number of gene copies.
                retval[i, j] <- pegas::theta.s(S, n)
            } else if (n >= 4L && S > 0L) {
                # Standard Tajima's D; undefined small/monomorphic cases stay NA.
                D <- pegas::tajima.test(ape::as.DNAbin(dna))$D
                if (is.finite(D)) retval[i, j] <- D
            }
        }
    }
    retval
}

landscape.theta.s <- function(rland) {
    .landscape.sequence.summary(rland)
}

# Internal: standard Tajima's D, retaining sampled haplotype multiplicities.
landscape.tajima.d <- function(rland) {
    .landscape.sequence.summary(rland, tajima = TRUE)
}
