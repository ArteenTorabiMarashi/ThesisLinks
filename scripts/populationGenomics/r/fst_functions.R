###making a function to do the numbering for our data 

chrNumbering <- function(dat, chrom) {
  
  datX <- subset(dat, chrom == "X")
  a <- dim(datX)[1]
  datX$number <- 1:a
  
  dat2L <- subset(dat, chrom == "2L")
  b <- dim(dat2L)[1]
  dat2L$number <- (a+1):(a+b)
  
  dat2R <- subset(dat, chrom == "2R")
  c <- dim(dat2R)[1]
  dat2R$number <- (a+b+1):(a+b+c)
  
  dat3L <- subset(dat, chrom == "3L")
  d <- dim(dat3L)[1]
  dat3L$number <- (a+b+c+1):(a+b+c+d)
  
  dat3R <- subset(dat, chrom == "3R")
  e <- dim(dat3R)[1]
  dat3R$number <- (a+b+c+d+1):(a+b+c+d+e)
  
  dat4 <- subset(dat, chrom == "4")
  f <- dim(dat4)[1]
  dat4$number <- (a+b+c+d+e+1):(a+b+c+d+e+f)
  
  return(rbind(datX, dat2L, dat2R, dat3L, dat3R, dat4))
  
}

middleChr <- function(dat, chr) {
  a <- dim(subset(dat, chrom == "X"))[1]
  b <- dim(subset(dat, chrom == "2L"))[1] 
  c <- dim(subset(dat, chrom == "2R"))[1] 
  d <- dim(subset(dat, chrom == "3L"))[1] 
  e <- dim(subset(dat, chrom == "3R"))[1] 
  f <- dim(subset(dat, chrom == "4"))[1]
  
  result <- as.vector(rep(NA, 6))
  result[1] <- a/2 
  result [2] <- a + (b/2)
  result[3] <- a + b + (c/2)
  result[4] <- a + b + c + (d/2)
  result[5] <- a + b + c + d + (e/2)
  result[6] <- a + b + c + d + e + (f/2)
  return(result)
  
}

#for cleaning up popoolation Fst calculation 

populaionFst_cleanup <- function(dat, x = "") {
  
  ccol <- ncol(dat)
  
  #this next for loop is getting rid of all of the comparison info that is in each column witht the fst values. The gsub command takes everything before the = (the .* means everything with the . being an escape from the wildcard) and replaces it with nothing.
  for (i in 6:ccol){
    dat[[i]] <- gsub(".*=","", dat[[i]])
  }
  
  #to change the columns to numeric so we can eventually get the mean
  for (i in 6:ccol){
    dat[[i]] <- as.numeric(dat[[i]])
  }
  
  
  
  colnames(dat) <- c('chr', 'window', "num", 'frac', 'meanCov', x)
  
  return(dat)
}

