read_sync <- function (file, gen, repl, polarization = c("minor", "rising", 
                                            "reference"), keepOnlyBiallelic = FALSE) 
{
  polarization <- match.arg(polarization)
  cat("Reading sync file ...\n")
  syncDt <- fread(file, sep = " ", header = FALSE, stringsAsFactors = FALSE)
  setDT(syncDt)
  if ((ncol(syncDt) - 3)%%length(gen) != 0 || (ncol(syncDt) - 
                                               3)%%length(repl) != 0) 
    stop("Either 'gen' (", length(gen), ") or 'repl' (", 
         length(repl), ") is not a multiple of the number of populations (", 
         ncol(syncDt) - 3, ") specified in the sync-file.")
  cat("Extracting biallelic counts ...\n")
  cppFunction(plugins = c("cpp11"), {
    "NumericMatrix Sync2Cnts(CharacterVector sync) {\n    NumericMatrix resMat(sync.length(), 4);\n\n    std::string current;\n    std::string token;\n    size_t posBeg = 0, posEnd = 0;\n    int cnt = 0;\n\n    for(int i=0; i<sync.length(); i++) {\n\n    cnt = 0;\n    posBeg = 0;\n    posEnd = 0;\n    current = Rcpp::as<std::string>(sync(i));\n    while ((posEnd = current.find(\":\", posBeg)) != std::string::npos && cnt <= 3) {\n    token = current.substr(posBeg, posEnd);\n    resMat(i, cnt) = std::stoi(token);\n    posBeg = posEnd+1;\n    cnt++;\n    }\n  }\n\n    return resMat;\n}"
  })
  syncCnts <- lapply(syncDt[, -1:-3, with = FALSE], function(col) {
    Sync2Cnts(col)
  })
  chr <- as.character(syncDt$V1)
  pos <- syncDt$V2
  ref <- syncDt$V3
  popCnt <- ncol(syncDt) - 3
  rm(syncDt)
  gc()
  sumCnts <- Reduce("+", syncCnts)
  ikeep = if (keepOnlyBiallelic) 
    which(rowSums(sumCnts > 0) == 2)
  else which(rowSums(sumCnts > 0) > 1)
  sumCnts <- sumCnts[ikeep, ]
  chr <- chr[ikeep]
  pos <- pos[ikeep]
  ref <- ref[ikeep]
  alleleRank <- rowRanks(sumCnts + rep(1:4/5, each = nrow(sumCnts)))
  rm(sumCnts)
  gc()
  alleleCnts <- lapply(syncCnts, function(pop) {
    cbind(major = t(pop[ikeep, ])[t(alleleRank == 4)], minor = t(pop[ikeep, 
    ])[t(alleleRank == 3)])
  })
  rm(syncCnts)
  gc()
  cat("Creating result object ...\n")
  chrNames <- unique(chr)
  chrID <- 1:length(chrNames)
  names(chrID) <- chrNames
  syncCntCol <- 1:4
  names(syncCntCol) <- c("A", "T", "C", "G")
  alleles <- data.table(chr = chr, pos = pos, ref = ref, major = names(syncCntCol)[which(t(alleleRank) == 
                                                                                           4) - 4 * seq(0, nrow(alleleRank) - 1)], minor = names(syncCntCol)[which(t(alleleRank) == 
                                                                                                                                                                     3) - 4 * seq(0, nrow(alleleRank) - 1)], rising = NA_character_)
  rm(alleleRank)
  gc()
  if (polarization == "reference" && any(alleles$ref != alleles$minor & 
                                         alleles$ref != alleles$major)) {
    warning("Cannot polarize for reference allele, because it is not among the two most common alleles for some SNPs. Changing polarization to 'minor'.")
    polarization <- "minor"
  }
  popInfo <- data.table(pop = 1:popCnt, gen = gen, repl = repl)
  for (r in unique(repl)) {
    for (i in seq(1:nrow(popInfo))[popInfo$repl == r]) {
      seqCov <- rowSums(alleleCnts[[i]])
      if (polarization == "minor" || polarization == "rising") {
        alleles[, `:=`(paste0("F", popInfo$gen[i], ".R", 
                              r, ".freq"), alleleCnts[[i]][, "minor"]/seqCov)]
      }
      else {
        alleles[, `:=`(paste0("F", popInfo$gen[i], ".R", 
                              r, ".freq"), ifelse(minor == ref, alleleCnts[[i]][, 
                                                                                "minor"]/seqCov, alleleCnts[[i]][, "major"]/seqCov))]
      }
      alleles[, `:=`(paste0("F", popInfo$gen[i], ".R", 
                            r, ".cov"), seqCov)]
    }
  }
  rm(alleleCnts)
  gc()
  if (polarization == "rising" && length(unique(popInfo$gen)) > 
      1) {
    minGen <- min(popInfo$gen)
    maxGen <- max(popInfo$gen)
    meanAF <- rowMeans(alleles[, grepl(paste0("F", maxGen, 
                                              "\\.R[0-9]+\\.freq"), colnames(alleles)), with = FALSE], 
                       na.rm = TRUE) - rowMeans(alleles[, grepl(paste0("F", 
                                                                       minGen, "\\.R[0-9]+\\.freq"), colnames(alleles)), 
                                                        with = FALSE], na.rm = TRUE)
    changePolarization <- meanAF < 0
    for (pop in grep(".freq", colnames(alleles), value = TRUE, 
                     fixed = TRUE)) {
      alleles[[pop]][changePolarization] <- 1 - alleles[[pop]][changePolarization]
    }
    alleles$rising <- ifelse(changePolarization, alleles$major, 
                             alleles$minor)
    rm(meanAF, changePolarization)
    gc()
  }
  return(new(Class = "sync", gen = as.numeric(sub("F([0-9]+)\\.R[0-9]+.*", 
                                                  "\\1", colnames(alleles)[-1:-6])), repl = as.numeric(sub(".*\\.R([0-9]+)\\..*", 
                                                                                                           "\\1", colnames(alleles)[-1:-6])), isAF = grepl(".*\\.freq$", 
                                                                                                                                                           colnames(alleles)[-1:-6]), polarization = polarization, 
             alleles = alleles))
}
