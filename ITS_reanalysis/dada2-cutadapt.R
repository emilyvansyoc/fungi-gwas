# remove primers if possible
# EVS 12/2024

suppressMessages(library(dada2))
suppressMessages(library(ShortRead))


# use external arguments to read in pathway to raw data and file of primers if available
args <- commandArgs(TRUE)
path_to_files <- args[[1]] # needs the path to raw files
primer_file <- args[[2]] # needs a fasta file of primers from the original publication


# get files in dada2 format
fnFs <- sort(list.files(path_to_files, pattern = "_1.fastq", full.names = TRUE))
fnRs <- sort(list.files(path_to_files, pattern = "_2.fastq", full.names = TRUE))

# remove N's 
fnFs.filtN <- file.path(path_to_files, "filtN", basename(fnFs))
fnRs.filtN <- file.path(path_to_files, "filtN", basename(fnRs))
filterAndTrim(fnFs, fnFs.filtN, fnRs, fnRs.filtN, maxN = 0, multithread = TRUE, verbose = TRUE)

# get primers; IF THERE IS NO PRIMER FILE, PRINT WARNING MESSAGE AND FINISH SCRIPT WITHOUT ERRORING
tryCatch({
  primers <- readDNAStringSet(primer_file)
  
  # get forward and reverse primers
  FWD <- unname(as.character(primers[names(primers) == "FWD"]))
  REV <- unname(as.character(primers[names(primers) == "REV"]))
  
  # dada2 function to get all primer orientations
  allOrients <- function(primer) {
    # Create all orientations of the input sequence
    #require(Biostrings)
    dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
    orients <- c(Forward = dna, 
                 Complement = Biostrings::complement(dna), 
                 Reverse = Biostrings::reverse(dna),
                 RevComp = Biostrings::reverseComplement(dna))
    return(sapply(orients, toString))  # Convert back to character vector
  }
  FWD.orients <- allOrients(FWD)
  REV.orients <- allOrients(REV)
  
  ### do standard ITS-friendly primer removal (RevComp orientation in FWD.ReverseReads and REV.ForwardReads)
  cutadapt <- "/Users/epb5360/micromamba/envs/itsenv/bin/cutadapt"
  system2(cutadapt, args = "--version") # sanity check to make sure shell can access cutadapt
  
  path.cut <- file.path(path_to_files, "cutadapt")
  if(!dir.exists(path.cut)) dir.create(path.cut)
  fnFs.cut <- file.path(path.cut, basename(fnFs.filtN))
  fnRs.cut <- file.path(path.cut, basename(fnRs.filtN))
  
  FWD.RC <- dada2:::rc(FWD)
  REV.RC <- dada2:::rc(REV)
  # Trim FWD and the reverse-complement of REV off of R1 (forward reads)
  R1.flags <- paste("-g", FWD, "-a", REV.RC) 
  # Trim REV and the reverse-complement of FWD off of R2 (reverse reads)
  R2.flags <- paste("-G", REV, "-A", FWD.RC) 
  # Run Cutadapt
  for(i in seq_along(fnFs.filtN)) {
    system2(cutadapt, args = c("--quiet", R1.flags, R2.flags, "-n", 2, # -n 2 required to remove FWD and REV from reads
                               "-o", fnFs.cut[i], "-p", fnRs.cut[i], # output files
                               fnFs.filtN[i], fnRs.filtN[i])) # input files
  }
  
  },
  
  error = function(e){cat("\n WARNING: primer file does not exist or otherwise fails primer removal \n")}
  
)


