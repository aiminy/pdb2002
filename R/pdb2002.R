# input.pdb.dir="~/Dropbox (BBSR)/Aimin_project/Noula_Shembade_protein_structure_prediction_models/cxcr2-k242R"
#
# input.pdb.dir="~/Dropbox (BBSR)/Aimin_project/Noula_Shembade_protein_structure_prediction_models/ModelCom"
# res <- pdb2002:::preparePdbFile(input.pdb.dir,"*.pdb$")

preparePdbFile <- function(input.pdb.dir,file.pattern){

  bam <- list.files(path = input.pdb.dir, pattern=file.pattern, all.files = TRUE,
                    full.names = TRUE, recursive = TRUE,
                    ignore.case = FALSE, include.dirs = TRUE, no.. = FALSE)
  #bam2 <- as.data.frame(cbind(basename(bam),bam))

  ID <- unlist(lapply(bam,function(u)
  {
    y <- basename(u)
    pos <- regexpr("\\.", y)
    pos <- pos - 1
    y <- substr(y, 1, pos)
    protein <- basename(dirname(u))
    id<-paste0(protein,"-",y)
    id
  }))

  pdb <- cbind(ID,bam)

  pdb.str <- lapply(pdb[,2],function(u){
    pdb <- read.pdb(u)
    pdb
  })

  names(pdb.str) <-  pdb[,1]

  pdb.xyz <- lapply(pdb.str,function(u){
    xyz <- u$xyz
    xyz
  })

  re <-list(pdb=pdb,pdb.str=pdb.str,pdb.xyz=pdb.xyz)
  re
}

# pdb2002:::analysis1(res)
#
analysis1 <- function(res) {
  # pdb<- res$pdb.str[[1]]
  # attributes(pdb)
  # pdb$atom[1:3, c("resno", "resid", "elety", "x", "y", "z")]
  # ca.inds <- atom.select(pdb, "calpha")
  # resnos <- pdb$atom[ca.inds$atom, "resno"]
  # bfacts <- pdb$atom[ca.inds$atom, "b"]
  # plot.bio3d(resnos, bfacts, sse = pdb, ylab = "B-factor", xlab = "Residue", typ = "l")
  #
  # loop <- pdb$sheet$end[3]:pdb$helix$start[1]
  # loop.inds <- atom.select(pdb, resno = loop, elety = "CA")
  #
  #
  # ids <- c("1TND_B", "1AGR_A", "1FQJ_A", "1TAG_A", "1GG2_A", "1KJY_A")
  # raw.files <- get.pdb(ids)
  # files <- pdbsplit(raw.files, ids)

  pdbs <- pdbaln(res$pdb.str)
  pdbs$id <- names(res$pdb.str)
  seqidentity(pdbs)

  #pdb.xyz <- res$pdb.xyz

  rd<-rmsd(pdbs,fit = TRUE)
  hist(rd, breaks = 40, xlab = "RMSD (Å)")
  hc.rd <- hclust(as.dist(rd))
  plot(hc.rd, labels = pdbs$id, xlab="",ylab = "RMSD", main = "RMSD Cluster Dendrogram between models")

  #pdbs2 <- res$pdb.str
  #rmsd(pdbs2$xyz, fit = TRUE)
}

analysis2 <- function() {
  myfiles <- list.files(mypath, ".pdb$", full.names=TRUE)
  xyz <- NULL
  for(i in 1:length(myfiles)) {
    xyz <- rbind(xyz, read.pdb(myfiles[i])$xyz)
  }
  pdb <- read.pdb(myfiles[1])
  ca.inds <- atom.select(pdb, "calpha")
  rmsd(xyz[,ca.inds$xyz], fit=TRUE)
}

# pdb2002:::calculateRmsd (res)

calculateRmsd <- function(res) {
pdbs <- pdbaln(res$pdb.str)
pdbs$id <- names(res$pdb.str)
seqidentity(pdbs)

#pdb.xyz <- res$pdb.xyz

rms<-rmsd(pdbs,fit = TRUE)
rms
}
