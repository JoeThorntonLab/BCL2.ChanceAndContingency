####################################
###Analysis of BCL2 PACE HTS data###
####################################

library(ggplot2)
library(reshape2)
library(MullerPlot)
library(pegas)
library(stringr)
library(ape)
library(DECIPHER)
library(aphid)
library(gridExtra)
library(ggpubr)
library(insect)
library(rgexf)
library(igraph)

setwd("")

###############################
###Sample to library mapping###
###############################
EXPERIMENTS    <- c(rep("B",20),rep("P",4),rep("C",4),rep("D",4),rep("J",4),rep("K12",2),
                 rep("K34",2),rep("L",20),rep("M",4),rep("Q",2),rep("W",2),rep("X",4))
EXPERIMENTS.WT <- c(rep("B",20),rep("B",4),rep("C",4),rep("D",4),rep("J",4),rep("K12",2),
                 rep("K34",2),rep("L",20),rep("M",4),rep("B",2),rep("L",2),rep("B",4))
WT.REF         <- c("B","C","D","J","K12","K34","L","M")

########################
###Alignment Cleaning###
########################
###Round 1###
#Remove sequences with Ns and low frequency insertions
FILES <- list.files("./Alignments/01.Start")
I <- 1:length(FILES)
for(i in I) {
  DATA <- read.FASTA(file = paste("./Alignments/01.Start/",FILES[i],sep=""))
  DATA <- as.matrix(DATA)
  
  #Remove sequences with Ns
  J <- 1:nrow(DATA) 
  ANY.N <- numeric(length(J))
  for(j in J) {
    ANY.N[j] <- base.freq(DATA[j,],freq=TRUE,all=TRUE)[15]
  }
  DATA <- DATA[!(ANY.N>0),]
  
  #Calculate allele frequency 
  K <- 1:ncol(DATA)
  ALLELE.FREQ <- matrix(0,nrow=5,ncol=length(K))
  rownames(ALLELE.FREQ) <- c("A","C","G","T","-")
  colnames(ALLELE.FREQ) <- K
  for(k in K) {
   ALLELE.FREQ[,k] <- base.freq(DATA[,k],freq=TRUE,all=TRUE)[c(1:4,16)]
  }
  
  #Remove positions that are mostly gaps. These either don't have sequence or are likely sequencing errors
  GAPS.REMOVE <- which(ALLELE.FREQ[5,] > colSums(ALLELE.FREQ)*.99)
  if(length(GAPS.REMOVE) > 0) {
    DATA <- DATA[,-GAPS.REMOVE]
    ALLELE.FREQ <- ALLELE.FREQ[,-GAPS.REMOVE]
  }
  #Remove positions at ends with no reference data 
  GAPS.REMOVE <- (max(which(as.character(DATA[1,]) != "-")) +1):ncol(DATA)
  if(length(GAPS.REMOVE) > 0) {
    DATA <- DATA[,-GAPS.REMOVE]
    ALLELE.FREQ <- ALLELE.FREQ[,-GAPS.REMOVE]
  }
  #Remove sequences with gaps at the ends
  J <- 1:nrow(DATA) 
  ANY.END.GAP <- numeric(length(J))
  for(j in J) {
    ANY.END.GAP[j] <- as.character(DATA[j,1]) == "-" | tail(as.vector(as.character(DATA[j,])),n=1) == "-" 
  }
  DATA <- DATA[!(ANY.END.GAP),]

  write.FASTA(DATA,paste("./Alignments/02.Edit/",FILES[i],sep=""))
}  

###Round 2###
#Standardize gaps within an experiment
FILES <- list.files("./Alignments/02.Edit")
H <- 1:length(names(table(EXPERIMENTS)))
for(h in H) {
  DATA <- list()
  FILES.H <- (3*(which(EXPERIMENTS==names(table(EXPERIMENTS))[h])[1]) - 2):(3*tail(which(EXPERIMENTS==names(table(EXPERIMENTS))[h]),1))
  I <- 1:(3*table(EXPERIMENTS)[h])
  GAPS <- numeric(length(I))
  #Identify which sequences have gaps 
  for(i in I) {
    DATA[[i]] <- as.matrix(read.FASTA(file = paste("./Alignments/02.Edit/",FILES[FILES.H[i]],sep="")))
    GAPS[i]   <- sum(as.character(DATA[[i]][1,]) == "-")
  }
  #Standardize
  J <- 1:3
  for(j in J) {
    if(j == 1) {
      K <- seq(1,3*table(EXPERIMENTS)[h],3)
      if(length(which(GAPS[K] != 0)) > 0) {
        K <- K
      } else {
        K <- NULL
      }
    } else if(j == 2) {
      K <- seq(2,3*table(EXPERIMENTS)[h],3)
      if(length(which(GAPS[K] != 0)) > 0) {
        K <- K
      } else {
        K <- NULL
      }
    } else if(j == 3) {
      K <- seq(3,3*table(EXPERIMENTS)[h],3)
      if(length(which(GAPS[K] != 0)) > 0) {
        K <- K
      } else {
        K <- NULL
      }
    }
    if(length(K) > 0) {
      for(k in K[-which(K == tail(K,n=1))]) {
        DATA1 <- DNAStringSet(unlist(lapply(as.character(as.list(DATA[[k]])),paste0,collapse="")))
        M <- K[(which(K == k)+1):length(K)]
        if(length(M) > 0) {
          for(m in M) {
            DATA2 <- DNAStringSet(unlist(lapply(as.character(as.list(DATA[[m]])),paste0,collapse="")))
            DATA.ALIGN <- AlignProfiles(DATA1,DATA2)
            DATA1 <- DATA.ALIGN[1:length(DATA1),]
            DATA2 <- DATA.ALIGN[(length(DATA1)+1):(length(DATA1) + length(DATA2)),]
            DATA[[m]] <- as.DNAbin(DATA2)
          }
        }
        DATA[[k]] <- as.DNAbin(DATA1)
      }
    }
  }
  for(i in I) {
    write.FASTA(DATA[[i]],paste("./Alignments/03.GapRefined/",FILES[FILES.H[i]],sep=""))
  }
}  

###Round 3###
#Trim ends
FILES <- list.files("./Alignments/03.GapRefined")
CUT.POINTS <- read.table("Cut.Points.txt",header=TRUE)
I <- 1:length(FILES)
for(i in I) {
  DATA <- as.matrix(read.FASTA(file = paste("./Alignments/03.GapRefined/",FILES[i],sep="")))[,CUT.POINTS$Start[i]:CUT.POINTS$Stop[i]]
  write.FASTA(DATA,paste("./Alignments/04.EndTrim/",FILES[[i]],sep=""))
}

###Round 4###
#Calcuate allele frequencies
FILES <- list.files("./Alignments/04.EndTrim")

I <- 1:length(FILES)
for(i in I) {
  DATA <- as.matrix(read.FASTA(file = paste("./Alignments/04.EndTrim/",FILES[i],sep="")))
  K <- 1:ncol(DATA)
  ALLELE.COUNT <- matrix(0,nrow=5,ncol=length(K))
  rownames(ALLELE.COUNT) <- c("a","c","g","t","-")
  colnames(ALLELE.COUNT) <- K
  COUNTS <- numeric(length(K))
  WT.COUNT <- matrix(0,nrow=5,ncol=length(K))
  rownames(WT.COUNT) <- c("a","c","g","t","-")
  colnames(WT.COUNT) <- K
  MUTANT.FREQ <- numeric(length(K))
  names(MUTANT.FREQ) <- K
    
  for(k in K) {
    ALLELE.COUNT[,k] <- base.freq(DATA[,k],freq=TRUE,all=TRUE)[c(1:4,16)]
    COUNTS[k] <- sum(ALLELE.COUNT[,k])
    WT.COUNT[which(as.character(DATA[1,k])[1] == rownames(ALLELE.COUNT)),k] <- COUNTS[k]
    MUTANT.FREQ[k] <- 1 - ALLELE.COUNT[which(as.character(DATA[1,k])[1] == rownames(ALLELE.COUNT)),k]/COUNTS[k]
  }
  write.table(ALLELE.COUNT,file=paste("./Alignments/Allele.Count/",substr(FILES[i],1,30),".txt",sep=""),quote=FALSE,sep="\t")
  write.table(WT.COUNT,    file=paste("./Alignments/WT.Count/",    substr(FILES[i],1,30),".txt",sep=""),quote=FALSE,sep="\t")
  write.table(MUTANT.FREQ, file=paste("./Alignments/Mutant.Freq/", substr(FILES[i],1,30),".txt",sep=""),quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE)
}


################
####Analysis####
################
#Read in previously analyzed data
FILES <- list.files("./Alignments/WT.Count")
CUT.POINTS <- read.table("Cut.Points.txt",header=TRUE)
I <- 1:length(FILES)
ALLELE.COUNT <- list()
WT.COUNT <- list()
MUTANT.FREQ <- list()
for(i in I) {
  ALLELE.COUNT[[i]] <- as.matrix(read.table(paste("./Alignments/Allele.Count/", FILES[i],sep=""),header=TRUE,row.names = 1))
  colnames(ALLELE.COUNT[[i]]) <- 1:ncol(ALLELE.COUNT[[i]])
  WT.COUNT[[i]]     <- as.matrix(read.table(paste("./Alignments/WT.Count/",FILES[i],sep=""),header=TRUE,row.names = 1))
  colnames(WT.COUNT[[i]]) <- 1:ncol(WT.COUNT[[i]])
  MUTANT.FREQ[[i]] <- read.table(paste("./Alignments/Mutant.Freq/",FILES[i],sep=""))
}

######################
####AA Frequencies####
######################
#Read in processed data
FILES <- list.files("./Alignments/04.EndTrim")
DATA <- list()
I <- 1:length(FILES)
for(i in I) {
  DATA[[i]] <- as.matrix(read.FASTA(file = paste("./Alignments/04.EndTrim/",FILES[i],sep="")))
}

#Amino acid mutant Freqs for specific comparison
MCL1.TS <- 45:60
BCL2.TS <- 1:20
MCL1.BB <- c(33:36,39:40,57:60,61:64)
BCL2.BB <- c(17:20,25:28,29:32)
BCL2X.EXPS <- c(17:20,21:24,69:72)
BCL2.RE <- 65:66
MCL1.NOXA <- 67:68
REP.SUBS <- list()
DEL.SUBS <- list()
INS.SUBS <- list()

#J <- MCL1.TS
#J <- BCL2.TS
#J <- MCL1.BB
#J <- BCL2.BB
#J <- BCL2X.EXPS
#J <- BCL2.RE
#J <- c(BCL2.BB,MCL1.BB)

for(j in J) {
  I <- 1:3
  N.VARIANTS <- numeric(length(I))
  PLOTS <- list()
  INDEX <- 1
  TRANSLATION.LENGTH <- numeric(length(I))
  
  for(i in I) {
    SET <- (j-1)*3+i
    DATA.SET <- apply(as.character(DATA[[SET]]),2,function(x) table(factor(x,levels=c("a","c","g","t","-"))))
    INSERTIONS <- which(as.character(DATA[[SET]][1,]) == "-")
    DELETIONS  <- which(DATA.SET[5,]/colSums(DATA.SET) > .01)
    DELETIONS  <- DELETIONS[-which(DELETIONS %in% INSERTIONS)]
    INDELS <- DATA[[SET]][,c(INSERTIONS,DELETIONS)]; colnames(INDELS) <- c(INSERTIONS,DELETIONS)
    INDEL.HAPLOTYPES <- haplotype(INDELS)
    HIGH.FREQ.INDELS <- which(sapply(attr(INDEL.HAPLOTYPES,"index"),length)/nrow(INDELS) > 0.01)
    N.VARIANTS[i] <- length( HIGH.FREQ.INDELS)

    if(!1 %in% HIGH.FREQ.INDELS) {
      HIGH.FREQ.INDELS <- c(1,HIGH.FREQ.INDELS)
    }
    
    for(k in HIGH.FREQ.INDELS) {
      STRUCTURAL.VARIANT <- k
      VARIANT <- attr(INDEL.HAPLOTYPES,"index")[[STRUCTURAL.VARIANT]]
      
      INSERT <- NULL
      INSERT <- INSERTIONS[which(as.character(DATA[[SET]][VARIANT[1],INSERTIONS]) != "-")]

      DELETE <- NULL
      DELETE <- DELETIONS[which(as.character(DATA[[SET]][VARIANT[1],DELETIONS]) == "-")]
      
      N.INSERT <- length(which(INSERTIONS %in% INSERT))
      N.DELETE <- length(which(DELETIONS  %in% DELETE))
    
      WT <- unlist(as.character(DATA[[SET]][1,]))
      if(N.INSERT == length(INSERTIONS) & N.DELETE == 0) {
        WT.INSERT <- as.matrix(as.DNAbin(WT))
      } else {
        WT.INSERT <- as.matrix(as.DNAbin(WT[-c(INSERTIONS[which(!(INSERTIONS %in% INSERT))],DELETIONS[which(DELETIONS %in% DELETE)])]))
      } 
      
      VARIANT.NO.GAPS  <- del.gaps(DATA[[SET]][VARIANT,])
      VARIANT.NO.GAPS2 <- sapply(VARIANT.NO.GAPS,length)
      VARIANT.MATRIX   <- as.matrix(VARIANT.NO.GAPS[VARIANT.NO.GAPS2 == length(WT.INSERT)])
      VARIANT.MATRIX   <- rbind(WT.INSERT,VARIANT.MATRIX)
      
      DNA.COUNT <- apply(as.character(VARIANT.MATRIX),    2,function(x) table(factor(x,levels=c("a","c","g","t","-"))))
      DNA.WT    <- apply(as.character(VARIANT.MATRIX[1,]),2,function(x) table(factor(x,levels=c("a","c","g","t","-"))))
      
      INSERT.POS <- NULL
      INSERT.POS <- which(as.character(VARIANT.MATRIX[1,]) == "-")
      DELETE.POS <- NULL
      DELETE.POS <- DELETE

      DNA.FREQ <- DNA.COUNT/colSums(DNA.COUNT)
      DNA.FREQ2 <- melt(DNA.FREQ); colnames(DNA.FREQ2) <- c("NUC","SITE","FREQ")
      DNA.FREQ2$WT <- rep(as.character(VARIANT.MATRIX[1,]),each=5)
      DNA.WT2 <- DNA.FREQ2[!duplicated(DNA.FREQ2[,c('WT','SITE')]),]
      DNA.WT2$WT <- factor(DNA.WT2$WT,levels=c("a","c","g","t","-"))
      DNA.WT2$FREQ <- 1; DNA.WT2$WT   <- as.integer(DNA.WT2$WT); DNA.WT2$SITE <- as.integer(DNA.WT2$SITE)
      
      INDEL.COUNT <- matrix(0,nrow=5,ncol=ncol(DNA.COUNT)); rownames(INDEL.COUNT) <- c("A","C","G","T","-")
      INDEL.COUNT[,which(as.character(WT.INSERT) == "-")] <- DNA.COUNT[,which(as.character(WT.INSERT) == "-")]
      INDEL.COUNT[,which(DNA.FREQ[5,] > 0.05)] <-  DNA.COUNT[,which(DNA.FREQ[5,] > 0.05)]
      
      TRANSLATION <- trans(VARIANT.MATRIX)
      TRANSLATION.COUNT <- apply(as.character(TRANSLATION),    2,function(x) table(factor(x,levels=c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y","*"))))
      TRANSLATION.WT    <- apply(as.character(TRANSLATION[1,]),2,function(x) table(factor(x,levels=c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y","*"))))
      TRANSLATION.STOP  <- which(TRANSLATION.COUNT[21,]/colSums(TRANSLATION.COUNT) > 0.05)[1]
      if(k == 1 | (length(INSERT.POS) == 0 & length(DELETE.POS) == 0)) {
        COUNT.STOP <- ncol(TRANSLATION.COUNT)
      } else {
        COUNT.STOP <- ceiling(min(c(INSERT.POS,DELETE.POS)-1)/3)
      }
      TRANSLATION.COUNT <- TRANSLATION.COUNT[,1:COUNT.STOP]
      TRANSLATION.WT    <- TRANSLATION.WT[,1:COUNT.STOP]
    
      if(k == 1) {
        POINT.SUBS <- TRANSLATION.COUNT
      } else {
        POINT.SUBS[,1:ncol(TRANSLATION.COUNT)] <- POINT.SUBS[,1:ncol(TRANSLATION.COUNT)] + TRANSLATION.COUNT
      }
    }
    if(length(INSERTIONS) > 0) {
      T <- 1:length(INSERTIONS)
      INSERTION.SITE <- numeric(length(T))
      for(t in T) {
        INSERTION.SITE[t] <- INSERTIONS[t] - (t-1) + 1
      }
      if(i == 2) {
        INSERTION.SITE <- INSERTION.SITE + (CUT.POINTS[SET-1,]$Stop - CUT.POINTS[SET-1,]$Start) 
      }
      if(i == 3){
        INSERTION.SITE <- INSERTION.SITE + (CUT.POINTS[SET-1,]$Stop - CUT.POINTS[SET-1,]$Start) + (CUT.POINTS[SET-2,]$Stop - CUT.POINTS[SET-2,]$Start)
      }
    }
    
    if(i == 1) {
       REP.SUBS[[which(J == j)]] <- POINT.SUBS
       if(length(INSERTIONS) > 0) {
         DEL.SUBS[[which(J == j)]] <- DATA.SET[,-INSERTIONS]
         INS.SUBS[[which(J == j)]] <- DATA.SET[,INSERTIONS]
         colnames(INS.SUBS[[which(J == j)]]) <- INSERTION.SITE
       } else {
         DEL.SUBS[[which(J == j)]] <- DATA.SET
         INS.SUBS[[which(J == j)]] <- DATA.SET[,-(1:ncol(DATA.SET))]
       }
    } else {
       REP.SUBS[[which(J == j)]] <- cbind(REP.SUBS[[which(J == j)]],POINT.SUBS)
       if(length(INSERTIONS) > 0) {
         DEL.SUBS[[which(J == j)]] <- cbind(DEL.SUBS[[which(J == j)]],DATA.SET[,-INSERTIONS])
         INS.SUBS[[which(J == j)]] <- cbind(INS.SUBS[[which(J == j)]],DATA.SET[,INSERTIONS])
         colnames(INS.SUBS[[which(J == j)]]) <- cbind(colnames(INS.SUBS[[which(J == j)]]),INSERTION.SITE)
       } else {
         DEL.SUBS[[which(J == j)]] <- cbind(DEL.SUBS[[which(J == j)]],DATA.SET)
       }
    }
  }
  colnames(REP.SUBS[[which(J == j)]]) <- 1:ncol(REP.SUBS[[which(J == j)]])
  colnames(DEL.SUBS[[which(J == j)]]) <- 1:ncol(DEL.SUBS[[which(J == j)]])

  if(j %in% c(33:38,41:60,67:68)) {
    COL.NAMES <- which(as.numeric(colnames(REP.SUBS[[which(J == j)]])) %/% (96+1) > 0)
    colnames(REP.SUBS[[which(J == j)]])[COL.NAMES] <- as.numeric(colnames(REP.SUBS[[which(J == j)]])[COL.NAMES]) + 51
    
    COL.NAMES <- which(as.numeric(colnames(DEL.SUBS[[which(J == j)]])) %/% (3*96+1) > 0)
    colnames(DEL.SUBS[[which(J == j)]])[COL.NAMES] <- as.numeric(colnames(DEL.SUBS[[which(J == j)]])[COL.NAMES]) + 3*51
    
    COL.NAMES <- which(as.numeric(colnames(INS.SUBS[[which(J == j)]])) %/% (3*96+1) > 0)
    colnames(INS.SUBS[[which(J == j)]])[COL.NAMES] <- as.numeric(colnames(INS.SUBS[[which(J == j)]])[COL.NAMES]) + 3*51
  }
  if(j %in% c(1:32,39:40,61:66,69:72)) {
    COL.NAMES <- which(as.numeric(colnames(REP.SUBS[[which(J == j)]])) %/% (6+1) > 0)
    colnames(REP.SUBS[[which(J == j)]])[COL.NAMES] <- as.numeric(colnames(REP.SUBS[[which(J == j)]])[COL.NAMES]) + 57
    
    COL.NAMES <- which(as.numeric(colnames(DEL.SUBS[[which(J == j)]])) %/% (3*6+1) > 0)
    colnames(DEL.SUBS[[which(J == j)]])[COL.NAMES] <- as.numeric(colnames(DEL.SUBS[[which(J == j)]])[COL.NAMES]) + 3*57
    
    COL.NAMES <- which(as.numeric(colnames(INS.SUBS[[which(J == j)]])) %/% (3*96+1) > 0)
    colnames(INS.SUBS[[which(J == j)]])[COL.NAMES] <- as.numeric(colnames(INS.SUBS[[which(J == j)]])[COL.NAMES]) + 3*57
    
    COL.NAMES <- which(as.numeric(colnames(REP.SUBS[[which(J == j)]])) %/% (198+1) > 0)
    colnames(REP.SUBS[[which(J == j)]])[COL.NAMES] <- as.numeric(colnames(REP.SUBS[[which(J == j)]])[COL.NAMES]) + 1
    
    COL.NAMES <- which(as.numeric(colnames(DEL.SUBS[[which(J == j)]])) %/% (3*198+1) > 0)
    colnames(DEL.SUBS[[which(J == j)]])[COL.NAMES] <- as.numeric(colnames(DEL.SUBS[[which(J == j)]])[COL.NAMES]) + 3*1
    
    COL.NAMES <- which(as.numeric(colnames(INS.SUBS[[which(J == j)]])) %/% (3*96+1) > 0)
    colnames(INS.SUBS[[which(J == j)]])[COL.NAMES] <- as.numeric(colnames(INS.SUBS[[which(J == j)]])[COL.NAMES]) + 3*1
  }
}
      
WT.AA <- trans(read.FASTA("wt_ref DNA.txt",type="DNA"))   

for(j in J) {
  AA.COUNT <- REP.SUBS[[which(J == j)]]
  AA.COUNT <- melt(AA.COUNT)
  colnames(AA.COUNT) <- c("AA","SITE","COUNT")
  WT <- which(WT.REF == EXPERIMENTS.WT[j])
  AA.COUNT$WT <- rep(unlist(as.character(WT.AA[WT])),each=21)
  AA.COUNT$EXP <- rep(j,nrow(AA.COUNT))
  AA.COUNT <- AA.COUNT[,c(5,2,4,1,3)]
  
  AA.FREQ <- t(t(REP.SUBS[[which(J == j)]])/colSums(REP.SUBS[[which(J == j)]]))
  AA.FREQ <- melt(AA.FREQ)
  colnames(AA.FREQ) <- c("AA","SITE","FREQ")
  WT <- which(WT.REF == EXPERIMENTS.WT[j])
  AA.FREQ$WT <- rep(unlist(as.character(WT.AA[WT])),each=21)
  AA.FREQ$EXP <- rep(j,nrow(AA.FREQ))
  AA.FREQ <- AA.FREQ[,c(5,2,4,1,3)]
  
  DEL.FREQ <- t(t(DEL.SUBS[[which(J == j)]])/colSums(DEL.SUBS[[which(J == j)]]))
  DEL.FREQ <- melt(DEL.FREQ)
  colnames(DEL.FREQ) <- c("DNA","SITE","FREQ")
  DEL.FREQ$EXP <- rep(j,nrow(DEL.FREQ))
  DEL.FREQ <- DEL.FREQ[,c(4,2,1,3)]
  
  INS.FREQ <- t(t(INS.SUBS[[which(J == j)]])/colSums(INS.SUBS[[which(J == j)]]))
  INS.FREQ <- melt(INS.FREQ)
  colnames(INS.FREQ) <- c("DNA","SITE","FREQ")
  INS.FREQ$EXP <- rep(j,nrow(INS.FREQ))
  INS.FREQ <- INS.FREQ[,c(4,2,1,3)]
  
  if(which(J == j) == 1) {
    ALL.AA.SITES.COUNT <- AA.COUNT
    ALL.AA.SITES       <- AA.FREQ
    AA.SITES           <- AA.FREQ[which(AA.FREQ$FREQ > 0.05 & AA.FREQ$AA != AA.FREQ$WT),]
    DEL.SITES          <- DEL.FREQ[which(DEL.FREQ$FREQ > 0.05 & DEL.FREQ$DNA == "-"),]
    INS.SITES          <- INS.FREQ[which(INS.FREQ$FREQ > 0.05 & INS.FREQ$DNA != "-"),]
  } else {
    ALL.AA.SITES.COUNT <- rbind(ALL.AA.SITES.COUNT, AA.COUNT)
    ALL.AA.SITES  <- rbind(ALL.AA.SITES,AA.FREQ)
    AA.SITES      <- rbind(AA.SITES,    AA.FREQ[which(AA.FREQ$FREQ > 0.05 & AA.FREQ$AA != AA.FREQ$WT),])
    DEL.SITES     <- rbind(DEL.SITES,   DEL.FREQ[which(DEL.FREQ$FREQ > 0.05 & DEL.FREQ$DNA == "-"),])
    INS.SITES     <- rbind(INS.SITES,   INS.FREQ[which(INS.FREQ$FREQ > 0.05 & INS.FREQ$DNA != "-"),])
  }
}

#write.table(ALL.AA.SITES.COUNT, file="ALL.AA.SITES.COUNT.txt",quote=FALSE)
ALL.AA.SITES.COUNT <- read.table("ALL.AA.SITES.COUNT.txt")

DNA.SITES <- rbind(DEL.SITES,INS.SITES)
DNA.SITES <- DNA.SITES[order(DNA.SITES$SITE),] 
DNA.SITES$SITE <- as.factor(DNA.SITES$SITE)
DNA.SITES$EXP  <- as.factor(DNA.SITES$EXP)

AA.SITES <- AA.SITES[order(AA.SITES$SITE),]
AA.SITES$SITE <- as.factor(AA.SITES$SITE)
AA.SITES$EXP  <- as.factor(AA.SITES$EXP)

#AA point sub Frequency
#MCL1 colors
ggplot(AA.SITES, aes(SITE:AA,EXP)) + coord_fixed(ratio=1) + geom_tile(aes(fill=FREQ)) +
  theme(legend.position = "top") +
  theme(axis.text.x = element_text(size=12)) +
  theme(panel.background = element_rect(fill = 'gray', colour = 'gray')) +
  scale_fill_gradientn(colors=c("#FFFFFF","#E066FF","#8D4896","#5D478B"),limits=c(0,1)) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(AA.SITES, aes(SITE,EXP)) + coord_fixed(ratio=1) + geom_tile(aes(fill=FREQ)) +
  theme(legend.position = "top") +
  theme(axis.text.x = element_text(size=8)) +
  theme(panel.background = element_rect(fill = 'gray', colour = 'gray')) +
  scale_fill_gradientn(colors=c("#FFFFFF","#E066FF","#8D4896","#5D478B"),limits=c(0,1)) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank())

#BCL2 colors
ggplot(AA.SITES, aes(SITE:AA,EXP)) + coord_fixed(ratio=1) + geom_tile(aes(fill=FREQ)) +
  theme(legend.position = "top") +
  theme(axis.text.x = element_text(size=12)) +
  theme(panel.background = element_rect(fill = 'gray', colour = 'gray')) +
  scale_fill_gradientn(colors=c("#FFFFFF","#00FF66","#228B22","#006400"),limits=c(0,1)) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(AA.SITES, aes(SITE,EXP)) + coord_fixed(ratio=1) + geom_tile(aes(fill=FREQ)) +
  theme(legend.position = "top") +
  theme(axis.text.x = element_text(size=8)) +
  theme(panel.background = element_rect(fill = 'gray', colour = 'gray')) +
  scale_fill_gradientn(colors=c("#FFFFFF","#00FF66","#228B22","#006400"),limits=c(0,1)) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank())


#Indel Frequency
#MCL1 colors
ggplot(DNA.SITES, aes(SITE:DNA,EXP)) + coord_fixed(ratio=1) + geom_tile(aes(fill=FREQ)) +
  theme(legend.position = "top") +
  theme(axis.text.x = element_text(size=8)) +
  theme(panel.background = element_rect(fill = 'gray', colour = 'gray')) +
  scale_fill_gradientn(colors=c("#FFFFFF","#E066FF","#8D4896","#5D478B"),limits=c(0,1)) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank())

#BCL2 colors
ggplot(DNA.SITES, aes(SITE:DNA,EXP)) + coord_fixed(ratio=1) + geom_tile(aes(fill=FREQ)) +
  theme(legend.position = "top") +
  theme(axis.text.x = element_text(size=8)) +
  theme(panel.background = element_rect(fill = 'gray', colour = 'gray')) +
  scale_fill_gradientn(colors=c("#FFFFFF","#00FF66","#228B22","#006400"),limits=c(0,1)) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank())


###################################################
###AA Frequency comparison for entire experiment###
###################################################
#All sites and reps
ALL.AA.SITES.COUNT <- ALL.AA.SITES.COUNT
AA.COUNT.FULL <- acast(ALL.AA.SITES.COUNT,ALL.AA.SITES.COUNT$EXP ~ ALL.AA.SITES.COUNT$SITE + ALL.AA.SITES.COUNT$AA)
ALL.AA.SITES.COUNT.WT <- ALL.AA.SITES.COUNT[seq(1,nrow(ALL.AA.SITES.COUNT),21),]
AA.WT.FULL <- acast(ALL.AA.SITES.COUNT.WT,ALL.AA.SITES.COUNT.WT$EXP ~ ALL.AA.SITES.COUNT.WT$SITE,value.var="WT")
AA.WT.FULL <- matrix(lapply(AA.WT.FULL, as.character), nrow=nrow(AA.WT.FULL),byrow = FALSE)
colnames(AA.WT.FULL) <- names(table(ALL.AA.SITES.COUNT$SITE))
rownames(AA.WT.FULL) <- names(table(ALL.AA.SITES.COUNT$EXP))

#Define set of sites and replicates used for main analyses
COMPARABLE.SITES <- c(1:6,64:96,148:198,200:269)
COMPARABLE.REPS  <- c(17:20,25:28,29:31,33:36,57:60,61:64)

REDUCED.REPS.SITES.COUNT <- subset(ALL.AA.SITES.COUNT, ALL.AA.SITES.COUNT$EXP %in% COMPARABLE.REPS & ALL.AA.SITES.COUNT$SITE %in% COMPARABLE.SITES)
AA.COUNT.COMPARE <- acast(REDUCED.REPS.SITES.COUNT,REDUCED.REPS.SITES.COUNT$EXP ~ REDUCED.REPS.SITES.COUNT$SITE + REDUCED.REPS.SITES.COUNT$AA)
REDUCED.REPS.SITES.COUNT <- REDUCED.REPS.SITES.COUNT[seq(1,nrow(REDUCED.REPS.SITES.COUNT),21),]
AA.WT.COMPARE <- acast(REDUCED.REPS.SITES.COUNT,REDUCED.REPS.SITES.COUNT$EXP ~ REDUCED.REPS.SITES.COUNT$SITE,value.var="WT")
AA.WT.COMPARE <- matrix(lapply(AA.WT.COMPARE, as.character), nrow=nrow(AA.WT.COMPARE),byrow = FALSE)
colnames(AA.WT.COMPARE) <- COMPARABLE.SITES
rownames(AA.WT.COMPARE) <- COMPARABLE.REPS

###################################################
###Fst like estimators of chance and contingency###
###################################################
#Based on Hivert et al. 2018. Genetics for Fst calculations of pool-seq data
#Levels are 1. Replicates, 2. Genotypes, 3. Total

N <- 10^9
COUNT.FILTER <- 1
G <- list("B"=c(17:20),"D"=c(29:31),"C"=c(25:28),"M"=c(61:64),"J"=c(33:36),"L"=c(57:60))

L <- c(1:6,64:96,148:198,200:269)

#Error bar function
error.bar <- function(x, y, upper, lower=upper, length=0.01,...){
  if(length(x) != length(y) | length(y) !=length(lower) | length(lower) != length(upper))
    stop("vectors must be same length")
  arrows(x,upper, x, lower, angle=90, code=3, length=length, ...)
}

#Fst calculation Function
MS.SITE.CALC <- function(SITEDATA = SITEDATA, SITE= l, TYPE="States") {

  if(!(TYPE %in% c("States","Sites"))) {
    print("Choose either States or Sites")
    break()
  }

  WT.STATE <- unlist(c(unique(AA.WT.COMPARE[,SITE])))
  if(TYPE == "States") {
    if(length(WT.STATE) == 1) {
        SITEDATA <- cbind(SITEDATA[,grep(WT.STATE[1],colnames(SITEDATA)),drop=FALSE],SITEDATA[,-grep(WT.STATE[1],colnames(SITEDATA)),drop=FALSE])
    } else {
      SITEDATA <- cbind(rowSums(SITEDATA[,grep(paste(WT.STATE,collapse="|"),colnames(SITEDATA))]),SITEDATA[,-grep(paste(WT.STATE,collapse="|"),colnames(SITEDATA))])
    }
    colnames(SITEDATA)[1] <- "WT"
  }
  
  if(TYPE == "Sites") {
   I <- 1:nrow(SITEDATA)
   SITEDATA.WT <- matrix(0,nrow=length(I),ncol=2)
   colnames(SITEDATA.WT) <- c("WT","MUT")
   for(i in I) {
     WT.COL <- grep(paste(WT.STATE,collapse="|"),colnames(SITEDATA))
     SITEDATA.WT[i,1] <- sum(SITEDATA[i, WT.COL])
     SITEDATA.WT[i,2] <- sum(SITEDATA[i,-WT.COL])
   }
   SITEDATA <- SITEDATA.WT
  }
  
  C1  <- sum(SITEDATA)
  C1i <- apply(SITEDATA,1,sum)
  C2  <- sum(C1i^2)
  D2  <- sum((C1i + N - 1)/N)
  D2S <- sum((C1i*(C1i + N - 1)/N))/C1
  Nc <- (C1 - C2/C1)/(D2 - D2S)
  
  SITEDATA.FREQ <- SITEDATA/(rowSums(SITEDATA))
  
  MSI <- (1/(C1 - D2))*sum(C1i*SITEDATA.FREQ*(1-SITEDATA.FREQ))
  MSP <- (1/(D2 - D2S))*sum(C1i*t(t(SITEDATA.FREQ) - colSums(SITEDATA)/C1)^2)
  
  if(is.finite(MSI) == FALSE | is.finite(MSP) == FALSE | Nc == 0) {
    MSI <- NA
    MSP <- NA
  }
  return(list("MSI"=MSI,"MSP"=MSP,"Nc"=Nc))
}

F.CALC <- function(DATA = AA.COUNT.COMPARE, SITES = L, GENOTYPES = G) {
  DATA <- DATA[sapply(unlist(GENOTYPES),FUN=function(x) {which(rownames(DATA) == x)}),]

  Frt_States <- 0
  Frg_States <- numeric(length(GENOTYPES))
  Frg_States_avg <- 0
  Fgt_States <- 0
  Frt_Sites  <- 0
  Frg_Sites  <- numeric(length(GENOTYPES))
  Frg_Sites_avg <- 0
  Fgt_Sites  <- 0
  
  Nc_Frt_States  <- numeric(length(SITES))
  Nc_Frg_States  <- list()
  Nc_Fgt_States  <- numeric(length(SITES))
  Nc_Frt_Sites   <- numeric(length(SITES))
  Nc_Frg_Sites   <- list()
  Nc_Fgt_Sites   <- numeric(length(SITES))
  
  MSI_Frt_States  <- numeric(length(SITES))
  MSI_Frg_States  <- list()
  MSI_Fgt_States  <- numeric(length(SITES))
  MSI_Frt_Sites   <- numeric(length(SITES))
  MSI_Frg_Sites   <- list()
  MSI_Fgt_Sites   <- numeric(length(SITES))
  
  MSP_Frt_States  <- numeric(length(SITES))
  MSP_Frg_States  <- list()
  MSP_Fgt_States  <- numeric(length(SITES))
  MSP_Frt_Sites   <- numeric(length(SITES))
  MSP_Frg_Sites   <- list()
  MSP_Fgt_Sites   <- numeric(length(SITES))
  
  GENOTYPEDATA <- matrix(0,nrow=length(GENOTYPES),ncol=ncol(DATA))
  for(g in 1:length(GENOTYPES)) {
    GENOTYPEDATA[g,] <- colSums(DATA[rownames(DATA) %in% GENOTYPES[[g]],])
  }
  colnames(GENOTYPEDATA) <- colnames(DATA)
  rownames(GENOTYPEDATA) <- names(GENOTYPES)
  
  #Frg
  for(g in 1:length(GENOTYPES)) {
    SUBDATA <- DATA[which(unlist(GENOTYPES) %in% GENOTYPES[[g]]),]
    
    Nc_Frg_States[[g]]  <- numeric(length(SITES))
    MSI_Frg_States[[g]] <- numeric(length(SITES))
    MSP_Frg_States[[g]] <- numeric(length(SITES))
    Nc_Frg_Sites[[g]]   <- numeric(length(SITES))
    MSI_Frg_Sites[[g]]  <- numeric(length(SITES))
    MSP_Frg_Sites[[g]]  <- numeric(length(SITES))
    
    for(l in SITES){
      SITEDATA <- SUBDATA[, (21*(which(L==l)-1) + 1):(which(L==l)*21)]
      STATE.CALC <- MS.SITE.CALC(SITEDATA,which(L==l),"States")
      MSI_Frg_States[[g]][which(SITES==l)] <- STATE.CALC[[1]]
      MSP_Frg_States[[g]][which(SITES==l)] <- STATE.CALC[[2]]
      Nc_Frg_States[[g]][which(SITES==l)]  <- STATE.CALC[[3]]
      
      SITE.CALC <- MS.SITE.CALC(SITEDATA,which(L==l),"Sites")
      MSI_Frg_Sites[[g]][which(SITES==l)] <- SITE.CALC[[1]]
      MSP_Frg_Sites[[g]][which(SITES==l)] <- SITE.CALC[[2]]
      Nc_Frg_Sites[[g]][which(SITES==l)]  <- SITE.CALC[[3]]
    }
  }
  
  MSP_Frg_States_C <- matrix(unlist(MSP_Frg_States),nrow=length(SITES))
  MSI_Frg_States_C <- matrix(unlist(MSI_Frg_States),nrow=length(SITES))
  Nc_Frg_States_C  <- matrix(unlist(Nc_Frg_States), nrow=length(SITES))
  MSP_Frg_Sites_C  <- matrix(unlist(MSP_Frg_Sites), nrow=length(SITES))
  MSI_Frg_Sites_C  <- matrix(unlist(MSI_Frg_Sites), nrow=length(SITES))
  Nc_Frg_Sites_C   <- matrix(unlist(Nc_Frg_Sites),  nrow=length(SITES))
  Ng <- sapply(GENOTYPES,length)
  
  Frg_States_MSP <- rowSums(t(t(MSP_Frg_States_C)*Ng))/sum(Ng,na.rm=TRUE)
  Frg_States_MSI <- rowSums(t(t(MSI_Frg_States_C)*Ng))/sum(Ng,na.rm=TRUE)
  Frg_States_Nc  <- rowSums(t(t(Nc_Frg_States_C)*Ng))/sum(Ng,na.rm=TRUE)
  Frg_Sites_MSP  <- rowSums(t(t(MSP_Frg_Sites_C)*Ng))/sum(Ng,na.rm=TRUE)
  Frg_Sites_MSI  <- rowSums(t(t(MSI_Frg_Sites_C)*Ng))/sum(Ng,na.rm=TRUE)
  Frg_Sites_Nc   <- rowSums(t(t(Nc_Frg_Sites_C)*Ng))/sum(Ng,na.rm=TRUE)
  
  Frg_States_Num <- rowSums(t(t(MSP_Frg_States_C - MSI_Frg_States_C)*Ng),na.rm=TRUE)/sum(Ng,na.rm=TRUE)
  Frg_States_Den <- rowSums(t(t(MSP_Frg_States_C + (Nc_Frg_States_C - 1)*MSI_Frg_States_C)*Ng),na.rm=TRUE)/sum(Ng,na.rm=TRUE)
  Frg_Sites_Num  <- rowSums(t(t(MSP_Frg_Sites_C - MSI_Frg_Sites_C)*Ng),na.rm=TRUE)/sum(Ng,na.rm=TRUE)
  Frg_Sites_Den  <- rowSums(t(t(MSP_Frg_Sites_C + (Nc_Frg_Sites_C - 1)*MSI_Frg_Sites_C)*Ng),na.rm=TRUE)/sum(Ng,na.rm=TRUE)
  Frg_States_avg <- Frg_States_Num/Frg_States_Den
  Frg_Sites_avg  <- Frg_Sites_Num/Frg_Sites_Den
  
  #Fgt
  for(l in SITES) {
    SITEDATA <- GENOTYPEDATA[, (21*(which(L==l)-1) + 1):(which(L==l)*21), drop=FALSE]
    STATE.CALC <- MS.SITE.CALC(SITEDATA,which(L==l),"States")
    MSI_Fgt_States[which(SITES==l)] <- STATE.CALC[[1]]
    MSP_Fgt_States[which(SITES==l)] <- STATE.CALC[[2]]
    Nc_Fgt_States[which(SITES==l)]  <- STATE.CALC[[3]]
    
    SITE.CALC <- MS.SITE.CALC(SITEDATA,which(L==l),"Sites")
    MSI_Fgt_Sites[which(SITES==l)] <- SITE.CALC[[1]]
    MSP_Fgt_Sites[which(SITES==l)] <- SITE.CALC[[2]]
    Nc_Fgt_Sites[which(SITES==l)]  <- SITE.CALC[[3]]
  }
  Fgt_States_Num <- (MSP_Fgt_States - MSI_Fgt_States)
  Fgt_Sites_Num  <- (MSP_Fgt_Sites -  MSI_Fgt_Sites)
  
  Fgt_States_Den <- (MSP_Fgt_States + (Nc_Fgt_States-1)*MSI_Fgt_States)
  Fgt_Sites_Den  <- (MSP_Fgt_Sites  + (Nc_Fgt_Sites-1)* MSI_Fgt_Sites)
  
  Fgt_States <- Fgt_States_Num/Fgt_States_Den
  Fgt_Sites  <- Fgt_Sites_Num/Fgt_Sites_Den
  
  #Frt
  for(l in SITES) {
    SITEDATA <- DATA[, (21*(which(L==l)-1) + 1):(which(L==l)*21), drop=FALSE]

    STATE.CALC <- MS.SITE.CALC(SITEDATA,which(L==l),"States")
    MSI_Frt_States[which(SITES==l)] <- STATE.CALC[[1]]
    MSP_Frt_States[which(SITES==l)] <- STATE.CALC[[2]]
    Nc_Frt_States[which(SITES==l)]  <- STATE.CALC[[3]]
    
    SITE.CALC <- MS.SITE.CALC(SITEDATA,which(L==l),"Sites")
    MSI_Frt_Sites[which(SITES==l)] <- SITE.CALC[[1]]
    MSP_Frt_Sites[which(SITES==l)] <- SITE.CALC[[2]]
    Nc_Frt_Sites[which(SITES==l)]  <- SITE.CALC[[3]]
  }
  Frt_States_Num <- (MSP_Frt_States - MSI_Frt_States)
  Frt_Sites_Num  <- (MSP_Frt_Sites -  MSI_Frt_Sites)
  
  Frt_States_Den <- (MSP_Frt_States + (Nc_Frt_States-1)*MSI_Frt_States)
  Frt_Sites_Den  <- (MSP_Frt_Sites  + (Nc_Frt_Sites-1)* MSI_Frt_Sites)
  
  Frt_States <- Frt_States_Num/Frt_States_Den
  Frt_Sites  <- Frt_Sites_Num/Frt_Sites_Den
  
  F <- rbind(Frg_States_MSP,Frg_States_MSI,Frg_States_Nc,
             Frg_Sites_MSP, Frg_Sites_MSI, Frg_Sites_Nc,
             MSP_Fgt_States,MSI_Fgt_States,Nc_Fgt_States,
             MSP_Fgt_Sites, MSI_Fgt_Sites, Nc_Fgt_Sites,
             MSP_Frt_States,MSI_Frt_States,Nc_Frt_States,
             MSP_Frt_Sites, MSI_Frt_Sites, Nc_Frt_Sites)

  return(F)
}
##Global D statistics
#Get MSP, MSI, and Nc for both states and sites
DST <- F.CALC(DATA = AA.COUNT.COMPARE, SITES = L, GENOTYPES = G)

#Estimate D statistics
D.STAT <- numeric(6)
D.STAT[1] <- 1/(1 - sum(DST[1,]  - DST[2,] )/sum(DST[1,]  + (DST[3,]-1 )*DST[2,] ))
D.STAT[2] <- 1/(1 - sum(DST[4,]  - DST[5,] )/sum(DST[4,]  + (DST[6,]-1 )*DST[5,] ))
D.STAT[3] <- sum(2*(DST[7,] - DST[8,]) )/sum(DST[9, ]*DST[8,] ) + 1
D.STAT[4] <- sum(2*(DST[10,]- DST[11,]))/sum(DST[12,]*DST[11,]) + 1
D.STAT[5] <- D.STAT[1]*D.STAT[3]
D.STAT[6] <- D.STAT[2]*D.STAT[4]

#Set up matrices for pairwise comparison of D statistics
CHANCE.STATES.MATRIX <- matrix(0,nrow=6,ncol=6,byrow=TRUE)
CHANCE.SITES.MATRIX  <- matrix(0,nrow=6,ncol=6,byrow=TRUE)
CONTIN.STATES.MATRIX <- matrix(0,nrow=6,ncol=6,byrow=TRUE)
CONTIN.SITES.MATRIX  <- matrix(0,nrow=6,ncol=6,byrow=TRUE)
TOTALD.STATES.MATRIX <- matrix(0,nrow=6,ncol=6,byrow=TRUE)
TOTALD.SITES.MATRIX  <- matrix(0,nrow=6,ncol=6,byrow=TRUE)

#Pairwise D statistic calculations
I <- 1:(nrow(CHANCE.STATES.MATRIX) -1)
for(i in I) {
  J <- (i+1):ncol(CHANCE.STATES.MATRIX)
  for(j in J){
    R.DATA <- F.CALC(DATA = AA.COUNT.COMPARE, SITES = L, GENOTYPES = G[c(i,j)])
    
    CHANCE.STATES.MATRIX[i,j]  <- 1/(1- sum(R.DATA[1,]  - R.DATA[2,] )/sum(R.DATA[1,]  + (R.DATA[3,]-1 )*R.DATA[2,] ))
    CHANCE.SITES.MATRIX[i,j]   <- 1/(1- sum(R.DATA[4,]  - R.DATA[5,] )/sum(R.DATA[4,]  + (R.DATA[6,]-1 )*R.DATA[5,] ))
    CONTIN.STATES.MATRIX[i,j]  <- sum(2*(R.DATA[7,] - R.DATA[8,]) )/sum(R.DATA[9, ]*R.DATA[8,] ) + 1
    CONTIN.SITES.MATRIX[i,j]   <- sum(2*(R.DATA[10,]- R.DATA[11,]))/sum(R.DATA[12,]*R.DATA[11,]) + 1
    TOTALD.STATES.MATRIX[i,j]  <- CHANCE.STATES.MATRIX[i,j]*CONTIN.STATES.MATRIX[i,j]
    TOTALD.SITES.MATRIX[i,j]   <- CHANCE.SITES.MATRIX[i,j ]*CONTIN.SITES.MATRIX[i,j]
  }
}

#Collate data based on statistic instead of comparison
D.DATA.PHYLO <- matrix(0,nrow=6,ncol=5)
D.DATA.PHYLO[1,] <- CHANCE.STATES.MATRIX[row(CHANCE.STATES.MATRIX) == col(CHANCE.STATES.MATRIX) - 1]
D.DATA.PHYLO[2,] <- CHANCE.SITES.MATRIX[ row(CHANCE.SITES.MATRIX)  == col(CHANCE.SITES.MATRIX)  - 1]
D.DATA.PHYLO[3,] <- CONTIN.STATES.MATRIX[row(CONTIN.STATES.MATRIX) == col(CONTIN.STATES.MATRIX) - 1]
D.DATA.PHYLO[4,] <- CONTIN.SITES.MATRIX[ row(CONTIN.SITES.MATRIX)  == col(CONTIN.SITES.MATRIX)  - 1]
D.DATA.PHYLO[5,] <- TOTALD.STATES.MATRIX[row(TOTALD.STATES.MATRIX) == col(TOTALD.STATES.MATRIX) - 1]
D.DATA.PHYLO[6,] <- TOTALD.SITES.MATRIX[ row(TOTALD.SITES.MATRIX)  == col(TOTALD.SITES.MATRIX)  - 1]

###Relationship of D statistics with phylogenetic distance
#Get phylogenetic distances between ancestors based on branch lengths
PHYLO.DIST <- as.matrix(read.table("Phylo.Dist.txt",header=TRUE))

#Remove AncMB1 as it only has two replicates
PHYLO.DIST <- PHYLO.DIST[-5,-5]

#Phylogenetically independent
par(mfrow=c(1,1))
plot(  PHYLO.DIST[row(PHYLO.DIST) == col(PHYLO.DIST) - 1],D.DATA.PHYLO[1,],pch=19,cex=0.8,xlim=c(0,2.5),ylim=c(1,3),col="#FF0000",xlab="Phylogenetic Distance",ylab="R")
points(PHYLO.DIST[row(PHYLO.DIST) == col(PHYLO.DIST) - 1],D.DATA.PHYLO[3,],pch=19,cex=0.8,col="#0000FF")
points(PHYLO.DIST[row(PHYLO.DIST) == col(PHYLO.DIST) - 1],D.DATA.PHYLO[5,],pch=19,cex=0.8,col="#00FF00")
abline(lm(D.DATA.PHYLO[1,] ~ PHYLO.DIST[row(PHYLO.DIST) == col(PHYLO.DIST) - 1]),col="#FF0000")
abline(lm(D.DATA.PHYLO[3,] ~ PHYLO.DIST[row(PHYLO.DIST) == col(PHYLO.DIST) - 1]),col="#0000FF")
abline(lm(D.DATA.PHYLO[5,] ~ PHYLO.DIST[row(PHYLO.DIST) == col(PHYLO.DIST) - 1]),col="#00FF00")
anova( lm(D.DATA.PHYLO[1,] ~ PHYLO.DIST[row(PHYLO.DIST) == col(PHYLO.DIST) - 1]))
anova( lm(D.DATA.PHYLO[3,] ~ PHYLO.DIST[row(PHYLO.DIST) == col(PHYLO.DIST) - 1]))
anova( lm(D.DATA.PHYLO[5,] ~ PHYLO.DIST[row(PHYLO.DIST) == col(PHYLO.DIST) - 1]))

#Phylogenetically non-independent
plot(  PHYLO.DIST[row(PHYLO.DIST) <= col(PHYLO.DIST) - 1],CHANCE.STATES.MATRIX[row(CHANCE.STATES.MATRIX) <= col(CHANCE.STATES.MATRIX) - 1],pch=19,cex=0.8,xlim=c(0,7),ylim=c(1,4),col=c("#FF0000"),xlab="Phylogenetic Distance",ylab="R")
points(PHYLO.DIST[row(PHYLO.DIST) <= col(PHYLO.DIST) - 1],CONTIN.STATES.MATRIX[row(CONTIN.STATES.MATRIX) <= col(CONTIN.STATES.MATRIX) - 1],pch=19,cex=0.8,col="#0000FF")
points(PHYLO.DIST[row(PHYLO.DIST) <= col(PHYLO.DIST) - 1],TOTALD.STATES.MATRIX[row(TOTALD.STATES.MATRIX) <= col(TOTALD.STATES.MATRIX) - 1],pch=19,cex=0.8,col="#00FF00")
abline(lm(CHANCE.STATES.MATRIX[row(CHANCE.STATES.MATRIX) <= col(CHANCE.STATES.MATRIX) - 1] ~  PHYLO.DIST[row(PHYLO.DIST) <= col(PHYLO.DIST) - 1]),col="#FF0000")
abline(lm(CONTIN.STATES.MATRIX[row(CONTIN.STATES.MATRIX) <= col(CONTIN.STATES.MATRIX) - 1] ~  PHYLO.DIST[row(PHYLO.DIST) <= col(PHYLO.DIST) - 1]),col="#0000FF")
abline(lm(TOTALD.STATES.MATRIX[row(TOTALD.STATES.MATRIX) <= col(TOTALD.STATES.MATRIX) - 1] ~  PHYLO.DIST[row(PHYLO.DIST) <= col(PHYLO.DIST) - 1]),col="#00FF00")
anova( lm(CHANCE.STATES.MATRIX[row(CHANCE.STATES.MATRIX) <= col(CHANCE.STATES.MATRIX) - 1]  ~ PHYLO.DIST[row(PHYLO.DIST) <= col(PHYLO.DIST) - 1]))
anova( lm(CONTIN.STATES.MATRIX[row(CONTIN.STATES.MATRIX) <= col(CONTIN.STATES.MATRIX) - 1]  ~ PHYLO.DIST[row(PHYLO.DIST) <= col(PHYLO.DIST) - 1]))
anova( lm(TOTALD.STATES.MATRIX[row(TOTALD.STATES.MATRIX) <= col(TOTALD.STATES.MATRIX) - 1]  ~ PHYLO.DIST[row(PHYLO.DIST) <= col(PHYLO.DIST) - 1]))
legend(5,4,c("Chance","Contingency","Total"),fill=c("red","blue","green"),cex=0.4,pt.cex=0.4)

points(PHYLO.DIST[row(PHYLO.DIST) == col(PHYLO.DIST) - 1],D.DATA.PHYLO[1,],pch=21,cex=1.4,col="#FF0000")
points(PHYLO.DIST[row(PHYLO.DIST) == col(PHYLO.DIST) - 1],D.DATA.PHYLO[3,],pch=21,cex=1.4,col="#0000FF")
points(PHYLO.DIST[row(PHYLO.DIST) == col(PHYLO.DIST) - 1],D.DATA.PHYLO[5,],pch=21,cex=1.4,col="#00FF00")
abline(lm(D.DATA.PHYLO[1,] ~ PHYLO.DIST[row(PHYLO.DIST) == col(PHYLO.DIST) - 1]),col="#FF0000",lty=2)
abline(lm(D.DATA.PHYLO[3,] ~ PHYLO.DIST[row(PHYLO.DIST) == col(PHYLO.DIST) - 1]),col="#0000FF",lty=2)
abline(lm(D.DATA.PHYLO[5,] ~ PHYLO.DIST[row(PHYLO.DIST) == col(PHYLO.DIST) - 1]),col="#00FF00",lty=2)
anova( lm(D.DATA.PHYLO[1,] ~ PHYLO.DIST[row(PHYLO.DIST) == col(PHYLO.DIST) - 1]))
anova( lm(D.DATA.PHYLO[3,] ~ PHYLO.DIST[row(PHYLO.DIST) == col(PHYLO.DIST) - 1]))
anova( lm(D.DATA.PHYLO[5,] ~ PHYLO.DIST[row(PHYLO.DIST) == col(PHYLO.DIST) - 1]))

##Site Frequency comparison
#Identify sites with a high frequency non-wt allele
COUNT.DATA <- acast(ALL.AA.SITES.COUNT, ALL.AA.SITES.COUNT$EXP ~ ALL.AA.SITES.COUNT$SITE ~ ALL.AA.SITES.COUNT$AA)
FREQ.DATA  <- apply(COUNT.DATA,c(1,2),function(x) { x/sum(x)})
FREQ.DATA  <- aperm(FREQ.DATA,c(2,3,1),resize = TRUE)

NO.WT.FREQ.DATA <- FREQ.DATA
for(l in L) {
  WT.STATE <- AA.WT.FULL[,l]
  I <- 1:dim(NO.WT.FREQ.DATA)[1]
  for(i in I) {
    NO.WT.FREQ.DATA[i,l,grep(WT.STATE[i],names(NO.WT.FREQ.DATA[i,l,]))] <- 0
  }
}
NO.WT.FREQ.DATA <- NO.WT.FREQ.DATA[dimnames(NO.WT.FREQ.DATA)[[1]] %in% COMPARABLE.REPS,,]

FREQS.NOZERO <- NO.WT.FREQ.DATA[which(NO.WT.FREQ.DATA != 0,arr.ind = TRUE)]
hist(log10(FREQS.NOZERO),breaks=100,xlab="Log10(Non-WT Frequency)",main="")
abline(v=mean(log10(FREQS.NOZERO)),col="black",lty=2)
abline(v=mean(log10(FREQS.NOZERO)) +   sd(log10(FREQS.NOZERO)),col="black")
abline(v=mean(log10(FREQS.NOZERO)) + 2*sd(log10(FREQS.NOZERO)),col="blue")
abline(v=mean(log10(FREQS.NOZERO)) + 3*sd(log10(FREQS.NOZERO)),col="green")
abline(v=mean(log10(FREQS.NOZERO)) + 4*sd(log10(FREQS.NOZERO)),col="orange")
abline(v=mean(log10(FREQS.NOZERO)) + 5*sd(log10(FREQS.NOZERO)),col="red")

FREQ.01 <- apply(NO.WT.FREQ.DATA,c(1,2),function(x) { sum(x > 0.01)})
FREQ.05 <- apply(NO.WT.FREQ.DATA,c(1,2),function(x) { sum(x > 0.05)})
FREQ.10 <- apply(NO.WT.FREQ.DATA,c(1,2),function(x) { sum(x > 0.10)})


#Repeated mutation counts
FREQS.MCL1 <- NO.WT.FREQ.DATA[dimnames(NO.WT.FREQ.DATA)[[1]] %in% c(57:60,33:36,61:64),,]
FREQS.BCL2 <- NO.WT.FREQ.DATA[dimnames(NO.WT.FREQ.DATA)[[1]] %in% c(17:20,29:31,25:28),,]

REPEAT.COUNT <- matrix(0,nrow=4,ncol=4)
rownames(REPEAT.COUNT) <- c("SR.SG","SR.MG","MR.SG","MR.MG")
colnames(REPEAT.COUNT) <- c("SITE.MCL1","SITE.BCL2","STATE.MCL1","STATE.BCL2")
 
#MCL1
SITE.REPLICATE.FREQ.COUNT.MCL1 <- apply(FREQS.MCL1,c(1,2),function(x) { sum(x > 0.05)})
SITE.GENOTYPES.FREQ.COUNT.MCL1 <- matrix(0,nrow=length(GENOTYPES),ncol=dim(FREQS.MCL1)[2])
for(g in 1:length(GENOTYPES)) {
  SITE.GENOTYPES.FREQ.COUNT.MCL1[g,] <- apply(FREQS.MCL1[dimnames(FREQS.MCL1)[[1]] %in% unlist(GENOTYPES[[g]]),,],2,function(x) { sum(apply(x,1,function(y) { sum(y > 0.05,na.rm=TRUE) }) > 0, na.rm=TRUE)})
}  
REPEAT.COUNT[1,1] <- length(which(colSums(SITE.GENOTYPES.FREQ.COUNT.MCL1 > 0, na.rm=TRUE) == 1 & colSums(SITE.GENOTYPES.FREQ.COUNT.MCL1 > 1, na.rm=TRUE) == 0))
REPEAT.COUNT[2,1] <- length(which(colSums(SITE.GENOTYPES.FREQ.COUNT.MCL1 > 0, na.rm=TRUE) >= 2 & colSums(SITE.GENOTYPES.FREQ.COUNT.MCL1 > 1, na.rm=TRUE) == 0))
REPEAT.COUNT[3,1] <- length(which(colSums(SITE.GENOTYPES.FREQ.COUNT.MCL1 > 0, na.rm=TRUE) == 1 & colSums(SITE.GENOTYPES.FREQ.COUNT.MCL1 > 1, na.rm=TRUE) >= 1))
REPEAT.COUNT[4,1] <- length(which(colSums(SITE.GENOTYPES.FREQ.COUNT.MCL1 > 0, na.rm=TRUE) >= 2 & colSums(SITE.GENOTYPES.FREQ.COUNT.MCL1 > 1, na.rm=TRUE) >= 1))
  
STATE.REPLICATE.FREQ.COUNT.MCL1 <- apply(FREQS.MCL1,c(1,2,3),function(x) { x > 0.05})
STATE.GENOTYPES.FREQ.COUNT.MCL1 <- array(0,c(length(GENOTYPES),dim(FREQS.MCL1)[2],dim(FREQS.MCL1)[3]))
for(g in 4:6) {
  STATE.GENOTYPES.FREQ.COUNT.MCL1[g,,] <- apply(FREQS.MCL1[dimnames(FREQS.MCL1)[[1]] %in% unlist(GENOTYPES[[g]]),,],2,function(x) { colSums(apply(x,2,function(y) { y > 0.05 }) > 0, na.rm=TRUE)})
}  
REPEAT.COUNT[1,3] <- length(which(apply(STATE.GENOTYPES.FREQ.COUNT.MCL1,c(2,3), function(x) { sum(x > 0,na.rm=TRUE) }) == 1 & apply(STATE.GENOTYPES.FREQ.COUNT.MCL1,c(2,3), function(x) { sum(x > 1,na.rm=TRUE) }) == 0))
REPEAT.COUNT[2,3] <- length(which(apply(STATE.GENOTYPES.FREQ.COUNT.MCL1,c(2,3), function(x) { sum(x > 0,na.rm=TRUE) }) >= 2 & apply(STATE.GENOTYPES.FREQ.COUNT.MCL1,c(2,3), function(x) { sum(x > 1,na.rm=TRUE) }) == 0))
REPEAT.COUNT[3,3] <- length(which(apply(STATE.GENOTYPES.FREQ.COUNT.MCL1,c(2,3), function(x) { sum(x > 0,na.rm=TRUE) }) == 1 & apply(STATE.GENOTYPES.FREQ.COUNT.MCL1,c(2,3), function(x) { sum(x > 1,na.rm=TRUE) }) >= 1))
REPEAT.COUNT[4,3] <- length(which(apply(STATE.GENOTYPES.FREQ.COUNT.MCL1,c(2,3), function(x) { sum(x > 0,na.rm=TRUE) }) >= 2 & apply(STATE.GENOTYPES.FREQ.COUNT.MCL1,c(2,3), function(x) { sum(x > 1,na.rm=TRUE) }) >= 1))

#BCL2
SITE.REPLICATE.FREQ.COUNT.BCL2 <- apply(FREQS.BCL2,c(1,2),function(x) { sum(x > 0.05)})
SITE.GENOTYPES.FREQ.COUNT.BCL2 <- matrix(0,nrow=length(GENOTYPES),ncol=dim(FREQS.BCL2)[2])
for(g in 1:length(GENOTYPES)) {
  SITE.GENOTYPES.FREQ.COUNT.BCL2[g,] <- apply(FREQS.BCL2[dimnames(FREQS.BCL2)[[1]] %in% unlist(GENOTYPES[[g]]),,],2,function(x) { sum(apply(x,1,function(y) { sum(y > 0.05,na.rm=TRUE) }) > 0, na.rm=TRUE)})
}  
REPEAT.COUNT[1,2] <- length(which(colSums(SITE.GENOTYPES.FREQ.COUNT.BCL2 > 0, na.rm=TRUE) == 1 & colSums(SITE.GENOTYPES.FREQ.COUNT.BCL2 > 1, na.rm=TRUE) == 0))
REPEAT.COUNT[2,2] <- length(which(colSums(SITE.GENOTYPES.FREQ.COUNT.BCL2 > 0, na.rm=TRUE) >= 2 & colSums(SITE.GENOTYPES.FREQ.COUNT.BCL2 > 1, na.rm=TRUE) == 0))
REPEAT.COUNT[3,2] <- length(which(colSums(SITE.GENOTYPES.FREQ.COUNT.BCL2 > 0, na.rm=TRUE) == 1 & colSums(SITE.GENOTYPES.FREQ.COUNT.BCL2 > 1, na.rm=TRUE) >= 1))
REPEAT.COUNT[4,2] <- length(which(colSums(SITE.GENOTYPES.FREQ.COUNT.BCL2 > 0, na.rm=TRUE) >= 2 & colSums(SITE.GENOTYPES.FREQ.COUNT.BCL2 > 1, na.rm=TRUE) >= 1))

STATE.REPLICATE.FREQ.COUNT.BCL2 <- apply(FREQS.BCL2,c(1,2,3),function(x) { x > 0.05})
STATE.GENOTYPES.FREQ.COUNT.BCL2 <- array(0,c(length(GENOTYPES),dim(FREQS.BCL2)[2],dim(FREQS.BCL2)[3]))
for(g in 1:3) {
  STATE.GENOTYPES.FREQ.COUNT.BCL2[g,,] <- apply(FREQS.BCL2[dimnames(FREQS.BCL2)[[1]] %in% unlist(GENOTYPES[[g]]),,],2,function(x) { colSums(apply(x,2,function(y) { y > 0.05 }) > 0, na.rm=TRUE)})
}  
REPEAT.COUNT[1,4] <- length(which(apply(STATE.GENOTYPES.FREQ.COUNT.BCL2,c(2,3), function(x) { sum(x > 0,na.rm=TRUE) }) == 1 & apply(STATE.GENOTYPES.FREQ.COUNT.BCL2,c(2,3), function(x) { sum(x > 1,na.rm=TRUE) }) == 0))
REPEAT.COUNT[2,4] <- length(which(apply(STATE.GENOTYPES.FREQ.COUNT.BCL2,c(2,3), function(x) { sum(x > 0,na.rm=TRUE) }) >= 2 & apply(STATE.GENOTYPES.FREQ.COUNT.BCL2,c(2,3), function(x) { sum(x > 1,na.rm=TRUE) }) == 0))
REPEAT.COUNT[3,4] <- length(which(apply(STATE.GENOTYPES.FREQ.COUNT.BCL2,c(2,3), function(x) { sum(x > 0,na.rm=TRUE) }) == 1 & apply(STATE.GENOTYPES.FREQ.COUNT.BCL2,c(2,3), function(x) { sum(x > 1,na.rm=TRUE) }) >= 1))
REPEAT.COUNT[4,4] <- length(which(apply(STATE.GENOTYPES.FREQ.COUNT.BCL2,c(2,3), function(x) { sum(x > 0,na.rm=TRUE) }) >= 2 & apply(STATE.GENOTYPES.FREQ.COUNT.BCL2,c(2,3), function(x) { sum(x > 1,na.rm=TRUE) }) >= 1))


#Permutation testing for overlap in use of sites
T <- 1:100000
PERMUTE.REPEAT.COUNT <- array(0,c(nrow(REPEAT.COUNT),ncol(REPEAT.COUNT),length(T)))
rownames(PERMUTE.REPEAT.COUNT) <- c("SR.SG","MR.SG","SR.MG","MR.MG")

for(t in T) {
  #MCL1
  PERMUTE.MCL1 <- apply(FREQS.MCL1,c(1,3),function(x) {sample(x,length(x),replace=FALSE)})
  PERMUTE.MCL1 <- aperm(PERMUTE.MCL1,c(2,1,3),resize = TRUE)
  dimnames(PERMUTE.MCL1)[[3]] <- dimnames(FREQS.MCL1)[[3]] 
  
  SITE.GENOTYPES.FREQ.COUNT.MCL1 <- matrix(0,nrow=length(GENOTYPES),ncol=dim(PERMUTE.MCL1)[2])
  for(g in 1:length(GENOTYPES)) {
    SITE.GENOTYPES.FREQ.COUNT.MCL1[g,] <- apply(PERMUTE.MCL1[dimnames(PERMUTE.MCL1)[[1]] %in% unlist(GENOTYPES[[g]]),,],2,function(x) { sum(apply(x,1,function(y) { sum(y > 0.05,na.rm=TRUE) }) > 0, na.rm=TRUE)})
  }  
  PERMUTE.REPEAT.COUNT[1,1,t] <- length(which(colSums(SITE.GENOTYPES.FREQ.COUNT.MCL1 > 0, na.rm=TRUE) == 1 & colSums(SITE.GENOTYPES.FREQ.COUNT.MCL1 > 1, na.rm=TRUE) == 0))
  PERMUTE.REPEAT.COUNT[2,1,t] <- length(which(colSums(SITE.GENOTYPES.FREQ.COUNT.MCL1 > 0, na.rm=TRUE) >= 2 & colSums(SITE.GENOTYPES.FREQ.COUNT.MCL1 > 1, na.rm=TRUE) == 0))
  PERMUTE.REPEAT.COUNT[3,1,t] <- length(which(colSums(SITE.GENOTYPES.FREQ.COUNT.MCL1 > 0, na.rm=TRUE) == 1 & colSums(SITE.GENOTYPES.FREQ.COUNT.MCL1 > 1, na.rm=TRUE) >= 1))
  PERMUTE.REPEAT.COUNT[4,1,t] <- length(which(colSums(SITE.GENOTYPES.FREQ.COUNT.MCL1 > 0, na.rm=TRUE) >= 2 & colSums(SITE.GENOTYPES.FREQ.COUNT.MCL1 > 1, na.rm=TRUE) >= 1))
  
  STATE.GENOTYPES.FREQ.COUNT.MCL1 <- array(0,c(length(GENOTYPES),dim(PERMUTE.MCL1)[2],dim(PERMUTE.MCL1)[3]))
  for(g in 4:6) {
    STATE.GENOTYPES.FREQ.COUNT.MCL1[g,,] <- apply(PERMUTE.MCL1[dimnames(PERMUTE.MCL1)[[1]] %in% unlist(GENOTYPES[[g]]),,],2,function(x) { colSums(apply(x,2,function(y) { y > 0.05 }) > 0, na.rm=TRUE)})
  }  
  PERMUTE.REPEAT.COUNT[1,3,t] <- length(which(apply(STATE.GENOTYPES.FREQ.COUNT.MCL1,c(2,3), function(x) { sum(x > 0,na.rm=TRUE) }) == 1 & apply(STATE.GENOTYPES.FREQ.COUNT.MCL1,c(2,3), function(x) { sum(x > 1,na.rm=TRUE) }) == 0))
  PERMUTE.REPEAT.COUNT[2,3,t] <- length(which(apply(STATE.GENOTYPES.FREQ.COUNT.MCL1,c(2,3), function(x) { sum(x > 0,na.rm=TRUE) }) >= 2 & apply(STATE.GENOTYPES.FREQ.COUNT.MCL1,c(2,3), function(x) { sum(x > 1,na.rm=TRUE) }) == 0))
  PERMUTE.REPEAT.COUNT[3,3,t] <- length(which(apply(STATE.GENOTYPES.FREQ.COUNT.MCL1,c(2,3), function(x) { sum(x > 0,na.rm=TRUE) }) == 1 & apply(STATE.GENOTYPES.FREQ.COUNT.MCL1,c(2,3), function(x) { sum(x > 1,na.rm=TRUE) }) >= 1))
  PERMUTE.REPEAT.COUNT[4,3,t] <- length(which(apply(STATE.GENOTYPES.FREQ.COUNT.MCL1,c(2,3), function(x) { sum(x > 0,na.rm=TRUE) }) >= 2 & apply(STATE.GENOTYPES.FREQ.COUNT.MCL1,c(2,3), function(x) { sum(x > 1,na.rm=TRUE) }) >= 1))
  
  #BCL2
  PERMUTE.BCL2 <- apply(FREQS.BCL2,c(1,3),function(x) {sample(x,length(x),replace=FALSE)})
  PERMUTE.BCL2 <- aperm(PERMUTE.BCL2,c(2,1,3),resize = TRUE)
  dimnames(PERMUTE.BCL2)[[3]] <- dimnames(FREQS.BCL2)[[3]] 
  
  SITE.GENOTYPES.FREQ.COUNT.BCL2 <- matrix(0,nrow=length(GENOTYPES),ncol=dim(PERMUTE.BCL2)[2])
  for(g in 1:length(GENOTYPES)) {
    SITE.GENOTYPES.FREQ.COUNT.BCL2[g,] <- apply(PERMUTE.BCL2[dimnames(PERMUTE.BCL2)[[1]] %in% unlist(GENOTYPES[[g]]),,],2,function(x) { sum(apply(x,1,function(y) { sum(y > 0.05,na.rm=TRUE) }) > 0, na.rm=TRUE)})
  }  
  PERMUTE.REPEAT.COUNT[1,2,t] <- length(which(colSums(SITE.GENOTYPES.FREQ.COUNT.BCL2 > 0, na.rm=TRUE) == 1 & colSums(SITE.GENOTYPES.FREQ.COUNT.BCL2 > 1, na.rm=TRUE) == 0))
  PERMUTE.REPEAT.COUNT[2,2,t] <- length(which(colSums(SITE.GENOTYPES.FREQ.COUNT.BCL2 > 0, na.rm=TRUE) >= 2 & colSums(SITE.GENOTYPES.FREQ.COUNT.BCL2 > 1, na.rm=TRUE) == 0))
  PERMUTE.REPEAT.COUNT[3,2,t] <- length(which(colSums(SITE.GENOTYPES.FREQ.COUNT.BCL2 > 0, na.rm=TRUE) == 1 & colSums(SITE.GENOTYPES.FREQ.COUNT.BCL2 > 1, na.rm=TRUE) >= 1))
  PERMUTE.REPEAT.COUNT[4,2,t] <- length(which(colSums(SITE.GENOTYPES.FREQ.COUNT.BCL2 > 0, na.rm=TRUE) >= 2 & colSums(SITE.GENOTYPES.FREQ.COUNT.BCL2 > 1, na.rm=TRUE) >= 1))
  
  STATE.GENOTYPES.FREQ.COUNT.BCL2 <- array(0,c(length(GENOTYPES),dim(PERMUTE.BCL2)[2],dim(PERMUTE.BCL2)[3]))
  for(g in 1:3) {
    STATE.GENOTYPES.FREQ.COUNT.BCL2[g,,] <- apply(PERMUTE.BCL2[dimnames(PERMUTE.BCL2)[[1]] %in% unlist(GENOTYPES[[g]]),,],2,function(x) { colSums(apply(x,2,function(y) { y > 0.05 }) > 0, na.rm=TRUE)})
  }  
  PERMUTE.REPEAT.COUNT[1,4,t] <- length(which(apply(STATE.GENOTYPES.FREQ.COUNT.BCL2,c(2,3), function(x) { sum(x > 0,na.rm=TRUE) }) == 1 & apply(STATE.GENOTYPES.FREQ.COUNT.BCL2,c(2,3), function(x) { sum(x > 1,na.rm=TRUE) }) == 0))
  PERMUTE.REPEAT.COUNT[2,4,t] <- length(which(apply(STATE.GENOTYPES.FREQ.COUNT.BCL2,c(2,3), function(x) { sum(x > 0,na.rm=TRUE) }) >= 2 & apply(STATE.GENOTYPES.FREQ.COUNT.BCL2,c(2,3), function(x) { sum(x > 1,na.rm=TRUE) }) == 0))
  PERMUTE.REPEAT.COUNT[3,4,t] <- length(which(apply(STATE.GENOTYPES.FREQ.COUNT.BCL2,c(2,3), function(x) { sum(x > 0,na.rm=TRUE) }) == 1 & apply(STATE.GENOTYPES.FREQ.COUNT.BCL2,c(2,3), function(x) { sum(x > 1,na.rm=TRUE) }) >= 1))
  PERMUTE.REPEAT.COUNT[4,4,t] <- length(which(apply(STATE.GENOTYPES.FREQ.COUNT.BCL2,c(2,3), function(x) { sum(x > 0,na.rm=TRUE) }) >= 2 & apply(STATE.GENOTYPES.FREQ.COUNT.BCL2,c(2,3), function(x) { sum(x > 1,na.rm=TRUE) }) >= 1))
}

#save(PERMUTE.MCL1, file="PERMUTE.MCL1.RData")
load("PERMUTE.MCL1.RData")

#save(PERMUTE.BCL2, file="PERMUTE.BCL2.RData")
load("PERMUTE.BCL2.RData")

#save(PERMUTE.REPEAT.COUNT, file="PERMUTE.REPEAT.COUNT.RData")
load("PERMUTE.REPEAT.COUNT.RData")



#P-values for all categories seperatley 
par(mfrow=c(4,4))
par(mar=c(2,2,2,2))
hist(PERMUTE.REPEAT.COUNT[1,1,],xlim=c(0,100),breaks=seq(0,100,2),xlab="PACE mutations in category",main="Sites-SR-SG-MCL1")
abline(v=REPEAT.COUNT[1,1],col="red")
hist(PERMUTE.REPEAT.COUNT[2,1,],xlim=c(0,100),breaks=seq(0,100,2),xlab="PACE mutations in category",main="Sites-SR-MG-MCL1")
abline(v=REPEAT.COUNT[2,1],col="red")
hist(PERMUTE.REPEAT.COUNT[3,1,],xlim=c(0,100),breaks=seq(0,100,2),xlab="PACE mutations in category",main="Sites-MR-SG-MCL1")
abline(v=REPEAT.COUNT[3,1],col="red")
hist(PERMUTE.REPEAT.COUNT[4,1,],xlim=c(0,100),breaks=seq(0,100,2),xlab="PACE mutations in category",main="Sites-MR-MG-MCL1")
abline(v=REPEAT.COUNT[4,1],col="red")

hist(PERMUTE.REPEAT.COUNT[1,2,],xlim=c(0,100),breaks=seq(0,100,2),xlab="PACE mutations in category",main="Sites-SR-SG-BCL2")
abline(v=REPEAT.COUNT[1,2],col="red")
hist(PERMUTE.REPEAT.COUNT[2,2,],xlim=c(0,100),breaks=seq(0,100,2),xlab="PACE mutations in category",main="Sites-SR-MG-BCL2")
abline(v=REPEAT.COUNT[2,2],col="red")
hist(PERMUTE.REPEAT.COUNT[3,2,],xlim=c(0,100),breaks=seq(0,100,2),xlab="PACE mutations in category",main="Sites-MR-SG-BCL2")
abline(v=REPEAT.COUNT[3,2],col="red")
hist(PERMUTE.REPEAT.COUNT[4,2,],xlim=c(0,100),breaks=seq(0,100,2),xlab="PACE mutations in category",main="Sites-MR-MG-BCL2")
abline(v=REPEAT.COUNT[4,2],col="red")

hist(PERMUTE.REPEAT.COUNT[1,3,],xlim=c(0,100),breaks=seq(0,100,2),xlab="PACE mutations in category",main="States-SR-SG-MCL1")
abline(v=REPEAT.COUNT[1,3],col="red")
hist(PERMUTE.REPEAT.COUNT[2,3,],xlim=c(0,100),breaks=seq(0,100,2),xlab="PACE mutations in category",main="States-SR-MG-MCL1")
abline(v=REPEAT.COUNT[2,3],col="red")
hist(PERMUTE.REPEAT.COUNT[3,3,],xlim=c(0,100),breaks=seq(0,100,2),xlab="PACE mutations in category",main="States-MR-SG-MCL1")
abline(v=REPEAT.COUNT[3,3],col="red")
hist(PERMUTE.REPEAT.COUNT[4,3,],xlim=c(0,100),breaks=seq(0,100,2),xlab="PACE mutations in category",main="States-MR-MG-MCL1")
abline(v=REPEAT.COUNT[4,3],col="red")

hist(PERMUTE.REPEAT.COUNT[1,4,],xlim=c(0,100),breaks=seq(0,100,2),xlab="PACE mutations in category",main="States-SR-SG-BCL2")
abline(v=REPEAT.COUNT[1,4],col="red")
hist(PERMUTE.REPEAT.COUNT[2,4,],xlim=c(0,100),breaks=seq(0,100,2),xlab="PACE mutations in category",main="States-SR-MG-BCL2")
abline(v=REPEAT.COUNT[2,4],col="red")
hist(PERMUTE.REPEAT.COUNT[3,4,],xlim=c(0,100),breaks=seq(0,100,2),xlab="PACE mutations in category",main="States-MR-SG-BCL2")
abline(v=REPEAT.COUNT[3,4],col="red")
hist(PERMUTE.REPEAT.COUNT[4,4,],xlim=c(0,100),breaks=seq(0,100,2),xlab="PACE mutations in category",main="States-MR-MG-BCL2")
abline(v=REPEAT.COUNT[4,4],col="red")

REPEAT.COUNT.P.LOW  <- matrix(0,nrow=4,ncol=4)
REPEAT.COUNT.P.HIGH <- matrix(0,nrow=4,ncol=4)
I <- 1:4
J <- 1:4
for(i in I) {
  for(j in J) {
    REPEAT.COUNT.P.LOW[i,j] <- sum(REPEAT.COUNT[i,j] < PERMUTE.REPEAT.COUNT[i,j,])
    REPEAT.COUNT.P.HIGH[i,j] <- sum(REPEAT.COUNT[i,j] > PERMUTE.REPEAT.COUNT[i,j,])
  }
}

1 - t(REPEAT.COUNT.P.LOW)/100000
1 - t(REPEAT.COUNT.P.HIGH)/100000

#P-values for MCL1/BCL2 combined
REPEAT.COUNT.COMBINED <- REPEAT.COUNT[,c(1,3)] + REPEAT.COUNT[,c(2,4)]
colnames(REPEAT.COUNT.COMBINED) <- c("SITE","STATE")

PERMUTE.REPEAT.COUNT.COMBINED <- PERMUTE.REPEAT.COUNT[,c(1,3),] + PERMUTE.REPEAT.COUNT[,c(2,4),]
dimnames(PERMUTE.REPEAT.COUNT.COMBINED)[[2]] <- c("SITE","STATE")

par(mfrow=c(2,4))
par(mar=c(2,2,2,2))
hist(PERMUTE.REPEAT.COUNT.COMBINED[1,1,],xlim=c(0,200),breaks=seq(0,200,2),xlab="PACE mutations in category",main="Sites-SR-SR")
abline(v=REPEAT.COUNT.COMBINED[1,1],col="red")
hist(PERMUTE.REPEAT.COUNT.COMBINED[2,1,],xlim=c(0,200),breaks=seq(0,200,2),xlab="PACE mutations in category",main="Sites-SR-MG")
abline(v=REPEAT.COUNT.COMBINED[2,1],col="red")
hist(PERMUTE.REPEAT.COUNT.COMBINED[3,1,],xlim=c(0,200),breaks=seq(0,200,2),xlab="PACE mutations in category",main="Sites-MR-SG")
abline(v=REPEAT.COUNT.COMBINED[3,1],col="red")
hist(PERMUTE.REPEAT.COUNT.COMBINED[4,1,],xlim=c(0,200),breaks=seq(0,200,2),xlab="PACE mutations in category",main="Sites-MR-MG")
abline(v=REPEAT.COUNT.COMBINED[4,1],col="red")

hist(PERMUTE.REPEAT.COUNT.COMBINED[1,2,],xlim=c(0,200),breaks=seq(0,200,2),xlab="PACE mutations in category",main="States-SR-SR")
abline(v=REPEAT.COUNT.COMBINED[1,2],col="red")
hist(PERMUTE.REPEAT.COUNT.COMBINED[2,2,],xlim=c(0,200),breaks=seq(0,200,2),xlab="PACE mutations in category",main="States-SR-MG")
abline(v=REPEAT.COUNT.COMBINED[2,2],col="red")
hist(PERMUTE.REPEAT.COUNT.COMBINED[3,2,],xlim=c(0,200),breaks=seq(0,200,2),xlab="PACE mutations in category",main="States-MR-SG")
abline(v=REPEAT.COUNT.COMBINED[3,2],col="red")
hist(PERMUTE.REPEAT.COUNT.COMBINED[4,2,],xlim=c(0,200),breaks=seq(0,200,2),xlab="PACE mutations in category",main="States-MR-MG")
abline(v=REPEAT.COUNT.COMBINED[4,2],col="red")

REPEAT.COUNT.COMBINED.P.LOW  <- matrix(0,nrow=4,ncol=2)
REPEAT.COUNT.COMBINED.P.HIGH <- matrix(0,nrow=4,ncol=2)
I <- 1:4
J <- 1:2
for(i in I) {
  for(j in J) {
    REPEAT.COUNT.COMBINED.P.LOW[i,j] <- sum(REPEAT.COUNT.COMBINED[i,j] < PERMUTE.REPEAT.COUNT.COMBINED[i,j,])
    REPEAT.COUNT.COMBINED.P.HIGH[i,j] <- sum(REPEAT.COUNT.COMBINED[i,j] > PERMUTE.REPEAT.COUNT.COMBINED[i,j,])
  }
}

1 - t(REPEAT.COUNT.COMBINED.P.LOW)/100000
1 - t(REPEAT.COUNT.COMBINED.P.HIGH)/100000

#P-values for repeated mutation categories combined
REPEAT.COUNT.ALL <- rbind(REPEAT.COUNT[1,],colSums(REPEAT.COUNT[2:4,]))
rownames(REPEAT.COUNT.ALL) <- c("SINGLE","REPEAT")
PERMUTE.REPEAT.COUNT.ALL <- array(apply(PERMUTE.REPEAT.COUNT,2,function(x) {rbind(x[1,],colSums(x[2:4,]))}), dim=c(2,100000,4))
PERMUTE.REPEAT.COUNT.ALL <- aperm(PERMUTE.REPEAT.COUNT.ALL,c(1,3,2),resize = TRUE)
dimnames(PERMUTE.REPEAT.COUNT.ALL)[[1]] <- c("SINGLE","REPEAT")
dimnames(PERMUTE.REPEAT.COUNT.ALL)[[2]] <- c("Sites-MCL1","Sites-BCL2","States-MCL1","States-BCL2")

par(mfrow=c(4,2))
hist(PERMUTE.REPEAT.COUNT.ALL[1,1,],xlim=c(0,100),breaks=seq(0,100,2),xlab="PACE mutations in category",main="Sites-Single-MCL1")
abline(v=REPEAT.COUNT.ALL[1,1],col="red")
hist(PERMUTE.REPEAT.COUNT.ALL[2,1,],xlim=c(0,100),breaks=seq(0,100,2),xlab="PACE mutations in category",main="Sites-Multi-MCL1")
abline(v=REPEAT.COUNT.ALL[2,1],col="red")
hist(PERMUTE.REPEAT.COUNT.ALL[1,2,],xlim=c(0,100),breaks=seq(0,100,2),xlab="PACE mutations in category",main="Sites-Single-BCL2")
abline(v=REPEAT.COUNT.ALL[1,2],col="red")
hist(PERMUTE.REPEAT.COUNT.ALL[2,2,],xlim=c(0,100),breaks=seq(0,100,2),xlab="PACE mutations in category",main="Sites-Multi-BCL2")
abline(v=REPEAT.COUNT.ALL[2,2],col="red")

hist(PERMUTE.REPEAT.COUNT.ALL[1,3,],xlim=c(0,100),breaks=seq(0,100,2),xlab="PACE mutations in category",main="States-Single-MCL1")
abline(v=REPEAT.COUNT.ALL[1,3],col="red")
hist(PERMUTE.REPEAT.COUNT.ALL[2,3,],xlim=c(0,100),breaks=seq(0,100,2),xlab="PACE mutations in category",main="States-Multi-MCL1")
abline(v=REPEAT.COUNT.ALL[2,3],col="red")
hist(PERMUTE.REPEAT.COUNT.ALL[1,4,],xlim=c(0,100),breaks=seq(0,100,2),xlab="PACE mutations in category",main="States-Single-BCL2")
abline(v=REPEAT.COUNT.ALL[1,4],col="red")
hist(PERMUTE.REPEAT.COUNT.ALL[2,4,],xlim=c(0,100),breaks=seq(0,100,2),xlab="PACE mutations in category",main="States-Multi-BCL2")
abline(v=REPEAT.COUNT.ALL[2,4],col="red")

REPEAT.COUNT.ALL.P.LOW  <- matrix(0,nrow=2,ncol=4)
REPEAT.COUNT.ALL.P.HIGH <- matrix(0,nrow=2,ncol=4)
I <- 1:2
J <- 1:4
for(i in I) {
  for(j in J) {
    REPEAT.COUNT.ALL.P.LOW[i,j] <- sum(REPEAT.COUNT.ALL[i,j] < PERMUTE.REPEAT.COUNT.ALL[i,j,])
    REPEAT.COUNT.ALL.P.HIGH[i,j] <- sum(REPEAT.COUNT.ALL[i,j] > PERMUTE.REPEAT.COUNT.ALL[i,j,])
  }
}
1 - t(REPEAT.COUNT.ALL.P.LOW)/100000
1 - t(REPEAT.COUNT.ALL.P.HIGH)/100000

#P-values for MCL1/BCL2 and repeated mutation categories combined
REPEAT.COUNT.ALL.COMBINED <- REPEAT.COUNT.ALL[,c(1,3)] + REPEAT.COUNT.ALL[,c(2,4)]
rownames(REPEAT.COUNT.ALL.COMBINED) <- c("SINGLE","REPEAT")
colnames(REPEAT.COUNT.ALL.COMBINED) <- c("SITES","STATES")
PERMUTE.REPEAT.COUNT.ALL.COMBINED <- PERMUTE.REPEAT.COUNT.ALL[,c(1,3),] + PERMUTE.REPEAT.COUNT.ALL[,c(2,4),]
dimnames(PERMUTE.REPEAT.COUNT.ALL.COMBINED)[[1]] <- c("SINGLE","REPEAT")
dimnames(PERMUTE.REPEAT.COUNT.ALL.COMBINED)[[2]] <- c("SITES","STATES")

par(mfrow=c(2,2))
hist(PERMUTE.REPEAT.COUNT.ALL.COMBINED[1,1,],xlim=c(0,200),breaks=seq(0,200,2),xlab="PACE mutations in category",main="Sites Single")
abline(v=REPEAT.COUNT.ALL.COMBINED[1,1],col="red")
hist(PERMUTE.REPEAT.COUNT.ALL.COMBINED[2,1,],xlim=c(0,200),breaks=seq(0,200,2),xlab="PACE mutations in category",main="Sites Multi")
abline(v=REPEAT.COUNT.ALL.COMBINED[2,1],col="red")
hist(PERMUTE.REPEAT.COUNT.ALL.COMBINED[1,2,],xlim=c(0,200),breaks=seq(0,200,2),xlab="PACE mutations in category",main="States Single")
abline(v=REPEAT.COUNT.ALL.COMBINED[1,2],col="red")
hist(PERMUTE.REPEAT.COUNT.ALL.COMBINED[2,2,],xlim=c(0,200),breaks=seq(0,200,2),xlab="PACE mutations in category",main="States Multi")
abline(v=REPEAT.COUNT.ALL.COMBINED[2,2],col="red")

REPEAT.COUNT.ALL.COMBINED.P.LOW  <- matrix(0,nrow=2,ncol=2)
REPEAT.COUNT.ALL.COMBINED.P.HIGH <- matrix(0,nrow=2,ncol=2)
I <- 1:2
J <- 1:2
for(i in I) {
  for(j in J) {
    REPEAT.COUNT.ALL.COMBINED.P.LOW[i,j] <- sum(REPEAT.COUNT.ALL.COMBINED[i,j] < PERMUTE.REPEAT.COUNT.ALL.COMBINED[i,j,])
    REPEAT.COUNT.ALL.COMBINED.P.HIGH[i,j] <- sum(REPEAT.COUNT.ALL.COMBINED[i,j] > PERMUTE.REPEAT.COUNT.ALL.COMBINED[i,j,])
  }
}
1 - t(REPEAT.COUNT.ALL.COMBINED.P.LOW)/100000
1 - t(REPEAT.COUNT.ALL.COMBINED.P.HIGH)/100000

#Bootstrap D statistics
T <- 1:1000
BOOTSTRAP <- list()
H <- G
H[[1]] <- integer(0)
for(t in T) {
  while(sum(sapply(H,length) < 2) > 0) {
    H <- G
    SAMPLE <- sample(unlist(G),length(unlist(G)),replace=TRUE)
    SAMPLE <- SAMPLE[order(SAMPLE)]
    for(g in 1:length(G)) {
      H[[g]] <- SAMPLE[which(SAMPLE %in% G[[g]])]
    }  
  }
  BOOTSTRAP[[t]] <- F.CALC(DATA = AA.COUNT.COMPARE, SITES = L, GENOTYPES = H)
  H[[1]] <- integer(0)
}

#save(BOOTSTRAP, file="R.BOOTSTRAP.RData")
load("R.BOOTSTRAP.RData")

T <- 1:1000
BOOTSTRAP.SUM <- matrix(0,nrow=length(T),ncol=6)
colnames(BOOTSTRAP.SUM) <- c("CHANCE_States","CHANCE_Sites",
                             "CONTIN_States","CONTIN_Sites",
                             "TOTALD_States","TOTALD_Sites")                       

for(t in T) {
  BOOTSTRAP.CHANCE.STATES <- 1/(1 - sum(BOOTSTRAP[[t]][1,]  - BOOTSTRAP[[t]][2,] )/sum(BOOTSTRAP[[t]][1,]  + (BOOTSTRAP[[t]][3,]-1 )*BOOTSTRAP[[t]][2,] ))
  BOOTSTRAP.CHANCE.SITES  <- 1/(1 - sum(BOOTSTRAP[[t]][4,]  - BOOTSTRAP[[t]][5,] )/sum(BOOTSTRAP[[t]][4,]  + (BOOTSTRAP[[t]][6,]-1 )*BOOTSTRAP[[t]][5,] ))
  BOOTSTRAP.CONTIN.STATES <- sum(2*(BOOTSTRAP[[t]][7,]- BOOTSTRAP[[t]][8,]))/sum( BOOTSTRAP[[t]][9,]* BOOTSTRAP[[t]][8,]) + 1
  BOOTSTRAP.CONTIN.SITES  <- sum(2*(BOOTSTRAP[[t]][10,]-BOOTSTRAP[[t]][11,]))/sum(BOOTSTRAP[[t]][12,]*BOOTSTRAP[[t]][11,]) + 1
  BOOTSTRAP.TOTALD.STATES <- BOOTSTRAP.CHANCE.STATES*BOOTSTRAP.CONTIN.STATES
  BOOTSTRAP.TOTALD.SITES  <- BOOTSTRAP.CHANCE.SITES*BOOTSTRAP.CONTIN.SITES

  BOOTSTRAP.SUM[t,] <- c(BOOTSTRAP.CHANCE.STATES, BOOTSTRAP.CHANCE.SITES,  
                         BOOTSTRAP.CONTIN.STATES, BOOTSTRAP.CONTIN.SITES,  
                         BOOTSTRAP.TOTALD.STATES, BOOTSTRAP.TOTALD.SITES)
}

BOOTSTRAP.SUM <- t(t(BOOTSTRAP.SUM) - D.STAT)

BOOTSTRAP.OUT <- matrix(0,nrow=6,ncol=4)
colnames(BOOTSTRAP.OUT) <- c("Mean","Var","Upper","Lower")
BOOTSTRAP.OUT[,1] <- apply(BOOTSTRAP.SUM,2,mean,na.rm=TRUE)
BOOTSTRAP.OUT[,2] <- apply(BOOTSTRAP.SUM,2,var,na.rm=TRUE)
BOOTSTRAP.OUT[,3] <- apply(BOOTSTRAP.SUM,2,mean,na.rm=TRUE) - apply(BOOTSTRAP.SUM,2,function(x) {quantile(x,0.025,na.rm=TRUE)})
BOOTSTRAP.OUT[,4] <- apply(BOOTSTRAP.SUM,2,mean,na.rm=TRUE) - apply(BOOTSTRAP.SUM,2,function(x) {quantile(x,0.975,na.rm=TRUE)})
  
barx <- barplot(t(D.STAT-1),axis.lty=1,ylim=c(0,2),ylab="D",xlab="",main="")
error.bar(barx,t(D.STAT-1),upper=t(D.STAT-1) + BOOTSTRAP.OUT[,3],lower=t(D.STAT-1) + BOOTSTRAP.OUT[,4],length=0.03)


#Checking Contingency variance due to finite sampling of populations
DST2 <- F.CALC(DATA = AA.COUNT.COMPARE, SITES = L, GENOTYPES = G[c(1,2)])
D.STAT2 <- numeric(6)
D.STAT2[1] <- sum(2*(DST2[1,]- DST2[2,]))/sum( DST2[3,]* DST2[2,]) + 1
D.STAT2[2] <- sum(2*(DST2[4,]- DST2[5,]))/sum( DST2[6,]* DST2[5,]) + 1
D.STAT2[3] <- sum(2*(DST2[7,]- DST2[8,]))/sum( DST2[9,]* DST2[8,]) + 1
D.STAT2[4] <- sum(2*(DST2[10,]-DST2[11,]))/sum(DST2[12,]*DST2[11,]) + 1
D.STAT2[5] <- sum(2*(DST2[13,]-DST2[14,]))/sum(DST2[15,]*DST2[14,]) + 1
D.STAT2[6] <- sum(2*(DST2[16,]-DST2[17,]))/sum(DST2[18,]*DST2[17,]) + 1

I <- 1:5
D.OUT <- array(0,dim=c(6,6,6))

for(i in I) {
  J <- (i+1):6
  for(j in J) {
    K <- 1:1000
    D.STAT2 <- matrix(0,nrow=6,ncol=length(K))
    for(k in K) {
      SAMPLE <- sample(unlist(G[c(i,j)]),length(G[[i]]) + length(G[[j]]),FALSE)
      TEST.LIST <- list("X"=unname(SAMPLE[1:length(G[[i]])]),"Y"=unname(SAMPLE[(length(G[[i]])+1):length(SAMPLE)]))
      DST2 <- F.CALC(DATA = AA.COUNT.COMPARE, SITES = L, GENOTYPES = TEST.LIST)
      D.STAT2[1,k] <- 1/(1 - sum(DST2[1,]  - DST2[2,] )/sum(DST2[1,]  + (DST2[3,]-1 )*DST2[2,] ))
      D.STAT2[2,k] <- 1/(1 - sum(DST2[4,]  - DST2[5,] )/sum(DST2[4,]  + (DST2[6,]-1 )*DST2[5,] ))
      D.STAT2[3,k] <- sum(2*(DST2[7,]- DST2[8,]))/sum( DST2[9,]* DST2[8,]) + 1
      D.STAT2[4,k] <- sum(2*(DST2[10,]-DST2[11,]))/sum(DST2[12,]*DST2[11,]) + 1
      D.STAT2[5,k] <- D.STAT2[1,k]*D.STAT2[3,k]
      D.STAT2[6,k] <- D.STAT2[2,k]*D.STAT2[4,k]
    }
    D.OUT[i,j,] <- apply(D.STAT2,1,mean)
  }
}
#save(D.OUT,file="D.OUT.RData")

plot(  PHYLO.DIST[row(PHYLO.DIST) <= col(PHYLO.DIST) - 1],CHANCE.STATES.MATRIX[row(CHANCE.STATES.MATRIX) <= col(CHANCE.STATES.MATRIX) - 1],pch=19,cex=0.8,xlim=c(0,7),ylim=c(1,4),col=c("#FF0000"),xlab="Phylogenetic Distance",ylab="R")
points(PHYLO.DIST[row(PHYLO.DIST) <= col(PHYLO.DIST) - 1],CONTIN.STATES.MATRIX[row(CONTIN.STATES.MATRIX) <= col(CONTIN.STATES.MATRIX) - 1],pch=19,cex=0.8,col="#0000FF")
points(PHYLO.DIST[row(PHYLO.DIST) <= col(PHYLO.DIST) - 1],TOTALD.STATES.MATRIX[row(TOTALD.STATES.MATRIX) <= col(TOTALD.STATES.MATRIX) - 1],pch=19,cex=0.8,col="#00FF00")

points(PHYLO.DIST[row(PHYLO.DIST) <= col(PHYLO.DIST) - 1],D.OUT[,,1][row(D.OUT[,,1]) <= col(D.OUT[,,1]) - 1],pch=18,cex=0.6,col="#FF0000")
points(PHYLO.DIST[row(PHYLO.DIST) <= col(PHYLO.DIST) - 1],D.OUT[,,3][row(D.OUT[,,3]) <= col(D.OUT[,,3]) - 1],pch=18,cex=0.6,col="#0000FF")
points(PHYLO.DIST[row(PHYLO.DIST) <= col(PHYLO.DIST) - 1],D.OUT[,,5][row(D.OUT[,,5]) <= col(D.OUT[,,5]) - 1],pch=18,cex=0.6,col="#00FF00")

