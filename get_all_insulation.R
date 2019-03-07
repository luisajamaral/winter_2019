library(ggplot2)

get_ratio <- function(file, pos) { 
  all_ints = read.csv(file,sep = "\t", header = F)
  all_ints = all_ints[,1:6]
  colnames(all_ints) = c("chr", "start1", "end1", "chr2", "start2", "end2")
  all_ints$dist = all_ints$start2- all_ints$start1
  all_ints = all_ints[-which(all_ints$dist < 10000), ]
  up = all_ints[which(all_ints$start1 < pos),]
  across = up[which(up$start2 > pos),]
 # cat("across insert location:", nrow(across), "\n")
  up_only = up[which(up$start2 < pos),]
 # cat("within upstream:", nrow(up_only), "\n")
  down =  all_ints[which(all_ints$start1 >pos),]
  down_only =  down[which(down$start2 >pos),]
#  cat("within downstream:", nrow(down_only), "\n")
  within = nrow(down_only) + nrow(up_only)
 # cat("Ratio of across insert vs on 100kb either side:", nrow(across)/within, "\n")
  ratioA = nrow(across)/within
  return(ratioA)
  #Aratios = c(Aratios, ratioA)  
}

get_insulation <- function(pgl_id,pgl_idB,insert_locs,chr, position, allele, status,rep1, rep2,clone) {
  bin = position/10000
  Aratios = c()
  Bratios = c()
  for (pos in seq((bin-10),(bin+10),1) ) {
    #cat(pos, "\n")
    comm = paste(chr, (pos-10)*10000, (pos+10)*10000 , chr, (pos-10)*10000, (pos+10)*10000, "1", "2", sep = "\t")
    sink(file = "t.pgl")
    cat(comm)
    sink()
    pgl_comm = paste("pgltools intersect -a t.pgl -b ", pgl_id," > A", sep = "")
    system(pgl_comm, intern = T)
    
    pgl_commB = paste("pgltools intersect -a t.pgl -b ",pgl_idB," > B", sep = "")
    system(pgl_commB, intern = T)
    pos = pos*10000
    info = file.info("A")
    empty = (is.na(info$size) || info$size == 0) 
    if (empty) {
      Aratios = c(Aratios, NA)
    }
    else {
      Aratios = c(Aratios, get_ratio("A",pos))
    }
    info = file.info("B")
    empty = (is.na(info$size) || info$size == 0) 
    if (empty) {
      Bratios = c(c(Bratios, NA))       
    } 
    else {
      Bratios = c(c(Bratios, get_ratio("B",pos)))       
    }
  }
  
  cat(Aratios, "\n")
  cat(Bratios, "\n")
  if(allele == status) {
    write.table(Aratios, paste("IN_ratios",status,position,rep1, sep = "_"), quote = F, row.names = F,col.names = F)
    write.table(Bratios, paste("IN_ratios",status, position, rep2,sep = "_"), quote = F , row.names = F,col.names = F)
    df = data.frame(rep1 = Aratios,rep2 = Bratios)
    final = rowMeans(df)
    write.table(final, paste("IN_avg_ratios",status, position, clone,sep = "_"), quote = F , row.names = F,col.names = F)
    
  }
  if(allele != status) {
    write.table(Bratios, paste("OUT_ratios",status,position, rep2,sep = "_"), quote = F , row.names = F,col.names = F)
    write.table(Aratios, paste("OUT_ratios",status,position, rep1,sep = "_"), quote = F , row.names = F,col.names = F)
    df = data.frame(rep1 = Aratios,rep2 = Bratios)
    final = rowMeans(df)
    write.table(final, paste("OUT_avg_ratios",status, position,clone, sep = "_"), quote = F , row.names = F,col.names = F)
  }
  return(list(Aratios,Bratios)) 
}


get_allele <- function(file_name1A,file_name1B,file_name2A,file_name2B, position) {
  counts = c()
  for(file in c(file_name1A,file_name1B,file_name2A,file_name2B)) {
    one_map = read.table(file, sep = " ")
    count = 0
    for (i in 1:nrow(one_map)) { 
      if (one_map$V2[i] < (position + 500000) && one_map$V2[i] > (position - 500000)) {
        count = count + 1
      }
    }
    counts = c(counts, count)
  }
  A_counts = counts[1] + counts[3]
  B_counts = counts[2] + counts[4]
  cat(counts, "\n")
  cat("A: ", A_counts, "| B: ", B_counts, "\n")
  if(A_counts == B_counts) {
    cat("nope: ", position, "\n")
    return(NA)
  }
  if(A_counts > B_counts) {
    return("A")
  }
  if(A_counts < B_counts) {
    return("B")
  }
}


#pieces = strsplit(system("pwd",intern=T),"/")[[1]]
#chr = pieces[length(pieces)]
#chr = strsplit(pgl_id, "_")[[1]][1]
insert_locs = read.table(file = "all_mid_ann.txt", header = T)
##insert_locs = insert_locs[which(insert_locs$chr == chr),]

for (p in c(27,28,29,43,25,42,23,24)) {
  chr = insert_locs$chr[p]
  position = insert_locs$Mid[p]
  prefix = paste(chr, "/", sep = "")
  clone = insert_locs$KI[p]
  if(clone == "KI1") {
    file_name1A = paste(prefix,chr, "_","RH942","_","Aallele.bam_one_mapped_herv.txt", sep= "")
    file_name1B = paste(prefix,chr, "_","RH942","_","Ballele.bam_one_mapped_herv.txt", sep= "")
    file_name2A = paste(prefix,chr, "_","RH943","_","Aallele.bam_one_mapped_herv.txt", sep= "")
    file_name2B = paste(prefix,chr, "_","RH943","_","Ballele.bam_one_mapped_herv.txt", sep= "")
    allele = get_allele(file_name1A,file_name1B,file_name2A,file_name2B, position)
    if(is.na(allele)) {
      next
    }
    rep1 = "RH942"
    rep2 = "RH943"
    
    pgl_id1A = paste(prefix,chr, "_",rep1,"_","Aallele.pgl", sep= "")
    pgl_id2A = paste(prefix,chr, "_",rep2,"_","Aallele.pgl", sep= "")
    A_allele = get_insulation(pgl_id1A,pgl_id2A,insert_locs,chr,position,allele, "A",rep1,rep2,clone) 
    pgl_id1B = paste(prefix,chr, "_",rep1,"_","Ballele.pgl", sep= "")
    pgl_id2B = paste(prefix,chr, "_",rep2,"_","Ballele.pgl", sep= "")
    B_allele = get_insulation(pgl_id1B,pgl_id2B,insert_locs,chr,position,allele,"B",rep1,rep2,clone) 
    
  }
  
  if(clone == "KI2") {
    file_name1A = paste(prefix,chr, "_","RH944","_","Aallele.bam_one_mapped_herv.txt", sep= "")
    file_name1B = paste(prefix,chr, "_","RH944","_","Ballele.bam_one_mapped_herv.txt", sep= "")
    file_name2A = paste(prefix,chr, "_","RH945","_","Aallele.bam_one_mapped_herv.txt", sep= "")
    file_name2B = paste(prefix,chr, "_","RH945","_","Ballele.bam_one_mapped_herv.txt", sep= "")
    allele = get_allele(file_name1A,file_name1B,file_name2A,file_name2B, position)
    if(is.na(allele)) {
      next
    }
    rep1 = "RH944"
    rep2 = "RH945"
    
    pgl_id1A = paste(prefix,chr, "_",rep1,"_","Aallele.pgl", sep= "")
    pgl_id2A = paste(prefix,chr, "_",rep2,"_","Aallele.pgl", sep= "")
    A_allele = get_insulation(pgl_id1A,pgl_id2A,insert_locs,chr,position,allele, "A",rep1,rep2,clone) 
    pgl_id1B = paste(prefix,chr, "_",rep1,"_","Ballele.pgl", sep= "")
    pgl_id2B = paste(prefix,chr, "_",rep2,"_","Ballele.pgl", sep= "")
    B_allele = get_insulation(pgl_id1B,pgl_id2B,insert_locs,chr,position,allele,"B",rep1,rep2,clone) 
  }
}
