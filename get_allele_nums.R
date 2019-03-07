

get_allele <- function(file_name1A,file_name1B,file_name2A,file_name2B, chr, position) {
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
  #cat(counts, "\n")
  cat(paste(chr, position , A_counts, B_counts, sep = "\t") ,"\n" )
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

for (p in 1:nrow(insert_locs)) {
  chr = insert_locs$chr[p]
  position = insert_locs$Mid[p]
  prefix = paste(chr, "/", sep = "")
  clone = insert_locs$KI[p]
  if(clone == "KI1") {
    file_name1A = paste(prefix,chr, "_","RH942","_","Aallele.bam_one_mapped_herv.txt", sep= "")
    file_name1B = paste(prefix,chr, "_","RH942","_","Ballele.bam_one_mapped_herv.txt", sep= "")
    file_name2A = paste(prefix,chr, "_","RH943","_","Aallele.bam_one_mapped_herv.txt", sep= "")
    file_name2B = paste(prefix,chr, "_","RH943","_","Ballele.bam_one_mapped_herv.txt", sep= "")
    allele = get_allele(file_name1A,file_name1B,file_name2A,file_name2B, chr, position)
    if(is.na(allele)) {
      next
    }
  }
  
  if(clone == "KI2") {
    file_name1A = paste(prefix,chr, "_","RH944","_","Aallele.bam_one_mapped_herv.txt", sep= "")
    file_name1B = paste(prefix,chr, "_","RH944","_","Ballele.bam_one_mapped_herv.txt", sep= "")
    file_name2A = paste(prefix,chr, "_","RH945","_","Aallele.bam_one_mapped_herv.txt", sep= "")
    file_name2B = paste(prefix,chr, "_","RH945","_","Ballele.bam_one_mapped_herv.txt", sep= "")
    allele = get_allele(file_name1A,file_name1B,file_name2A,file_name2B,chr,position)
    if(is.na(allele)) {
      next
    }
  }
}
