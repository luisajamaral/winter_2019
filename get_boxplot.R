

diff = read.csv("diffOutputCondensed_edgeR.txt", sep = "\t")
diff = diff[-which(rowMeans(diff[,9:12])<4),]
locs = read.table("KI2_KI1_strand_mid.txt", sep = "\t", header = T)
locs$st_mid = locs$Mid-100000
locs$end_mid = locs$Mid+100000
locs$start= locs$Mid-1000000
locs$end= locs$Mid+1000000

KI1_logfcs_500up = c()
KI1_logfcs_500dwn = c()
KI1_insert_region = c()

KI2_logfcs_500up = c()
KI2_logfcs_500dwn = c()
KI2_insert_region = c()
not_near = c()
for(i in 1:nrow(diff)) {
  curr = locs[which(locs$chr == paste(diff$chr[i])), ]
  d_start = diff$start[i]
  d_end = diff$end[i]
  if (nrow(curr)!=0 && !is.na(diff$KI1.vs..KI2.Log2.Fold.Change[i])) {
    for (j in 1:nrow(curr)) {
      if (d_start < curr$st_mid[j] && d_end < curr$st_mid[j]) {
        if (d_start > curr$start[j]) {
          #this gene is 500kb upstream
          if (curr$X[j] == "KI1") {
            if (curr$strand[j] == "+") {
              KI1_logfcs_500up = c(KI1_logfcs_500up, diff$KI1.vs..KI2.Log2.Fold.Change[i])
            }
            if (curr$strand[j] == "-") {
              KI1_logfcs_500dwn = c(KI1_logfcs_500dwn, diff$KI1.vs..KI2.Log2.Fold.Change[i])
            }
          }
          if (curr$X[j] == "KI2") {
            if (curr$strand[j] == "+") {
              KI2_logfcs_500up = c(KI2_logfcs_500up, diff$KI1.vs..KI2.Log2.Fold.Change[i])
            }
            if (curr$strand[j] == "-") {
              KI2_logfcs_500dwn = c(KI2_logfcs_500dwn, diff$KI1.vs..KI2.Log2.Fold.Change[i])
            }
          }
        }
        else {
          #this gene is not near an insert
          not_near = c(not_near, diff$KI1.vs..KI2.Log2.Fold.Change[i])
          
        }
      }
      else if (d_start > curr$end_mid[j] && d_end > curr$end_mid[j]) {
        if (d_start < curr$end[j] ) {
          #this gene is 500kb downstream
          if (curr$X[j] == "KI1") {
            if (curr$strand[j] == "-") {
              KI1_logfcs_500up = c(KI1_logfcs_500up, diff$KI1.vs..KI2.Log2.Fold.Change[i])
            }
            if (curr$strand[j] == "+") {
              KI1_logfcs_500dwn = c(KI1_logfcs_500dwn, diff$KI1.vs..KI2.Log2.Fold.Change[i])
            }
          }
          if (curr$X[j] == "KI2") {
            if (curr$strand[j] == "-") {
              KI2_logfcs_500up = c(KI2_logfcs_500up, diff$KI1.vs..KI2.Log2.Fold.Change[i])
            }
            if (curr$strand[j] == "+") {
              KI2_logfcs_500dwn = c(KI2_logfcs_500dwn, diff$KI1.vs..KI2.Log2.Fold.Change[i])
            }
          }
        }
        else {
          #this gene is not near an insert
          not_near = c(not_near, diff$KI1.vs..KI2.Log2.Fold.Change[i])
        }
      }
      else {
        #this gene is within 100kb of the insert/ could be chimeric
        if (curr$X[j] == "KI1") {
          KI1_insert_region = c(KI1_insert_region, diff$KI1.vs..KI2.Log2.Fold.Change[i])
        }
        if (curr$X[j] == "KI2") {
          KI2_insert_region = c(KI2_insert_region, diff$KI1.vs..KI2.Log2.Fold.Change[i])
        }
        cat(paste(diff$chr[i]),
            diff$start[i],
            diff$end[i],
            diff$KI1.vs..KI2.Log2.Fold.Change[i],
            "\n")
      }
    }
  }
}

KI1_logfcs_500up = KI1_logfcs_500up*-1
KI1_logfcs_500dwn = KI1_logfcs_500dwn*-1
KI1_insert_region = KI1_insert_region*-1

#KI1_insert_region = KI1_insert_region[-1]


plot_ly(y = KI1_logfcs_500up, type = "box", name = "1Mb Upstream", jitter = 0.3, pointpos = -1.8, boxpoints = 'all') %>%
  add_trace(y = KI1_logfcs_500dwn, name = "1Mb Downstream") %>%
  add_trace(y = KI1_insert_region, name = "Within 200kb") %>%
  #add_trace(y = not_near, name = "Not near") %>%
  
  layout(title = "KI1",yaxis = list(title = "logFC in KI1 vs KI2", range = c(-2, 5)) )


#KI2_insert_region = KI2_insert_region[-1]

plot_ly(y = KI2_logfcs_500up, type = "box", name = "1Mb Upstream", jitter = 0.3, pointpos = -1.8, boxpoints = 'all') %>%
  add_trace(y = KI2_logfcs_500dwn, name = "1Mb Downstream") %>%
  add_trace(y = KI2_insert_region, name = "Within 200kb") %>%
 # add_trace(y = diff$KI1.vs..KI2.Log2.Fold.Change, name = "Not near") %>%
  
  layout(title = "KI2",yaxis = list(title = "logFC in KI2 vs KI1", range = c(-2, 5)) )





#plot_ly(y = logfcs_500up, type = "box", name = "100-900kb Upstream") %>%
#  add_trace(y = logfcs_500dwn, name = "100-900kb Downstream") %>%
#  add_trace(y = insert_region, name = "Within 100kb")
  #add_trace(y = not_near, name = "not near insert")

