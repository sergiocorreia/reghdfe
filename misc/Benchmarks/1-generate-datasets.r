library(data.table)
K <- 100
set.seed(1)
for (file in c("2e6", "1e7", "1e8")){
	N <- as.integer(file)
	DT <- data.table(
	  id1 = sample(sprintf("id%03d",1:K), N, TRUE),      # large groups (char)
	  id2 = sample(sprintf("id%03d",1:K), N, TRUE),      # large groups (char)
	  id3 = sample(sprintf("id%010d",1:(N/K)), N, TRUE), # small groups (char)
	  id4 = sample(K, N, TRUE),                          # large groups (int)
	  id5 = sample(K, N, TRUE),                          # large groups (int)
	  id6 = sample(N/K, N, TRUE),                        # small groups (int)
	  v1 =  sample(5, N, TRUE),                          # int in range [1,5]
	  v2 =  sample(1e6, N, TRUE),                        # int in range [1,1e6]
	  v3 =  sample(round(runif(100,max=100),4), N, TRUE) # numeric e.g. 23.5749
	)
	write.table(DT,paste0(file,".csv"),row.names=F, sep="\t")

	if (file == "2e6"){
		DT_merge <- unique(DT, by = c("id1", "id3"))
		write.table(DT_merge,"merge.csv",row.names=F, sep="\t")
	}
}