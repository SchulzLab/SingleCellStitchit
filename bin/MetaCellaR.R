### To summerize single cells based on k-medoids ###
### The clustering is done per cell type ###
### parameter k is defined based on the number of potential groups that have at least 30 cells inside (k <- ncol(CT_cluster) / 30)

set.seed(0)
library(cluster)
library(knn.covertree)
library(irlba)
library(data.table)
library(ggplot2)
library(ggdendro)
library(cowplot)
library(Matrix)
########################
summary_method <- "kmed_means" #kmed
iter_flag <- F
csv_flag <- F
merge_flag <- T
#args <- commandArgs(trailingOnly= T)
#file_name <- args[1] # Path to where the Seurat object is stored
#RNA_count_expression <- args[2]
#celltype_expression <- args[3]
########################
## Parse arguments

argsexpr <- commandArgs(trailingOnly= T)
#argsexpr <- c("-file /projects/expregulation/work/singleCell/heart_data/Metacell_Rda_files/3005_.Rds", "-RNA assays$RNA@counts", "-celltype meta.data$celltype", "-output MetaCellar_3005", "-umap T")
#argsexpr <- c("-file /projects/expregulation/work/singleCell/heart_data/Metacell_Rda_files/3008_.Rds", "-RNA assays$RNA@counts", "-celltype meta.data$celltype", "-output MetaCellar_in_GAZE_3008", "-umap T", "-e 30", "-d 2")
#argsexpr <- c("-file /projects/triangulate/archive/rna_atac_merged_coembedded_harmonized.rds", "-RNA assays$RNA@counts", "-celltype meta.data$Celltypes_refined", "-output testUMAPfullKoptimized_100_UMAP2", "-umap T", "-assay meta.data$datasets2", "-ATAC assays$peaks@counts", "-e 100", "-d 20")
#-RNA 'assays$RNA@data' -celltype 'meta.data$Celltypes_refined'

defined_args <- c("-file", "-RNA", "-celltype", "-output", "-k", "-assay", "-umap", "-ATAC", "-e", "-t", "-d", "-reduction")
arg_tokens <- unlist(strsplit(argsexpr, split= " "))
file_hit <- which(arg_tokens == defined_args[1])
RNA_hit <- which(arg_tokens == defined_args[2])
celltype_hit <- which(arg_tokens == defined_args[3])
output_hit <- which(arg_tokens == defined_args[4])
k_hit <- which(arg_tokens == defined_args[5])
assay_hit <- which(arg_tokens == defined_args[6])
umap_hit <- which(arg_tokens == defined_args[7])
ATAC_hit <- which(arg_tokens == defined_args[8])
expCnt_hit <- which(arg_tokens == defined_args[9])
t_hit <- which(arg_tokens == defined_args[10])
d_hit <- which(arg_tokens == defined_args[11])
reduction_hit <- which(arg_tokens == defined_args[12])
if(length(reduction_hit)){
	reduction <- arg_tokens[reduction_hit + 1]
}else{
	reduction <- "umap" ## default UMAP that should be in the Seurat integration
}
if(length(umap_hit)){
	umap_flag <- as.logical(arg_tokens[umap_hit + 1])
}else{
	umap_flag <- F
}
if(length(file_hit)){
  if(length(RNA_hit) && length(celltype_hit)){
    file_name <- arg_tokens[file_hit + 1]
    RNA_count_expression <- arg_tokens[RNA_hit + 1]
    ATAC_count_expression <- arg_tokens[ATAC_hit + 1]
    celltype_expression <- arg_tokens[celltype_hit + 1]
    assay_expression <- arg_tokens[assay_hit + 1]
  }
  else{# CSV file
    if(!length(celltype_hit)){
      stop("For the csv input file the -celltype options must be present! Please provide a two column file with headers, linking the cell names (first column) to cell types (second columns).")
    }
    print(paste("Reading", arg_tokens[file_hit + 1], "..."))
    csv_flag <- T
    csv_data <- read.csv(arg_tokens[file_hit + 1], row.names= 1)
    csv_cells <- read.table(arg_tokens[celltype_hit + 1], header= T)
		general_data <- csv_data
    print("Done!")
  }
}else{
  stop("Argument -file is missing. Please provide the path to your input file using the -file option followed by the path to your file!")
}
if(!length(output_hit)){
  output_file <- getwd()
}else{
  output_file <- arg_tokens[output_hit + 1]
}
if(length(expCnt_hit)){
	expected_cells <- as.numeric(arg_tokens[expCnt_hit + 1])
}else{
	expected_cells <- 30
}
if(length(t_hit)){
	threshold <- as.numeric(arg_tokens[t_hit + 1])
}else{
	threshold <- 3 * expected_cells
}
if(length(d_hit)){
	umap_dim <- as.numeric(arg_tokens[d_hit + 1])
}else{
	umap_dim <- 20
}
########################
if(!csv_flag){
  library(Seurat)
  if(strsplit(file_name, "\\.")[[1]][2] == "rds"){
    print("Reading the Seurat object")
    Sub <- readRDS(file_name)
    Rds_name <- "Sub"
  }else{
    print("Loading the Seurat object")
    Rds_name <- load(file_name)
    Sub <- get(Rds_name) ##To make it generalizable for Suerat object stored with any name
  }
  print("Done loading the Seurat object")
  print("Done loading Rds!")
  celltypes <- eval(parse(text= paste0(Rds_name, "@", celltype_expression)))
	original_celltypes <- celltypes
	## Identify NA cell annotation and remove them from the object
	na_idx <- is.na(celltypes)
	if(sum(na_idx)){
		## TODO:
		Sub <- subset(Sub, !is.na(celltypes))## I know that it doesn't work as I've already tested it on a Seurat obj... gotta find another way
	}
	print(paste("RNA_count_expression:", RNA_count_expression))
  RNAcounts <- eval(parse(text= paste0(Rds_name, "@", RNA_count_expression))) #Sub@assays$RNA@counts
	gene_names <- rownames(RNAcounts)
	RNAcounts <- as.matrix(RNAcounts)
  cell_types <- unique(as.character(celltypes))
  if(length(assay_hit)){
  	ATACcounts <- eval(parse(text= paste0(Rds_name, "@", ATAC_count_expression))) #Sub@assays$RNA@counts
	  assays <- eval(parse(text= paste0(Rds_name, "@", assay_expression)))

	  rna_hits <- which(tolower(assays) == "scrna-seq" | tolower(assays) == "scrna" | tolower(assays) == "scrna" | tolower(assays) == "rna")
	  ATACcounts <- ATACcounts[, -rna_hits]
	  RNAcounts <- RNAcounts[, rna_hits]
		general_data <- RNAcounts
	  RNA_celltypes <- celltypes[rna_hits]
	  ATAC_celltypes <- celltypes[-rna_hits]
		RNA_umap <- Sub[[reduction]]@cell.embeddings[rna_hits, ]
		ATAC_umap <- Sub[[reduction]]@cell.embeddings[-rna_hits, ]
		print("ATAC_umap")
		print(head(ATAC_umap))
	  celltypes <- RNA_celltypes
	  cell_types <- unique(celltypes)
  }
}else{
  cell_types <- unique(csv_cells[, 2])
}
if(umap_flag){
	print("Computing UMAP")
	ptm <- proc.time()
	if(csv_flag){
		## Florian's way ##
		#pca_res <- prcomp_irlba(csv_data, n= 30, scale.= T)
		#umap_res <- umap::umap(t(csv_data), n_components= 20)
		#umap_res <- umap::umap(pca_res$rotation, n_components= 20, n_threads= 10)
		##Seurat way###
    ptm_umap <- proc.time(); 
    umap_layout <- uwot::umap(t(as.matrix(csv_data)), pca= 30, pca_center= T, scale= T, n_components= umap_dim)
		rownames(umap_layout) <- colnames(csv_data)
		print(paste("UMAP computation time:", proc.time() - ptm_umap))
		celltypes <- csv_cells[, 2]
	}else{
		## Florian's way ##
		#print("Running PCA")
		#pca_res <- prcomp_irlba(RNAcounts, n=30, scale.= T)
		#print("Done running PCA")
		#print(paste("PCA computation time:", (proc.time() - ptm)))
		#umap_res <- umap::umap(t(as.matrix(RNAcounts)), n_components= 20)
		#ptm2 <- proc.time()
		#print("Starting UMAP calculation...")
		#umap_res <- umap::umap(pca_res$rotation, n_components= 20, n_threads= 10)
		#print("UMAP calculation done!")
		#print(paste("UMAP computation time:", (proc.time() - ptm)))
		##Seurat way###
		ptm_umap <- proc.time();
		umap_layout <- uwot::umap(t(as.matrix(RNAcounts)), pca= 30, pca_center= T, scale= T, n_components= umap_dim)
		print(paste("UMAP computation time:", proc.time() - ptm_umap))
		rownames(umap_layout) <- colnames(RNAcounts)
	}
	RNA_umap <- umap_layout
}
clusters <- list()
all_mediods <- list()
RNA_metacell_umap <- NULL
cell2metacell_info <- NULL
cluster_data <- list()
print(dim(RNAcounts))

####### Begin functions #######
################################
my_knn <- function(train, test, k= 1, threshold= 3 * expected_cells){
	dist_mat <- matrix(NA, nrow= nrow(test), ncol= nrow(train)) ## Matrix holding the distances of any ATAC point to any RNA metacell
	MC_tracking <-seq(ncol(dist_mat)) * 0
	for(i in seq(nrow(test))){
		## Compute the distance between each test point and all the training points
		rep_test <- matrix(rep(test[i, ], nrow(train)), byrow= T, nrow= nrow(train), ncol= ncol(train))
		dist <- sqrt(rowSums((rep_test - train)^2)) ## Euclidean distance
		dist_mat[i, ] <- dist		
	}

	### Assign the labels while keeping track of the label satuaration for any specific metacell
	i <- 1;
	labels <- vector("character", length= nrow(test))
	while(i <= nrow(dist_mat)){
		hit <- which.min(dist_mat[i, ]);
		if(MC_tracking[hit] < threshold){
			MC_tracking[hit] <- MC_tracking[hit] + 1;
			labels[i] <- rownames(train)[hit];
			i <- i + 1
		}else{
			dist_mat[i, hit] <- Inf;
		}
	}
	return(labels)

}

#####

invoke_clara <- function(CT_cluster, original_CT_cluster, iter_flag, clusters, RNA_metacell_umap, ct, k){
	class(CT_cluster)
	if(iter_flag){
		for(iter in seq(10)){
			clara_res <- clara(t(as.matrix(CT_cluster)), k, metric = "euclidean", stand = FALSE, samples= 30, pamLike = FALSE)
			dim(clara_res$medoids)
			dim(CT_cluster)

			clusters[[ct]] <- clara_res
			all_mediods[[ct]] <- rbind(all_mediods[[ct]], clara_res$medoids)
		}
	}
	else{
		clusters[[ct]] <- clara(t(as.matrix(CT_cluster)), k, metric = "euclidean", stand = FALSE, samples= 30, pamLike = FALSE)
	}
	if(length(assay_hit)){
		RNA_metacell_umap_ct <- NULL
		for(i in unique(clusters[[ct]]$clustering)){
			data_subset <- RNA_umap[which(clusters[[ct]]$clustering == i), ];
			if(is.null(dim(data_subset))){
				RNA_metacell_umap_ct <- rbind(RNA_metacell_umap_ct, data_subset);
			}else{
				RNA_metacell_umap_ct <- rbind(RNA_metacell_umap_ct, colMeans(data_subset))
			}
		}
		rownames(RNA_metacell_umap_ct) <- paste(ct, seq(nrow(RNA_metacell_umap_ct)), sep= "_")
		print(paste("Done clustering", ct))
		RNA_metacell_umap <- rbind(RNA_metacell_umap, RNA_metacell_umap_ct)
	}
	#return(list(clusters= clusters, RNA_metacell_umap= RNA_metacell_umap))
	return(clusters= clusters)
}

#####

invoke_clara_simplified <- function(CT_cluster, k){
	class(CT_cluster)
	cluster_data <- t(as.matrix(CT_cluster))
	clusters <- clara(t(as.matrix(CT_cluster)), k, metric = "euclidean", stand = FALSE, samples= 30, pamLike = FALSE)
	return(list(clusters= clusters))
}
######
get_k_grid <- function(ncells, expected_cells){
	mid_k <- floor(ncells / expected_cells)
	k <- seq(floor(mid_k * .75), mid_k * 1.5, by= 2)
	## If k is an even number, make them odd:
	if(k[1] %% 2 == 0){
		k <- k + 1
	}
	return(k)
}

metacell_density <- function(clara_out, k_search){
	densities <- NULL
	for(i in seq(ncol(clara_out))){
		expected_cell_num <- length(clara_out[1, 1][[1]][[1]]$clustering) / k_search[i]
		metacell_counts <- table(clara_out[1, i][[1]][[1]]$clustering)
		densities <- c(densities, sum(metacell_counts > expected_cell_num) / length(unique(clara_out[1, i][[1]][[1]]$clustering)))
	}
	return(densities)
}

########
metacell_density1 <- function(clara_out, k_search){
	densities <- NULL
	for(i in seq(length(k_search))){
		expected_cell_num <- length(clara_out[[i]]$clustering) / k_search[i]
		metacell_counts <- table(clara_out[[i]]$clustering)
		#densities <- c(densities, sum(metacell_counts > expected_cell_num) / length(unique(clara_out[[i]]$clustering)))
		densities <- c(densities, sum(metacell_counts > expected_cell_num))
	}
	return(densities)
}
########

metacell_quality <- function(clara_out, expected_cell_num= 30, slack_ratio= .15){
	cell_num_slack <- expected_cell_num - floor(expected_cell_num * slack_ratio)
	metacell_counts <- table(clara_out$clustering)
	qualities <- metacell_counts > cell_num_slack
	return(qualities)
}

get_outliers <- function(clara_out, outlier_cell_num= 10){
	metacell_counts <- table(clara_out$clustering)
	return(which(metacell_counts < outlier_cell_num))
}

########
cluster_means <- function(clusters){
	temp_cl <- NULL
	for(clst in unique(clusters$clustering)){
		if("numeric" %in% class(clusters$data[clusters$cluster == clst, ])){
			temp_cl <- rbind(temp_cl, clusters$data[clusters$cluster == clst, ])
		}
		else{
			temp_cl <- rbind(temp_cl, apply(clusters$data[clusters$cluster == clst, ], 2, FUN= mean))
		}
	}
	return(temp_cl)
}
########
merge_small_mc <- function(clustering, thresh= 30){
	flag <- T
	while(flag){
		metacell_counts <- table(clustering)
		mc_table_sorted <- sort(metacell_counts)
		low_mc_idx <- which(metacell_counts < thresh)
		if(!length(low_mc_idx)){
			return(clustering)
		}
		idx <- order(metacell_counts[low_mc_idx], decreasing= F)
		new_mc <- NULL;
		mc_sum <- metacell_counts[low_mc_idx[idx][1]]
		new_mc <- names(metacell_counts[low_mc_idx[idx][1]])
		max_cluster_id <- max(clustering) + 1
		pointer <- 2
		while(mc_sum < thresh){
			if(length(idx) == 1){
				all_idx <- order(metacell_counts, decreasing= F)
				new_mc <- c(new_mc, names(metacell_counts[all_idx][2]))
				clustering[which(clustering %in% as.numeric(new_mc))] <- max_cluster_id
				flag <- F
				break;
			}
			mc_sum <- mc_sum + metacell_counts[low_mc_idx[idx][pointer]]
			new_mc <- c(new_mc, names(metacell_counts[low_mc_idx[idx][pointer]]))
			pointer <- pointer + 1
			if(pointer > length(idx)){
				break;
			}
		}
		clustering[which(clustering %in% as.numeric(new_mc))] <- max_cluster_id
		if(length(which(table(clustering) < thresh)) == 0){
			flag <- F
		}
	}
	return(clustering)
}
########
####### End of functions #######
################################
## run with another dataset
## catalogue the first and second pass
## optimize with 30 or different slack and see if it makes sense to do the second pass
mc_quality_info <- list()
mc_outlier_info <- list()
RNA_barcodes_ct <- NULL
mc_distr <- list()
ks_info <- list()
library(pheatmap)
pdf(paste0(output_file, "_mc_dendrogram.pdf"))
#pdf(paste0(output_file, "_cluster_pass_pheatmap.pdf"))
for(ct in cell_types){
  print(ct)
  if(!csv_flag){
	  if(!length(assay_hit)){
	    CT_cluster <- RNAcounts[, rownames(eval(parse(text= paste0(Rds_name, "@meta.data")))[celltypes == ct, ])]
		}
    CT_cluster <- RNAcounts[, celltypes == ct]
  }else{
    CT_cluster <- csv_data[, csv_cells[csv_cells[, 2] == ct, 1]]
  }
	original_CT_cluster <- CT_cluster
	if(umap_flag){
		CT_cluster <- t(umap_layout[celltypes == ct, ]) # I need to transform it bcuz the data gets transformed in the clara call
		original_CT_cluster <- RNAcounts[, celltypes == ct]
	}
  print(c("dim(CT_cluster):", dim(CT_cluster)))
  if(length(k_hit)){
    k <- as.integer(arg_tokens[k_hit + 1])
  }else{
    k <- floor(ncol(CT_cluster) / expected_cells)
		if(!k){
			k <- ncol(CT_cluster)
			print(paste("Setting k to", k))
		}else{
	    print(paste("Setting k to", k))
		}
  }
  if(summary_method == "kmed" || summary_method == "kmed_means"){
    print(paste("k=", k, "ncol(CT_cluster)=", ncol(CT_cluster)))
	#########################%%%%%%%%%%%%%%%%%%%^^^^^^^^^^^^^^^%%%%%%%%%%%%%%%%%%%###############
	#########################%%%%%%%%%%%%%%%%%%%^^^^^^^^^^^^^^^%%%%%%%%%%%%%%%%%%%###############
	#########################%%%%%%%%%%%%%%%%%%%^^^^^^^^^^^^^^^%%%%%%%%%%%%%%%%%%%###############
		if(k >= ifelse(is.null(ncol(original_CT_cluster)), 1, ncol(original_CT_cluster))){
			print("too small... gotta merge all cells into one metacell")
			clusters[[ct]]$clustering <- rep(1, ifelse(is.null(ncol(original_CT_cluster)), 1, ncol(original_CT_cluster)))
			cluster_data[[ct]] <- t(as.matrix(original_CT_cluster))
			## When original_CT_cluster is a vector, it loses its colname. That's why I had to add the line below to extract the correct annotation, especially for the RNA_barcodes_ct later. However, this fix will only apply to Seurat input data, and not the csv data. N.B. I'm using RNAcounts to infer the names.
			rownames(cluster_data[[ct]]) <- colnames(RNAcounts)[which(celltypes == ct)]
			RNA_barcodes_ct <- c(RNA_barcodes_ct, rownames(cluster_data[[ct]]))
      cell2metacell_info <- c(cell2metacell_info, paste(ct, clusters[[ct]]$clustering, sep= "_"))

			if(length(assay_hit)){
				RNA_metacell_umap_ct <- NULL
				for(i in unique(clusters[[ct]]$clustering)){
					data_subset <- RNA_umap[which(clusters[[ct]]$clustering == i), ];
					if(is.null(dim(data_subset))){
						RNA_metacell_umap_ct <- rbind(RNA_metacell_umap_ct, data_subset);
					}else{
						RNA_metacell_umap_ct <- rbind(RNA_metacell_umap_ct, colMeans(data_subset))
					}
				}
				#rownames(RNA_metacell_umap_ct) <- paste(ct, seq(nrow(RNA_metacell_umap_ct)), sep= "_")
				rownames(RNA_metacell_umap_ct) <- paste(ct, unique(clusters[[ct]]$clustering), sep= "_")
				print(paste("Done clustering", ct))
				RNA_metacell_umap <- rbind(RNA_metacell_umap, RNA_metacell_umap_ct)
			}
			### merge all cells into one metacell
		#}else if(length(k) >0 && k > 3 && k < ncol(CT_cluster)){
		}else if(length(k) >0 && k < ncol(CT_cluster)){
			library(parallel)
			#clust <- makeCluster(ceiling(detectCores()/2))
			pass_num <- 1
			pass_max <- 1 #3
			while(pass_num <= pass_max){
				#clara_out <- invoke_clara(CT_cluster, original_CT_cluster, iter_flag, clusters, RNA_metacell_umap, ct, k)
				clara_out <- invoke_clara_simplified(CT_cluster, k)
				mc_qual <- metacell_quality(clara_out$clusters, expected_cell_num= expected_cells, slack_ratio= .15)
				#outlier_idx <- which(as.character(clara_out$clustering) %in% names(mc_qual[which(mc_qual == F)]))
				outlier_idx <- get_outliers(clara_out$clusters, 10)
				print(c("Pass:", pass_num, "valid metacells:", length(which(mc_qual == T)), "outliers:", length(outlier_idx)))
				mc_quality_info[[ct]][[pass_num]] <- length(which(mc_qual == T))
				mc_outlier_info[[ct]][[pass_num]] <- length(outlier_idx)
				mc_distr[[ct]][[pass_num]] <- clara_out$clusters$clustering
				## Plot heatmap
				ann_col <- data.frame(metacells= factor(clara_out$clusters$clustering))
				rownames(ann_col) <- seq(ncol((CT_cluster)))
				colnames(CT_cluster) <- seq(ncol((CT_cluster)))
				#pheatmap(CT_cluster, annotation_col= ann_col, fontsize= 6, show_colnames= F, main= paste(ct, pass_num))
				###############
				cluster_data[[ct]] <- t(as.matrix(original_CT_cluster))
				if(pass_max > 1){
					if(length(outlier_idx)){
						CT_cluster <- CT_cluster[, -outlier_idx]
						cluster_data[[ct]] <- t(as.matrix(original_CT_cluster[, -outlier_idx]))
					}else{
						cluster_data[[ct]] <- t(as.matrix(original_CT_cluster))
						pass_num <- pass_max
					}
				}
				## plot the metacell count table with a dendrogram obtained from mean clusters
				cluster_means_res <- cluster_means(clara_out$clusters)
				print(paste("dim(cluster_means_res):", dim(cluster_means_res)))
				tryCatch(
					 {
						hc_res <- hclust(d = dist(x = cluster_means_res))
						df <- data.frame(table(clara_out$clusters$clustering));
						df$y <- as.factor(rep(1, nrow(df)))
						df_ordered <- df[hc_res$order, ]
						df_ordered$Var1 <- factor(df_ordered$Var1, levels= df_ordered$Var1)
						p1 <- ggplot(df_ordered, aes(y= Var1, x= y)) + geom_tile(show.legend= F, color= "gray", fill= "white") + theme_classic() + geom_text(aes(label= Freq), size= 1) + theme(axis.text = element_text(size = 3)) + theme(axis.text.y = element_text(size = rel(1), hjust = 1, angle = 0), 
							# margin: top, right, bottom, and left ->  plot.margin = unit(c(1, 0.2, 0.2, -0.5), "cm"), 
							panel.grid.minor = element_blank()) + ggtitle(ct)
						p2 <- ggdendrogram(data = as.dendrogram(hc_res), rotate= T, segments= T, leaf_labels= F) + geom_text(size= 1)
						print(plot_grid(p1, p2, rel_widths=c(.25, 1), ncol= 2))#align= "h"
					},
					error= function(e){
						message(paste("dim(cluster_means_res):", dim(cluster_means_res)))
						message(e)
					}
				)
				## update k after removing outliars
				ks_info[[ct]][[pass_num]] <- k
				k <- floor(ncol(CT_cluster) / expected_cells)
				pass_num <- pass_num + 1
			}

			## Tried parallelization of the k_search, but it still was pretty slow and we decided to do one k (the heuristic one)
			## and apply the outlier (metacells with few cells) removal, then cluster again as the second pass (k needs to be accordingly).
			#k_search <- get_k_grid(ncol(CT_cluster), 30)
			#clara_out <- mclapply(k_search, function(k){return(invoke_clara(CT_cluster, original_CT_cluster, iter_flag, clusters, RNA_metacell_umap, ct, k))}, mc.cores= min(ceiling(detectCores()/2), length(k_search)))

			#density_res <- metacell_density(clara_out, k_search)

			clusters[[ct]] <- clara_out$clusters
			if(merge_flag){
				clusters[[ct]]$clustering <- merge_small_mc(clusters[[ct]]$clustering, thresh= expected_cells)
			}
			print(table(clusters[[ct]]$clustering))
			if(length(assay_hit)){
				RNA_metacell_umap_ct <- NULL
				for(i in unique(clusters[[ct]]$clustering)){
					data_subset <- RNA_umap[which(clusters[[ct]]$clustering == i), ];
					if(is.null(dim(data_subset))){
						RNA_metacell_umap_ct <- rbind(RNA_metacell_umap_ct, data_subset);
					}else{
						RNA_metacell_umap_ct <- rbind(RNA_metacell_umap_ct, colMeans(data_subset))
					}
				}
				#rownames(RNA_metacell_umap_ct) <- paste(ct, seq(nrow(RNA_metacell_umap_ct)), sep= "_")
				rownames(RNA_metacell_umap_ct) <- paste(ct, unique(clusters[[ct]]$clustering), sep= "_")
				print(paste("Done clustering", ct))
				RNA_metacell_umap <- rbind(RNA_metacell_umap, RNA_metacell_umap_ct)
			}
			#RNA_barcodes_ct <- c(RNA_barcodes_ct, colnames(original_CT_cluster))
			RNA_barcodes_ct <- c(RNA_barcodes_ct, rownames(cluster_data[[ct]]))
			cell2metacell_info <- c(cell2metacell_info, paste(ct, clusters[[ct]]$clustering, sep= "_"))
	}	
	#########################%%%%%%%%%%%%%%%%%%%^^^^^^^^^^^^^^^%%%%%%%%%%%%%%%%%%%###############
	#########################%%%%%%%%%%%%%%%%%%%^^^^^^^^^^^^^^^%%%%%%%%%%%%%%%%%%%###############
	#########################%%%%%%%%%%%%%%%%%%%^^^^^^^^^^^^^^^%%%%%%%%%%%%%%%%%%%###############
  }else if(summary_method == "kmeans"){
    if(length(k) > 0 && k > 3){
      print(dim(CT_cluster))
      if(class(CT_cluster) == "numeric"){
        next;
      }
      clusters[[ct]] <- kmeans(t(as.matrix(CT_cluster)), k)
			cell2metacell_info <- c(cell2metacell_info, paste(ct, clusters[[ct]]$clustering, sep= "_"))
    }
  }
  else{
    error("Undefined method of summarization. Please pick either kmed, kmeans, or kmed_means!")
  }
}
dev.off()
save(ks_info, mc_distr, mc_quality_info, mc_outlier_info, cluster_data, clusters, file= paste0(output_file, "_clusters_heart_debug.Rdata"))
if(pass_max > 1){
	mc_qual <- as.data.table(sapply(seq(length(mc_quality_info)), function(i) unlist(mc_quality_info[[i]])))
	colnames(mc_qual) <- names(mc_quality_info)
	mc_qual[, pass := seq(1, 3)]

	mc_outlier <- as.data.table(sapply(seq(length(mc_outlier_info)), function(i) unlist(mc_outlier_info[[i]])))
	colnames(mc_outlier) <- names(mc_outlier_info)
	mc_outlier[, pass := seq(1, 3)]


	ks_dt <- as.data.table(sapply(seq(length(ks_info)), function(i) unlist(ks_info[[i]])))
	colnames(ks_dt) <- names(ks_info)
	ks_dt[, pass := seq(1, 3)]



	dt <- NULL;
	for(i in seq(length(mc_distr)))
		for(j in seq(length(mc_distr[[i]]))){
			dt <- rbind(dt, data.table(clusters= mc_distr[[i]][[j]], pass= as.character(j), celltype= names(mc_distr)[i]))
		}

	pdf(paste0(output_file, "_mc_qual_freq.pdf")); ggplot(melt(mc_qual, id.vars= "pass", variable.name= "cell_type", value.name= "quality"), aes(x= pass, y= quality)) + geom_point() + geom_line() + facet_wrap(~cell_type) + theme(axis.text.x= element_text(angle= 45)); ggplot(melt(mc_outlier, id.vars= "pass", variable.name= "cell_type", value.name= "outlier"), aes(x= pass, y= outlier)) + geom_point() + geom_line() + facet_wrap(~cell_type) + theme(axis.text.x= element_text(angle= 45));
	ggplot(melt(ks_dt, id.vars= "pass", variable.name= "cell_type", value.name= "k"), aes(x= pass, y= k)) + geom_point() + geom_line() + facet_wrap(~cell_type, scales= "free") + theme(axis.text.x= element_text(angle= 45));
#for(i in seq(length(unique(dt$celltype))))
#ggplot(dt) + geom_histogram(aes(x= clusters, fill= pass)) + facet_wrap(~celltype, scales= "free");
	dt$clusters <- factor(dt$clusters)
	ggplot(dt) + geom_bar(astat = "identity", position = 'dodge', aes(x= clusters, fill= pass)) + geom_hline(yintercept= 10, linetype= "dashed", color= "black")+ facet_wrap(~celltype, scales= "free") + theme(axis.text.x= element_text(angle= 45, size= 3));
	dev.off()
}

print("Done clustering!")
mat <- NULL
mat_sum <- NULL
mc_names <- NULL
for(i in seq(length(clusters))){
  if(summary_method == "kmed"){
    mat <- cbind(mat, t(clusters[[i]]$medoids))
  }else if(summary_method == "kmeans"){
    mat <- cbind(mat, t(clusters[[i]]$centers))
  }
  else{##kmed_means
    temp_cl <- NULL
    temp_cl_sum <- NULL
		if(umap_flag){
			for(clst in unique(clusters[[i]]$clustering)){
				idx <- which(clusters[[i]]$cluster == clst)
				mc_names <- c(mc_names, paste0(names(clusters)[i], "_", clst))
				if("numeric" %in% class(cluster_data[[i]][idx, ])){
					temp_cl <- rbind(temp_cl, cluster_data[[i]][idx, ])
					temp_cl_sum <- rbind(temp_cl_sum, cluster_data[[i]][idx, ])
				}else{
					temp_cl <- rbind(temp_cl, apply(cluster_data[[i]][idx, ], 2, FUN= mean))
					temp_cl_sum <- rbind(temp_cl_sum, apply(cluster_data[[i]][idx, ], 2, FUN= sum))
				}
			}
		}else{
    	for(clst in unique(clusters[[i]]$clustering)){
				idx <- which(clusters[[i]]$cluster == clst)
      	if("numeric" %in% class(clusters[[i]]$data[idx, ])){
        	temp_cl <- rbind(temp_cl, clusters[[i]]$data[idx, ])
        	temp_cl_sum <- rbind(temp_cl_sum, clusters[[i]]$data[idx, ])
      	}
      	else{
					idx <- which(clusters[[i]]$cluster == clst)
        	temp_cl <- rbind(temp_cl, apply(clusters[[i]]$data[idx, ], 2, FUN= mean))
        	temp_cl_sum <- rbind(temp_cl_sum, apply(clusters[[i]]$data[idx, ], 2, FUN= sum))
      	}
    	}
		}
		print(dim(temp_cl))
		print(dim(temp_cl_sum))
    mat <- cbind(mat, t(temp_cl))
    mat_sum <- cbind(mat_sum, t(temp_cl_sum))
  }
}

print("done making mat")

colnames(mat) <- colnames(mat_sum) <- mc_names

##############################
##############################
final_umap_res <- uwot::umap(t(mat), pca= 30, pca_center= T, scale= T, n_components= umap_dim)
rownames(final_umap_res) <- colnames(mat)
colnames(final_umap_res) <- paste0("UMAP", seq(umap_dim))

celltypes <- sapply(colnames(mat), function(i) strsplit(i, "_")[[1]][1])
df <- data.frame(UMAP1= final_umap_res[, 1], UMAP2= final_umap_res[, 2], celltype= celltypes)
print("creating directory!")
dir.create(output_file)
library(ggplot2)
pdf(paste0(output_file, "/umap_", summary_method, ".pdf"))
print(ggplot(df, aes(x= UMAP1, y= UMAP2)) + geom_point(aes(color= celltype)) + theme_classic() + geom_text(aes(label= celltype),hjust=0, vjust=0, size= 3, check_overlap = T))
dev.off()
##############################
##############################

if(length(assay_hit)){
	#kk <- 5L
	#knn_res <- class::knn(train= RNA_metacell_umap, test= ATAC_umap, cl= rownames(RNA_metacell_umap), k= kk)
	knn_res <- my_knn(train= RNA_metacell_umap, test= ATAC_umap, threshold= threshold)
	atac2metacell_info <- data.frame(barcode= rownames(ATAC_umap), metacell= knn_res)
	write.csv(atac2metacell_info, paste0(output_file, "/ATAC_cell2metacell_info_", summary_method, ".csv"), row.names= F)
}

rna2metacell_info <- data.frame(barcode= RNA_barcodes_ct, metacell= cell2metacell_info)
write.csv(rna2metacell_info, paste0(output_file, "/RNA_cell2metacell_info_", summary_method, ".csv"), row.names= F)
write.csv(mat, paste0(output_file, "/cellSummarized_", summary_method, ".csv"))
write.csv(mat_sum, paste0(output_file, "/cellSummarized_", summary_method, "_sum.csv"))
write.csv(final_umap_res, paste0(output_file, "/RNA_metacell_umap_", summary_method, ".csv"))
if(length(assay_hit)){
	save(atac2metacell_info, ATACcounts, clusters, RNA_metacell_umap, ATAC_umap, mc_names, file= paste0(output_file, "/", summary_method, "_clustered.RData"))

	uniq_mc <- unique(atac2metacell_info$metacell)
	atac_metacell <- NULL;
	atac_metacell_sum <- NULL;
	for(i in seq(length(uniq_mc))){
		hits <- atac2metacell_info$barcode[which(atac2metacell_info$metacell == uniq_mc[i])];
		if(length(hits) > 1){
			atac_metacell <- cbind(atac_metacell, Matrix::rowMeans(ATACcounts[, hits]))
			atac_metacell_sum <- cbind(atac_metacell_sum, Matrix::rowSums(ATACcounts[, hits]))
		}else{
			atac_metacell <- cbind(atac_metacell, ATACcounts[, hits])
			atac_metacell_sum <- cbind(atac_metacell_sum, ATACcounts[, hits])
		}
	}
	colnames(atac_metacell) <- uniq_mc
	colnames(atac_metacell_sum) <- uniq_mc
	write.csv(atac_metacell, paste0(output_file, "/cellSummarized_ATAC_", summary_method, ".csv"))
	write.csv(atac_metacell_sum, paste0(output_file, "/cellSummarized_ATAC_", summary_method, "_sum.csv"))
}
print("Done!")
Rtnse_plot <- F
if(Rtnse_plot){
  library(Rtsne)
  Rtsne_whole_res <- Rtsne(as.matrix(RNAcounts), check_duplicates= F)
  pdf(paste0(output_file, "/tSNE_", summary_method, "_Rtsne.pdf"))
  plot(Rtsne_whole_res$Y)
  dev.off()
}


if(length(assay_hit)){
#################################################
############# BEGIN VISUALIZATION ###############

	library(ggplot2)
	addSmallLegend <- function(myPlot, pointSize = 0.5, textSize = 3, spaceLegend = 0.1) {
		myPlot +
			guides(shape = guide_legend(override.aes = list(size = pointSize)),
						 color = guide_legend(override.aes = list(size = pointSize))) +
	theme(legend.title = element_text(size = textSize), legend.text  = element_text(size = textSize), legend.key.size = unit(spaceLegend, "lines"))
	}
####
###
##
	ATAC_info <- atac2metacell_info
	RNA_info <- rna2metacell_info
	df <- data.frame(umap1= Sub[[reduction]]@cell.embeddings[, 1], umap2= Sub[[reduction]]@cell.embeddings[, 2], assays= assays)
	df$Seurat_celltype <- original_celltypes
#df$Seurat_celltype <- as.character(Sub@meta.data$seurat_clusters)
	df <- cbind(df[c(ATAC_info$barcode, RNA_info$barcode), ], rbind(ATAC_info, RNA_info))
	celltypes <- sapply(as.character(df$metacell), function(i) strsplit(i, "_")[[1]][1])
	uniq_celltypes <- unique(celltypes)
	df$celltype <- celltypes

	pdf(paste0(output_file, "/plot_umap_mc_", expected_cells, ".pdf"));
	print(ggplot(df, aes(x= umap1, y= umap2, shape= assays, colour= Seurat_celltype)) + geom_point(alpha= .6) + theme_classic() + geom_text(aes(label= Seurat_celltype),hjust=0, vjust=0, size= 3, check_overlap = T));
	for(i in seq(length(unique(celltypes)))) {
		myPlot <- ggplot(subset(df, celltype == uniq_celltypes[i]), aes(x= umap1, y= umap2, shape= assays)) + geom_point(aes(colour= metacell)) + ggtitle(uniq_celltypes[i]) + theme_classic(); print(addSmallLegend(myPlot))
	};
	dev.off()

###############
	df2 <- data.frame(umap1= Sub[[reduction]]@cell.embeddings[, 1], umap2= Sub[[reduction]]@cell.embeddings[, 2], assays= assays, Seurat_celltype= original_celltypes)
	pdf(paste0("plot_umap_seurat_", expected_cells, ".pdf")); for(i in seq(length(unique(celltypes)))) {print(ggplot(subset(df2, Seurat_celltype == uniq_celltypes[i]), aes(x= umap1, y= umap2, shape= Seurat_celltype, colour= assays)) + geom_point(alpha= .6) + ggtitle(uniq_celltypes[i]) + theme_classic())}; dev.off()
###############

	rna_cnt <- table(RNA_info$metacell)
	atac_cnt <- table(ATAC_info$metacell)

	df_main <- data.frame(RNA_mc_cnt= as.numeric(rna_cnt), row.names= names(rna_cnt))
	df_main$ATAC_mc_cnt <- 0
	df_main[names(atac_cnt), ]$ATAC_mc_cnt <- as.numeric(atac_cnt)
	df_main$celltype <- sapply(rownames(df_main), function(i) strsplit(i, "_")[[1]][1])
	print(paste("expected_cells=", expected_cells))
	pdf(paste0(output_file, "/MC_scatterplots_", expected_cells, ".pdf"))
	print(ggplot(df_main, aes(x= RNA_mc_cnt, y= ATAC_mc_cnt)) + geom_point(colour= "blue") + theme_classic() + ggtitle("all cell types"))
	for(ct in unique(df_main$celltype))
		print(ggplot(subset(df_main, celltype == ct), aes(x= RNA_mc_cnt, y= ATAC_mc_cnt)) + geom_point(colour= "blue") + theme_classic() + ggtitle(ct))

	library(data.table)
	dt_main <- as.data.table(df_main)
	dt <- melt(dt_main, id.vars= "celltype", variable.name= "assay", value.name= "mc_count")


	print(ggplot(dt, aes(x= celltype, y= mc_count)) + geom_boxplot(aes(fill= assay)) + geom_hline(yintercept= expected_cells, linetype="dashed", color = "green") + theme_classic() + geom_text(aes(-1, expected_cells, label = "exp. num. cells", hjust = 0, vjust= 1), color= "green") + theme(axis.text.x = element_text(angle = 45)))
#######
	rna_rc <- mat_sum
	atac_rc <- atac_metacell_sum 
	rna_rc_sum <- colSums(rna_rc)
	atac_rc_sum <- colSums(atac_rc)
#######

	rc_dt <- data.table(read_counts= c(rna_rc_sum, atac_rc_sum), assay= c(rep("RNA", length(rna_rc_sum)), rep("ATAC", length(atac_rc_sum))))
	rc_dt_log <- rc_dt
	rc_dt_log$read_counts <- log2(1 + rc_dt_log$read_counts)

	print(ggplot(rc_dt_log, aes(x= assay, y= read_counts)) + geom_boxplot() + theme_classic() + ylab("log2(1 + sum of read counts per metacell)"))

	dev.off()

	if(F){## I have to add the UMAP of ATAC counts and then make the plots... right now it's only for the RNAs.
		pdf(paste0(output_file, "/MC_umaps_", expected_cells, ".pdf"))
		umap_df <- data.frame(UMAP1= umap_res[, 1], UMAP2= umap_res[, 2])
		umap_df_reordered <- umap_df[RNA_info$barcode,]
		umap_df_reordered$MC <- RNA_info$metacell
		ggplot(umap_df_reordered, aes(x= UMAP1, y= UMAP2)) + geom_point(aes(colour= MC), alpha= .5) + theme_classic() + theme(legend.position='none')
		umap_celltypes <- sapply(umap_df_reordered$MC, function(i) strsplit(i, "_")[[1]][1])
		for(ct in unique(umap_celltypes)){
			hits <- which(umap_celltypes == ct)
			print(ggplot(umap_df_reordered[hits, ], aes(x= UMAP1, y= UMAP2)) + geom_point(aes(colour= MC), alpha= .5) + theme_classic() + theme(legend.position='none') + ggtitle(ct))
		}
		dev.off()
	}
}
