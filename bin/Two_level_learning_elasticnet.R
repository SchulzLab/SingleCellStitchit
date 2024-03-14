args <- commandArgs(TRUE)
library("methods")
library("glmnet")
library("doMC")
library("parallel")
if(length(args) < 1) {
	args <- c("--help")
}
 
## Help 
if("--help" %in% args) {
	cat("
	INVOKE offers linear regression with Lasso, Ridge, and Elastic Net regularisation.
	Arguments:
	--outDir Output directory (will be created if it does not exist)
	--dataDir Directory containing the data
	--response Name of the response variable
	--cores Number of cores to be use (default 1)
	--fixedAlpha Use a fixed value for the alpha parameter in elastic net regulatisation, do not perform a grid search
	--alpha Stepsize to optimise the alpha parameter in elastic net regularisation (default 0.05)
	--testsize Size of test data[%] (default 0.2)
	--regularisation L for Lasso, R for Ridge, and E for Elastic net (default E)
	--innerCV Number of folds for inner cross-validation (default 6)
	--outerCV Number of iterations of outer cross-validation to determine test error (default 3)
	--constraint Specifies a constraint on the coefficent sign, enter N for negative and P for positive constraint
	--performance Flag indiciating whether the performance of the model should be assessed (default TRUE)
	--seed Random seed used for random number generation (default random)
	--leaveOneOutCV Flag indicating whether a leave one out cross-validation should be used (default FALSE)
	--asRData Store feature coefficients as RData files (default FALSE)
	--randomise Randomise the feature matrix (default FALSE) 
	--logResponse Flag indicating whether the response variable should be log transformed (default TRUE)
	--ftest Flag indicating whether partial F-test should be computed to assess the significance of each feature (default FALSE)
	--coefP p-value threshold for model coefficient (default 1, all OLS coefs will be returned)
	--help=print this text
")
	q(save="no")
}

# Process command arguments
parseArgs <- function(x) strsplit(sub("^--", "", x), "=")
argsDF <- as.data.frame(do.call("rbind", parseArgs(args)))
print(argsDF)
argsL <- as.list(as.character(argsDF$V2))
print(argsL)
names(argsL) <- argsDF$V1
print(names(argsL))

if(is.null(argsL$outDir)) {
	cat("No output directory specified. Use the --outDir option to specify an output directory.")
	q(save="no")
}
argsL$outDir<-paste0(argsL$outDir,"/")

if(is.null(argsL$dataDir)) {
	cat("No data directory specified. Use the --dataDir option to specify a data directory.")
	q(save="no")
}
argsL$dataDir<-paste0(argsL$dataDir,"/")
Data_Directory <- argsL$dataDir

if(is.null(argsL$response)) {
	cat("No response variable name specified. Use the --response option to specify a response variable.")
	q(save="no")
}

if(is.null(argsL$ftest)){
	argsL$ftest<-FALSE
}

if(is.null(argsL$testsize)){
	argsL$testsize <- 0.2
}

if(is.null(argsL$innerCV)){
	argsL$innerCV <- 6
}


if(is.null(argsL$outerCV)){
	argsL$outerCV<- 3
}

if (as.numeric(argsL$outerCV) < 2){
	cat("Number of outer cross validation folds must be at least 2")
	q(save="no")
}

if(is.null(argsL$alpha)) {
	argsL$alpha <- 0.05
}

if(is.null(argsL$cores)) {
	argsL$cores <- 1
}

if(is.null(argsL$coefP)) {
	argsL$coefP <- 1
}

if(is.null(argsL$regularisation)){
	argsL$regularisation<-c("E")
}

if(is.null(argsL$constraint)){
	lower_bound <- NULL
	upper_bound <- NULL
}else if(argsL$constraint=="P"){
	lower_bound <- 0
}else if(argsL$constraint=="N"){
	upper_bound <- 0
}

if(is.null(argsL$fixedAlpha)){
	argsL$fixedAlpha <- -1
}

if (is.null(argsL$performance)){
	argsL$performance <- TRUE
}

if (! is.null(argsL$seed)){
	set.seed(as.numeric(argsL$seed))
}

if(is.null(argsL$leaveOneOutCV)){
     argsL$leaveOneOutCV <- FALSE
}

if (is.null(argsL$asRData)){
	argsL$asRData <- FALSE
}

if (is.null(argsL$randomise)){
	argsL$randomise <- FALSE
}

if (is.null(argsL$logResponse)){
	argsL$logResponse <- TRUE
}

if ((argsL$leaveOneOutCV==TRUE) & (argsL$ftest==TRUE)){
     cat("The F-Test can not be combined with leave one out cross validation.")
     q(save="no")
	}

registerDoMC(cores = argsL$cores)

permute<-function(x,resPos){
s<-sample(length(x))
s<-s[which(s != resPos)]
c(x[s],x[resPos])
}



#Check output directory, create it if necessary
dir.create(argsL$outDir,showWarning=FALSE)

#lrumpf
#error log file
error_log_file <- paste0(argsL$outDir,"failed_genes.log")


#Initilaise lists for storage of intermediate results
FileList<-list.files(path=Data_Directory)
numFiles=length(FileList)
pearson_correlation<-vector("list",numFiles)
spearman_correlation<-vector("list",numFiles)
test_error<-vector("list",numFiles)
rss_error<-vector("list",numFiles)
ftest_result<-vector("list",numFiles)
coefficients<-vector("list",numFiles)
coefficientsF<-vector("list",numFiles)
Sample_View<-vector("list",numFiles)
validSamples<-vector("logical",numFiles)
spearmanPassed<-vector("numeric",numFiles)

#lrumpf
#store test/training partition in inner CV
#innerCV_partition <- vector("list",numFiles)
#store training/test partition in outer CV 
outerCV_partition <- vector("list",numFiles)

#Print sample names
print(paste("Total number of samples:",as.character(length(FileList)),sep=" "))
if (length(FileList)==0){
	print("No samples available! Aborting")
	exit()
}
counter<-0
for(Sample in FileList){
	counter<-counter+1
	print(paste(as.character(counter),unlist(unlist(strsplit(Sample,".txt")))))
	}

#Loop through sample files
i<-0
for(Sample in FileList){
	i<-i+1
	print(paste("Learning sample ",as.character(i),sep=" "))	
	#Loading and preprocessing data
	print("Processing sample matrix. This can take a few minutes. Please wait.")
	print(paste(Data_Directory,Sample,sep=""))
	M<-read.table(paste(Data_Directory,Sample,sep=""),header=TRUE,sep="",row.names=1)
	
	M<-unique(M)
	M<-data.frame(M)
	FeatureNames_temp<-colnames(M)
	Response_Variable_location_temp <- grep(argsL$response,FeatureNames_temp)
	if (min(M[,Response_Variable_location_temp]) >= 0){
		if (argsL$logResponse == TRUE){
			M<-log2(M+1)
		}
	}else{
		print("Applying log2 transformation only to features")
		Response_Variable_location_temp <- grep(argsL$response,FeatureNames_temp)
		M[,-Response_Variable_location_temp]<-log2(M[,-Response_Variable_location_temp]+1)
	}
	SD<-apply(M,2,sd)
	Feature_zero_SD<-as.vector(which(SD==0))
	if(length(Feature_zero_SD)>0){
		print("Warning, there are constant features. These are not considered for further analysis.")
		if (Response_Variable_location_temp %in% Feature_zero_SD){
			print("Warning, response is constant, this sample is excluded")
			validSamples[i]=FALSE
			next;
			}
		M<-M[,-c(Feature_zero_SD)]
	}
	if (is.null(dim(M))){
		validSamples[i]=FALSE
		spearmanPassed[i]=1
		print("Warning, sample matrix is null, this sample is excluded")
		next;
	}
	if (dim(M)[2] < 2){
		print("Warning, no data included")
		validSamples[i]=FALSE
		spearmanPassed[i]=1
		next;
	}
	#print(length(which(M==0)))
	if (length(which(M==0))>(dim(M)[1]*dim(M)[2]*0.8)){
		validSamples[i]=FALSE
		spearmanPassed[i]=1
		print("Warning, insufficient data coverage")
		next;
	}

	FeatureNames<-colnames(M)
	#lrumpf: scaling after partition into train and test set
	#M<-data.frame(scale(M,center=TRUE, scale=TRUE))
     if (dim(M)[1] < 25){
          print("Warning, less then 30 samples available. This file is not processed.")
          validSamples[i]=FALSE
		spearmanPassed[i]=1
          next;
     }else{
          validSamples[i]=TRUE
     }

	if (argsL$leaveOneOutCV==FALSE){
		vectorLength=as.numeric(argsL$outerCV)
		pearson_correlation[[i]]<-vector("list",vectorLength)
		spearman_correlation[[i]]<-vector("list",vectorLength)
		test_error[[i]]<-vector("list",)
		coefficients[[i]]<-vector("list",vectorLength)
		#lrumpf
		outerCV_partition[[i]]<-vector("list",vectorLength) 
	}else{
		vectorLength<-nrow(M)    
		pearson_correlation[[i]]<-cbind(vector("list",vectorLength),vector("list",vectorLength))
		spearman_correlation[[i]]<-cbind(vector("list",vectorLength),vector("list",vectorLength))
		test_error[[i]]<-vector("list",vectorLength)
		coefficients[[i]]<-vector("list",vectorLength)
		#lrumpf
		outerCV_partition[[i]]<-vector("list",vectorLength)
	}
	if (argsL$ftest ==TRUE){
		ftest_result[[i]]<-vector("list",dim(M)[2]-1)
	}

	name<-unlist(unlist(strsplit(Sample, ".txt")))
	Response_Variable_location<- grep(argsL$response,FeatureNames)
	#Randomise the data
	if (argsL$randomise == TRUE){
		MP<-t(apply(M,1,permute,Response_Variable_location))
		colnames(MP)<-colnames(M)
		M<-data.frame(scale(MP,center=TRUE, scale=TRUE))  
	}

	if (argsL$performance == TRUE){
		if (argsL$leaveOneOutCV==FALSE){
			predictedAll<-c()
			measuredAll<-c()
			
		#lrumpf	
		#store best model for each fold in outer CV		
		best_models <- vector("list",vectorLength)
		#inner_fold_ids <- vector("list",as.numeric(argsL$outerCV))
		
		#error handling
		#error.count <- 0
		
		#Looping through the outer fold
		for (k in 1:argsL$outerCV){
			print(paste("Outer cross validation fold: ",as.character(k),sep=" "))
			
			# Partition data into test and training data sets
			Test_size<-round(nrow(M)/(1/as.numeric(argsL$testsize)))
			#TODO: save partitions for best model
			rndselect<-sample(x=nrow(M), size=Test_size)
			Test_Data<-M[rndselect,]
			Train_Data<-M[-rndselect,]
			
			#lrumpf: scale training data
			Train_Data<-scale(Train_Data,center=TRUE, scale=TRUE)
			#scaling parameters for test data 
			training_means <- attr(Train_Data,"scaled:center")
			training_sds <- attr(Train_Data,"scaled:scale")
			Train_Data<-data.frame(Train_Data)
			
			#TODO:
			#print(ncol(Train_Data))
			#print(any(is.na(df)))
		  #print(Train_Data)
	
			#lrumpf: scale test data with training parameters
			Test_Data <- data.frame(scale(Test_Data, center = training_means, scale = training_sds))

			# Split the features from response
			x_train<-as.matrix(Train_Data[,-Response_Variable_location])
			x_test<-as.matrix(Test_Data[,-Response_Variable_location])
			y_train<-as.vector(unlist(Train_Data[,Response_Variable_location,drop=FALSE]))
			y_test<-as.vector(unlist(Test_Data[,Response_Variable_location]))
			
			#Creating alpha vector
		 	A<-c()
			if(argsL$regularisation=="L"){
				alphaslist <- c(1.0)
				print("The value of alpha is set to 1.0 (Lasso penalty)")
			}else{
				if(argsL$regularisation=="R"){
				alphaslist <- c(0.0)
				print("The value of alpha is set to 0.0 (Ridge penalty)")
				}else{
					alphaslist<-seq(0,1,by=as.numeric(argsL$alpha))
				}
			}
				
			#Learning model on training data
		 	#TODO: keep=TRUE to keep test/train partitions (fold ids for each observation)
		 	if(argsL$regularisation=="E"){   
				if(argsL$fixedAlpha==-1){
					if(is.null(argsL$constraint)){
					  elasticnet<-mclapply(alphaslist, function(x){;cv.glmnet(x_train, y_train,alpha=x,nfolds=as.numeric(argsL$innerCV),keep=TRUE)}, mc.cores=argsL$cores)
            
					 # cv_glmnet_with_fit <- function(alpha_value) {
					 # model <- cv.glmnet(x_train, y_train, alpha = alpha_value, nfolds = as.numeric(argsL$innerCV))
					 #   return(model$glmnet.fit) # Extract glmnet.fit object
					 # }
					 # elasticnet <- mclapply(alphaslist, cv_glmnet_with_fit, mc.cores = argsL$cores)  
										
						#print(elasticnet[[]])
						#print(elasticnet$foldid)
						#v <- elasticnet$fit.preval
						#print(v)
					}else{ 
						if(argsL$constraint=="P"){
							elasticnet<-mclapply(alphaslist, function(x){cv.glmnet(x_train, y_train,alpha=x,lower=0,nfolds=as.numeric(argsL$innerCV),keep=TRUE)}, mc.cores=argsL$cores)
							}else{
								if(argsL$constraint=="N"){
								elasticnet<-mclapply(alphaslist, function(x){cv.glmnet(x_train, y_train,alpha=x,upper=0,nfolds=as.numeric(argsL$innerCV),keep=TRUE)}, mc.cores=argsL$cores)
								}
								
							}
			      		}
				}else{
				x=argsL$fixedAlpha
				if(is.null(argsL$constraint)){
					elasticnet<-cv.glmnet(x_train, y_train,alpha=x,nfolds=as.numeric(argsL$innerCV),keep=TRUE, parallel=TRUE)
				}else{ 
					if(argsL$constraint=="P"){
							elasticnet<-cv.glmnet(x_train, y_train,alpha=x,lower=0,nfolds=as.numeric(argsL$innerCV),keep=TRUE, parallel=TRUE)
					}else{
						if(argsL$constraint=="N"){
							elasticnet<-cv.glmnet(x_train, y_train,alpha=x,upper=0,nfolds=as.numeric(argsL$innerCV),keep=TRUE, parallel=TRUE)
								}
						}
		      	}
				}   
			}else{
			x=alphaslist[1]
			if(is.null(argsL$constraint)){
					elasticnet<-cv.glmnet(x_train, y_train,alpha=x,nfolds=as.numeric(argsL$innerCV),keep=TRUE,parallel=TRUE)
					}else{ 
					if(argsL$constraint=="P"){
						elasticnet<-cv.glmnet(x_train, y_train,alpha=x,lower=0,nfolds=as.numeric(argsL$innerCV),keep=TRUE,parallel=TRUE)
					}else{
						if(argsL$constraint=="N"){
							elasticnet<-cv.glmnet(x_train, y_train,alpha=x,upper=0,nfolds=as.numeric(argsL$innerCV),keep=TRUE, parallel=TRUE)
							}
						}
	     			}
				}

		 #constant.response <- FALSE
		 #tryCatch( { standardization <- log("$ operator is invalid for atomic vectors"); print(standardization) },
		 #error = function(e) {constant.response <- TRUE})
			elasticnet.fail <- FALSE
			if(length(elasticnet[[1]]) > 1){ #catches missing data in input for glmnet?
				if (argsL$regularisation=="E"){
						if(argsL$fixedAlpha==-1){
								for (j in 1:length(alphaslist)) {
								#print(is.atomic(elasticnet[[j]]))								  #print(elasticnet[[j]]$foldid)
								#v <- elasticnet[[j]]$fit.preval
								#print(v)
								#print(elasticnet[[j]])
								
								#print(elasticnet[[j]])
								#A[j]<-min(elasticnet[[j]]$cvm)
								#print(alphaslist[j])
								#tryCatch(A[j]<-min(elasticnet[[j]]$cvm),
								#error = function(e) {break})
								#if(constant.response) break;

								#if (constant.response){
									#print(elasticnet[[1]])
									#print(length(elasticnet[[1]]))
									#elasticnet[[1]] <- elasticnet[[-c(1:length(elasticnet[[1]]))]]

									#Can never be best alpha value
								#	A[j] <- 100000000
								
								#}

								#ERROR HANDLING: lrumpf
							  elasticnet.catch <- tryCatch( A[j]<-min(elasticnet[[j]]$cvm),error=function(e) {
							    msg <- conditionMessage(e)
							    err.msg <- attr(elasticnet[[j]],"condition")
							    #print(msg)
							    #if(msg == "$ operator is invalid for atomic vectors"){
					        #print exact error message from elastic net model
							    print(err.msg); 
							    write(toString(Sample), error_log_file, append=TRUE);
							    write(toString(err.msg), error_log_file, append=TRUE);
							    elasticnet.fail <<-TRUE;
							    #print(e); constant.response <<-TRUE
							    #}
							    })
							  #break inner CV evaluation
								if(elasticnet.fail) {
								  #set gene as not valid (according to Florians error handling)
								  #length(elasticnet[[1]]) = 0
								  #error.count <- error.count+1
								  #print(error.count)
									break;	
								}
							  #print("elastic net")
							  #print(j)  
							  #print(elasticnet[[j]]$cvm)
								}
					    
						  #break outer CV: allow one failed inner CV
						  #if (elasticnet.fail && error.count < 2 ){
						  #  next
						  #}
							#if(elasticnet.fail && error.count == 2) {
							#  validSamples[i]=FALSE
							#  break
							#}
						  
						  if(elasticnet.fail) {
						    validSamples[i]=FALSE
						    break
						  }
						  
						  
						  #print(alphaslist)
						  #print(A)
							#Determine best alpha value from training data
							#lrumpf: multiple minimum error values possible -> choose first index
							index<-which(A==min(A), arr.ind=TRUE)
							if (length(index > 1)) {index <- index[1]}
							#TODO: store model (best lambda param obtained from inner CV )
							model<-elasticnet[[index]]
					}
				}else{
				  #TODO: store model
					model<-elasticnet
					}
			}

			#lrumpf: debug
			#print(alphaslist)
			#print(A)
			#print(elasticnet)
			#print(model$foldid)
			if (length(elasticnet[[1]]) > 1){ 
			  
				######## lrumpf ########
				#store best model (alpha, lambda and coefficients) -> all lambda values 
				#store training and test partition outer for best alpha and for inner folds (best lambda) 
				#cross-validation error inner and outer folds (performance)
				#print(model$foldid)
			  
				#for each outer fold k store elastic net models for min. or fixed alpha
				#store glmnet object for reconstructing the best model (lambda.min) (glmnet object needed for SHAP values)
				#write to per gene file without collecting results for all genes 
			  
				inner_fold_ids <- data.frame(
				  "sample_id" = rownames(Train_Data),
				  "fold_id" = model$foldid
				)
				
				min.alpha <- alphaslist[[index]]
				#min.alpha <- model$call$alpha
				best_models[[k]] <- list("model" = model, "alpha" = min.alpha, 
				                         "lambda" = model$lambda.min, "inner_folds" = inner_fold_ids) #model$glmnet.fit 
			
				#print(best_models[[k]])
				#samples used for validation in outer CV-fold k for gene i 
				#TODO: in performance overview validation fold index for each sample 
				#test_data
				outerCV_partition[[i]][[k]] <- data.frame("sample_ids"=rownames(Test_Data),"outer_folds"= rep(k,nrow(Test_Data)))
				
				#store test/training partition inner folds (lambda-tuning)
				#get fold id for each observation -> keep=TRUE in cv.glmnet() function
				#TODO: add sample id 
				#innerCV_partition[[i]][k] <- model$foldid
				########################
        
				#Determine error of the best alpha model on hold out data and on training data (outer CV)
				predict_fit<-predict(model, x_test, s="lambda.min")
				predict_fit_train<-predict(model, x_train, s="lambda.min")
				coefficients[[i]][[k]]<-coef(model, s = "lambda.min")
				
				#TODO: adapt to new model performance format  
				if (var(predict_fit,y_test) > 0){
					pearson_correlation[[i]][k]<-cor(predict_fit,y_test)
					spearman_correlation[[i]][k]<-cor(predict_fit,y_test,method='spearman')
					predictedAll<-c(predictedAll,predict_fit)
					measuredAll<-c(measuredAll,y_test)
				}else{
					pearson_correlation[[i]][k]<-0.0
					spearman_correlation[[i]][k]<-0.0
				}
				#performance: 
				test_error[[i]][k]<-sum((y_test-predict_fit)^2)/length(y_test)
				rss_error[k]<-sum((y_test-predict_fit)^2)
			}else{
			  #print error message
			  err.msg <- elasticnet[[1]]
			  print(err.msg) #missing input -> why??
			  write(toString(Sample), error_log_file, append=TRUE);
			  write(toString(err.msg), error_log_file, append=TRUE);
				#coefficients[[i]][k]<-c()
				#pearson_correlation[[i]][k]<-0.0
				#spearman_correlation[[i]][k]<-0.0
				#test_error[[i]][k]<-1.0
				#rss_error[k]<-1.0
				validSamples[i]=FALSE
				break
				#spearmanPassed[i]=1
			}
		}
	################# outer CV end ###################

		#lrumpf
		#catch sparse data: skip sample (gene) when inner CV results twice in constant predictors 
		if (!validSamples[i]){
		  coefficients[[i]][k]<-c()
		  pearson_correlation[[i]][k]<-0.0
		  spearman_correlation[[i]][k]<-0.0
		  test_error[[i]][k]<-1.0
		  rss_error[k]<-1.0
		  spearmanPassed[i]=1
		  print("Skip sample...")
		  next
		}
			
		if (length(predictedAll) > 0){
			if (var(predictedAll,measuredAll) > 0){
				spearmanP<-cor.test(predictedAll,measuredAll,method='spearman')$p.value
				spearmanPassed[i]=spearmanP
				}else{
				spearmanPassed[i]=1
				}
			}
		else{
			spearmanPassed[i]=1
		}
		if (argsL$ftest==TRUE){
			numFeatures=dim(M)[2]
			if (numFeatures-1 > 1){
				featureMatrixF<-c()
		          for (j in 1:length(coefficients[[i]])){
		               if (length(coefficients[[i]][[j]]>1)){
		                    featureMatrixF<-rbind(featureMatrixF,coefficients[[i]][[j]][,1])
		               }
		          }
				featureMatrixF<-featureMatrixF[,-1]
		          if (length(featureMatrixF > 1)){
						ftest_result[[i]]<-rep(1.0,numFeatures-1)
		        meanFeature<-apply(featureMatrixF,2,median)
						nonZeroFeatures<-which(meanFeature != 0)
						print(nonZeroFeatures)
						print(paste0("Running F-test using ",length(nonZeroFeatures)," features"))
						featureIndex<-nonZeroFeatures
						for (m in featureIndex){
							print(paste0("Excluding feature ",m," ",FeatureNames[m]))
							MF<-M[,-m]
							error<-c(1:argsL$outerCV)
							for (k in 1:argsL$outerCV){
									Test_size<-round(nrow(MF)/(1/as.numeric(argsL$testsize)))
									rndselect<-sample(x=nrow(MF), size=Test_size)
									Test_Data<-MF[rndselect,]
									Train_Data<-MF[-rndselect,]
									Response_Variable_location_MF<- grep(argsL$response,colnames(MF))
									x_train<-as.matrix(Train_Data[,-Response_Variable_location_MF])
									x_test<-as.matrix(Test_Data[,-Response_Variable_location_MF])
									y_train<-as.vector(unlist(Train_Data[,Response_Variable_location_MF,drop=FALSE]))
									y_test<-as.vector(unlist(Test_Data[,Response_Variable_location_MF]))
							 		A<-c()
									if(argsL$regularisation=="L"){
										alphaslist <- c(1.0)
									}else{
										if(argsL$regularisation=="R"){
											alphaslist <- c(0.0)
										}else{
											alphaslist<-seq(0,1,by=as.numeric(argsL$alpha))
										}
									}
								if(argsL$regularisation=="E"){   
									if(argsL$fixedAlpha==-1){
										if(is.null(argsL$constraint)){
											elasticnet<-mclapply(alphaslist, function(x){cv.glmnet(x_train, y_train,alpha=x,nfolds=as.numeric(argsL$innerCV),keep=TRUE)}, mc.cores=argsL$cores)
										}else{ 
											if(argsL$constraint=="P"){
												elasticnet<-mclapply(alphaslist, function(x){cv.glmnet(x_train, y_train,alpha=x,lower=0,nfolds=as.numeric(argsL$innerCV),keep=TRUE)}, mc.cores=argsL$cores)
											}else{
												if(argsL$constraint=="N"){
													elasticnet<-mclapply(alphaslist, function(x){cv.glmnet(x_train, y_train,alpha=x,upper=0,nfolds=as.numeric(argsL$innerCV),keep=TRUE)}, mc.cores=argsL$cores)
												}
											}
							      		}
									}else{
										x=argsL$fixedAlpha
										if(is.null(argsL$constraint)){
											elasticnet<-cv.glmnet(x_train, y_train,alpha=x,nfolds=as.numeric(argsL$innerCV),keep=TRUE,parallel=TRUE)
										}else{ 
										if(argsL$constraint=="P"){
											elasticnet<-cv.glmnet(x_train, y_train,alpha=x,lower=0,nfolds=as.numeric(argsL$innerCV),keep=TRUE,parallel=TRUE)
										}else{
											if(argsL$constraint=="N"){
												elasticnet<-cv.glmnet(x_train, y_train,alpha=x,upper=0,nfolds=as.numeric(argsL$innerCV),keep=TRUE,parallel=TRUE)
													}
												}
								     	 	}
										}			   
								}else{
								x=alphaslist[1]
								if(is.null(argsL$constraint)){
									elasticnet<-cv.glmnet(x_train, y_train,alpha=x,nfolds=as.numeric(argsL$innerCV),keep=TRUE,parallel=TRUE)
								}else{ 
									if(argsL$constraint=="P"){
										elasticnet<-cv.glmnet(x_train, y_train,alpha=x,lower=0,nfolds=as.numeric(argsL$innerCV),keep=TRUE,parallel=TRUE)
									}else{
										if(argsL$constraint=="N"){
											elasticnet<-cv.glmnet(x_train, y_train,alpha=x,upper=0,nfolds=as.numeric(argsL$innerCV),keep=TRUE,parallel=TRUE)
											}
										}
					     			}
								}
							elasticnet.fail <- FALSE 		
							if(length(elasticnet[[1]]) > 1){
								if (argsL$regularisation=="E"){
									if(argsL$fixedAlpha==-1){
										for (j in 1:length(alphaslist)) {
										    A[j]<-min(elasticnet[[j]]$cvm)
										    #elasticnet.catch <- tryCatch( A[j]<-min(elasticnet[[j]]$cvm),error=function(e) {
										      # msg <- conditionMessage(e)
										      # print(attr(elasticnet[[j]],"condition")); 
										      # elasticnet.fail <<-TRUE;
										      # })
										    #break inner CV evaluation
										    #if(elasticnet.fail) {
										      #  break;	
										    #}
										}
									  #if(elasticnet.fail) {
									  # validSamples[i]=FALSE
									  # break
									  #}
									index<-which(A==min(A), arr.ind=TRUE)
									if (length(index > 1)) {index <- index[1]}
									model<-elasticnet[[index]]
									}
								}else{
									model<-elasticnet
									}
								}
							if (length(elasticnet[[1]]) > 1){ 
								predict_fit<-predict(model, x_test, s="lambda.min")
								predict_fit_train<-predict(model, x_train, s="lambda.min")
								error[k]<-sum((y_test-predict_fit)^2)
							}else{
								error[k]<-1.0
							}
						}
						mRSS<-mean(error)
						mORSS<-mean(unlist(rss_error[[i]]),na.rm=TRUE)
						MSe<-mORSS/(dim(MF)[1]-(length(nonZeroFeatures)))
						partialRSS<-mRSS-mORSS
						fvalue<-partialRSS/MSe
						pValue<-1.0-pf(as.numeric(fvalue),length(nonZeroFeatures)-1,(dim(MF)[1]-length(nonZeroFeatures)))
						ftest_result[[i]][m]<-pValue
						}
					}
				###Generate output
			}else{
				print(paste0("Computing significance for feature ",FeatureNames[1]))
				mORSS<-mean(unlist(rss_error[[i]]),na.rm=TRUE)
				MSe<-mORSS/(dim(MF)[1]-2)
				fvalue<-mORSS/MSe
				pValue<-1.0-pf(as.numeric(fvalue),1,dim(MF)[1]-2)
				ftest_result[[i]][1]<-pValue
			}
	}
	}else{
	#Leave one out cross validation
	for (k in 1:nrow(M)){                                                                                                                                                                                                                                                                                                                      
		#print(paste("Leave one out cross validation fold: ",as.character(k),sep=" "))
		# Partition data into test and training data sets
		Test_Data<-M[k,]
		Train_Data<-M[-k,]
		
		#lrumpf: scale training data
		Train_Data<-scale(Train_Data,center=TRUE, scale=TRUE)
		#scaling parameters for test data 
		training_means <- attr(Train_Data,"scaled:center")
		training_sds <- attr(Train_Data,"scaled:scale")
		Train_Data<-data.frame(Train_Data)
		
		#lrumpf: scale test data with training parameters
		Test_Data <- data.frame(scale(Test_Data, center = training_means, scale = training_sds))
		
		# Split the features from response
		x_train<-as.matrix(Train_Data[,-Response_Variable_location])
		x_test<-as.matrix(Test_Data[,-Response_Variable_location])
		y_train<-as.vector(unlist(Train_Data[,Response_Variable_location,drop=FALSE]))
		y_test<-as.vector(unlist(Test_Data[,Response_Variable_location]))
		#Creating alpha vector
		A<-c()
		if(argsL$regularisation=="L"){
			print("The value of alpha is set to 1.0 (Lasso penalty)")
			alphaslist <- c(1.0)
		}else{
		if(argsL$regularisation=="R"){
			print("The value of alpha is set to 0.0 (Ridge penalty)")
			alphaslist <- c(0.0)
			}else{
				alphaslist<-seq(0,1,by=as.numeric(argsL$alpha))
			}
		}
	
		#Learning model on training data
		if(argsL$regularisation=="E"){   
			if(argsL$fixedAlpha==-1){
		if(is.null(argsL$constraint)){
			elasticnet<-mclapply(alphaslist, function(x){cv.glmnet(x_train, y_train,alpha=x,nfolds=as.numeric(argsL$innerCV),keep=TRUE)}, mc.cores=argsL$cores)
		}else{ 
			if(argsL$constraint=="P"){
			elasticnet<-mclapply(alphaslist, function(x){cv.glmnet(x_train, y_train,alpha=x,lower=0,nfolds=as.numeric(argsL$innerCV),keep=TRUE)}, mc.cores=argsL$cores)
			}else{
				if(argsL$constraint=="N"){
					elasticnet<-mclapply(alphaslist, function(x){cv.glmnet(x_train, y_train,alpha=x,upper=0,nfolds=as.numeric(argsL$innerCV),keep=TRUE)}, mc.cores=argsL$cores)
					}
				}
			}
		}else{
			x=argsL$fixedAlpha
			if(is.null(argsL$constraint)){
				elasticnet<-cv.glmnet(x_train, y_train,alpha=x,nfolds=as.numeric(argsL$innerCV),keep=TRUE, parallel=TRUE)
				}else{ 
					if(argsL$constraint=="P"){
						elasticnet<-cv.glmnet(x_train, y_train,alpha=x,lower=0,nfolds=as.numeric(argsL$innerCV),keep=TRUE, parallel=TRUE)
					}else{
						if(argsL$constraint=="N"){
							elasticnet<-cv.glmnet(x_train, y_train,alpha=x,upper=0,nfolds=as.numeric(argsL$innerCV),keep=TRUE, parallel=TRUE)
						}
					}
				}
			}   
	}else{
	x=alphaslist[1]
	if(is.null(argsL$constraint)){
		elasticnet<-cv.glmnet(x_train, y_train,alpha=x,nfolds=as.numeric(argsL$innerCV),keep=TRUE,parallel=TRUE)
		}else{ 
			if(argsL$constraint=="P"){
				elasticnet<-cv.glmnet(x_train, y_train,alpha=x,lower=0,nfolds=as.numeric(argsL$innerCV),keep=TRUE,parallel=TRUE)
				}else{
				if(argsL$constraint=="N"){
					elasticnet<-cv.glmnet(x_train, y_train,alpha=x,upper=0,nfolds=as.numeric(argsL$innerCV),keep=TRUE,parallel=TRUE)
					}
				}
			}
		}
	elasticnet.fail <- FALSE
	if(length(elasticnet[[1]]) > 1){
		if (argsL$regularisation=="E"){
			if(argsL$fixedAlpha==-1){
				for (j in 1:length(alphaslist)) {
				  #ERROR HANDLING: lrumpf
				  elasticnet.catch <- tryCatch( A[j]<-min(elasticnet[[j]]$cvm),error=function(e) {
				    #msg <- conditionMessage(e)
				    err.msg <- attr(elasticnet[[j]],"condition");
				    print(err.msg)
				    write(toString(Sample), error_log_file, append=TRUE);
				    write(toString(err.msg), error_log_file, append=TRUE);
				    elasticnet.fail <<-TRUE;
				  })
				  #break inner CV evaluation
				  if(elasticnet.fail) {
				    break;	
				  }

				}
			  if(elasticnet.fail) {
			    validSamples[i]=FALSE
			    break
			  }
				#Determine best alpha value from training data
				index<-which(A==min(A), arr.ind=TRUE)
				if (length(index > 1)) {index <- index[1]}
				model<-elasticnet[[index]]
				}
			}else{
			model<-elasticnet
		}
	}  
	if (length(elasticnet[[1]]) > 1){ 
	
	  #lrumpf
	  #store elastic net model
	  inner_fold_ids <- data.frame(
	    "sample_id" = rownames(Train_Data),
	    "fold_id" = model$foldid
	  )
	  
	  min.alpha <- alphaslist[[index]]
	  best_models[[k]] <- list("model" = model, "alpha" = min.alpha, 
	                           "lambda" = lambda.min, "inner_folds" = inner_fold_ids) #model$glmnet.fit 
	  
	  #samples used for validation in outer CV-fold k for gene i 
	  #TODO: in performance overview validation fold index for each sample 
	  #test_data
	  outerCV_partition[[i]][[k]] <- list("sample_id" = rownames(Test_Data), "outer_fold"= k)
	  
	  #Determine error of the best alpha model on hold out data and on training data
		predict_fit<-predict(model, x_test, s="lambda.min")
		predict_fit_train<-predict(model, x_train, s="lambda.min")
		coefficients[[i]][[k]]<-coef(model, s = "lambda.min")
		pearson_correlation[[i]][k,1]<-y_test
		pearson_correlation[[i]][k,2]<-predict_fit
		spearman_correlation[[i]][k,1]<-y_test
		spearman_correlation[[i]][k,2]<-predict_fit
		test_error[[i]][k]<-(y_test-predict_fit)^2

	}else{
	  err.msg <- elasticnet[[1]]
	  print(err.msg) #missing input -> why??
	  write(toString(Sample), error_log_file, append=TRUE);
	  write(toString(err.msg), error_log_file, append=TRUE);
	  validSamples[i]=FALSE
	  break
		#coefficients[[i]][k]<-c()
		#pearson_correlation[[i]][k]<-0.0
		#spearman_correlation[[i]][k]<-0.0
		#test_error[[i]][k]<-1.0
		#validSamples[i]<-FALSE
		}
	}		  
	  ### outer CV Leave-One-Out ###
	  if (!validSamples[i]){
	    coefficients[[i]][k]<-c()
	    pearson_correlation[[i]][k]<-0.0
	    spearman_correlation[[i]][k]<-0.0
	    test_error[[i]][k]<-1.0
	    print("Skip sample...")
	    next
	  }
		}
	}

	#Learning the model once on the full data set
	if (! is.null(argsL$seed)){
		set.seed(as.numeric(argsL$seed))
		}
	print(paste0("Learning OLS model on reduced feature space for sample ",i))
	######################
	#Determine nonzero model coefficients
	#modelCoefMatrix<-c()
	#for (j in 1:length(coefficients[[i]])){
	#	if (length(coefficients[[i]][[j]]>1)){
	#		modelCoefMatrix<-rbind(modelCoefMatrix,coefficients[[i]][[j]][,1])
	#		}
	#	}
	# if (length(modelCoefMatrix) != 0){
	# medianModelCoefMatrix<-apply(modelCoefMatrix,2,median)[-1]
	# nObs<-dim(M)[1]
	#Partition data into test and training data sets
	#if (length(which(medianModelCoefMatrix!=0))){
	#	if (length(which(medianModelCoefMatrix!=0))>=nObs){
	#		ols_Data<-M[,c(order(abs(medianModelCoefMatrix),decreasing=T)[1:(nObs-2)],Response_Variable_location)]
	#	}else{
	#		ols_Data<-M[,c(which(medianModelCoefMatrix!=0),Response_Variable_location)]
  # } 
	 #######################
	  #lrumpf
	  #change to non-zero coefficients from best model -> pick best model with minimum cv-error across outer folds
	  #print(unlist(test_error[[i]]))
	  min.err <- min(unlist(test_error[[i]]))
	  min.cv.err.ind <- which(unlist(test_error[[i]]) == min.err)
	  if (length(min.cv.err.ind) > 1){
	    min.cv.err.ind <- min.cv.err.ind[1]
	  }
	  #print(min.cv.err.ind)
	  #print(best_models[[min.cv.err.ind]])
	  #print(best_models[[min.cv.err.ind]]$model)
	  best.outer.model <- best_models[[min.cv.err.ind]][[1]]
	  #print(best.outer.model)
	  model.coefficients <- coef(best.outer.model, s = "lambda.min")[,1]
	  #print(model.coefficients)
	  #extract non-zero coefficients without intercept
	  #model.coefficients[which(model.coefficients!=0)][-1] 
	  ##### save best elastic net model #####
	  #list "model", "alpha", "lambda", "inner_folds", "outer_fold_ids"
	  elasticnet_model <- c(best_models[[min.cv.err.ind]],outerCV_partition[[i]][[min.cv.err.ind]])
	  save(elasticnet_model, file = paste0(argsL$outDir,"Elasticnet_Regression_Model_",unlist(unlist(strsplit(Sample,".txt"))),".RData"))
	  #######################################
    if (length(model.coefficients) != 0){
      nObs<-dim(M)[1]
      print(nObs)
      print(length(model.coefficients))
      #Partition data into test and training data sets (--> OLS on complete dataset)
      #if number of features is higher than number of observations choose coefficients with highest absolute value
        if (length(model.coefficients)>=nObs){
          #why -2? intercept
      		ols_Data<-M[,c(order(abs(model.coefficients!=0),decreasing=T)[1:(nObs-2)],Response_Variable_location)]
      		#print(ols_Data)
      	}else{
      		ols_Data<-M[,c(which(model.coefficients!=0),Response_Variable_location)]
      		#print(ols_Data)
      	} 
	 ########################
    #OLS  
		model<-lm(Expression~.,ols_Data)	
		model.coefs<-summary(model)$coefficients[,c(1,4)]
		signif.coefs<-which(model.coefs[,2]<=as.numeric(argsL$coefP))
		model.coefs.signif<-model.coefs[signif.coefs,]
		if (length(signif.coefs > 0)){
			for (j in  1:length(row.names(model.coefs.signif))){
				row.names(model.coefs.signif)[j]<-gsub(".","\t",row.names(model.coefs.signif)[j],fixed=T)
			}
			if (length(signif.coefs)>1){
				write.table(model.coefs.signif,file=paste0(argsL$outDir,"Selected_Regions_",unlist(unlist(strsplit(Sample,".txt"))),".bed"),quote=F,row.names=T,col.names=F,sep="\t")
			}else{
			cat(paste0(gsub(".","\t",row.names(model.coefs)[signif.coefs],fixed=T),"\t",model.coefs.signif[1],"\t",model.coefs.signif[2],"\n"),file=paste0(argsL$outDir,"Selected_Regions_",unlist(unlist(strsplit(Sample,".txt"))),".bed"))
				}
			}
		}
	
}
	
###############################
###Writing model performance###
###############################
#TODO: adapt -> store elastic net model
if (argsL$performance == TRUE){
	if (argsL$leaveOneOutCV == FALSE){
	     for (i in 1:length(FileList)){
	          if (validSamples[i]==FALSE){
     	          next;
	          }
	          cm<-mean(unlist(pearson_correlation[[i]]),na.rm=TRUE)
	          csd<-var(unlist(pearson_correlation[[i]]),na.rm=TRUE)
	          cms<-mean(unlist(spearman_correlation[[i]]),na.rm=TRUE)
     	     csds<-var(unlist(spearman_correlation[[i]]),na.rm=TRUE)
	          erm<-mean(unlist(test_error[[i]]),na.rm=TRUE)
	          ersd<-var(unlist(test_error[[i]]),na.rm=TRUE)
	          Sample_View[[i]]<-data.frame(Sample_Name=FileList[i],Pearson=cm,PearsonVar=csd,Spearman=cms,SpearmanVar=csds,MSE=erm,MSEVar=ersd,pVal=spearmanPassed[i],qVal=p.adjust(spearmanPassed,method="BY")[i])
	     }
	}else{
	     for (i in 1:length(FileList)){
	          if (validSamples[i]==FALSE){
	               next;
     	     }
	          cm<-cor(unlist(pearson_correlation[[i]][,1]),unlist(pearson_correlation[[i]][,2]))
	          csd<-0.0
	          cms<-cor(unlist(spearman_correlation[[i]][,1]),unlist(spearman_correlation[[i]][,2]),method="spearman")
	          csds<-0.0
	          erm<-mean(unlist(test_error[[i]]),na.rm=TRUE)
	          ersd<-var(unlist(test_error[[i]]),na.rm=TRUE)
	          Sample_View[[i]]<-data.frame(Sample_Name=FileList[i],Pearson=cm,PearsonVar=csd,Spearman=cms,SpearmanVar=csds,MSE=erm,MSEVar=ersd)
			write.table(cbind(pearson_correlation[[i]][,1],pearson_correlation[[i]][,2]),paste0(argsL$outDir,"Leave_One_Out_Predicitions_",FileList[i]),quote=FALSE,sep="\t",row.names=F)
     	}
	}
	Sample_ViewF<-do.call("rbind",Sample_View)
	if (argsL$asRData==TRUE){
		save(coefficients,file=paste(argsL$outDir,"Coefficients.RData",sep=""))
	}
	write.table(Sample_ViewF,paste(argsL$outDir,"Performance_Overview.txt",sep=""),quote=FALSE,sep="\t",row.names=F)
}

