source("https://bioconductor.org/biocLite.R")
biocLite()
library(DNAshapeR)
library(caret)

path <- "D:/GitHub/BISC577/gcPBM/"
fa.files <- paste0(path, dir(path, pattern="\\.txt.fa$"))
txt.files <- paste0(path, dir(path, pattern="\\.txt$"))
seq.names <- substr(basename(txt.files), 1, 3)

all.1mer.r2 <- vector("numeric")
all.1mer.shape.r2 <- vector("numeric")
for (i in 1:length(fa.files)) {
	## Question 4
	## Predict DNA shapes
	pred <- getShape(fa.files[i])
	exp.data <- read.table(txt.files[i])
	label <- seq.names[i]
	print(i)
	print(label)
	
	## Encode feature vectors
	featureType1 <- c("1-mer")
	featureType2 <- c("1-mer", "1-shape")
	featureVector1 <- encodeSeqShape(fa.files[i], pred, featureType1)
	featureVector2 <- encodeSeqShape(fa.files[i], pred, featureType2)
	
	## Build MLR model
	trainControl <- trainControl(method="cv", number=10, savePredictions=TRUE)
	df.1mer <- data.frame(affinity=exp.data$V2, featureVector1)
	df.1mer.shape <- data.frame(affinity=exp.data$V2, featureVector2)
	model.1mer <- train(affinity ~ ., data=df.1mer, trControl=trainControl, 
						 method="glmnet", 
						 tuneGrid=data.frame(alpha=0, lambda=c(2^c(-15:15))))
	model.1mer.shape <- train(affinity ~ ., data=df.1mer.shape, 
							   trControl=trainControl,  method="glmnet",
							   tuneGrid=data.frame(alpha=0, 
												   lambda=c(2^c(-15:15))))
	
	## Print R^2 values
	model.1mer.r2 <- mean(model.1mer$results$Rsquared, na.rm=TRUE)
	model.1mer.shape.r2 <- mean(model.1mer.shape$results$Rsquared, na.rm=TRUE)
	all.1mer.r2 <- c(all.1mer.r2, model.1mer.r2)
	all.1mer.shape.r2 <- c(all.1mer.shape.r2, model.1mer.shape.r2)
	print(paste(label, "1-mer:", model.1mer.r2))
	print(paste(label, "1-mer+shape:", model.1mer.shape.r2))	
}

