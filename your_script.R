library(DNAshapeR)
library(caret)
library(ggplot2)
library(grid)

path <- "D:/GitHub/BISC577/gcPBM/"
fa.files <- paste0(path, dir(path, pattern="\\.txt.fa$"))
txt.files <- paste0(path, dir(path, pattern="\\.txt$"))
seq.names <- substr(basename(txt.files), 1, 3)

all.1mer.r2 <- vector("numeric")
all.1mer.shape.r2 <- vector("numeric")
trainControl <- trainControl(method="cv", number=10, savePredictions=TRUE)
for (i in 1:length(fa.files)) {
    # Question 4
    # Predict DNA shapes
    pred <- getShape(fa.files[i])
    exp.data <- read.table(txt.files[i])
    label <- seq.names[i]
    
    # Encode feature vectors
    featureVector1 <- encodeSeqShape(fa.files[i], pred, "1-mer")
    featureVector2 <- encodeSeqShape(fa.files[i], pred, c("1-mer", "1-shape"))
    
    # Build MLR model
    df.1mer <- data.frame(affinity=exp.data$V2, featureVector1)
    df.1mer.shape <- data.frame(affinity=exp.data$V2, featureVector2)
    model.1mer <- train(affinity ~ ., data=df.1mer, trControl=trainControl, 
                         method="glmnet", 
                         tuneGrid=data.frame(alpha=0, lambda=c(2^c(-15:15))))
    model.1mer.shape <- train(affinity ~ ., data=df.1mer.shape, 
                               trControl=trainControl,  method="glmnet",
                               tuneGrid=data.frame(alpha=0, 
                                                   lambda=c(2^c(-15:15))))
    
    # Print R^2 values
    model.1mer.r2 <- mean(model.1mer$results$Rsquared, na.rm=TRUE)
    model.1mer.shape.r2 <- mean(model.1mer.shape$results$Rsquared, na.rm=TRUE)
    print(paste(label, "1-mer:", model.1mer.r2))
    print(paste(label, "1-mer+shape:", model.1mer.shape.r2))    
    all.1mer.r2 <- c(all.1mer.r2, max(model.1mer$results$Rsquared, na.rm=TRUE))
    all.1mer.shape.r2 <- c(all.1mer.shape.r2, 
                           max(model.1mer.shape$results$Rsquared, na.rm=TRUE))
}

# Question 5
# Plot R^2 for 1mer vs. 1mer+shape
r2 <- data.frame(x=all.1mer.r2, y=all.1mer.shape.r2, type=seq.names)
my.theme <- theme(
    plot.margin=unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
    axis.text=element_text(colour="black", size=12),
    axis.title.x=element_text(colour="black", size=12),
    axis.title.y=element_text(colour="black", size=12),
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),
    panel.background=element_blank(),
    axis.line=element_line(colour = "black"),
    axis.ticks=element_line(colour = "black"),
    legend.position="right",
    legend.background=element_blank(),
    legend.title=element_blank(),
    legend.key=element_blank()
)
p <- ggplot(r2)
p + geom_point(aes(x=x, y=y, color=type)) + 
    geom_abline(slope=1) + geom_vline(xintercept=0) + geom_hline(yintercept=0) + 
    coord_fixed(ratio = 1, xlim = c(0,1), ylim = c(0,1)) +
    scale_x_continuous(expand = c(0, 0)) + 
    scale_y_continuous(expand = c(0, 0)) + 
    xlab(expression(1-mer~R^{2})) + ylab(expression(1-mer+shape~R^{2})) + 
    my.theme

# Question 7
path <- "D:/GitHub/BISC577/CTCF/"
fa.files <- paste0(path, dir(path, pattern="\\_500.fa$"))
for (i in 1:length(fa.files)) {
}