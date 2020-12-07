# Set up
library(ggplot2)
library(corrr)
library(purrr)
library(tidyr)
library(vegan)
library(chemometrics)
library(BiodiversityR)
library(dtree)
library(rpart.plot)
library(rpart.utils)

# Import the data
Uranium <- read.delim("Data Final.txt", header =T)

##############################
#####EXPLORATORY ANALYSIS#####
##############################

# Plot parameter density distributions and incldue units for all labels
var.labs <- list("Al" = "Al (g/kg)", "BD"= expression("BD (g/cm"^3*")"), "Carbonate"="Carbonate (%)", 
                 "Clay"="Clay (%)", "Corg"="Corg (%)", "Fe"="Fe (g/kg)", 
                 "FineSand"="FineSand (%)","MedSand"="MedSand (%)", "Mn"="Mn (g/kg)",  
                 "pH"= "pH", "Silt"="Silt (%)", "Turbidity"="Turbidity (NTU)","U"="U (ppm)")
labels <- function(variable,value){return(var.labs[value])}

Uranium[7:length(Uranium)] %>% keep(is.numeric) %>% gather() %>% 
ggplot(aes(value)) +
  facet_wrap(~ key, scales = "free", labeller = labels)+
  geom_histogram()+
  theme_classic()+
  theme(text=element_text(size=14))

# Plot Uranium contents distribution along locations and slopes
  # Rename Base1 and Base2 as Base to plot them together as a single slope position at each location
Uranium$Zone <- as.character(Uranium$Zone)
Uranium[Uranium == "BBase1"] <- "BBase"
Uranium[Uranium == "BBase2"] <- "BBase"
Uranium[Uranium == "EBase1"] <- "EBase"
Uranium[Uranium == "EBase2"] <- "EBase"
Uranium[Uranium == "IBase1"] <- "IBase"
Uranium[Uranium == "IBase2"] <- "IBase"
Uranium$Zone <- as.factor(Uranium$Zone)

# Final Plot
data.frame(Zone = Uranium$Zone, Position = Uranium$Position, U = as.numeric(Uranium$U)) %>%
ggplot(aes(x=Zone, y=U, fill=Position)) + 
  geom_violin(trim = F)+
  geom_boxplot(width=0.1, fill = "white")+
  geom_vline(xintercept = 3.4, linetype = "longdash")+
  geom_vline(xintercept = 6.5, linetype = "longdash")+
  labs(x = "Location and Slope Position", y = "U (ppm)")+
  theme_classic()+
  theme(text=element_text(size=14))

# Re-import original data frame to separate Base as Base1 and Base2
Uranium <- read.delim("Data Final.txt", header =T)

# Tukey HSD - Compare U-content among all sites (location + slope position)
aovU<-summary(Uaov<-aov(Uranium$U~Uranium$Zone))
U1<-cbind("Df"=Uaov[[1]][[1]],"SumSq"=Uaov[[1]][[2]],
          "MeanSq"=Uaov[[1]][[3]],"Fvalue"=Uaov[[1]][[4]],
          "pValue"=Uaov[[1]][[5]])
Utuk<-TukeyHSD(Uaov)
Utuk<-data.frame(Utuk[1][1])
colnames(Utuk)<-c("diff","lwr","upr","padj")
View(Utuk)
which(Utuk$padj<0.05) 

###################
##Correlation map##
###################

# Create a correlation matrix and plot the correlation map
cor(Uranium [,c(6:length(Uranium))], method = "spearman") %>%
network_plot(min_cor = 0.5,repel = T)
  
# Scatter plots for U against variables (examples)
plot(Uranium$U~Uranium$Carbonate, xlab = "Carbonate (%)", ylab = "Uranium (ppm)")
plot(Uranium$U~Uranium$pH, xlab = "pH", ylab = "Uranium (ppm)")
plot(Uranium$U~Uranium$Al, xlab = "Aluminum (g/kg)", ylab = "Uranium (ppm)")
plot(Uranium$U~Uranium$Fe, xlab = "Iron (g/kg)", ylab = "")

############
####PCA#####
############

### Check requirements
## Multivariate normality
# Plot a QQ-plot
chisqplot.multi <- function(m, main="QQ plot", ylab=expression(paste(chi^2, " Quantile")),
                                                                     xlab=expression("Ordered Mahalanobis D"^2*"")){
  # n x p numeric matrix
  x <- as.matrix(m)
  # centroid
  center <- colMeans(x) 
  n <- ncol(x)
  cov <- cov(x)
  # distances 
  d <- mahalanobis(x,center,cov) 
  s <- sort(d, index=TRUE)
  q <- (0.5:length(d))/length(d)
  par(las=1, cex=1.2)
  plot(s$x, qchisq(q,df=n), main=main, xlab=expression("Ordered Mahalanobis D"^2*""), ylab= ylab)
  abline(a=0,b=1)
}
chisqplot.multi(Uranium[6:length(Uranium)])

## Look for outliers
## Maybe use the methodology described here: https://towardsdatascience.com/multivariate-outlier-detection-in-high-dimensional-spectral-data-45878fd0ccb8
# Plot score distances and othogonal distances
par(mfrow=c(1,2),cex=2)
pcaDiagplot(Uranium[6:length(Uranium)], 
            princomp(Uranium[6:length(Uranium)],cor = TRUE),
            a =2, quantile=0.975)
#### Prepare data for PCA
# Select variables and add vectors for color and symbol shape for plotting purposes
legend <- c("BTop", "BMiddle", "BBase1", "BBase2","ETop", "EMiddle", "EBase1", "EBase2",
            "ITop", "IMiddle", "IBase1", "IBase2") #legend text for "Location"
symbol <- c(17,19,14,12,17,19,14,12,17,19,14,12) # symbols for "Location"
color  <- c("orange","orange","orange","orange","dodgerblue","dodgerblue","dodgerblue","dodgerblue", 
            "green", "green","green","green") # color for "Positions"
Uranium_PCA <- data.frame(Uranium[c(3,6:length(Uranium))],
                          color = color[match(Uranium$Zone,legend)],
                          symbol = symbol[match(Uranium$Zone,legend)])

### Perform PCA
dat.pca <- rda(Uranium_PCA[2:15], scale = T)   
spec <- summary(dat.pca)$species

### Plot biplot
par(mfrow=c(1,1),cex=1.2)
biplot(dat.pca,display="species",col="black",axes=TRUE,
       xlab= paste0("PC1 (",round(summary(dat.pca)$cont$importance[2,1]*100,1)," %)"), 
       ylab= paste0("PC2 (",round(summary(dat.pca)$cont$importance[2,2]*100,1)," %)"),
       ylim = c(-1.6,1.6),
       cex=3)
points(dat.pca,"sites", pch=Uranium_PCA$symbol, col=as.character(Uranium_PCA$color),cex=1)
legend("topright",legend=legend,
       pch=symbol,col=color,cex=0.6)

#### Determining the number of components with significant effects
# Apply Kaiser-Guttman criterion: select axes where eigenvalues are larger than the mean eigenvalue
dat.pca$CA$eig[dat.pca$CA$eig > mean(dat.pca$CA$eig)]
# Apply screeplot
screeplot(dat.pca, type="lines")
# Apply broken stick criterion
PCAsignificance(dat.pca)

########################
### CLUSTER ANALYSIS ###
########################
# Simple hierarchical clustering
# Scale data
Uranium_Clust <- as.data.frame(scale(Uranium[7:length(Uranium)]))
rownames(Uranium_Clust) <- paste0(Uranium[[3]],"-",c(1:57))

# Building of the divisive clusters (average linkage)
dist_mat <- dist(Uranium_Clust, method = 'euclidean')
hclust <- hclust(dist_mat, "average")
# Plot the dendrogram 
factoextra::fviz_dend(hclust, k = 3,
                      type = "phylogenic",
                      phylo_layout = "layout.gem", 
                      repel = T,
                      sub = "", cex = 1.0, lwd = 0.7,
                      ggtheme = theme_void(base_size = 0) 
)

# Plot Euclidean distances against Cophenetic distances
plot(cophenetic(hclust), dist_mat, xlab = "Cophenetic Distance", ylab= "Euclidean Distance")
abline(0,1)
lines(lowess(cophenetic(hclust), dist_mat), col = "green")

# Correlation of matrices (distance and cluster)
cor(cophenetic(hclust), dist_mat)

# Compute Stress1 (analogous to NMDS)
sqrt(sum((dist_mat - cophenetic(hclust))^2) / sum(dist_mat^2) )

##################################
#####     ELASTIC NET       ######
##################################
### Elastic net function (using gaussian family)
elastic_net <- function(X, seed, output_column, number_iterations){
  #X <- Uranium[6:length(Uranium)]
  #seed <- 0
  #output_column <- length(Uranium[6:length(Uranium)])
  #number_iterations <- 1 
  library(tidyverse)  
  library(caret)
  library(glmnet)
  library(robustHD)
  library(e1071)
  library(tidyr)
  Results <- list()
  Selected_masses <- list()
  Stand_DATA <- standardize(X)
  f <- paste(paste(names(X)[output_column]),"~.")
  BOOT <- trainControl(method = "boot", number = 50)
  for (i in 1:number_iterations){
    #Data preparation
    set.seed(seed+i)
    Sel <- createDataPartition(Stand_DATA[[output_column]], p = 0.8, list = F)
    Training <- Stand_DATA[Sel,]
    Training[[output_column]] <- as.numeric(Training[[output_column]])
    Test <- Stand_DATA[-Sel,]
    #Build the model - tunelength: numbers of different values for alpha and beta to be tested
    set.seed(seed+i)
    model <- train(as.formula(f), 
                   data = Training, 
                   method = "glmnet", 
                   family = "gaussian", 
                   trControl = BOOT, 
                   tuneLength = 5)
    #Best tuning parameters
    Masses <- data.frame(coef(model$finalModel,
                              model$bestTune$lambda)@Dimnames[[1]][coef(model$finalModel,model$bestTune$lambda)@i+1],
                         coef(model$finalModel,model$bestTune$lambda)@x
    )
    Selected_masses[[i]] <- coef(model$finalModel,
                                 model$bestTune$lambda)@Dimnames[[1]][coef(model$finalModel,model$bestTune$lambda)@i+1]
    #prediction for test
    x.test <- model.matrix(as.formula(f), Test)
    predictions <- as.numeric(as.character(predict(model,x.test)))
    Results[[i]] <-  list(Masses, 
                          data.frame(RMSE = RMSE(predictions, Test[[output_column]]), 
                                     Rsquare = cor(predictions, Test[[output_column]])^2
                          )
    )
  }
  names(Results) <- c(1:i)
  Results[[i+1]] <- Selected_masses
  return(Results)
}

### Models building
iterations <- 50
UElN_net <- elastic_net(Uranium[6:length(Uranium)],
                        0,
                        length(Uranium[6:length(Uranium)]),
                        iterations)

### Plot results R2 and variables used in the models
Graph <- function(X){
  #R2
  C <- c()
  for (i in 1:(length(X)-1)){C <- c(C,X[[i]][[2]][[2]])}
  G1 <- ggplot()+ 
    geom_boxplot(aes(C))+
    coord_flip()+
    labs(x = expression("R"^2*""), y = "")+
    theme_classic()+
    theme(axis.text.x = element_blank(), axis.ticks.x=element_blank(),text=element_text(size=15))
  
  #Regression coefficients
  A <- data.frame()
  Variables <- names(Uranium[6:(length(Uranium)-1)])
  for (i in 1:(length(X)-1)){
    A <- rbind(A,
               data.frame(var = X[[i]][[1]][[1]], coef = X[[i]][[1]][[2]])
    )
    if (any(!Variables %in% X[[i]][[1]][[1]] == TRUE)){
      A <- rbind(A,
                 data.frame(var = Variables[!Variables %in% X[[i]][[1]][[1]]], coef = 0)
      )
    }
  }
  A <- A[-which(A[[1]] == "(Intercept)"),]
  G2 <- ggplot(A, aes(x=var, y=coef))+ 
    geom_boxplot(outlier.colour="white")+
    #geom_jitter(color="red", size=0.8, alpha=0.9)+
    labs(x = "", y = "Regression coefficients")+
    theme_classic()+
    theme(axis.text.x= element_text(angle = 45, hjust = 1), text=element_text(size=15))
  #Combine graphs
  gridExtra::grid.arrange(G1, G2, ncol=2, heights = c(5,1))
  }
Graph(UElN_net)

#######################################
### Partial least square regression ###
#######################################
PLS <- function(X, seed, output_column, number_iterations){
  library(pls)
  Results <- list()
  for (i in 1:number_iterations){
    # Preparing training and test sets
    set.seed(seed+i)
    Sel <- X[[output_column]] %>% createDataPartition(p = 0.8, list = FALSE)
    Train <- X[Sel,]
    Test <- X[-Sel,]
    # Building the model
    set.seed(seed+i)
    model <- train( U~., data = Train, method = "pls",
                    scale = TRUE,
                    trControl = trainControl("cv", number = 20),
                    tuneLength = 10)
    # Test the model
    Predictions <- model %>% predict(Test)
    # Save model data and statistics
    Results[[i]] <- list(nb_comp = model$bestTune$ncomp,
                         RMSE = caret::RMSE(Predictions, Test[[output_column]]),
                         Rsquare = caret::R2(Predictions, Test[[output_column]])
    )
  }
  names(Results) <- c(1:i)
  return(Results)
}
Plot_pls <- function(X){
  C <- c()
  D <- c()
  for (i in 1:length(X)){
    C <- c(C,X[[i]][[3]])
    D <- c(D,X[[i]][[1]])
  }
  boxplot(C, ylab = expression("R"^2*""))
  hist(D, xlab =  "Number of components", main = "", 
       breaks = seq(from = min(D)-1, to = max(D), by = 1))
}

### Models building
PLS_results <- PLS(Uranium[6:length(Uranium)],
                   0,
                   length(Uranium[6:length(Uranium)]),
                   50)

### Plot results R2 and number of components used in the models
par(mfrow=c(1,2))
Plot_pls(PLS_results)

#####################
### Decision Tree ###
#####################
Arbre <- dtree(U~., data = Uranium[6:length(Uranium)], methods = c("bump"), bump.rep = 100)
summary(Arbre)

#Inspect the bumps 
Arbre$bump.list
Arbre$bump.matrix
which.min(Arbre$bump.matrix[,4])
#Chose a particular bump (use the one with the lowest RMSE)
par(mfrow=c(1,1))
rpart.plot(Arbre$bump.list[[which.min(Arbre$bump.matrix[,4])]])
