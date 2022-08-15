library(SFSI)

pheno <- read.table("Pheno1.txt", header = FALSE) ####### phenotypic file with 
b<-pheno[,1]
Y<- pheno[,2]
X <- pheno[,c(3:7,9:10)]
#X <-pheno[,3]
geno<- read.table("geno1.txt",header=FALSE)
c<-rownames(geno)

b<-pheno[,1]
b<-as.data.frame(b)
colnames(b)<-"V1"
geno1<- merge(b,geno,by='V1')
c<-geno1[,1]




geno1 <- data.frame(geno1[,-1], row.names = geno1[,1])
Adjusted_genotypes= snpQC(gen= geno1, MAF=0.05, impute=TRUE, remove = TRUE)   
geno1<- as.matrix(Adjusted_genotypes)
c<-rownames(geno1)



## generating the genetic relatedness matrix 
K = NAM::GRM(IMP(geno1),T) ##

##############
y <- as.vector(Y)

# Save file
save(y,K,X,file="prepared_data.RData")

##tatility in variance components###################################################################################################
library(SFSI)

# Load data 
load("prepared_data.RData")
evdK <- eigen(K)                 # Decomposition of K

# Fit model
y <- as.vector(scale(y))
fm0 <- fitBLUP(y,U=evdK$vectors,d=evdK$values,BLUP=FALSE) 
c(fm0$varU,fm0$varE,fm0$h2)

save(fm0,evdK,file="varComps.RData")
load("varComps.RData")

###################
load("prepared_data.RData") # load data

#---------- parameters ------------#
nPart <- 5            # Number of TRN-TST partitions to perform
 # Time-points to analyze
#----------------------------------#
nTST <- ceiling(0.3*length(y))    # Number of elements in TST set
partitions <- matrix(1,nrow=length(y),ncol=nPart)   # Matrix to store partitions
seeds <- round(seq(1E3, .Machine$integer.max/10, length=nPart))

for(k in 1:nPart)
{   set.seed(seeds[k])
  partitions[sample(1:length(y), nTST),k] <- 2
}
save(partitions,nPart, file="partitions.RData")  # Save partitions

#Genetic Covariances estimation

load("prepared_data.RData"); load("partitions.RData")

gencov <- phencov <- c()   # Matrices to store covariances
  for(k in 1:nPart)
  { cat("  partition = ",k,"of",nPart,"\n")
    indexTRN <- which(partitions[,k]==1)
    
    # Training set
    xTRN <- scale(X[indexTRN,])
    yTRN <- as.vector(scale(y[indexTRN]))
    KTRN <- K[indexTRN,indexTRN]   # Relationship matrix (given by replicates)
    
    # Get genetic variances and covariances
    fm <- getGenCov(y1=yTRN,y2=xTRN,K=KTRN,scale=FALSE,warn=FALSE,mc.cores=1)
    
    gencov <- cbind(gencov,fm$covU)
    phencov <- cbind(phencov,fm$covU + fm$covE)
  }
  save(gencov,phencov,file="covariances.RData")
  
  #########regression coefficients
  
  #load("prepared_data.RData"); load("parameters.RData")
  
  load("covariances.RData")
    
    # Objects to store regression coefficients
    bSI <- bPCSI <- bL1PSI <- vector("list",nPart)
    for(k in 1:nPart)
    { cat("  partition = ",k,"of",nPart,"\n")
      indexTRN <- which(partitions[,k]==1)
      
      # Training set
      xTRN <- scale(X[indexTRN,])
      VARx <- var(xTRN)
      EVDx <- eigen(VARx)
      
      # Standard SI
      VARinv <- EVDx$vectors %*% diag(1/EVDx$values) %*% t(EVDx$vectors) ##phenotypic covariance matrix of xi  ##gencov[,k] genetic covariance matrix of xi 
      bSI[[k]] <- VARinv %*%  gencov[,k]    
      # PC-based SI
      VARinv <- diag(1/EVDx$values)
      gamma <- as.vector(VARinv %*% t(EVDx$vectors) %*%  gencov[,k])
      beta <- apply(EVDx$vectors %*% diag(gamma),1,cumsum)
      bPCSI[[k]] <- data.frame(I(beta),df=1:nrow(beta),lambda=NA)
      
      # L1-PSI
      fm <- solveEN(VARx, gencov[,k],nLambda=100)
      # fm <- LARS(VARx, gencov[,k])  # Second option
      beta <- t(fm$beta)[-1,]
      bL1PSI[[k]] <- data.frame(I(beta),df=fm$df[-1],lambda=fm$lambda[-1])
    }
    save(bSI,bPCSI,bL1PSI,file=("coefficients.RData"))
    
 
###Accuracy of the index 
      load(("coefficients.RData"))
      accSI <- c()    # Object to store accuracy components
      
      for(k in 1:nPart)
      { cat("  partition = ",k,"of",nPart,"\n")
        indexTST <- which(partitions[,k]==2)
        
        # Testing set
        xTST <- scale(X[indexTST,])
        yTST <- as.vector(scale(y[indexTST]))
        KTST <- K[indexTST,indexTST]   # Connection given by replicates
        
        # Calculate the indices
        SI <- xTST %*% bSI[[k]]
        PCSI <- tcrossprod(xTST, bPCSI[[k]]$beta)
        L1PSI <- tcrossprod(xTST, bL1PSI[[k]]$beta)
        
        # Fit genetic models   
        fitSI <- as.matrix(data.frame(I(SI),I(PCSI),I(L1PSI)))
        fm <- getGenCov(y1=yTST,y2=fitSI,K=KTST,warn=FALSE,mc.cores=1)
        fm$varU2 <- ifelse(fm$varU2<0.05,0.05,fm$varU2)
        
        # Retrieve accuracy components
        h <- sqrt(fm$varU2/(fm$varU2 + fm$varE2))
        gencor <- fm$covU/sqrt(fm$varU1*fm$varU2)
        acc2 <- fm$covU/sqrt(fm$varU1)
        accuracy <- abs(gencor)*h
        
        df <- c(nrow(bSI[[k]]),bPCSI[[k]]$df,bL1PSI[[k]]$df)
        lambda <- c(min(bL1PSI[[k]]$lambda),bPCSI[[k]]$lambda,bL1PSI[[k]]$lambda)
        accSI <- rbind(accSI,data.frame(rep=k,SI=colnames(fitSI),h,gencor,accuracy,df,lambda))
      }
      save(accSI,file=("accuracy_tp.RData"))
   

####### plot accuracy
load("accuracy_tp.RData")

accSI = split(accSI,as.character(accSI$SI))
accSI=data.frame(do.call(rbind,lapply(accSI,function(x)apply(x[,-c(1,2)],2,mean,na.rm=T))))
accSI$SI = unlist(lapply(strsplit(rownames(accSI),"\\."),function(x)x[1]))
           
                    
# Plot of PC-SI
requireNamespace("reshape")  
  requireNamespace("ggplot2")
  
  png("PCSI.png", width = 700, height = 600)
    dat <- reshape::melt(accSI[accSI$SI=="PCSI",],id=c("df","lambda","SI"))
    plot1 <- ggplot2::ggplot(dat,ggplot2::aes(df,value,group=variable,color=variable)) +
   ggplot2::theme_bw(base_size = 20) + ggplot2::geom_line(size=1.5) + 
      ggplot2::labs(title="Principal components-based selection index",y="value",x="Number of PCs")+
      ggeasy::easy_center_title()+
    ggplot2::scale_color_brewer(palette = "Dark2")
    plot1
    dev.off()
             
    # Plot of L1-PSI
    png("L1-PSI.png", width = 700, height = 600)
   dat <- reshape::melt(accSI[accSI$SI=="L1PSI",],id=c("df","lambda","SI"))
  plot2 <- ggplot2::ggplot(dat,ggplot2::aes(df,value,group=variable,color=variable)) +
    ggplot2::theme_bw(base_size = 20) + ggplot2::geom_line(size=1.5) + 
    ggplot2::labs(title="Lasso penalized selection index",y="value",x="Number of active predictors")+
     ggeasy::easy_center_title()+
     ggplot2::scale_color_brewer(palette = "Dark2")          
  plot2
  dev.off()
  
  

#  Optimal regularized vs standard selection indices 
  load("accuracy_tp.RData")
accSI$Model <- unlist(lapply(strsplit(as.character(accSI$SI),"\\."),function(x)x[1]))
accSI <- split(accSI,paste(accSI$rep,"_",accSI$Model))
accSI <- do.call(rbind,lapply(accSI,function(x)x[which.max(x$accuracy),]))
h<-0.607312358
accSI<-as.data.frame(accSI)
accSI$RE <- (accSI$accuracy)/h
write.csv(accSI, file="highestacc.csv")

           dat = aggregate(accuracy ~ Model, mean, data=accSI)
           dat$sd = aggregate(accuracy ~ Model,sd,data=accSI)$accuracy
           dat$n = aggregate(accuracy ~ Model,length,data=accSI)$accuracy
           dat$se = qnorm(0.975)*dat$sd/sqrt(dat$n)
           write.csv(dat, "optimalSI.csv")
           png("optimal_model.png", width = 700, height = 500)
              ggplot2::ggplot(dat,ggplot2::aes(Model,accuracy, fill=dat$Model)) +
              ggplot2::geom_bar(stat="identity",width=0.6) + 
              ggplot2::scale_fill_brewer(palette = "Dark2") +
              ggplot2::geom_errorbar(ggplot2::aes(ymin=accuracy-se,ymax=accuracy+se),width=0.2) +
              ggplot2::geom_text(ggplot2::aes(label=sprintf("%.2f",accuracy),y=accuracy*0.5)) +
              ggplot2::theme_bw(base_size = 20) + ggplot2::geom_line(size=1.5) +
              ggplot2::labs(title="Optimal regularized selection index",y="Accuracy",x="Model", fill = "Model")+
              ggeasy::easy_center_title()
              dev.off()
    
              
#Estimating efficiency of indirect selection reltive to mass phenotypic selection (RE)            
h<-0.607312358
accSI<-as.data.frame(accSI)
accSI$RE <- (accSI$accuracy)/h
write.csv(accSI, file="highestacc.csv")            
    
##selecting the highest accuracy models 
              
indexTST <- which(partitions[,5]==2)

# Testing set
xTST <- scale(X[indexTST,])
yTST <- as.vector(scale(y[indexTST]))
KTST <- K[indexTST,indexTST]   # Connection given by replicates

# Calculate the indices
x<- c(-1,-1,1,1,1,1,1,1)
x1<-c(-2,-2,2,2,2,2,2,2)
x2<-c(-3,-3,3,3,3,3,3,3)
x3<-c(-4,-4,4,4,4,4,4,4)
X_all<-cbind(x,x1,x2,x3)


# Calculate the indices
SI <- t(X_all) %*% bSI[[5]]              

a<-bL1PSI[[5]][99,]
a<-as.matrix(a)
a<- a[,1:8]
L1PSI <- tcrossprod(t(X_all),t(a))

PCSI <- tcrossprod(t(X_all)[,1:4],bPCSI[[2]]$beta[,1:4])
               
              
              
              
# Calculate the indices
library(fmsb)
SIplot <- cbind(t(X_all),SI)
SIplot<- as.data.frame(SIplot)
colnames(SIplot)<-c('Raffinose','stachyose','seed size','flowering','maturity date','oil seed','seed protein','plant height','sucrose')
max<-rep(4,9)
min<-rep(-4,9)
SIplot<-rbind(SIplot,max,min)

# Define colors and titles
colors <- c("#00AFBB", "#E7B800", "#FC4E07", "blue")
titles <- c("x1", "x2", "x3","x4")

# Reduce plot margin using par()
# Split the screen in 3 parts
op <- par(mar = c(1, 1, 1, 1))

  create_beautiful_radarchart(
    data = SIplot, caxislabels = c(-4,-3,-2,-1,0,1,2,3,4),
    color = colors[i], title = titles[i]
  )




PCSIplot <- tcrossprod(xTST, bPCSI[[k]]$beta)
L1PSIplot<- tcrossprod(xTST, bL1PSI[[k]]$beta)            
              
              
              
              
              
              
              
              
              
              
              
              
              
              
              
              
              
              
              
              
              
              
              
              
              
              
              
              
              
              
              
              
  