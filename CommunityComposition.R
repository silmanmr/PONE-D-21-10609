#clear workspace
rm(list=ls())

#set working directory
setwd("") #where csv of data is (add your own file path)

#load libraries
library(vegan)
library(reshape)
library(dplyr)
library(rstatix)

#read in csv of species at each site (includes both diver obs and 3 images/site)
dat <- read.csv("CamValObsSp.csv")

#add new type_site_trt column
dat$Type_Site_Trt <- paste(dat$Type, dat$Site_Trt, sep='_')



######################################## NMDS for cams and diver observations #################################

#create a site-by-species matrix with values
site.sp.all<-as.data.frame(cast(dat, Type_Site_Trt ~ Species, value="Amount")) #records presence/absence, not counts

#replace NAs with 0s
site.sp.all[is.na(site.sp.all)] <- 0

#replace 2s (from multiple camera timepoints at the same location recording that species) with 1s
site.sp.all[site.sp.all == "2"] <- "1"

#remove non-numeric columns
site.sp.subA <- subset(site.sp.all, select=-c(Type_Site_Trt))

#remove ALGAE, FLOODED, EMPTY columns, etc.
site.sp.subA <- site.sp.subA[,-c(1,15,16,18,46)]

#write data to csv and reload to change data type to numeric (not the most efficient method, but it works)
write.csv(site.sp.subA, "AllObs.csv", row.names = FALSE)
site.sp.A <- read.csv("AllObs.csv")
site.sp.A <- as.data.frame(site.sp.A)

#assign zeros to the NA values
site.sp.A[is.na(site.sp.A)==TRUE] <- 0

#remove rows with no values in them
site.sp.A <- site.sp.A[-c(1,7,11,13,15),]


#plot a dendrogram and look at clustering

      #standardize/"scale" the data (generally good practice, not particularly important here)
      site.sp.sA <- scale(site.sp.A)
    
      #create a distance matrix
      site.sp.dA <- dist(site.sp.sA, method = "euclidean") #use site.sp.sA for cam-by-species matrix later

      #use Ward's Hierarchical Clustering
      site.sp.fitA <- hclust(site.sp.dA, method = "ward.D")
      
      #display the dendogram
      plot(site.sp.fitA)
      
      #cut tree into 4 clusters
      groupsA <- cutree(site.sp.fitA, k=4)
      
      #draw dendogram with red borders around the 4 clusters
      rect.hclust(site.sp.fitA, k=4, border="red")


#create a distance matrix using the Bray-Curtis distance (the default for vegdist)
all.bray <- vegdist(site.sp.A, binary=TRUE)

#run the NMDS on the Bray-Curtis distance matrix with k=2 axes
mds.all <- metaMDS(all.bray, k=2)

#look at the results
stressplot(mds.all)

#these are all the various data attributes produced by the MDS results
attributes(mds.all)

#defined so that sum of squared values is equal to squared stress;
#large values indicate poorer fit
goodness(mds.all)

#summary of the MDS analysis
mds.all

#plot of axis 1 vs. axis 2
plot(mds.all)

#draw NMDS ordination diagram with sites
plot (mds.all, display = 'sites', type = 't', main = 'Goodness of fit')

#add the points with size reflecting goodness of fit (bigger = worse fit)
points (mds.all, display = 'sites', type = 'p', cex = goodness(mds.all)*200)

#extract the data for MDS1 and MDS2 (i.e. the axis data)
tempA <- as.data.frame(mds.all$points[,1:2])

#make dataframe with sorting variables
sort <- as.data.frame(site.sp.all)
sort$Type <- c("cam","cam","cam","cam","cam","cam","cam","cam","cam","cam","cam","cam","cam","cam","cam","cam","cam","cam","dive","dive","dive","dive","dive","dive","dive","dive","dive","dive","dive","dive","dive","dive","dive","dive","dive","dive","dive","dive")
sort$Habitat <- c("grass", "halo", "grass","halo","halo","halo","grass","grass","halo","halo","grass","halo","grass","halo","grass","halo","grass","halo","grass","halo","grass","halo","halo","grass","halo","grass","grass","halo","grass","halo","grass","halo","grass","halo","grass","halo","grass","halo") 
sort <- sort[-c(1,7,11,13,15),]
     
#now put them into a dataframe with the group data from the cluster analysis
all.mat <- data.frame(Type=sort$Type, Habitat=sort$Habitat, groups=groupsA, tempA)

#set up parameters for plotting
par(mar=c(5.1,5,4.1,2.1))

#plot the NMDS results
plot(all.mat$MDS1,all.mat$MDS2, xlab="NMDS AXIS 1", ylab="NMDS AXIS 2",
     cex.lab=1, cex.axis=1, xlim=(c(-0.7,0.7)), ylim=(c(-0.5,0.5)), family="A")

#color points by Halo/Grass

      #grass, cams
      points(all.mat$MDS1[all.mat$Type=="cam"&all.mat$Habitat=="grass"],all.mat$MDS2[all.mat$Type=="cam"&all.mat$Habitat=="grass"],
             pch=17, col="seagreen4",cex=2.5) #making these points slightly larger so that overlapping points are not obscured
      
      #grass, divers
      points(all.mat$MDS1[all.mat$Type=="dive"&all.mat$Habitat=="grass"],all.mat$MDS2[all.mat$Type=="dive"&all.mat$Habitat=="grass"],
             pch=19, col="seagreen4",cex=2)
      
      #halo, divers
      points(all.mat$MDS1[all.mat$Type=="dive"&all.mat$Habitat=="halo"],all.mat$MDS2[all.mat$Type=="dive"&all.mat$Habitat=="halo"],
             pch=19, col="mediumorchid4",cex=2) 
      
      #halo, cams
      points(all.mat$MDS1[all.mat$Type=="cam"&all.mat$Habitat=="halo"],all.mat$MDS2[all.mat$Type=="cam"& all.mat$Habitat=="halo"],
             pch=17, col="mediumorchid4",cex=2) 
      
      
      legend("topright", legend = c("Halo, Diver", "Halo, Cam", "Grass, Diver", "Grass, Cam"),
             lty=NA, bty = "o", lwd = 2, pt.cex=1.5, cex = 1,
             col = c("mediumorchid4", "mediumorchid4","seagreen4","seagreen4"), pch = c(19, 17, 19, 17), xjust=0)


      

################################ ANOSIM analysis ##########################################

ano <- anosim(all.bray, all.mat$Habitat)
summary(ano)

#An R value close to "1.0" suggests dissimilarity between groups, 
#while an R value close to "0" suggests an even distribution of high and low ranks within and between groups. 
#R values below "0" suggest that dissimilarities are greater within groups than between groups. 
#See Clarke and Gorley (2001) for a guide to interpreting ANOSIM R values.

#Significance of the R statistic is determined by permuting group membership a large number of times 
#to obtain the null distribution of the R statistic. 
#Comparing the position of the observed R value to the null distribution allows an assessment 
#of statistical significance.




############################### Measure of diver effects ###################################

#clean data to remove rows with dat$Class == NA (i.e., dat$Species == "EMPTY" and others)
data <- na.omit(dat)


#subset by dat$Time (=="pre", =="on") and sum rows for each unique value in dat$Site_Trt
pre <- data[which(data$Time=="pre"),]
now <- data[which(data$Time=="on"),]
post <- data[which(data$Time=="post"),]

#create a vector of unique values for dat$Site_Trt
sites <- unique(data$Site_Trt)
      
fish_pre <- vector(length=length(sites))
      
for (i in 1:length(sites)){
        fish_pre[i] <- sum(data[which(pre$Site_Trt==sites[i]),2])
}
      
fish_now <- vector(length=length(sites))
      
for (i in 1:length(sites)){
  fish_now[i] <- sum(data[which(now$Site_Trt==sites[i]),2])
}
      
      
fish_post <- vector(length=length(sites))
    
for (i in 1:length(sites)){
  fish_post[i] <- sum(data[which(post$Site_Trt==sites[i]),2])
}

#combine fish counts from pre and now (while divers are present) into one dataframe with sites
comp <- data.frame(matrix(nrow=20, ncol=0))
comp$pre <- fish_pre
comp$now <- fish_now
comp$post <- fish_post
comp$site <- sites


#check for outliers, normality, etc.
mshapiro_test(comp[ , 1:3]) #multiple versions of the Shapiro-Wilk normality test suggest the data is not normal;
                            #this particular one is a multivariate normality test


#friedman test is the non-parametric equivalent to the one-way repeated measures anova

      mdat <- melt(comp) #convert data to long format for stats
      
      fried <- mdat %>% friedman_test(value~variable|site)
      
      fried
      
      mdat%>% friedman_effsize(value~variable|site) #relatively small effect size with Kendall's W

      
      
      
################################ PERMANOVA to compare methods ##############################

#create dataframe that corresponds to the 38 rows of data in site.sp.A (the basis for all.bray)
pdat <- unique(dat[,c(5,7,9,11)])

#remove rows with no values in them from site.sp.A above
pdat <- pdat[-c(1,7,11,13,15),]

#run the permanova on method, habitat type, and the interaction between the two
permanova <- adonis(all.bray~Type + Habitat + Type:Habitat, data=pdat, permutations=10000)

permanova



