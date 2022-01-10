#clear workspace
rm(list=ls())

#set working directory
setwd("") #where csv of data is (add your own file path)

#load libraries
library(reshape)
library(dplyr)
library(lubridate)
library(ggplot2)

#read in csv of model probabilities for each site
dat <- read.csv("PhotoProbabilitiesFromModel.csv")

#create new columns for better sorting of data, based on existing treatment designations and probabilities
dat <- dat %>%
  mutate(MPA = case_when(
    startsWith(Treatment, "MP") ~ "IN",
    startsWith(Treatment, "AR") ~ "IN",
    startsWith(Treatment, "OU") ~ "OUT",
    startsWith(Treatment, "BH") ~ "OUT"
    )
  )
  
dat <- dat %>%
  mutate(Habitat = case_when(
    endsWith(Treatment, "H") ~ "halo",
    endsWith(Treatment, "R") ~ "grass",
    endsWith(Treatment, "C") ~ "grass"
    )
  )

dat <- dat %>%
  mutate(Cover = case_when(
    startsWith(Treatment, "MPAa") ~ "algae",
    startsWith(Treatment, "MPAg") ~ "seagrass",
    startsWith(Treatment, "MPATr") ~ "seagrass",
    startsWith(Treatment, "OUTa") ~ "algae",
    startsWith(Treatment, "OUTg") ~ "seagrass",
    startsWith(Treatment, "AR") ~ "seagrass",
    startsWith(Treatment, "BH") ~ "seagrass"
    )
  )

dat <- dat %>%
  mutate(Fish = ifelse(
    Fish_proba>=0.5, "Yes", "No"
    )
  )

#make sure dates and times are stored as the correct format
dat$Date <- as_date(dat$Date)

dat$Date_Time <- as_datetime(dat$Date_Time)

      TimeNew <- as.POSIXlt(dat$Date_Time) #replacing day, month, year of all time items with the same values
  
      TimeNew$mday <- TimeNew[1]$mday
  
      TimeNew$mon <- TimeNew[1]$mon
  
      TimeNew$year <- TimeNew[1]$year
  
      TimeNew <- as.POSIXct(TimeNew)

dat$Time <- TimeNew #replaces dates for this column with the date of the first entry for time-of-day comparisons only


#differentiate multiple deployments at the same site
int1 <- interval(ymd("2018-12-01"), ymd("2020-01-01")) #all secondary deployments began after December, 2018

dat$Deployment <- vector(length=length(dat))

for(i in 1:nrow(dat)){
  
    if(dat[i,1]=="OUTg1H"){
      if(dat[i,3] %within% int1){
        dat[i,11] <- "2"
    } else {
      dat[i,11] <- "1"
    }
      } else if(dat[i,1]=="OUTg1C"){
        if(dat[i,3] %within% int1){
          dat[i,11] <- "2"
        } else {
          dat[i,11] <- "1"
        }
      } else if(dat[i,1]=="MPAg3H"){
        if(dat[i,3] %within% int1){
          dat[i,11] <- "2"
        } else {
          dat[i,11] <- "1"
        }
      } else if(dat[i,1]=="MPAg3C"){
        if(dat[i,3] %within% int1){
          dat[i,11] <- "2"
        } else {
          dat[i,11] <- "1"
        }
      } else {
      dat[i,11] <- "1"
      }
    }

dat$Treatment <- paste0(dat$Treatment, dat$Deployment)


#subset data to only include the first month of each deployment

      dat <- dat[with(dat, order(Date_Time, Treatment)),] #order dataframe by Date_Time and Treatment
      
      sdate <- dat[match(unique(dat$Treatment), dat$Treatment), 1:2] #create a dataframe with just the first date for each Treatment
      
      sdate$StartDate <- sdate$Date_Time
      
      sdate <- sdate[,-2]
      
      dat <- left_join(dat, sdate) #add a column with start dates for each Treatment (including separate deployments)
      
      dat$OneMonth <- dat$StartDate %m+% months(1) #add a column for one month from the start date
      
      dat1m <- dat[which(dat$Date_Time<=dat$OneMonth),] #create a subset with only the first month of each deployment


#compare probability of Fish photos ("detections") IN/OUT of MPA (no difference detected by divers)
MPAx <- c(length(dat1m[which(dat1m$Fish=="Yes" & dat1m$MPA=="IN"),5]), length(dat1m[which(dat1m$Fish=="Yes" & dat1m$MPA=="OUT"),5]))
MPAn <- c(length(dat1m[which(dat1m$MPA=="IN"),5]), length(dat1m[which(dat1m$MPA=="OUT"),5]))

prop.test(MPAx, MPAn, p=NULL, alternative="two.sided") #included in base R stats



#compare probability of Fish photos ("detections") between Habitat types, i.e. reef/halo vs. seagrass/algae (difference detected by divers/cams)
Habitatx <- c(length(dat1m[which(dat1m$Fish=="Yes" & dat1m$Habitat=="halo"),5]), length(dat1m[which(dat1m$Fish=="Yes" & dat1m$Habitat=="grass"),5]))
Habitatn <- c(length(dat1m[which(dat1m$Habitat=="halo"),5]), length(dat1m[which(dat1m$Habitat=="grass"),5]))

prop.test(Habitatx, Habitatn, p=NULL, alternative="two.sided") #included in base R stats


#compare probability of Fish photos ("detections") between Cover types, i.e. algae or seagrass
Coverx <- c(length(dat1m[which(dat1m$Fish=="Yes" & dat1m$Cover=="algae"),5]), length(dat1m[which(dat1m$Fish=="Yes" & dat1m$Cover=="seagrass"),5]))
Covern <- c(length(dat1m[which(dat1m$Cover=="algae"),5]), length(dat1m[which(dat1m$Cover=="seagrass"),5]))

prop.test(Coverx, Covern, p=NULL, alternative="two.sided") #included in base R stats



#graph proportion of images with fish in different habitats and protection states

      barsbase <- c((length(dat1m[which(dat1m$MPA=="IN"&dat1m$Habitat=="halo"&dat1m$Fish=="Yes"),5])/length(dat1m[which(dat1m$MPA=="IN"&dat1m$Habitat=="halo"),5])), 
                (length(dat1m[which(dat1m$MPA=="IN"&dat1m$Habitat=="grass"&dat1m$Fish=="Yes"),5])/length(dat1m[which(dat1m$MPA=="IN"&dat1m$Habitat=="grass"),5])),
                (length(dat1m[which(dat1m$MPA=="OUT"&dat1m$Habitat=="halo"&dat1m$Fish=="Yes"),5])/length(dat1m[which(dat1m$MPA=="OUT"&dat1m$Habitat=="halo"),5])),
                (length(dat1m[which(dat1m$MPA=="OUT"&dat1m$Habitat=="grass"&dat1m$Fish=="Yes"),5])/length(dat1m[which(dat1m$MPA=="OUT"&dat1m$Habitat=="grass"),5]))
                )
      
          Crossx <- c((length(dat1m[which(dat1m$MPA=="IN"&dat1m$Habitat=="halo"&dat1m$Fish=="Yes"),5])), 
                  (length(dat1m[which(dat1m$MPA=="IN"&dat1m$Habitat=="grass"&dat1m$Fish=="Yes"),5])),
                  (length(dat1m[which(dat1m$MPA=="OUT"&dat1m$Habitat=="halo"&dat1m$Fish=="Yes"),5])),
                  (length(dat1m[which(dat1m$MPA=="OUT"&dat1m$Habitat=="grass"&dat1m$Fish=="Yes"),5]))
          )
        
          Crossn <- c(length(dat1m[which(dat1m$MPA=="IN"&dat1m$Habitat=="halo"),5]), 
                  length(dat1m[which(dat1m$MPA=="IN"&dat1m$Habitat=="grass"),5]),
                  length(dat1m[which(dat1m$MPA=="OUT"&dat1m$Habitat=="halo"),5]),
                  length(dat1m[which(dat1m$MPA=="OUT"&dat1m$Habitat=="grass"),5])
          )
        
          prop.test(Crossx, Crossn, p=NULL, alternative="two.sided") #included in base R stats
          
          Crossp <- Crossx/Crossn
      
      sdbars <- sqrt(Crossn*Crossp*(1-Crossp))
      
      sebars <- sdbars/sqrt(Crossn)
          
      bars <- data.frame(
        Location=c("MPA Halo", "MPA Seagrass", "Outside Halo", "Outside Seagrass"),
        Proportion=barsbase,
        sd=(sdbars/Crossn)
      )
      
      barcomp <- ggplot(bars) +
          geom_bar(aes(x=Location, y=Proportion), stat="identity", fill=c("mediumorchid4", "seagreen4", "mediumorchid2", "seagreen2")) +
          geom_errorbar(aes(x=Location, ymin=Proportion-sd, ymax=Proportion+sd), width=0.3, color="black") +
          labs(y="Proportion of Images with Fish", x="Location") +
          theme(text=element_text(size=12,  family="Ariel"), 
                panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.background=element_blank(), 
                axis.ticks.x=element_blank(), axis.text.x=element_text(vjust=10), axis.title.y=element_text(margin=margin(r=25))) +
          geom_segment(aes(y=0, yend=0.5, x=-Inf, xend=-Inf))
      
      barcomp
  
  
#create a new column with the four treatment groups from this analysis
dat1m$HMPA <- paste0(dat1m$Habitat, dat1m$MPA) 


#graph fish detections (probality over time for MPA and Habitat by hour and month)

      colors <- c("seagreen4", "seagreen2", "mediumorchid4", "mediumorchid2") #create color pallet
      
      int2 <- interval(ymd_hms('2018-07-15 04:15:00'), ymd_hms('2018-07-15 19:15:00'))
      
      times <- dat1m[which(dat1m$Time %within% int2),] #subset out only observations falling during daylight hours (when most camears were deployed)
      
      time_plot <- ggplot(times, aes(x=Time, y=Fish_proba)) +
        theme(text=element_text(size=12,  family="Ariel"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.background = element_blank(), axis.line = element_line(color = "gray"),
              legend.key = element_blank(), legend.justification=c(1,1), legend.position=c(0.95,0.95), axis.title.y=element_text(margin=margin(r=20)),
              axis.title.x=element_text(margin=margin(b=25), vjust=-5)) +
        geom_smooth(aes(color = HMPA), method="loess", span=0.1, se=TRUE, size=1, alpha=0.2) +
        xlab("Time of Day") +
        ylab("Probability of Fish in Photo") +
        ylim(0,1) +
        scale_color_manual(values=colors, labels=c("MPA Seagrass", "Outside Seagrass", "MPA Halo", "Outside Halo")) +
        guides(color=guide_legend("Location"))
      
      time_plot
      
      
      season_plot <- ggplot(dat1m, aes(x=Date, y=Fish_proba)) +
        theme(text=element_text(size=12,  family="Ariel"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.background = element_blank(), axis.line = element_line(color = "gray"),
              legend.key = element_blank(), legend.position="none",
              axis.title.x=element_text(margin=margin(b=25), vjust=-5), axis.title.y=element_blank()) +
        geom_smooth(aes(color = HMPA), method="loess", span=.4, se=TRUE, size=1, alpha=0.2) +
        xlab("Month") +
        scale_color_manual(values=colors) +
        ylim(0,1) +
        scale_x_date(date_breaks = "2 months" , date_labels = "%b")
      
      season_plot




  


