.libPaths(c(.libPaths(), "C:/Users/zuazo-p/Documents/R/win-library/3.0", "C:/Users/zuazo-p/Documents/R/win-library/3.2.2"))
require(data.table)
require(ggplot2)
####################################################################################
# read files
####################################################################################
### csv
path <- "Z:/Group_workspace/Week_2_2016/Birds"
myfile <- paste0(path, "/hobson_birds_2H_feath_precip.csv")
# data <- fread("Hobson_bird_feathers.csv")
data <- fread(myfile)

### raster
data(wrld_simpl)

# 1 kriging
ebk.mean <- raster(paste0(path, "/ebk_2h_feat"))
crs(ebk.mean) <- projection(wrld_simpl)

ebk.se <- raster(paste0(path, "/ebk_2h_se"))
crs(ebk.se) <- projection(wrld_simpl)

# 2 precip isoscape (it needs to be transform using the model defined later)
precip <- raster(paste0(path, "/precip"))
crs(precip) <- projection(wrld_simpl)


# 
# ebk.mean <- raster("ebk_2h_feat")
# crs(ebk.mean) <- projection(wrld_simpl) # change coordinate reference system to latlon
# 
# ebk.se <- raster("ebk_2h_se")
# crs(ebk.se) <- projection(wrld_simpl)
# 
# precip <- raster("precip")
# crs(precip) <- projection(wrld_simpl)
# 
# rangemap <- raster("./btbwrange.asc")
# crs(rangemap) <- projection(wrld_simpl) # project to latlon

####################################################################################
# Define new functions (from mike)
####################################################################################

## function to generate posterior probabilities for raster cells
calcCellProb <- function(x,isoscape,stdv){
        m <- getValues(isoscape)
        if(class(stdv)[[1]]=="RasterLayer"){
                stdv <- getValues(stdv)
                m <- dnorm(x,mean=m,sd=stdv) 
        }else{
                m <- dnorm(x,mean=m,sd=stdv)   
        }

        cell.dens <- setValues(isoscape,m)
        return(cell.dens)
}

## function to get posterior probability distribution
calcPostProb <- function(x){
        pp <- x/cellStats(x,sum)
        return(pp)
}

## function to get normalized cell probabilites
calcNormProb <- function(x){
        np <- x/cellStats(x,max)
        return(np)
}

## function to make raster of odds ratios
calcOR<-function(x){ # x is a probability density raster
        oddsratio<-(maxValue(x)/(1-maxValue(x)))/(x/(1-x))
        # odds ratio for max of posterior density to cell of interest
        return(oddsratio)
}

## function to convert degrees (DD) to radians
deg2rad <- function(deg) return(deg*pi/180)

## function to get geodesic distance between two points using Spherical Law of Cosines
## xy specified by radian latitude/longitude 
calcDist <- function(long1, lat1, long2, lat2) {
        long1 <- deg2rad(long1)
        lat1 <- deg2rad(lat1)
        long2 <- deg2rad(long2)
        lat2 <- deg2rad(lat2)
        R <- 6371 # earth mean radius [km]
        d <- acos(sin(lat1)*sin(lat2) + cos(lat1)*cos(lat2) * cos(long2-long1)) * R
        return(d) # distance in km
}
####################################################################################
# Define transfer functions from precip isoscape to feather
####################################################################################
### Linear model
data[, d2H_precip:= as.numeric(d2H_precip)]
data[, Delta2H:= as.numeric(Delta2H)]

fit <- lm(Delta2H ~ d2H_precip, data=data)
summary(fit)
# hist(fit$residuals)
# fit.sd <- sd(fit$residuals)

# explore model
# 
par(mfrow=c(2,2))
plot(fit)
title(paste(fit$call[2]), outer=T, line = -2)
# 
par(mfrow=c(1,1))
qqPlot(fit, labels=rownames(data), simulate=TRUE, main=paste("Q-Q Plot \n",fit$cal[2]), id.method="y", id.n = 2)
# 
residplot <- function(fit, nbreaks=18) {
        z <- rstudent(fit)
        hist(z, breaks=nbreaks, freq=FALSE,
             xlab="Studentized Residual",
             main=paste("Distribution of Errors \n",fit$cal[2]))
        rug(jitter(z), col="brown")
        curve(dnorm(x, mean=mean(z), sd=sd(z)),
              add=TRUE, col="blue", lwd=2)
        lines(density(z)$x, density(z)$y,
              col="red", lwd=2, lty=2)
        legend("topright",
               legend = c( "Normal Curve", "Kernel Density Curve"),
               lty=1:2, col=c("blue","red"), cex=.7)
}
residplot(fit)
# 
hist(fit$residuals)
fit.sd <- sd(fit$residuals)
# 
# Cook s D values greater than 4/(n – k – 1), where n is the sample size and k is the number of predictor variables, indicate influential observations. D > 1 is more commonly used

cutoff <- 4/(nrow(data)-length(fit$coefficients)-2)
plot(fit, which=4, cook.levels=cutoff)
abline(h=cutoff, lty=2, col="red")
title("influential observations", outer=T, line = -2)
####################################################################################

## function to rescale precip for feathers
slope <- fit $ coefficient[[2]]
int <- fit $ coefficient[[1]]
rescale <- function(x){slope*x-int}

# choose sd for feathers
stdv<-c(11, 8, 10, 10, 8, 13, 4, 6, 12, 9, 10, 13, 5, 8, 10, 14, 9, 
        9, 11, 10, 9, 13, 10, 11, 8, 12, 8, 6, 18, 3, 9, 8, 5, 2, 6, 
        5, 5, 6, 5, 4, 3, 14, 17, 8, 7, 11, 17, 8, 11, 6, 4, 9, 6, 1, 
        2, 5, 9, 3, 4, 5, 6, 6, 13, 3, 11, 5, 4, 6, 5, 9, 5, 4, 7, 8,
        11, 4, 8, 8, 8, 13, 2, 10, 6, 5, 5, 7, 9, 14, 14, 14, 7, 10, 
        8, 18, 11, 12, 6, 6, 17, 16, 19, 18, 14, 9, 21, 15, 12, 12, 13, 
        8, 5, 4) # use this one for bootstrapping
stdv <- as.vector( sd(fit$residuals) ) # from residuals

precip.rescale <- rescale(precip) # rescale precip model for feathers

#################################################################################
## Fit data
## Get summary results
## extract odds value for the capture location
## extract distance and bearing to cell with maximum odds value
#################################################################################
#distance.list <- vector("list", length(birds[,1])) 
bird.assign.precip <- stack()
top10.precip <- stack()

bird.assign.krieg <- stack()
top10.krieg <- stack()

data[, nprob.precip:= as.numeric(NA)]
data[, nprob.krieg:= as.numeric(NA)]

# precip assignment
for (i in 1:nrow(data)){
        pp <- calcCellProb(data$Delta2H[i],precip.rescale, sample(c(stdv, stdv),size=1, replace=T)) # compute probabilities
        top <- pp >= quantile(getValues(pp),na.rm=T,0.9)[[1]]  # find top 10% of pp
        np <- calcNormProb(pp)

        data$nprob.precip[i] <- extract(np,cbind(data$Longitude[i],data$Latitude[i]))
        # bird.assign.precip <- raster::stack(bird.assign.precip, np)
        # top10.precip <- raster::stack(top10.precip, top)
        
        cat(i," ")
}

# krieg assignment
for (i in 1:nrow(data)){
        pp <- calcCellProb(data$Delta2H[i],ebk.mean, stdv = ebk.se) # compute probabilities
        top <- pp >= quantile(getValues(pp),na.rm=T,0.9)[[1]]  # find top 10% of pp
        np <- calcNormProb(pp)
        
        data$nprob.krieg[i] <- extract(np,cbind(data$Longitude[i],data$Latitude[i]))
        # bird.assign.krieg <- raster::stack(bird5assign.krieg, np)
        # top10.krieg <- raster::stack(top10.krieg, top)
        
        cat(i," ")
}

rm(i,pp,top, np)
data[, nprob.diff.krig_precp:= nprob.krieg - nprob.precip]
write.csv(data,"birds.csv")


hist(data$nprob.diff.krig_precp, na.rm=T, breaks=100, freq =F)

diff.results <- PToneway(nprob.precip ~ nprob.krieg, data, n=1)
summary.oneway(diff.results)

diff.results.W.MW <- W.MW.oneway(nprob.precip ~ nprob.krieg, data, n=1)
summary(diff.results.W.MW)

kruskal.test(nprob.precip ~ nprob.krieg, data)
t.test(data$nprob.precip, data$nprob.krieg, paired=T)

quantile(sort(data$nprob.diff.krig_precp), probs=0.95)
library(boot)
set.seed(1234)

# function to obtain R-Squared from the data 
mean2 <- function(data, indices) {
        d <- data[indices,] # allows boot to select sample 
        mymean <- mean(d, na.rm=T)
        return(mymean)
}

d <- data[!is.na(nprob.diff.krig_precp),nprob.diff.krig_precp]
boot <- as.vector(NA)
for(i in 1:10000){
        myd <- sample(d, size = nrow(data), replace = TRUE)
        boot[i] <- mean(myd, na.rm=T)
        # boot[i] <- mymean
}
quantile(sort(boot), probs=c(.025, .5, 0.975))

results <- boot(data=data[!is.na(nprob.diff.krig_precp), nprob.diff.krig_precp], statistic=mean2, R=1000)



results <- boot(data=data[!is.na(nprob.diff.krig_precp)], statistic=mean, formula=nprob.diff.krig_precp,
                R=1000)

results <- boot(data=data$nprob.diff.krig_precp, statistic=mean, R=1000)

print(results)

kruskal.test(data$nprob.precip, data$nprob.krieg,paired=TRUE)

require(multcomp)
require(car)
require(HH)
require(gvlma)
require(bootstrap)
require(MASS)
require(leaps)
require(coin)
require(lmPerm)
require(rmarkdown)
library(npar)


plot(density(data$nprob.krieg))
plot(density(data$nprob.krieg, na.rm=T, from =0, to =1))
plot(density(data$nprob.precip, na.rm=T, from =0, to =1))
lines(density(data$nprob.krieg, na.rm=T, from =0, to =1), col="red")

hist(data$nprob.precip, na.rm=T, breaks=100, ,
     freq =F)
par(new=TRUE)
hist(data$nprob.krieg, na.rm=T, breaks=100, freq =F)

hist(data$nprob.krieg, breaks=100, , freq =F)

hist(data$nprob.krieg, freq =F)
#################################################################################
## summary plots
#################################################################################
spplot(bird.assign.precip,col.regions=rev(terrain.colors(20)))
spplot(bird.assign.precip[[1]],col.regions=rev(terrain.colors(20)))
spplot(bird.assign.precip[[2]],col.regions=rev(terrain.colors(20)))

## Show maps of top 10% of prob densities and distances
pdf("top10maps.pdf")
opar<-par()
par(mfrow=c(2,2),mar=c(2,3,2,2)+0.01)
xl<-max(unlist(distance.list))
for (i in seq(along=birds[,1])){
        #  plot(bird.oddrast[[i]])
        plot(bird.top[[i]])
        plot(wrld_simpl,add=T)
        points(birds$lon[i],birds$lat[i],col="black",cex=1.5, pch=18)
        text(-65,35,paste("Bird",birds$sID[i]))
        text(-65,33,birds$age[i])
        text(-65,31,"Odds Ratio (Max:Cap)")
        text(-65,29,paste(round(birds$or4cap[i]),":1",sep=""))
        hist(distance.list[[i]],xlim=c(0,xl),xlab="",ylab="",main="")
}
par(opar)
dev.off()
rm(i)