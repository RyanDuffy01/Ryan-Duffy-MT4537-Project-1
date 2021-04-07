library(spatstat)

####SIMULATION####

#Create a function that creates a thomas cliuster process that takes in the desired intensity of parent points, the standard deviation of the multivariate normals, the desired mean number of points per cluster,and 2 lists of the length and breadth of the square window  
thomasprocess <- function(intensitypp,standarddev,meanpointcluster,xrange,yrange){
  #Assigns variables that are the maximum and minumumx values of the square window being created  
  xupper <- xrange[2]
  xlower <- xrange[1]
  yupper <- yrange[2]
  ylower <- yrange[1]
  #creates a square owin object that gives the window that the process is in
  windowprocess <- owin(xrange,yrange)
  #assigns a variable that stores the area of the window
  areaofwindow <- area(windowprocess)
  #assigns a variable that stores the number of parent points which is one random draw from a poisson distribution with lambda= the desired intensity per 1x1 square times the area overall of the window
  ppnumberofpoints <- rpois(1,intensitypp*areaofwindow)
  #assigns variables that store the randomly selected x and y coordinates of the parent points which are random draws from a random uniform distribution with minimum and max values being the max and min x and y values of the window  
  pplocationsx <- runif(ppnumberofpoints,min=xlower,max=xupper)
  pplocationsy <- runif(ppnumberofpoints,min=ylower,max=yupper)
  #create two empty lists in order to store the locations of the cluster points 
  locationsxcluster<- c()
  locationsycluster <- c()
  #creates a loop that loops the same number as the number of points in the parent point process
  for (i in c(1:length(pplocationsx))){
    #assigns a variable that stores the number of points in one cluster which is one random draw from a poisson distribution with lambda= mean number of points per cluster
    numberofpointsincluster <- rpois(1,meanpointcluster)
    #assigns variables that store the randomly selcted x and y coordinates of the parent points which are random draws from a normal distribution with mu= the location of a parent point and sd=the inputted standard deviation of the clusters 
    locationsxcluster <- c(locationsxcluster,rnorm(numberofpointsincluster,mean=pplocationsx[i],sd=standarddev))
    locationsycluster <- c(locationsycluster,rnorm(numberofpointsincluster,mean=pplocationsy[i],sd=standarddev))
  }
  #checks if any points of the clusters are outwith the window and if they are removes them
  pointstoomitxupper <- xrange[2] < locationsxcluster 
  pointstoomitxuppernumber <- which(pointstoomitxupper==TRUE)
  pointstoomityupper <- yrange[2] < locationsycluster 
  pointstoomityuppernumber <- which(pointstoomityupper==TRUE)
  pointstoomitxlower <- xrange[1] > locationsxcluster 
  pointstoomitxlowernumber <- which(pointstoomitxlower==TRUE)
  pointstoomitylower <- yrange[1] > locationsycluster 
  pointstoomitylowernumber <- which(pointstoomitylower==TRUE)
  pointstoomit <- unique(c(pointstoomitxuppernumber,pointstoomitxlowernumber,pointstoomityuppernumber,pointstoomitylowernumber))  
  locationsxcluster <- locationsxcluster[-pointstoomit]
  locationsycluster <- locationsycluster[-pointstoomit]
  #creates a point process from the new list of the cluster point x and y locations 
  clusteredpointprocess <- ppp(locationsxcluster,locationsycluster,window=windowprocess)
  #gets the function to return this point process
  return(clusteredpointprocess)
}


par(mfrow=c(1,2))

#defines a window for the 2 simulations
xrange1 <- c(0,10)
yrange1 <- c(0,10)
windowpp1 <- owin(xrange1,yrange1)


xrange2 <- c(0,40)
yrange2 <- c(0,60)
windowpp2 <- owin(xrange2,yrange2)

#defines the parameters of the process
process1parameters <- c(1,0.3,10)
process2parameters <- c(0.5,0.06,3)


#simulates two thomas processes using the function defined above 
thomasprocesssample1 <- thomasprocess(process1parameters[1],process1parameters[2],process1parameters[3],xrange1,yrange1)
plot(thomasprocesssample1)

thomasprocesssample2 <- thomasprocess(process2parameters[1],process2parameters[2],process2parameters[3],xrange2,yrange2)
plot(thomasprocesssample2)

#empirically finds the pair correlation function of the two simulations and plots them
sample1pcf <- pcf(thomasprocesssample1,correction="translate")
sample2pcf <- pcf(thomasprocesssample2,correction="translate")

plot(sample1pcf)
plot(sample2pcf)

#finds two fits of the simulations one as a thomas process and the other as a matern by using minimum contrast using the pcf
thomasfit1 <- thomas.estpcf(sample1pcf,method="SANN")
maternfit1 <- matclust.estpcf(sample1pcf,method="SANN")

thomasfit2 <- thomas.estpcf(sample1pcf,method="SANN")
maternfit2 <- matclust.estpcf(sample2pcf,method="SANN")

#extracts the fitted parameters and finds mu, the average number of points per cluster 
thomasfit1parameters <- c(thomasfit1$par[1],thomasfit1$par[2],thomasprocesssample1$n/(thomasfit1$par[1]*area(windowpp1)))
maternfit1parameters <- c(maternfit1$par[1],maternfit1$par[2],thomasprocesssample1$n/(maternfit1$par[1]*area(windowpp1)))

thomasfit2parameters <- c(thomasfit2$par[1],thomasfit2$par[2],thomasprocesssample2$n/(thomasfit2$par[1]*area(windowpp2)))
maternfit2parameters <- c(maternfit2$par[1],maternfit2$par[2],thomasprocesssample2$n/(maternfit2$par[1]*area(windowpp2)))

names(thomasfit1parameters)[3] <- "mu"
names(maternfit1parameters)[3] <- "mu"
names(thomasfit2parameters)[3] <- "mu"
names(maternfit2parameters)[3] <- "mu"

#creates a point process of the two simulated point processes based on the fitted parameters and plots them
thomasfit1pp <- rThomas(kappa=thomasfit1parameters[1],scale=thomasfit1parameters[2],mu=thomasfit1parameters[3],win=windowpp1)
maternfit1pp <- rThomas(kappa=maternfit1parameters[1],scale=maternfit1parameters[2],mu=maternfit1parameters[3],win=windowpp1)

thomasfit2pp <- rThomas(kappa=thomasfit2parameters[1],scale=thomasfit2parameters[2],mu=thomasfit2parameters[3],win=windowpp2)
maternfit2pp <- rThomas(kappa=maternfit2parameters[1],scale=maternfit2parameters[2],mu=maternfit2parameters[3],win=windowpp2)

#shows the absolute differences between the fitted parameters and the actual values of the parameters 
abs(thomasfit1parameters-process1parameters)
abs(maternfit1parameters-process1parameters)
abs(thomasfit2parameters-process2parameters)
abs(maternfit2parameters-process2parameters)

#plots  the two simulated point processes 
plot(thomasfit1pp)
plot(maternfit1pp)
plot(thomasfit2pp)
plot(maternfit2pp)






####FITTING TO DATA####

#finds an estimate of the intensity function using two different bandwidth selection methods
intensityfit1 <- density(bei,sigma=bw.ppl)
intensityfit2 <- density(bei,sigma=50)

#plots these intensity functions
plot(intensityfit1,main="Initial Intensity Estimate") 
plot(intensityfit2,main="Intensity Estimate with Sigma=50")

#fits intensity models (by finding coefficients) using the x and y coordinates 
fitbei1 <- ppm(bei,trend=~x+y) 
fitbei2 <- ppm(bei,trend=~poly(x,2)+poly(y,2))
fitbei3 <- ppm(bei,trend=~x+y+x:y)
fitbei4 <- ppm(bei,trend=~poly(x,5)+poly(y,5)+x:y)

#finds the AIC of these models
AICfitbei1 <- AIC(fitbei1)
AICfitbei2 <- AIC(fitbei2)
AICfitbei3 <- AIC(fitbei3)
AICfitbei4 <- AIC(fitbei4)

#plots the fitted intensity function of fit 4 which has the preffered fit
plot(fitbei4)


#fits intensity models (by finding coefficients) using the gradient and alt at certain points
fitbeiextra1 <- ppm(bei,trend=~elev+grad,data=bei.extra)
fitbeiextra2 <- ppm(bei,trend=~elev+grad+elev:grad,data=bei.extra)
fitbeiextra3 <- ppm(bei,trend=~poly(elev,2)+poly(grad,2)+elev:grad,data=bei.extra)
fitbeiextra4 <- ppm(bei,trend=~poly(elev,5)+poly(grad,5)+elev:grad,data=bei.extra)


#finds the AIC of these models
AICfitbeiextra1 <- AIC(fitbeiextra1)
AICfitbeiextra2 <- AIC(fitbeiextra2)
AICfitbeiextra3 <- AIC(fitbeiextra3)
AICfitbeiextra4 <- AIC(fitbeiextra4)

#plots the fitted intensity function of fit 4 which has the preffered fit
plot(fitbeiextra4)

#simulates a point process based on the preffered fitted model 
simulationprefmod <- simulate(fitbeiextra4,drop=T)

#finds the K-function, the spherical contact distribution function and pair correlation function for this simulation from the preffered model and plots them against the empirically found functions from the data 
kfunctionsim <- Kest(simulationprefmod)
plot(Kest(bei),main="K-function")
plot(kfunctionsim,col=6,add=T)
sphericalcontactsim <- Hest(simulationprefmod)
plot(Hest(bei),main="Spherical Contact")
plot(sphericalcontactsim,col=6,add=T)
pcfunctionsim <- pcf(simulationprefmod)
plot(pcf(bei),main="Pair Correlation Function")
plot(pcfunctionsim,col=6,add=T)

#finds the pair correlation function of bei
pcfbei <- pcf(bei)

#finds a fit of the parameters from the perspective that the processs is matern by using minimum contrast using the pcf
fitofbeimatern <- thomas.estpcf(pcfbei)

#creates a variable that stores the fit parameters
maternfitbeiparameters <- c(fitofbeimatern$par[1],fitofbeimatern$par[2],bei$n/(fitofbeimatern$par[1]*area(as.owin(bei))))

#creates and plots two simulations using this fit
thomassim1bei <- rThomas(kappa=maternfitbeiparameters[1],scale=maternfitbeiparameters[2],mu=maternfitbeiparameters[3],win=as.owin(bei)) 
thomassim2bei <- rThomas(kappa=maternfitbeiparameters[1],scale=maternfitbeiparameters[2],mu=maternfitbeiparameters[3],win=as.owin(bei))
                            
plot(thomassim1bei)
plot(thomassim2bei)

#plots the fitted intensity function of these two simulations 
plot(density(thomassim1bei))
plot(density(thomasssim2bei))

#plots the K-function, L-function and pair correlation function of the two simulations
plot(Kest(thomassim1bei))
plot(Kest(thomassim2bei))

plot(Lest(thomassim1bei))
plot(Lest(thomassim2bei))

plot(pcf(thomassim1bei))
plot(pcf(thomassim2bei))








