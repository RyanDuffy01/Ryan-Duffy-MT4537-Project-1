library(spatstat)
library(rlist)


maternprocess <- function(lambdapp,lambdabaseline,range,radius){
  xupper <- range[2]
  xlower <- range[1]
  yupper <- range[2]
  ylower <- range[1]
  windowprocess <- owin(range,range)
  areaofwindow <- area(windowprocess)
  numberofpoints <- rpois(1,lambdabaseline*areaofwindow)
  locationsx <- runif(numberofpoints,min=xlower,max=xupper)
  locationsy <- runif(numberofpoints,min=ylower,max=yupper)
  pointprocess <- ppp(locationsx,locationsy,window=windowprocess)
  numberofparentpoints <- rpois(1,lambdapp*areaofwindow) 
  parentpoints <- sample(1:length(locationsx),numberofparentpoints)
  pplocationsx <- locationsx[parentpoints] 
  pplocationsy <- locationsy[parentpoints]
  parentpointprocess <- ppp(pplocationsx,pplocationsy,window=windowprocess)
  pointstoretain <- c()
  for (i in c(1:length(parentpoints))){
    circlepp <- disc(radius=radius,centre=c(locationsx[parentpoints[i]],locationsy[parentpoints[i]]))
    retained <- inside.owin(locationsx,locationsy,circlepp)
    retainedpoints <- which(retained==TRUE)
    pointstoretain <- c(pointstoretain,retainedpoints)
  }
  retainedpoints2 <- unique(pointstoretain)
  newprocesslocationsx <- locationsx[retainedpoints2]
  newprocesslocationsy <- locationsy[retainedpoints2]
  newpointprocess <- ppp(newprocesslocationsx,newprocesslocationsy,window=windowprocess)
  return(newpointprocess)
}

thomasprocess <- function(intensitypp,standarddev,meanpointcluster,range){
  xupper <- range[2]
  xlower <- range[1]
  yupper <- range[2]
  ylower <- range[1]
  windowprocess <- owin(range,range)
  areaofwindow <- area(windowprocess)
  ppnumberofpoints <- rpois(1,intensitypp*areaofwindow)
  pplocationsx <- runif(ppnumberofpoints,min=xlower,max=xupper)
  pplocationsy <- runif(ppnumberofpoints,min=ylower,max=yupper)
  pppointprocess <- ppp(pplocationsx,pplocationsy,window=windowprocess)
  locationsxcluster<- c()
  locationsycluster <- c()
  for (i in c(1:length(pplocationsx))){
    numberofpointsincluster <- rpois(1,meanpointcluster)
    locationsxcluster <- c(locationsxcluster,rnorm(numberofpointsincluster,mean=pplocationsx[i],sd=standarddev))
    locationsycluster <- c(locationsycluster,rnorm(numberofpointsincluster,mean=pplocationsy[i],sd=standarddev))
  }
  
  pointstoomitxupper <- range[2] < locationsxcluster 
  pointstoomitxuppernumber <- which(pointstoomitxupper==TRUE)
  pointstoomityupper <- range[2] < locationsycluster 
  pointstoomityuppernumber <- which(pointstoomityupper==TRUE)
  pointstoomitxlower <- range[1] > locationsxcluster 
  pointstoomitxlowernumber <- which(pointstoomitxlower==TRUE)
  pointstoomitylower <- range[1] > locationsycluster 
  pointstoomitylowernumber <- which(pointstoomitylower==TRUE)
  pointstoomit <- unique(c(pointstoomitxuppernumber,pointstoomitxlowernumber,pointstoomityuppernumber,pointstoomitylowernumber))  
  
  locationsxcluster <- locationsxcluster[-pointstoomit]
  locationsycluster <- locationsycluster[-pointstoomit]
  clusteredpointprocess <- ppp(locationsxcluster,locationsycluster,window=windowprocess)
  return(clusteredpointprocess)
}




par(mfrow=c(1,2))

range1 <- c(0,1)
windowpp1 <- owin(range1,range1)
range2 <- c(0,4)
windowpp2 <- owin(range2,range2)


thomasprocesssample1 <- thomasprocess(50,0.01,10,range1)
plot(thomasprocesssample)

thomasprocesssample2 <- thomasprocess(25,0.06,5,range2)
plot(thomasprocesssample2)


thomasfit1 <- thomas.estK(thomasprocesssample)
maternfit1 <- matclust.estK(thomasprocesssample)

thomasfit2 <- thomas.estK(thomasprocesssample2)
maternfit2 <- matclust.estK(thomasprocesssample2)

thomasfit1parameters <- c(thomasfit1$par[1],thomasfit1$par[2],thomasprocesssample1$n/thomasfit1$par[1])
maternfit1parameters <- c(maternfit1$par[1],maternfit1$par[2],thomasprocesssample1$n/maternfit1$par[1])

thomasfit2parameters <- c(thomasfit2$par[1],thomasfit2$par[2],thomasprocesssample2$n/thomasfit2$par[1])
maternfit2parameters <- c(maternfit2$par[1],maternfit2$par[2],thomasprocesssample2$n/maternfit2$par[1])

thomasfit1pp <- rThomas(kappa=thomasfit1parameters[1],scale=thomasfit1parameters[2],mu=thomasfit1parameters[3],win=windowpp1)
maternfit1pp <- rMatClust(kappa=maternfit1parameters[1],scale=maternfit1parameters[2],mu=maternfit1parameters[3],win=windowpp1)

thomasfit2pp <- rThomas(kappa=thomasfit2parameters[1],scale=thomasfit2parameters[2],mu=thomasfit2parameters[3],win=windowpp2)
maternfit2pp <- rMatClust(kappa=maternfit2parameters[1],scale=maternfit2parameters[2],mu=maternfit2parameters[3],win=windowpp2)

plot(thomasfit1pp)
plot(maternfit1pp)


plot(thomasfit2pp)
plot(maternfit2pp)

plot(bei)

densityfunctionbei <- density(bei,at="points")
densityfunctionbei.extra <- density()

bei.extradataframeelev <- as.data.frame(bei.extra$elev)
bei.extradataframegrad <- as.data.frame(bei.extra$grad)



beidataframe <- data.frame(
  x=bei$x,
  y=bei$y,
  density=densityfunction
)

beiextradataframe <- data.frame(
  x=bei.extradataframeelev$x,
  y=bei.extradataframeelev$y,
  elevation=bei.extradataframeelev$value,
  gradient=bei.extradataframegrad$value
)


fitbei.gaussian1 <- glm(density~x+y+x:y,family=gaussian,data=beidataframe)
fitbeiextra.gaussian1 <- glm()


summary(fitbei.poisson1)

bei$n









