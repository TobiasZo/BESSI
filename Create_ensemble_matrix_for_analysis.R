#get sensitivity code
library(randtoolbox)


n <-ensemblesize/11# number of simulations to run

# List described the distribution of each variable
vals <- list(list(var="albedo_snow_new_spam",
                  dist="unif",
                  params=c(0.65,0.9)),
             list(var="albedo_snow_wet_spam",
                  dist="unif",
                  params=c(0.45,0.7)),
             list(var="albedo_ice_spam", 
                  dist="unif",
                  params=c(0.3,0.4)),
             list(var="D_sf_spam", 
                  dist="unif",
                  params=c(5,25)),
             list(var="ratio", #factor to calculate d_lf out of D-sf #zaker 0.65-0.97, #kondo 2.0 vs 2.1 #bintanja 0.36 #friehe CH=5*CE #In practice, the specific humidity term (q¯z ⫺ q¯s) is often replaced by (0.622/p)(e¯z ⫺ e¯0)
                  dist="unif",
                  params=c(0.5,1.2)),
             list(var="eps_air_spam", 
                  dist="unif",
                  params=c(0.6,0.9)),
             list(var="latent_heat_switch", #new part for different albedo parameters
                  dist="unif",
                  params=c(0,1)),
             list(var="albedo_module", #new part for different albedo parameters
                  dist="unif",
                  params=c(0,4)),
             list(var="max_lwc_spam",#Snow research in Svalbard an overview 10,8-20%; investigation of wet snow metamorphism in respect of lwc up to 14%)
                  dist="unif",
                  params=c(0.05,0.15)))

# Generate the quasi-random numbers [0,1]
Xall <- sobol(n, 2*length(vals))
Zall <- sobol(n, 2*length(vals),scrambling = 1)
X1 <- data.frame(Xall[1:n,1:length(vals)])
X2 <- data.frame(Xall[1:n,(length(vals)+1):(2*length(vals))])

#--------------------------------------------
# SOBOL MODEL
#--------------------------------------------
detach("package:randtoolbox", unload=TRUE)
library(sensitivity)
xx <- sobol2007(model = NULL, X1, X2, nboot = 10000)
#--------------------------------------------


#--------------------------------------------
# ADAPT DISTRIBUTIONS
#--------------------------------------------
samp <- matrix(rep(0,dim(xx$X)[1]*(length(vals))), nrow=dim(xx$X)[1])
for (i in 1:length(vals)) {
  l <- vals[[i]]
  print(i)
  dist <- l$dist
  params <- l$params
  samp[,i] <- eval(call(paste("q",dist,sep=""),xx$X[,i],params[1],params[2]))
}


# for (i in 1:dim(x$X)[1]) {
#   if (x$X[i,8]>3) {
#     x$X[i,8]=5
#   } else if (x$X[i,8]<=3 & x$X[i,8]>2) {
#     x$X[i,8]=4
#   } else if (x$X[i,8]<=2 & x$X[i,8]>1) {
#     x$X[i,8]=3
#   }  else {
#     x$X[i,8]=1
#   }
# } 


for (i in 1:dim(samp)[1]) {
  if (samp[i,7]>0.5) {
    samp[i,7]=1
  } else {
    samp[i,7]=0
  } 
  
}

for (i in 1:dim(samp)[1]) {
  if (samp[i,8]>3) {
    samp[i,8]=5
  } else if (samp[i,8]<=3 & samp[i,8]>2) {
    samp[i,8]=4
  } else if (samp[i,8]<=2 & samp[i,8]>1) {
    samp[i,8]=3
  }  else {
    samp[i,8]=1
  }
} 
#vector<-array(0,2000)


write.table(signif(samp,digits=4),'Sobol.sample',sep='\t',quote=FALSE,row.names=FALSE,col.names=FALSE)
write.table(signif(samp[,1],digits=4),'albedo_new.sample',sep='\t',quote=FALSE,row.names=FALSE,col.names=FALSE)
write.table(signif(samp[,2],digits=4),'albedo_wet.sample',sep='\t',quote=FALSE,row.names=FALSE,col.names=FALSE)
write.table(signif(samp[,3],digits=4),'albedo_ice.sample',sep='\t',quote=FALSE,row.names=FALSE,col.names=FALSE)
write.table(signif(samp[,4],digits=4),'D_sf.sample',sep='\t',quote=FALSE,row.names=FALSE,col.names=FALSE)
write.table(signif(samp[,5],digits=4),'D_lf.sample',sep='\t',quote=FALSE,row.names=FALSE,col.names=FALSE)
write.table(signif(samp[,6],digits=4),'e_air.sample',sep='\t',quote=FALSE,row.names=FALSE,col.names=FALSE)
write.table(signif(samp[,7],digits=4),'latent_heat_switch.sample',sep='\t',quote=FALSE,row.names=FALSE,col.names=FALSE)
write.table(signif(samp[,8],digits=4),'albedo_module.sample',sep='\t',quote=FALSE,row.names=FALSE,col.names=FALSE)
write.table(signif(samp[,9],digits=4),'max_lwc.sample',sep='\t',quote=FALSE,row.names=FALSE,col.names=FALSE)
# library(R.matlab)
