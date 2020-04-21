#' Simulates the COVID-19 dynamics of two inderacting populations
#' All variables and parameters come with an "o" and "y" subscript to indicate the "Over 65" and "under 65" ("Young") populations. 
#' @param bo.mult.init 
#' @param by.mult.init
#' @param base.mortality.ratio ratio of mortality to recovery rate
#' @param age.mortality.ratio ratio of the mortality rate for over_65 and under_65. The mortality ratio for under_65 is base.mortality.ratio*age.mortality.ratio  
#' @param N total population size. Default is 300e6, roughly the population of the US
#' @param sf Susceptible fraction. sf*N is the total number of susceptible individuals a time 0
#' @param step.up The time point in days at which social moditifaction  begins to linearly return to normal levels
#' @param overcrowdcost The multiplicative factor by which hospital overcrowding increases mortality
#' @param overcrowdthresh The threshold above which hospitals are overcrowded
#' @param nstep. The number of steps to run the model
#' @param step.size The size of the step in days
#' @param plot Make a plot
#' @return 
#' @examples


simulateSIR2pop=function(bo.mult.init, by.mult.init, base.mortality.ratio=1/30, age.mortality.ratio=1/50, N=300e6, sf=.25, step.up=180,overcrowdcost=2, overcrowdthresh=.5e6 , nstep=1000, step.size=1, plot=F){
  
  #hard coded parameters"
  over=.15 #fraction of population over-65, 
  Io=1000 #number of infected individuals
  Iy=1000
  year=365 #year, time point at which transmission levels return to normal values
  revertR=2.8 #COVID-19 R0 (the average number of individuals infected by a single individual in a completely susceptible population)"
  recoverytime=14 #recovery time in days. Taken to be 14 day
  hmratio=10 #the hospitalization/mortality ratio
  
  
  
  #population variables
  No0=N*over
  Ny0=N-No0

  #setting up initial conditions
  So=No0*sf #number susceptible
  Sy=Ny0*sf
  Ro=No0-So #number recovered
  Ry=Ny0-Sy
  My=Hy=0 #number of mortalities and hospitalizations is 0 at the beginning
  Mo=Ho=0

  
  

  

  
  
  #default values for gamma and delta variables that specify the rates of leaving the infected state
  go=gy=1/recoverytime #recovery rate
  do=go*base.mortality.ratio #mortality rate for the "over 65" population
  dy=(gy*base.mortality.ratio)*age.mortality.ratio #mortality rate for the "under 65" population,
  


  #define alpha, the total rate of leaving the infected state
  ao=go+do
  ay=gy+dy
  
 
  
  #set the actual beta parameters
  boo=bo.mult.init*(ao*No0/N+ay*Ny0/N)/sf
  byy=by.mult.init*(ao*No0/N+ay*Ny0/N)/sf
  byo=min(boo,byy) #we assume that the cross age group transmission is the minimum of the two
#  show("alphas")
#  show(c(ao,ay))
#show("betas")
 #   show(c(boo,byy,byo))
  
  out=matrix(nrow=nstep+1, ncol=10)
  colnames(out)=c("So", "Io", "Ro", "Ho","Mo", "Sy", "Iy", "Ry", "Hy","My")
  
  
  
  #compute the final betas's. These are the betas that correspond to natural transmission levels without intervention that the model reverts to at 365 days
  boo.final=revertR*(ao*No0/N+ay*Ny0/N)/sf
  byy.final=revertR*(ao*No0/N+ay*Ny0/N)/sf
  byo.final=min(boo.final,byy.final)
  
  #save the difference between final and initial values so we can step up
  boo.diff=boo.final-boo
  byy.diff=byy.final-byy
  byo.diff=byo.final-byo
  
 
  for(s in 1:nstep){
    totalhospital=Ho+Hy
    if (totalhospital>overcrowdthresh){
      overfrac=(totalhospital-overcrowdthresh)/totalhospital #the ratio of hospitalized individuals that are over capacity
      multiplier=overfrac*overcrowdcost+(1-overfrac)*1 # a weighted average of mortality for the below capacity and above capacity fraction of hospitalized individuals 
      thisdo=multiplier*do
        thisdy=multiplier*dy
  
    }
    else{
      thisdo=do
      thisdy=dy
    }
    thisao=go+thisdo
    thisay=gy+thisdy
    #save the current step
   
    out[s,]=c(So, Io, Ro, Ho,Mo, Sy, Iy, Ry, Hy,My)
    
    
    #compute the next step according to the SIR model
   
    
    
    newIy=Iy+(byy*Sy*Iy/N+byo*Sy*Io/N-thisay*Iy)*step.size
    Sy=Sy-(byy*Sy*Iy/N+byo*Sy*Io/N)*step.size
    Ry=Ry+gy*Iy*step.size
    My=My+thisdy*Iy*step.size
    Hy=hmratio*recoverytime*dy*Iy #dy instead of thisdy is used so that the excess mortalities (due to overcrowding) doesn't creat extra hospitalizations
    
    newIo=Io+(byo*So*Iy/N+boo*So*Io/N-thisao*Io)*step.size
    So=So-(byo*So*Iy/N+boo*So*Io/N)*step.size
    Ro=Ro+go*Io*step.size
    Mo=Mo+thisdo*Io*step.size
    Ho=hmratio*recoverytime*do*Io 
    
 
   
    
    Io=newIo
    Iy=newIy
    
    #assumption number 2: mitigation strategies must survive reintroduction
    Io=max(Io,1)
    Iy=max(Iy,1)
    
    #step up infection rates to normal levels by one year starting at step.up
    if (s*step.size>step.up && s*step.size<=year){
      boo=boo+boo.diff/((year-step.up)/step.size)
      byy=byy+byy.diff/((year-step.up)/step.size)
      byo=byo+byo.diff/((year-step.up)/step.size)
     
    }

    
  }
  #save the output of the last step
  out[nstep+1,]=c(So, Io, Ro, Ho,Mo, Sy, Iy, Ry, Hy,My)
  #print(c(My, Mo, My+Mo))

  if(plot){
    print( plotResults(out, step.size = step.size, oct=overcrowdthresh))
  }
  
#transmission values should be this at the end.
  #      boo=revertR*(ao*No0/N+ay*Ny0/N)/sf
   #     byy=revertR*(ao*No0/N+ay*Ny0/N)/sf
  #      byo=min(boo,byy)

  
  return(out)
  
}

#generates the simple plot with output from simulateSIR2pop
plotResults=function(out, step.size, oct=5e5){
  require(ggpubr) 
  require(reshape2)
  out=as.data.frame(out)
  
  out$time=1:nrow(out)*step.size
  outDF=melt(out,id.vars = "time")
  
  outDF$age=rep("over 65", nrow(outDF))
  outDF$age[substr(outDF$variable,2,2)=="y"]="under 65"
  outDF$variable=as.factor(substr(as.character(outDF$variable), 1,1))
  levels(outDF$variable) <- c('Hospitalizations', 'Infections', 'Mortalities', 'Recovered', 'Susceptibles')
  outDF=outDF[outDF$variable %in% c("Hospitalizations", "Infections", "Mortalities"),]
  ggline(outDF, x="time",
         y="value",
         color="variable",
         palette = c('green', 'red', 'black'), 
         facet.by="age", 
         plot_type = "l", 
         scales="free", 
         numeric.x.axis = T,
         ylab = "Number of individuals",
         xlab="Days") + 
    geom_hline(yintercept = oct, linetype=2) + 
    geom_vline(xintercept = 365, linetype=2) +
    labs(caption = 'Vertical lines show the duration of the social modification window which is 365 days. \n Horizontal lines are the hospital overcrowding threshold parameter.')
    #geom_text(aes(x=350, label="social modification window", y=1e6), angle=90)
    
  
  
}