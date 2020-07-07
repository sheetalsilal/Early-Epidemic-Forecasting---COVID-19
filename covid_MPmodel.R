#Code developed by Modelling and Simulation Hub, Africa (MASHA) www.masha.uct.ac.za
#
#Granularity: Provincial  
#Disease: Covid-19
#Time Frame (Data): 2020-03-05 - 2020-03-14
#edit

#Notes 


# INSTRUCTIONS
# Set initial set up parameters (lines 25-32)
# Load in data set and define sheets accordingly
# Set names for RData image and top parameter sets object.
# Set number of cores for parallelisation (default=24)
# Set number of parameter sets to generate. 

library(deSolve)
library(adaptivetau)
library(ggplot2)
library(readxl)
library(openxlsx)

#set up to run in parallel
library(doParallel)
registerDoParallel(24)


#function to read in all sheets in a workbook
read_excel_allsheets <- function(filename, tibble = FALSE) {
  sheets <- readxl::excel_sheets(filename)
  x <- lapply(sheets, function(X) readxl::read_excel(filename, sheet = X))
  if(!tibble) x <- lapply(x, as.data.frame)
  names(x) <- sheets
  x
}
# ************************************************************************************* #
#Set working directory, read in C files, load packages, start up model dimensions

setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # set wd to source file
system("R CMD SHLIB eq0.c")     # compiles c code
dyn.load("eq0.so")              # loads compiled code to memory (.so - OS X, .dll - Windows)

N<-9   # number of patches
B<-10   # number of variables per patch
A<-22  # number of transitions per patch
V<-N*B # total number of variables
L<-N*A #total number of transitions
startyear=0 # starting year of simulation 0=2020-03-01
tyears<-(10*7)/365 # total years of simulation 1
dtout<-1/365 # output timestep
tsteps<-round(tyears/dtout) # number of time steps
time<-startyear+seq(0,tyears,dtout) # time vector

# Call C function from eq0.so in R
EQ<-function(L, N, oldeq, transit,transitionsiu1,transitionsiu2,transitionsiv1,transitionsiv2){
  len<-length(oldeq)
  .C("EQ",
     as.integer(L),
     as.integer(N),
     as.double(oldeq),
     as.double(transit),
     as.integer(transitionsiu1),
     as.integer(transitionsiu2),
     as.integer(transitionsiv1),
     as.integer(transitionsiv2),
     as.double(vector("double",length = len))
  )[[9]]
  
}

# ************************************************************************************* #
# import data
# ************************************************************************************* #
alldata = read_excel_allsheets("za_data.xlsx") 
parameters = alldata$parameters[,c(1,3,4,5)]

# population
pvxy<-alldata$dem  #demographic data  
imp<-alldata$imp; imp[is.na(imp)]=0 # initial imported cases
loc<-alldata$local; loc[is.na(loc)]=0 # initial local cases
cases<-alldata$cases; cases[is.na(cases)]=0 # cumulative cases

#Seed infection 
seedst<-alldata$seed; seedst[is.na(seedst)]=0 # initial seed cases

# ************************************************************************************* #
# define variables
# ************************************************************************************* #
# Covid- 19
# 1=S: uninfected non-immune
# 2=E: infected & exposed
# 3=Ip: exposed and infectious 
# 4=Im: mild and infectious
# 5=Is: severe and infectious
# 6=Ic: critical and infectious 
# 7=Imt: mild and confirmed
# 8=Ist: severe and confirmed
# 9=Ict: critical and confirmed
# 10=R: Recovered

covpop<-1:10

# ************************************************************************************* #
# define indices
# ************************************************************************************* #
varind<-matrix(0,nrow=B,ncol=N)
traind<-matrix(0,nrow=A,ncol=N)
for (n in 1:N){
  for (b in 1:B){
    varind[b,n]<-(n-1)*B+b
  }
  for (a in 1:A){
    traind[a,n]<-(n-1)*A+a
  }
}

# ************************************************************************************* #
# define transitions
# ************************************************************************************* #
# first transition is given without index
transitions =ssa.maketrans(V,rbind(varind[4,1], +1)) 
for (n in 1:N){
  transitions[traind[1,n]]<-ssa.maketrans(V,rbind(varind[1,n], 0,varind[7,n], +1)) # imported cases -> Im=4 
  transitions[traind[2,n]]<-ssa.maketrans(V,rbind(varind[1,n], -1,varind[2,n],+1)) # incidence S to E 
  transitions[traind[3,n]]<-ssa.maketrans(V,rbind(varind[2,n], -1,varind[3,n],+1)) # incubation E to Ip 
  transitions[traind[4,n]]<-ssa.maketrans(V,rbind(varind[3,n], -1,varind[4,n],+1)) # mild infection Ip to Im
  transitions[traind[5,n]]<-ssa.maketrans(V,rbind(varind[3,n], -1,varind[5,n],+1)) # severe infection Ip to Is
  transitions[traind[6,n]]<-ssa.maketrans(V,rbind(varind[3,n], -1,varind[6,n],+1)) # critical infection Ip to Ic  
  transitions[traind[7,n]]<-ssa.maketrans(V,rbind(varind[4,n], -1,varind[10,n],+1)) # natural recovery Im to R 
  transitions[traind[8,n]]<-ssa.maketrans(V,rbind(varind[5,n], -1,varind[10,n],+1)) # natural recovery Is to R
  transitions[traind[9,n]]<-ssa.maketrans(V,rbind(varind[6,n], -1,varind[10,n],+1)) # natural recovery Ic to R
  transitions[traind[10,n]]<-ssa.maketrans(V,rbind(varind[4,n], -1,varind[7,n],+1)) # trt seeking Im to Imt
  transitions[traind[11,n]]<-ssa.maketrans(V,rbind(varind[5,n], -1,varind[8,n],+1)) # trt seeking Is to Ist
  transitions[traind[12,n]]<-ssa.maketrans(V,rbind(varind[6,n], -1,varind[9,n],+1)) # trt seeking Ic to Ict
  transitions[traind[13,n]]<-ssa.maketrans(V,rbind(varind[7,n], -1,varind[10,n],+1)) # recovery Imt to R
  transitions[traind[14,n]]<-ssa.maketrans(V,rbind(varind[8,n], -1,varind[10,n],+1)) # recovery Ist to R 
  transitions[traind[15,n]]<-ssa.maketrans(V,rbind(varind[9,n], -1,varind[10,n],+1)) # recovery Ict to R
  transitions[traind[16,n]]<-ssa.maketrans(V,rbind(varind[6,n], -1,varind[7,n],0)) # death Ic -> 
  transitions[traind[17,n]]<-ssa.maketrans(V,rbind(varind[9,n], -1,varind[8,n],0)) # death Ict -> 
  transitions[traind[18,n]]<-ssa.maketrans(V,rbind(varind[3,n], -1,varind[10,n],+1)) # natural recovery Ip to R
  transitions[traind[19,n]]<-ssa.maketrans(V,rbind(varind[8,n], -1,varind[9,n],+1)) # progression to critical disease Ist to Ict 
  transitions[traind[20,n]]<-ssa.maketrans(V,rbind(varind[8,n], 0,varind[4,n],+1)) #  seed -> Imt=7  
  transitions[traind[21,n]]<-ssa.maketrans(V,rbind(varind[1,n], 0,varind[8,n], +1)) # imported cases -> Is=4 
  transitions[traind[22,n]]<-ssa.maketrans(V,rbind(varind[1,n], 0,varind[9,n], +1)) # imported cases -> Ic=4 
  
  
}
#Alternate formulation of transitions matrix (used in epimodel function)
transitions2<-NULL
for (i in 1: length(transitions)){
  transitions2<-rbind(transitions2,cbind(as.integer(names(transitions[[i]]))[1],as.integer(names(transitions[[i]]))[2], transitions[[i]][1], transitions[[i]][2]))
}
row.names(transitions2)<-NULL
transitionsiu1<-transitions2[,1]
transitionsiu2<-transitions2[,2]
transitionsiv1<-transitions2[,3]
transitionsiv2<-transitions2[,4]


# ************************************************************************************* #
# Set the parameters (vivax-specific parameters prefixed with v)
# ************************************************************************************* #
pars = {list(
         gamma1 = parameters[1,2],  #1/non-infectious incubation duration
         gamma2 = parameters[2,2],  #1/infectious incubation duration
         pm=parameters[3,2], #proportion of mild cases
         ps=parameters[4,2], #proportion of severe cases
         taum=parameters[5,2], # 1/trt seeking mild
         ptm=parameters[6,2], # prob trt seeking mild
         pts=parameters[7,2], # prob trt seeking severe
         ptc=parameters[8,2], # prob trt seeking critical 
         r1=parameters[9,2], # 1/dur infectiousness
         r2=parameters[10,2], # 1/dur infectiousness
         mu = parameters[11,2], #death rate
         pd1 = parameters[12,2], # death rate (critical untreated)
         pd2 = parameters[13,2], # death rate (critical treated)
         alpha1 = parameters[14,2], # prop sus due to non symp social distancing
         alpha2 = parameters[15,2], # prop sus due to SI without trt seeking
         alpha3 = parameters[16,2], # prop sus due to SI with trt seeking
         alpha4 = 0*parameters[17,2], # prop sus due to hospital quarantine
         beta0 = rep(parameters[18,2],N), # no. of contacts
         R0 = parameters[19,2], # basic reproductive number
         ptrans = parameters[20,2], # probability of transmission
         indexhet = rep(parameters[21,2], N), # connectivity between patches
         taus=parameters[22,2], # 1/trt seeking severe
         tauc=parameters[23,2], # 1/trt seeking critical
         pa = parameters[24,2], # asymptomatic proportion of cases
         pprogc = parameters[25,2], # proportion of severe cases who progress to critical
         tauprog = parameters[26,2], # 1/ duration of progression of severe cases to critical 
         zeta1 = parameters[27,2], #relative infectiousness of asymptomatic population
         pmu = parameters[28,2], #proportion of UNDETECTED mild cases
         psu = parameters[29,2], #proportion of UNDETECTED severe case        
         psd=parameters[30,2], #proportion of DETECTED severe cases
         pcd=parameters[31,2], #proportion of DETECTED critical cases
         r3=parameters[31,2], #1/dur stay hospital (severe)
         r4=parameters[31,2] #1/dur stay ICU (critical)
         
);
} 


#Compute B0 from R0
R0estimator<-function(pars){
  with(as.list(c( pars)),
    {
      A1 = (1-pa)*pm*gamma2/(taum*ptm+(1-ptm)*r1)
      A2 = (1-pa)*ps*gamma2/(taus*pts+(1-pts)*r1)
      A3 = (1-pm-ps)*(1-pa)*gamma2/(tauc*ptc+(1-ptc)*r1 +mu*pd1)
      R0e = mean(beta0)*ptrans*(alpha1*zeta1*+alpha2*(A1+A2+A3))/((1-pa)*gamma2+pa*r1)
      betae = R0*((1-pa)*gamma2+pa*r1)/(ptrans*(alpha1*zeta1*+alpha2*(A1+A2+A3)))
      
    return(c(R0e, betae))
  })
}
R0e<-R0estimator(pars)[1]
betae<-R0estimator(pars)[2] #<--Use this beta estimated from R0 specified in pars vector

# ************************************************************************************* #
# Function to calculate inputs to transition rates
# ************************************************************************************* #

inputs<-function(parmal, scenario){
  
 
  return(list())
  }
  
# ************************************************************************************* #
#Set up matrices for rate function

# ************************************************************************************* #
# Function to calculate transition rates, given variables and parameters
# ************************************************************************************* #
covrates <- function(x, input, parmal, t,ti, scenario) {
  with(as.list(c( parmal, scenario)),
       {
  t_internal<-(ti-1)*dtout+t+startyear
  #Set up matrices
  seas<-c(rep(1,N)) # seasonality
  pop<-popc<-c(rep(0,N))    # population sizes  
  foi<-c(rep(0,N))     # forces of infection 
  import<-seed<-c(rep(0,N))
  
  betae<-R0estimator(parmal)[2] #<--Use this beta estimated from R0 specified in pars vector
  
  for (n in 1:N){popc[n]<-sum(x[varind[covpop,n]])} # list of the variable indices for the human population
  pop[n]<-sum(x[varind[covpop,n]])
  tranrate<-matrix(0,nrow=N,ncol=A)   # transition rate matrix
   for (n in 1:N){

    foi<-seas*(beta0[n]*ptrans*(alpha1*(zeta1*x[varind[3,n]])+alpha2*(x[varind[4,n]]+x[varind[5,n]]+x[varind[6,n]])+alpha3*(x[varind[7,n]])+alpha4*(x[varind[8,n]]+x[varind[9,n]]))/popc[n])

    
    if (t_internal<0.5) {
      import[n]<-approx(imp$Step, imp[,n+3], t_internal, rule=2)$y
      seed[n]<-approx(seedst$Step, seedst[,n+3], t_internal, rule=2)$y
    }
    
    tranrate[n,]<-c(      
      pm*import[n], #  importation to Imt 1
      foi[n]*alpha1*x[varind[1,n]],       # incidence S to E     2
      gamma1*x[varind[2,n]],       # incubation E to Ip     3
      (1-pa)*pm*gamma2*x[varind[3,n]],       # mild infection Ip to Im       4
      (1-pa)*ps*gamma2*x[varind[3,n]],       #  severe infection Ip to Is       5
      (1-pa)*(1-pm-ps)*gamma2*x[varind[3,n]],       #  critical infection Ip to Ic         6
      r1*(1-ptm)*x[varind[4,n]],       # natural recovery Im to R         7
      r1*(1-pts)*x[varind[5,n]],       # natural recovery Is to R         8
      r1*(1-pd1)*(1-ptc)*x[varind[6,n]],       # natural recovery Ic to R         9
      taum*ptm*x[varind[4,n]],       # trt seeking Im to Imt              10
      taus*pts*x[varind[5,n]],       # trt seeking Is to Ist             11
      tauc*ptc*x[varind[6,n]],       # trt seeking Ic to Ict                12
      r2*x[varind[7,n]],  #      recovery Imt to R   13
      (1-pprogc)*r3*x[varind[8,n]],  #      recovery Ist to R  14
      r4*(1-pd2)*x[varind[9,n]],  #      recovery Ict to R   15
      mu*pd1*(1-ptc)*x[varind[6,n]],      #     death Ic ->    16
      mu*pd2*x[varind[9,n]],      #     death Ict ->    17
      r1*pa*x[varind[3,n]],       # natural recovery Ip to R   18
      pprogc*tauprog*x[varind[8,n]], # progression from severe disease to critical disease 19
      seed[n], # initial seed - at least 50% undetected  20
      psd*import[n], #importation to Ist 1 21
      pcd*import[n] #importation to Ict 22
      )
  }
  return(c(t(tranrate)))
       })
}

# POST PROCESSING function
postproc <- function(parpro,out,tran) {
  with(as.list(c( parpro)),
       {
         # ************************************************************************************* #
         # for outputting the  time series for each patch
         # ************************************************************************************* #
         
         # Case outputs
         locd_pred<-impd_pred<-prev_crit<-prev_sev<-prev_pred<-mildb_pred<-sevb_pred<-critb_pred<-deathb_pred<-mild_pred<-sev_pred<-crit_pred<-death_pred<-matrix(0,nrow=length(out[,1]),ncol=N)
        for (n in 1:N){
           mildb_pred[,n]<-tran[,c(traind[4,n])]/365 + tran[,c(traind[20,n])]
           sevb_pred[,n]<-(tran[,c(traind[5,n])])/365
           critb_pred[,n]<-(tran[,c(traind[c(6),n])])/365
           mild_pred[,n]<-(tran[,c(traind[10,n])])/365 + tran[,c(traind[1,n])]
           sev_pred[,n]<-(tran[,c(traind[11,n])])/365 + tran[,c(traind[21,n])]
           crit_pred[,n]<-rowSums(tran[,c(traind[c(12,19),n])])/365 + tran[,c(traind[22,n])]
           deathb_pred[,n]<-(tran[,c(traind[16,n])])/365
           death_pred[,n]<-(tran[,c(traind[17,n])])/365
           prev_pred[,n]<-rowSums(out[,c(varind[2:9,n])+1])/rowSums(out[,(varind[covpop,n]+1)])
           prev_sev[,n]<-out[,c(varind[8,n])+1]
           prev_crit[,n]<-out[,c(varind[9,n])+1]
           impd_pred[,n]<-rowSums(tran[,c(traind[c(1,21,22),n])])
           locd_pred[,n]<-rowSums(tran[,c(traind[10:12,n])]/365)
        }
         
         
         return(cbind(mildb_pred,    #1
                      sevb_pred,     #2
                      critb_pred,  #3
                      mild_pred,   #4
                      sev_pred,    #5
                      crit_pred,    #6
                      deathb_pred,    #7
                      death_pred,  #8
                      prev_pred, #9
                      prev_sev, #10
                      prev_crit,  #11
                      impd_pred, #12
                      locd_pred #13
         ))
         
       })
}

ti<-1

# ************************************************************************************* #
# ************************************************************************************* #
# ************************************************************************************* #
# ODE SOLVER
# ************************************************************************************* #
# ************************************************************************************* #

epiModel<-function(t,state, parode,input, scenario) 
{ 
  
  with(as.list(c(state, parode)),
       {
         
         #   # ************************************************************************************* #
         #   # define variables
         #   # ************************************************************************************* #
         
         Z=state
         
         # rates of change
         ti<-1
         transit<-covrates(Z[1:V],input,parode,Z[V+1],ti, scenario)  
         
         if (sum(is.na(transit))>0)  {
           stop("transit NA   ",Z[V+1], "                                      ", 
                as.data.frame(transit))
         }
         
         # if (min(transit)< -1){
         #   which(sapply(transit, function(x){x < -1.07e+03}))
         #   browser()
         # }
         
         #    transit2<-apply(cbind((rpois(L, transit)),rep(0,L)), 1, FUN=max, na.rm=T)
         transit2<-transit#apply(cbind((rnbinom(L, mu=transit, size=4)),rep(0,L)), 1, FUN=max, na.rm=T)
         
         # if (sum(is.na(transit2))> 0){
         #   browser()
         # }
         
         eq<-rep(0.0, V)
         
         eq<-EQ(L, N, eq, transit2,transitionsiu1,transitionsiu2,transitionsiv1,transitionsiv2)
         
         eq[V+1]<-1
         
         dZ<-eq
         
         # return the rate of change
         list(c(dZ))
       }
  ) 
  # end with(as.list ...
}

#*************************************************************************
#*************************************************************************
scenario<-NULL

# MODEL RUN FUNCTION
run <- function(parrun, scenario){
  
# ************************************************************************************* #
  # define initial conditions
  initcondrun<-NULL
  for (n in 1:N) { initcondrun<-c(initcondrun, c(pvxy[n,2],rep(0,9)))}
  # initcondrun[varind[7,]]<-as.numeric(cases[18,4:12])  #starting cumulative cases at 18 March 2020
  # initcondrun[varind[7,]][initcondrun[varind[7,]]==0]<-1 # starting infection on 18 March in uninfected provinces

  # all initial conditions must be integers
  initoderun<-initcondrun
  staterun <- c(initoderun,0)
  ti<-1
  inp<-inputs(parrun)
  transitrun <- covrates(initoderun,inp,parrun,0,ti,scenario )
  
  #
  # # SOLVE THE ODEs and get output
  timesrun <- seq(0, tyears, by = dtout) # Model run time
  #Solve ODE
  outoderun <- ode(y = staterun, times = timesrun, func = epiModel, parms = parrun, method  = "vode", input=inp,scenario=scenario )
  # Compute transitions at each time step
  tranoderun<-matrix(0,nrow=length(outoderun[,1]),ncol=length(transitions))
  for (ti in 1:(tsteps+1)){
    tranoderun[ti,]<-t(covrates(outoderun[ti,2:(1+V)],inp,parrun,0,ti,scenario ))
  }
  #Compute outputs
  ppout<-postproc(parrun,outoderun,tranoderun)
  modeltimes<-outoderun[,1]+startyear
  
  mildb_pred_ode<-ppout[,1:N]
  sevb_pred_ode<-ppout[,(N+1):(2*N)]
  critb_pred_ode<-ppout[,(2*N+1):(3*N)]
  mild_pred_ode<-ppout[,(3*N+1):(4*N)]
  sev_pred_ode<-ppout[,(4*N+1):(5*N)]
  crit_pred_ode<-ppout[,(5*N+1):(6*N)]
  deathb_pred_ode<-ppout[,(6*N+1):(7*N)]
  death_pred_ode<-ppout[,(7*N+1):(8*N)]
  prev_pred_ode<-ppout[,(8*N+1):(9*N)]
  prev_sev_ode<-ppout[,(9*N+1):(10*N)]
  prev_crit_ode<-ppout[,(10*N+1):(11*N)]
  impd_pred_ode<-ppout[,(11*N+1):(12*N)]
  locd_pred_ode<-ppout[,(12*N+1):(13*N)]
  incb_pred_ode<- mildb_pred_ode+sevb_pred_ode+critb_pred_ode
  inc_pred_ode<-mild_pred_ode+sev_pred_ode+crit_pred_ode
  
  #Weekly sums
  # Counts for data fitting
  modelweeks<-1:(tyears*365/7)
  dweeks<-dtimes[seq(1,length(dtimes),7)]
  
  mildbW<-sevbW<-critbW<-mildW<-sevW<-critW<-deathbW<-deathW<-incbW<-incW<-impdW<-locdW<-array(0,c(tyears*365/7,N))
  for (n in 1:N){
  mildbW[,n]<-tapply(head(mildb_pred_ode[,n],-1), cut(1:length(head(mildb_pred_ode[,n],-1)),tyears*365/7),FUN=sum)
  sevbW[,n]<-tapply(head(sevb_pred_ode[,n],-1), cut(1:length(head(sevb_pred_ode[,n],-1)),tyears*365/7),FUN=sum)
  critbW[,n]<-tapply(head(critb_pred_ode[,n],-1), cut(1:length(head(critb_pred_ode[,n],-1)),tyears*365/7),FUN=sum)
  mildW[,n]<-tapply(head(mild_pred_ode[,n],-1), cut(1:length(head(mild_pred_ode[,n],-1)),tyears*365/7),FUN=sum)
  sevW[,n]<-tapply(head(sev_pred_ode[,n],-1), cut(1:length(head(sev_pred_ode[,n],-1)),tyears*365/7),FUN=sum)
  critW[,n]<-tapply(head(crit_pred_ode[,n],-1), cut(1:length(head(crit_pred_ode[,n],-1)),tyears*365/7),FUN=sum)
  deathbW[,n]<-tapply(head(deathb_pred_ode[,n],-1), cut(1:length(head(deathb_pred_ode[,n],-1)),tyears*365/7),FUN=sum)
  deathW[,n]<-tapply(head(death_pred_ode[,n],-1), cut(1:length(head(death_pred_ode[,n],-1)),tyears*365/7),FUN=sum)
  incbW[,n]<-tapply(head(incb_pred_ode[,n],-1), cut(1:length(head(incb_pred_ode[,n],-1)),tyears*365/7),FUN=sum)
  incW[,n]<-tapply(head(inc_pred_ode[,n],-1), cut(1:length(head(inc_pred_ode[,n],-1)),tyears*365/7),FUN=sum)
  impdW[,n]<-tapply(head(impd_pred_ode[,n],-1), cut(1:length(head(impd_pred_ode[,n],-1)),tyears*365/7),FUN=sum)
  locdW[,n]<-tapply(head(locd_pred_ode[,n],-1), cut(1:length(head(locd_pred_ode[,n],-1)),tyears*365/7),FUN=sum)
  }
  
  COVout<-list(modeltimes, #1
          mildb_pred_ode,#2
          sevb_pred_ode,#3
          critb_pred_ode,#4
          mild_pred_ode,#5
          sev_pred_ode,#6
          crit_pred_ode,#7
          deathb_pred_ode,#8
          death_pred_ode,#9
          prev_pred_ode,#10
          incb_pred_ode,#11
          inc_pred_ode,#12
          prev_sev_ode,#13
          prev_crit_ode,#14
          inc_pred_ode,#15
          impd_pred_ode, #16
          locd_pred_ode, #17
          outoderun, #18
          tranoderun #19
  )
  
  COVoutW<-list(modelweeks, #1
                mildbW,#2
                sevbW,#3
                critbW,#4
                mildW,#5
                sevW,#6
                critW,#7
                deathbW,#8
                deathW,#9
                incbW,#10
                incW,#11
                impdW, #12
                locdW, #13
                dweeks #14
  )
  
   return(list(COVout, COVoutW))
}

dtimes<-format(as.Date(time*365,origin="2020-03-01"), "%d/%m/%Y")
dweeks<-dtimes[seq(1,length(dtimes),7)]


#### ODE  ####

st <- Sys.time()
outode <- run(pars, scenario)
en <- Sys.time()
en - st

