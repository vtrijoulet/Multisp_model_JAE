###############################################################################
######      Script that fits the multispecies model to the data sets     ######
######      Author: Vanessa Trijoulet, NEFSC, NOAA, Woods Hole, USA      ######
######      Current contact information:                                 ######
######      DTU Aqua, Lyngby, Denmark                                    ######
######      vtri@aqua.dtu.dk                                             ######
###############################################################################


require(TMB)


#### To change depending on simulation choice ##############

args <-c(1,1,1,1,1) # see below for corresponding arguments to change

############################################################

process_rec<-as.numeric(args[1]) # 1 = process error on recruitment, 0 = no process error
process_survival<-as.numeric(args[2]) # 1 = process error on numbers at age, 0 = no process error
predation_on <- as.numeric(args[3]) # 0 = predation mortality=0, 1 = predation on
loop1<-as.numeric(args[4]) # number of simulated data set to start from (from 1 to 1000)
n_loop<-as.numeric(args[5]) # number of simulated data set to stop at so n_loop-loop1+1=total number of simulation


#### Necessary scripts for running the optimization ####
compile("MS_SSM.cpp")
dyn.load(dynlib("MS_SSM"))
source("MS_SSM_make.map.fn.R") # file that creates the map argument for objective function
source("MS_SSM_make.random.fn.R") # file that creates the random argument for the objective function


#### Load data files needed ####
load("init.RData") # initial values for the optimization
load("sim_data_all.RData") # full list of 1000 simulated data sets


#### Extract data needed depending on args[4:5] ####
data_sim_loop<-sim_data_all[loop1:n_loop] # extract data set that are going to be used if you don't want to fit the 1000 data sets at once
rm(sim_data_all) # remove initial large list that takes too much memory
gc() # free memory from that list


#### Prepare necessary lists ####
name_lists<-names(data_sim_loop)
opt_all<-list()
rep_all<-list()
non_opt<-list()
non_rep<-list()
obj_all<-list()
se_all<-list()
gc_all<-list()


#### Optimization and saving the results loop ####
for (it in 1:(n_loop-loop1+1)){
  
  cat("ITERATION",it,"_", process_rec, process_survival, predation_on, loop1, n_loop, name_lists[it])
  
  # Get the specific simulated data set
  data_sim<-data_sim_loop[[name_lists[it]]]
  
  # Read in options for EM (no change needed)
  data_sim$process_rec<-process_rec
  data_sim$process_survival<-process_survival
  data_sim$predation_on<-predation_on
  data_sim$flag_nll_diet <- 3 # diet model option used in the paper
  data_sim$scale_M_upper[]<-10 # upper limit on M scaler

  
  ######## Options to change after reading data set for different EMs ##############
  
  #data_sim$recruit_model<-1 # only needed if want to estimate random walk on recruitment (EM1b)
  data_sim$M_model[]<-4 # 4 = M is estimated as single constant over time and ages (EM1-EM3), 3 = Lorenzen M (EM4)
  
  ##################################################################################
  
  
  # Prepare map, random and init arguments
  map = make.map.fn(data_sim)
  random = make.random.fn(data_sim)

  # Write the objective function and optimize it
  obj_sim<-MakeADFun(data=data_sim, parameters=init, DLL="MS_SSM", random=random, map=map)
  opt_sim<-nlminb(obj_sim$par*0.8, obj_sim$fn, obj_sim$gr, control = list(eval.max = 5000, iter.max = 5000,rel.tol=1e-7))
  
  cat("ITERATION",it,"EXTRA NEWTON STEPS")
  
  test <- try(
    for(i in seq_len(10)) { # Take 10 extra newton steps
      g <- as.numeric( obj_sim$gr(opt_sim$par) )
      h <- optimHess(opt_sim$par, obj_sim$fn, obj_sim$gr)
      opt_sim$par <- opt_sim$par - solve(h, g)
      opt_sim$objective <- obj_sim$fn(opt_sim$par)
    })
  
  if (length(test)==0){ # no error when inversing the hessian
    
  cat("ITERATION",it,"SDREPORT")
  
  sd <- sdreport(obj_sim)
  se <- as.list(sd,"Std. Error")
  se_all[[name_lists[it]]]<-se[which(names(se)%in%names(sd$par.fixed))]
  gc_all[[name_lists[it]]]<-sd$gradient.fixed

  # Save outputs into list
  
  if (sd$pdHess==TRUE && max(abs(sd$gradient.fixed))<=0.0001) { # convergence criteria
    opt_all[[name_lists[it]]]<-opt_sim
    rep_all[[name_lists[it]]]<-obj_sim$report()
    obj_all[[name_lists[it]]]<-obj_sim
  } else {
    opt_all[[name_lists[it]]]<-NA
    rep_all[[name_lists[it]]]<-NA
    obj_all[[name_lists[it]]]<-obj_sim
    non_opt[[name_lists[it]]]<-opt_sim
    non_rep[[name_lists[it]]]<-obj_sim$report()
  }
  
  } else {
    opt_all[[name_lists[it]]]<-NA
    rep_all[[name_lists[it]]]<-NA
    obj_all[[name_lists[it]]]<-obj_sim
    non_opt[[name_lists[it]]]<-opt_sim
    non_rep[[name_lists[it]]]<-obj_sim$report()
    se_all[[name_lists[it]]]<-NA
    gc_all[[name_lists[it]]]<-NA
    
  }
  
  

}


#### Save results ####
name.simul<-paste0("rec-m",data_sim$recruit_model,"_rec",data_sim$process_rec,"_surv",data_sim$process_survival,"_pred",data_sim$predation_on,"_M",data_sim$M_model[1],"_data",loop1,"-",n_loop)
save(opt_all,file=paste("sim_opt_",name.simul,".Rdata",sep=""))
save(rep_all,file=paste("sim_rep_",name.simul,".Rdata",sep=""))
save(non_opt,file=paste("non_opt_",name.simul,".Rdata",sep=""))
save(non_rep,file=paste("non_rep_",name.simul,".Rdata",sep=""))
save(obj_all,file=paste("sim_obj_",name.simul,".Rdata",sep=""))
save(se_all,file=paste("sim_se_",name.simul,".Rdata",sep=""))
save(gc_all,file=paste("sim_gc_",name.simul,".Rdata",sep=""))

