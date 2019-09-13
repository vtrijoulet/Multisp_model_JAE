make.random.fn = function(data)
{

  if (data$process_rec==1 & data$process_survival==1){
    random<-c("log_rec","log_NAA")
  }
  if (data$process_rec==1 & data$process_survival==0){
    random<-c("log_rec")
  }
  if (data$process_rec==0 & data$process_survival==1){
    random<-c("log_NAA")
  }
  if (data$process_rec==0 & data$process_survival==0){
    random<-NULL
  }
  if (data$process_M==1){
    if (data$M_model[1]==2)
      random<-c(random,"log_MAA")
    if (data$M_model[1]==3)
      random<-c(random,"log_lorenzen1","lorenzen2")
  }
  
  return(random)
}

