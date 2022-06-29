Var_2_change<-data.frame(
  Name=c("fates_fire_fdi_a",
         "fates_fire_fdi_b",
         "fates_fire_nignitions",
         "fates_fire_SAV",
         "fates_CWD_frac",
         "fates_fire_FBD" ),
  Value=c(17.62,
          243.12,
          2.25,
          list(13,3.58,0.98 ,0.2,100,66),
          list( 0.09,0.15,0.42,0.34),
          list(15.4,16.8,19.6,999,6,4)),
  pftl=c(0,0,0,0,0),
  pfth=c(0,0,0,0,0),
  tisl=c(0,0,0,0,0),
  tish=c(0,0,0,0,0)
)



Drive='C:/Users/345578/Documents/GitHub/Transfer_HPC/'
file='fates_params_CA_StuntRanch.nc'
newname="fates_params_CA_SR_6_28.nc"
file.copy(paste0(Drive,file),paste0(Drive,newname),overwrite = T)
fates_para=nc_open(paste0(Drive,newname),write=T)
ncvar_put(fates_para,"fates_fire_fdi_a",17.62) 
ncvar_put(fates_para,"fates_fire_fdi_b",243.12) 
ncvar_put(fates_para,"fates_fire_nignitions",2.25) 
ncvar_put(fates_para,"fates_fire_SAV",c(13,3.58,0.98 ,0.2,100,66)) 
ncvar_put(fates_para,"fates_CWD_frac",c(0.09,0.15,0.42,0.34)) 
ncvar_put(fates_para,"fates_fire_FBD",c(15.4,16.8,19.6,999,6,4)) 
nc_close(fates_para)


