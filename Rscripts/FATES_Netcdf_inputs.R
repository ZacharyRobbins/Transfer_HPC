### FATES Variable changing

library(ncdf4)
####################### Functions #################################
assign_fates_para<-function(fates.para.values,var.name,fast.val,pft.dim.lower,pft.dim.upper,tissue.dim.lower,tissue.dim.upper){

    if(length(dim(fates.para.values))>1){ # pft & tissue
    if(pft.dim.lower==0 | tissue.dim.lower==0){ 
      stop (paste("The pft/tissue dimsion should be more than 0 for ",var.name,"!")) 
    }
    fates.para.values [pft.dim.lower:pft.dim.upper,tissue.dim.lower:tissue.dim.upper]<-fast.val 
  }else{ #only pft dimension or tissue dimension
    if(pft.dim.lower>0 & tissue.dim.lower>0){ 
      stop (paste("Error: Both pft and tissue dimsion are more than 0 for ",var.name,"!\n",'Set the pft indices (pft.start	and pft.end) to zero if only one pft is present!')) 
    }
    if(pft.dim.lower>0){ #more than 1 pft
      fates.para.values [pft.dim.lower:pft.dim.upper]<-fast.val
    }
    if(tissue.dim.lower>0){ #more than 1 tissue
      fates.para.values [tissue.dim.lower:tissue.dim.upper]<-fast.val 
    }
    if(pft.dim.lower==0 & tissue.dim.lower==0){ # 1 pft and 1 tissue
      fates.para.values<-fast.val
    }
  }    
  return (fates.para.values)
}
readwritenetcdf<-function(file,var_name,update,pftl,pfth,tisl,tish){
  fates_para=nc_open(file,write=T)
  para_values <- ncvar_get(fates_para,var_name)
  print(paste0("Original value: ",para_values))
  fates_para_values_updated<-assign_fates_para(para_values,var_name,update,pftl,pfth,tisl,tish) 
  ncvar_put(fates_para,var_name,fates_para_values_updated) 
  nc_close(fates_para)
}
createnewcdf<-function(Drive,file,newname,varDF){
  file.copy(paste0(Drive,file),paste0(Drive,newname),overwrite = T)
  ### Updating the values 
  for(i in 1:length(Var_2_change$Name)){
    ####Writing in values 
    readwritenetcdf(paste0(Drive,newname),Var_2_change$Name[i],Var_2_change$Value[i],Var_2_change$pftl[i],
                    Var_2_change$pfth[i],Var_2_change$tisl[i],Var_2_change$tish[i])
    ####Checking 
    var_name=Var_2_change$Name[i]
    fates_para=nc_open(paste0(Drive,newname),write=F)
    para_values <- ncvar_get(fates_para,var_name)
    print(paste0("Updated value: ",para_values))
  }
}
###################################################################
#file='C:/Users/345578/Documents/GitHub/Transfer_HPC/fates_params_CA_StuntRanch.nc'

#"fates_fire_fdi_alpha","fates_fire_fdi_b"
fates_para=nc_open(paste0(Drive,file))
#fates_para=file
names(fates_para$var)
var_name="fates_fire_nignitions"
fates_para$var[[var_name]]$unit
fates_para$var[[var_name]]$longname
"fates_pftname" 


ncvar_get(fates_para,"fates_fire_fdi_a")
# 0.0082 from data
####### Building update dataframe
Var_2_change<-data.frame(
  Name=c("fates_fire_fdi_a",
         "fates_fire_fdi_b"),
  Value=c(12,
          86),
  pftl=c(0,
         0),
  pfth=c(0,
         0),
  tisl=c(0,
         0),
  tish=c(0,
         0)
)



#filename1<-nfname.in
#filename2<-paste(nfname.out.base,i,".nc",sep="")


#### File that is the baseline
Drive='C:/Users/345578/Documents/GitHub/Transfer_HPC/'
file='fates_params_CA_StuntRanch.nc'
newname="fates_params_CA_SR_update2.nc"


createnewcdf(Drive,file,newname,Var_2_change)
