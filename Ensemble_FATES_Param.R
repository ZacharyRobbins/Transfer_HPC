library(ncdf4)
######################################################
#conversion functions
conversion<-function(y,l,u){
  x<-l+y*(u-l)
}
#================================================================
#input/output section
#setwd("C:/Users/242621/Documents/FAST")
input.para.name = "params.input.ens.csv"
input.para = read.table(input.para.name,header=T,sep=",")
nfname.in = "./params/fates_params_default_c20200213.saltmarsh.update.nc" #input netcdf file nmae
nfname.out.base <- "./params/fates_params_default_c20200213.saltmarsh.wuma"
n.ens = 10
#================================================
#test the parameter values
test.flag=F
if(test.flag){
   fates_para <-nc_open(nfname.in, write=F)
   list_var <- names(fates_para$var)
   #write.table(file="varname.txt",list_var)
   var.name<-"fates_hydr_avuln_gs"
   fates_para$var[[var.name]]$unit
   fates_para$var[[var.name]]$longname
   para_values <- ncvar_get(fates_para,var.name)
   para_values
   nc_close(fates_para )
   tt<-grepl('d2ca',list_var) #found the right variable
   list_var[tt]
}
#================================================
#set the dimension to update
n.par<-nrow(input.para)
for (i in 1:n.ens){
 #set the parameter value
 #================================================
 #create the files
 nfname.out = paste(nfname.out.base,".ens.",i, ".nc", sep="")#output netcdf file name
 file.copy(nfname.in,nfname.out,overwrite = T)
 fates_para <-nc_open(nfname.out, write=T)
 #------------------------------------------------- 
 for(j in 1:n.par){ 
    para.OI.fates<-input.para[,'parameter']
    pft.dim.lower<-input.para[,'pft.start']
    pft.dim.upper<-input.para[,'pft.end']
    tissue.dim.lower<-input.para[,'org.start']
    tissue.dim.upper<-input.para[,'org.end']   
    var.name<-as.character(para.OI.fates[j])
    if(var.name!="NA"){
      valuemin<-input.para[j,'valuemin']
      valuemax<-input.para[j,'valuemax']
      para.val<-runif(1,valuemin,valuemax)
      fates.para.values <- ncvar_get(fates_para,var.name) #read trait values from the netcdf file
      if(pft.dim.lower[j]>0){
        if (tissue.dim.lower[j]==0){
          fates.para.values [pft.dim.lower[j]:pft.dim.upper[j]]<-para.val
        }else{
          fates.para.values [pft.dim.lower[j]:pft.dim.upper[j],tissue.dim.lower[j]:tissue.dim.upper[j]]<-para.val 
        }
      }else{
        fates.para.values<-para.val
      }
      ncvar_put(fates_para,var.name,fates.para.values ) 
     } # if(var.name!="NA"){
  } #j
  nc_close(fates_para)
}

