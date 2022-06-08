######################################################################
#this script generates FATES parameter files from the FAST sampling
#######################################################################
#Special notices:
#set the pft indices (pft.start	and pft.end) to zero if only one pft is present
#set the organ indicees if no organ indices is present
# install.packages('ncdf4',repos="http://cran.us.r-project.org")
library(ncdf4)
#a conversion function
conversion<-function(y,l,u){
  x<-l+y*(u-l)
}
#assign values for the ncfile
assign.fates.para<-function(fates.para.values,fast.val,pft.dim.lower,pft.dim.upper,tissue.dim.lower,tissue.dim.upper){
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
########################################################
#user inputs
#set the working space, with ./FASTPara for store fates parameter files
#setwd("/turquoise/usr/projects/veg/cxu/ACME_cases/BCI/params/FASTParams")
setwd("/usr/projects/veg/zjrobbins/ACME_cases/BCI/params")
# name of your parameter files from FAST
fast.para.name ="/usr/projects/veg/zjrobbins/Fates_Sen/new_HYDRAULICTRAITS.FAST.Sam" 
#the input fates parameter file name
#nfname.in = "fates_params_default_056193a_mod1PFTs_pureEDbareground_exp4.nc"  
nfname.in = "/usr/projects/veg/zjrobbins/ACME_cases/BCI/params/fates_params_default_1trop_bciopt224.c201022.nc"  
#the start of sample index
sample.start.index = 1
#the end of sample index
sample.end.index = 100
#=======================================================
# read the fast parameter file
fast.para = read.table(fast.para.name,header=T,sep="\t")
fast.para <- fast.para[,1:(ncol(fast.para)-1)] #remove the last column of NAs
#================================================
#read the paramater settings
input.para.name = "/usr/projects/veg/zjrobbins/ACME_cases/BCI/params/ParameterSettings.csv"
input.para = read.table(input.para.name,header=T,sep=",")
n.par<-nrow(input.para)
para.OI.fates<-input.para[,'parameter']
pft.dim.lower<-input.para[,'pft.start']
pft.dim.upper<-input.para[,'pft.end']
tissue.dim.lower<-input.para[,'org.start']
tissue.dim.upper<-input.para[,'org.end']
reverse.flag<-input.para[,'reverse.flag']
reciprocal.flag<-input.para[,'reciprocal.flag']
conversion.flag<-input.para[,'conversion.flag']
conversion.lower<-input.para[,'conversion.lower']
conversion.upper<-input.para[,'conversion.upper']
scale<-input.para[,'scale']
shift<-input.para[,'shift']
#========================================================
#check the validity of parameter setting
if(n.par!=ncol(fast.para)){
  stop ("The FATES parameter setting does not match the FATES input file") 
}
#================================================
#test the parameter values
#setwd("./FASTParams")
test.flag=T
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
#create the files
n.sam<-nrow(fast.para)
sample.end.index = min(sample.end.index,n.sam)
g <- gregexpr(".", nfname.in, fixed=TRUE)
loc <- g[[1]]
file.ext.loc<-loc[length(loc)]
nfname.out.base<-substr(nfname.in,1,file.ext.loc)
for(i in sample.start.index:sample.end.index){
  filename1<-nfname.in
  filename2<-paste(nfname.out.base,i,".nc",sep="")
  file.copy(filename1,filename2,overwrite = T)
}
#-------------------------------------------------
#set the parameter value
pb <- txtProgressBar(min =sample.start.index, max = sample.end.index, style = 3)
for(i in sample.start.index:sample.end.index){
  #filename2<-paste("fates_params_2troppftclones.c171018.hydro",i,".nc",sep="")
  setTxtProgressBar(pb, i)
  filename2<-paste(nfname.out.base,i,".nc",sep="")
  print(filename2)
  fates_para <-nc_open(filename2, write=T)
  for(j in 1:n.par){   
    var.name<-as.character(para.OI.fates[j])
    if(!is.na(var.name)){
      fast.val<-(fast.para[i,j]+shift[j])*scale[j]
      if(reverse.flag[j])fast.val=-fast.val;
      if(reciprocal.flag[j])fast.val=1.0/fast.val;
      if(conversion.flag[j])fast.val=conversion(fast.val,conversion.lower[j],conversion.upper[j])
      fates.para.values <- ncvar_get(fates_para,var.name)
      fates.para.values.updated<-assign.fates.para(fates.para.values,fast.val,pft.dim.lower[j],pft.dim.upper[j],tissue.dim.lower[j],tissue.dim.upper[j])
      ncvar_put(fates_para,var.name,fates.para.values.updated) 
      if(var.name == "fates_hydr_p50_gs"){
        var.name = "fates_hydr_avuln_gs"
        fast.val.avuln = -4*60.15*(-fast.val)^(-1.25)/100*fast.val
        fates.para.values <- ncvar_get(fates_para,var.name)
        fates.para.values.updated<-assign.fates.para(fates.para.values,fast.val.avuln, pft.dim.lower[j],pft.dim.upper[j],tissue.dim.lower[j],tissue.dim.upper[j])
        ncvar_put(fates_para,var.name,fates.para.values.updated) 
      } 
      if(var.name == "fates_hydr_pitlp_node"){  
        var.name = "fates_hydr_epsil_node"
        fates.epsil.values <- ncvar_get(fates_para,var.name)
        if(length(dim(fates.para.values))>1){
          epsil = mean(fates.epsil.values[pft.dim.lower[j]:pft.dim.upper[j],tissue.dim.lower[j]:tissue.dim.upper[j]])
        }else{
          epsil = mean(fates.epsil.values[tissue.dim.lower[j]:tissue.dim.upper[j]])
        }
        fast.val.pinot =  fast.val *epsil/(epsil-fast.val)
        var.name = "fates_hydr_pinot_node"
        fates.para.values <- ncvar_get(fates_para,var.name)
        fates.para.values.updated<-assign.fates.para(fates.para.values,fast.val.pinot, pft.dim.lower[j],pft.dim.upper[j],tissue.dim.lower[j],tissue.dim.upper[j])
        ncvar_put(fates_para,var.name,fates.para.values.updated) 
      }
      if(var.name == "fates_hydr_resid_node"){  
        var.name.thetas = "fates_hydr_thetas_node"
        fates.thetas.values <- ncvar_get(fates_para,var.name.thetas)
        if(length(dim(fates.para.values))>1){
            thetas = mean(fates.thetas.values[pft.dim.lower[j]:pft.dim.upper[j],tissue.dim.lower[j]:tissue.dim.upper[j]])
        }else{
            thetas = mean(fates.thetas.values[tissue.dim.lower[j]:tissue.dim.upper[j]])
        }
        fast.val.resid =  fast.val * thetas
        fates.para.values.updated<-assign.fates.para(fates.para.values,fast.val.resid, pft.dim.lower[j],pft.dim.upper[j],tissue.dim.lower[j],tissue.dim.upper[j])
        ncvar_put(fates_para,var.name,fates.para.values.updated) 
      }
    } ## var.name!="NA"
  } ##j    
  nc_close(fates_para)
} #i
