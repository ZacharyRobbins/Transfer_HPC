# -*- coding: utf-8 -*-
"""
Created on Thu Jun  9 12:57:04 2022

@author: 345578
"""


import numpy as np
from numpy import *
import sys
import getopt
import code  # For development: code.interact(local=locals())
from datetime import datetime
from scipy.io import netcdf
import os
#import matplotlib.pyplot as plt
os.chdir("C:/Users/345578/Documents/GitHub/Transfer_HPC/")
def main(input_fname,output_fname,donor_pft_indices):

    # Interpret the arguments to the script
  
    num_pft_out = len(donor_pft_indices)

    # Open the netcdf files
    fp_out = netcdf.netcdf_file(output_fname, 'w')

    fp_in  = netcdf.netcdf_file(input_fname, 'r')

    for key, value in sorted(fp_in.dimensions.items()):
        if(key==pft_dim_name):
            fp_out.createDimension(key,int(num_pft_out))
            print('Creating Dimension: {}={}'.format(key,num_pft_out))
        else:
            fp_out.createDimension(key,int(value))
            print('Creating Dimension: {}={}'.format(key,value))

    for key, value in sorted(fp_in.variables.items()):
        print('Creating Variable: ',key)
        #   code.interact(local=locals())


        in_var  = fp_in.variables.get(key)


        # Idenfity if this variable has pft dimension
        pft_dim_found = -1
        prt_dim_found = -1
        hydro_dim_found = -1
        litt_dim_found  = -1
        string_dim_found = -1
        pft_dim_len   = len(fp_in.variables.get(key).dimensions)

        for idim, name in enumerate(fp_in.variables.get(key).dimensions):

            # Manipulate data
            if(name==pft_dim_name):
                pft_dim_found = idim
            if(name==prt_dim_name):
                prt_dim_found = idim
            if(name==litt_dim_name):
                litt_dim_found = idim
            if(name==hydro_dim_name):
                hydro_dim_found = idim
            if(name==string_dim_name):
                string_dim_found = idim

        # Copy over the input data
        # Tedious, but I have to permute through all combinations of dimension position
        if( pft_dim_len == 0 ):
            out_var = fp_out.createVariable(key,'d',(fp_in.variables.get(key).dimensions))
            out_var.assignValue(float(fp_in.variables.get(key).data))
        elif( (pft_dim_found==-1) & (prt_dim_found==-1) & (litt_dim_found==-1) & (hydro_dim_found==-1)  ):
            out_var = fp_out.createVariable(key,'d',(fp_in.variables.get(key).dimensions))
            out_var[:] = in_var[:]
        elif( (pft_dim_found==0) & (pft_dim_len==1) ):           # 1D fates_pft
            out_var = fp_out.createVariable(key,'d',(fp_in.variables.get(key).dimensions))
            tmp_out  = np.zeros([num_pft_out])
            for id,ipft in enumerate(donor_pft_indices):
                tmp_out[id] = fp_in.variables.get(key).data[ipft-1]
            out_var[:] = tmp_out

        # 2D   hydro_organ - fates_pft
        # or.. prt_organ - fates_pft
        elif( (pft_dim_found==1) & (pft_dim_len==2) ):
            out_var = fp_out.createVariable(key,'d',(fp_in.variables.get(key).dimensions))
            dim2_len = fp_in.dimensions.get(fp_in.variables.get(key).dimensions[0])
            tmp_out  = np.zeros([dim2_len,num_pft_out])
            for id,ipft in enumerate(donor_pft_indices):
                for idim in range(0,dim2_len):
                    tmp_out[idim,id] = fp_in.variables.get(key).data[idim,ipft-1]
            out_var[:] = tmp_out

        elif( (pft_dim_found==0) & (pft_dim_len==2) ):          # fates_pft - string_length
            out_var = fp_out.createVariable(key,'c',(fp_in.variables.get(key).dimensions))
            dim2_len = fp_in.dimensions.get(fp_in.variables.get(key).dimensions[1])
            out_var[:] = np.empty([num_pft_out,dim2_len], dtype="S{}".format(dim2_len))
            for id,ipft in enumerate(donor_pft_indices):
                out_var[id] = fp_in.variables.get(key).data[ipft-1]
                
        elif( (prt_dim_found==0) & (pft_dim_len==2) ):         
            out_var = fp_out.createVariable(key,'c',(fp_in.variables.get(key).dimensions))
            out_var[:] = in_var[:]

        elif( (hydro_dim_found==0) & (string_dim_found>=0) ):      
            out_var = fp_out.createVariable(key,'c',(fp_in.variables.get(key).dimensions))
            out_var[:] = in_var[:]
        
        elif( (litt_dim_found==0) & (string_dim_found>=0) ):       
            out_var = fp_out.createVariable(key,'c',(fp_in.variables.get(key).dimensions))
            out_var[:] = in_var[:]   

        elif( prt_dim_found==0 ): # fates_prt_organs - indices
            out_var = fp_out.createVariable(key,'d',(fp_in.variables.get(key).dimensions))
            out_var[:] = in_var[:]

        elif( litt_dim_found==0 ):
            out_var = fp_out.createVariable(key,'d',(fp_in.variables.get(key).dimensions))
            out_var[:] = in_var[:]
        elif( hydro_dim_found==0):
            out_var = fp_out.createVariable(key,'d',(fp_in.variables.get(key).dimensions))
            out_var[:] = in_var[:]
        else:
            print('This variable has a dimensioning that we have not considered yet.')
            print('Please add this condition to the logic above this statement.')
            print('Aborting')
            for idim, name in enumerate(fp_in.variables.get(key).dimensions):
               print("idim: {}, name: {}".format(idim,name))
            exit(2)

        out_var.units     = in_var.units
        out_var.long_name = in_var.long_name

    fp_out.history = "This file was made from FatesPFTIndexSwapper.py \n Input File = {} \n Indices = {}"\
          .format(input_fname,donor_pft_indices)

    #var_out.mode = var.mode
    #fp.flush()

    fp_in.close()
    fp_out.close()

    print('Cloneing complete!')



pft_dim_name = 'fates_pft'
prt_dim_name = 'fates_prt_organs'
hydro_dim_name = 'fates_hydr_organs'
litt_dim_name = 'fates_litterclass'
string_dim_name = 'fates_string_length'

main('fates_params_CA_StuntRanch.nc','testonepft2.nc',[1,1,1,1,1,1])
