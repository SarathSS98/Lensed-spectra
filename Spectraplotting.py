import numpy as np
import matplotlib.pyplot as plt
import os
from astropy.visualization import astropy_mpl_style
from numpy.polynomial.polynomial import polyfit
from scipy.optimize import curve_fit
from astropy.utils.data import get_pkg_data_filename
from astropy.io import fits
from scipy import optimize
import pandas as pd
import math
import statistics
from scipy.ndimage import gaussian_filter1d 
import csv
import sys


def gauss4(x,c,a,ga,gb,gc,ga1,gb1,gc1,ga2,gb2,gc2,ga3,gb3,gc3):
    return ( c+a*x+ga*np.exp(-((x-gb)/gc)**2 )+ga1*np.exp(-((x-gb1)/gc1)**2)+ga2*np.exp(-((x-gb2)/gc2)**2 )
           +ga3*np.exp(-((x-gb3)/gc3)**2 ))
def gauss3(x,c,a,ga,gb,gc,ga1,gb1,gc1,ga2,gb2,gc2):
    return ( c+a*x+ga*np.exp(-((x-gb)/gc)**2 )+ga1*np.exp(-((x-gb1)/gc1)**2)+ga2*np.exp(-((x-gb2)/gc2)**2 ))
def gauss2(x,c,a,ga,gb,gc,ga1,gb1,gc1):
    return ( c+a*x+ga*np.exp(-((x-gb)/gc)**2 )+ga1*np.exp(-((x-gb1)/gc1)**2))
def gauss1(x,c,a,ga,gb,gc):
     return ( c+a*x+ga*np.exp(-((x-gb)/gc)**2 ))


path='/media/sarath/DATA/UNAB-PHD/Semester-1/Project-lensed Quasars/new_lenses/Spectra for plotting/'













DF=pd.read_csv('NTT.csv',header=0)
Spectra=DF.Spectra
z=DF.z
name=DF.Name
for i in range(0,len(DF)): 

    source=np.load(path+Spectra[i],allow_pickle=True)  
    print(Spectra[i])
    z=DF.z[i]
    print(z)
    wavelength=source[2]
    flux=source[0]
    bright_flux=flux1=flux[0]
    faint_flux=flux2=flux[1]
    fig, ax=plt.subplots(2,1,figsize=(10,10))
###BRIGHT SPECTRUM---------------------------------------------------------------------------------------
    x=wavelength
    pos=np.where((x >3760) & (x <8900)) 
    X=wavelength[pos]
    
    x=X/(z+1)
    y1=flux1[pos]
    x= np.nan_to_num(x)
    y1= np.nan_to_num(y1)
    #print(len(y1))
    n=len(x)
    mean=sum(x*y1)/n
    gc=sum(y1*(x-mean)**2)/n
    gc1=gc2=gc3=gc

#y1=gaussian_filter1d(y1,0.1)## This value should be adjusted
    med=statistics.median(y1)
    if ((z>=2) and (z<=2.178)):
    
        p0=[med,0,10,1215,1,10,1549,1,10,1908,1,10,2800,1] #change this when : Optimal parameters not found:
                                             #Number of calls to function has reached maxfev = 1800.
        popt1,pcov1= curve_fit(gauss4,x,y1,p0=p0)
        FWHM_LyA=round(2*np.sqrt(2*np.log(2))*popt1[4])
        FWHM_CIV=round(2*np.sqrt(2*np.log(2))*popt1[7])
        FWHM_CIII=round(2*np.sqrt(2*np.log(2))*popt1[10])
        FWHM_MgII=round(2*np.sqrt(2*np.log(2))*popt1[13])


        LyA_low=(popt1[0]+popt1[1]*popt1[3])
        LyA_high=(gauss4(popt1[3],*popt1))    
        LyA_dif=LyA_high-LyA_low

        CIV_low=(popt1[0]+popt1[1]*popt1[6])
        CIV_high=(gauss4(popt1[6],*popt1))    
        CIV_dif=CIV_high-CIV_low
        
        CIII_low=(popt1[0]+popt1[1]*popt1[9])
        CIII_high=(gauss4(popt1[9],*popt1))    
        CIII_dif=CIII_high-CIII_low
    
        MgII_low=(popt1[0]+popt1[1]*popt1[12])
        MgII_high=(gauss4(popt1[12],*popt1))    
        MgII_dif=MgII_high-MgII_low
                   
       # print('name[i],LyA_dif,CIV_dif,CIII_dif,MGII_dif,FWHM_LyA,FWHM_CIV,FWHM_CIII,FWHM_MGII')
        value=[name[i],LyA_dif,CIV_dif,CIII_dif,MGII_dif,FWHM_LyA,FWHM_CIV,FWHM_CIII,FWHM_MGII]
   # print(value)
        ax[0].plot(x,y1,label='redshifted data')
        ax[0].plot(x, gauss4(x, *popt1),'r', label='guassianfit')
        ax[0].set_xlabel('Wavelength')
        ax[0].set_ylabel('Flux')       
        ax[0].legend()
        
        
        
        
        
    elif ((z>=1.42) and (z<=2.178)):
        p0=[med,0,10,1549,1,10,1908,1,10,2800,1]
        popt1,pcov1= curve_fit(gauss3,x,y1,p0=p0)
        
        
        FWHM_CIV=round(2*np.sqrt(2*np.log(2))*popt1[4])
        FWHM_CIII=round(2*np.sqrt(2*np.log(2))*popt1[7])
        FWHM_MgII=round(2*np.sqrt(2*np.log(2))*popt1[10])



        CIV_low=(popt1[0]+popt1[1]*popt1[3])
        CIV_high=(gauss3(popt1[3],*popt1))    
        CIV_dif=CIV_high-CIV_low
    
        CIII_low=(popt1[0]+popt1[1]*popt1[6])
        CIII_high=(gauss3(popt1[6],*popt1))    
        CIII_dif=CIII_high-CIII_low
    
        MgII_low=(popt1[0]+popt1[1]*popt1[9])
        MgII_high=(gauss3(popt1[9],*popt1))    
        MgII_dif=MgII_high-MgII_low
               
        print('name[i],LyA_dif,CIV_dif,CIII_dif,MGII_dif,FWHM_LyA,FWHM_CIV,FWHM_CIII,FWHM_MGII')
        value=[name[i],0,CIV_dif,CIII_dif,MGII_dif,0,FWHM_CIV,FWHM_CIII,FWHM_MGII]
   # print(value)
        ax[0].plot(x,y1,label='redshifted data')
        ax[0].plot(x, gauss3(x, *popt1),'r', label='guassianfit')
        ax[0].set_xlabel('Wavelength')
        ax[0].set_ylabel('Flux')       
        ax[0].legend()
        
        
    
    elif ((z>=2.092) and (z<=3.664)):
        p0=[med,0,10,1215,1,10,1549,1,10,1908,1]
        popt1,pcov1= curve_fit(gauss3,x,y1,p0=p0)
        
        FWHM_LyA=round(2*np.sqrt(2*np.log(2))*popt1[4])
        FWHM_CIV=round(2*np.sqrt(2*np.log(2))*popt1[7])
        FWHM_CIII=round(2*np.sqrt(2*np.log(2))*popt1[10])


        LyA_low=(popt1[0]+popt1[1]*popt1[3])
        LyA_high=(gauss3(popt1[3],*popt1))    
        LyA_dif=LyA_high-LyA_low

        CIV_low=(popt1[0]+popt1[1]*popt1[6])
        CIV_high=(gauss3(popt1[6],*popt1))    
        CIV_dif=CIV_high-CIV_low
    
        CIII_low=(popt1[0]+popt1[1]*popt1[9])
        CIII_high=(gauss3(popt1[9],*popt1))    
        CIII_dif=CIII_high-CIII_low
    
        print('name[i],LyA_dif,CIV_dif,CIII_dif,MGII_dif,FWHM_LyA,FWHM_CIV,FWHM_CIII,FWHM_MGII')
        value=[name[i],LyA_dif,CIV_dif,CIII_dif,0,FWHM_LyA,FWHM_CIV,FWHM_CIII,0]
   # print(value)
        ax[0].plot(x,y1,label='redshifted data')
        ax[0].plot(x, gauss3(x, *popt1),'r', label='guassianfit')
        ax[0].set_xlabel('Wavelength')
        ax[0].set_ylabel('Flux')       
        ax[0].legend()
        
    elif ((z>=0.9) and (z<=2.178)):
        p0=[med,0,10,1908,1,10,2800,1]
        popt1,pcov1= curve_fit(gauss2,x,y1,p0=p0)
        
        
        FWHM_CIII=round(2*np.sqrt(2*np.log(2))*popt1[4])
        FWHM_MgII=round(2*np.sqrt(2*np.log(2))*popt1[7])


        CIII_low=(popt1[0]+popt1[1]*popt1[3])
        CIII_high=(gauss2(popt1[3],*popt1))    
        CIII_dif=CIII_high-CIII_low

        MgII_low=(popt1[0]+popt1[1]*popt1[6])
        MgII_high=(gauss2(popt1[6],*popt1))    
        MgII_dif=MgII_high-MgII_low
    
        
    
        print('name[i],LyA_dif,CIV_dif,CIII_dif,MGII_dif,FWHM_LyA,FWHM_CIV,FWHM_CIII,FWHM_MGII')
        value=[name[i],0,0,CIII_dif,MgII_dif,0,0,FWHM_CIII,FWHM_MgII]
    #print(value)
        ax[0].plot(x,y1,label='redshifted data')
        ax[0].plot(x, gauss2(x, *popt1),'r', label='guassianfit')
        ax[0].set_xlabel('Wavelength')
        ax[0].set_ylabel('Flux')       
        ax[0].legend()
        
        
        
        
    elif ((z>=1.42) and (z<=3.6)):
        p0=[med,0,10,1549,1,10,1908,1]
        popt1,pcov1= curve_fit(gauss2,x,y1,p0=p0)
        
        FWHM_CIV=round(2*np.sqrt(2*np.log(2))*popt1[4])
        FWHM_CIII=round(2*np.sqrt(2*np.log(2))*popt1[7])


        CIV_low=(popt1[0]+popt1[1]*popt1[3])
        CIV_high=(gauss2(popt1[3],*popt1))    
        CIV_dif=CIV_high-CIV_low

        CIII_low=(popt1[0]+popt1[1]*popt1[6])
        CIII_high=(gauss2(popt1[6],*popt1))    
        CIII_dif=CIII_high-CIII_low
    
        
    
        print('name[i],LyA_dif,CIV_dif,CIII_dif,MGII_dif,FWHM_LyA,FWHM_CIV,FWHM_CIII,FWHM_MGII')
        value=[name[i],0,CIV_dif,CIII_dif,0,0,FWHM_CIV,FWHM_CIII,0]
   # print(value)
        ax[0].plot(x,y1,label='redshifted data')
        ax[0].plot(x, gauss2(x, *popt1),'r', label='guassianfit')
        ax[0].set_xlabel('Wavelength')
        ax[0].set_ylabel('Flux')       
        ax[0].legend()
    
    elif ((z>=0.34) and (z<=2.178)):
        p0=[med,0,10,2800,1]
        popt1,pcov1= curve_fit(gauss1,x,y1,p0=p0)
        
          
        FWHM_MgII=round(2*np.sqrt(2*np.log(2))*popt1[4])


        MgII_low=(popt1[0]+popt1[1]*popt1[3])
        MgII_high=(gauss1(popt1[3],*popt1))    
        MgII_dif=MgII_high-MgII_low
    
    
        print('name[i],LyA_dif,CIV_dif,CIII_dif,MGII_dif,FWHM_LyA,FWHM_CIV,FWHM_CIII,FWHM_MGII')
        value=[name[i],0,0,0,MgII_dif,0,0,0,FWHM_MgII]
       # print(value)
        ax[0].plot(x,y1,label='redshifted data')
        ax[0].plot(x, gauss1(x, *popt1),'r', label='guassianfit')
        ax[0].set_xlabel('Wavelength')
        ax[0].set_ylabel('Flux')       
        ax[0].legend()

# with open('CSVFILE_bright.csv', 'a', newline='',) as f_object:
#     writer_object = csv.writer(f_object)
#     writer_object.writerow(value)  
#     f_object.close()
    
#perr1 = np.sqrt(np.diag(pcov1))
#ax[0].plot(x,y1,label='redshifted data')
#ax[0].plot(x, gauss(x, *popt1),'r', label='guassianfit')
#plt.ylim(-500,600)    

# ax[0].set_ylim(-500,500)

##FAINT FLUX---------------------------------------------------------------------

    y2=flux2[pos]
    x= np.nan_to_num(x)
    y2= np.nan_to_num(y2)

    n=len(x)
    mean=sum(x*y2)/n
    gc=sum(y2*(x-mean)**2)/n
    gc1=gc2=gc3=gc
#y2=gaussian_filter1d(y2,0.1)## This value should be adjusted
    med=statistics.median(y2)

    if ((z>=2) and (z<=2.178)):
    
        p0=[med,0,10,1215,1,10,1549,1,10,1908,1,10,2800,1] #change this when : Optimal parameters not found:
                                                 #Number of calls to function has reached maxfev = 1800.
        popt1,pcov1= curve_fit(gauss4,x,y2,p0=p0)
        FWHM_LyA=round(2*np.sqrt(2*np.log(2))*popt1[4])
        FWHM_CIV=round(2*np.sqrt(2*np.log(2))*popt1[7])
        FWHM_CIII=round(2*np.sqrt(2*np.log(2))*popt1[10])
        FWHM_MgII=round(2*np.sqrt(2*np.log(2))*popt1[13])


        LyA_low=(popt1[0]+popt1[1]*popt1[3])
        LyA_high=(gauss4(popt1[3],*popt1))    
        LyA_dif=LyA_high-LyA_low

        CIV_low=(popt1[0]+popt1[1]*popt1[6])
        CIV_high=(gauss4(popt1[6],*popt1))    
        CIV_dif=CIV_high-CIV_low
    
        CIII_low=(popt1[0]+popt1[1]*popt1[9])
        CIII_high=(gauss4(popt1[9],*popt1))    
        CIII_dif=CIII_high-CIII_low
    
        MgII_low=(popt1[0]+popt1[1]*popt1[12])
        MgII_high=(gauss4(popt1[12],*popt1))    
        MgII_dif=MgII_high-MgII_low
               
        print('name[i],LyA_dif,CIV_dif,CIII_dif,MGII_dif,FWHM_LyA,FWHM_CIV,FWHM_CIII,FWHM_MGII')
        value1=[name[i],LyA_dif,CIV_dif,CIII_dif,MGII_dif,FWHM_LyA,FWHM_CIV,FWHM_CIII,FWHM_MGII]
        #print(value1)
        ax[1].plot(x,y2,label='redshifted data')
        ax[1].plot(x, gauss4(x, *popt1),'r', label='guassianfit')
        ax[1].set_xlabel('Wavelength')
        ax[1].set_ylabel('Flux')       
        ax[1].legend()
        
        
        
        
        
    elif ((z>=1.42) and (z<=2.178)):
        p0=[med,0,10,1549,1,10,1908,1,10,2800,1]
        popt1,pcov1= curve_fit(gauss3,x,y2,p0=p0)
        
        
        FWHM_CIV=round(2*np.sqrt(2*np.log(2))*popt1[4])
        FWHM_CIII=round(2*np.sqrt(2*np.log(2))*popt1[7])
        FWHM_MgII=round(2*np.sqrt(2*np.log(2))*popt1[10])



        CIV_low=(popt1[0]+popt1[1]*popt1[3])
        CIV_high=(gauss3(popt1[3],*popt1))    
        CIV_dif=CIV_high-CIV_low
    
        CIII_low=(popt1[0]+popt1[1]*popt1[6])
        CIII_high=(gauss3(popt1[6],*popt1))    
        CIII_dif=CIII_high-CIII_low
    
        MgII_low=(popt1[0]+popt1[1]*popt1[9])
        MgII_high=(gauss3(popt1[9],*popt1))    
        MgII_dif=MgII_high-MgII_low
               
        print('name[i],LyA_dif,CIV_dif,CIII_dif,MGII_dif,FWHM_LyA,FWHM_CIV,FWHM_CIII,FWHM_MGII')
        value1=[name[i],0,CIV_dif,CIII_dif,MGII_dif,0,FWHM_CIV,FWHM_CIII,FWHM_MGII]
       # print(value1)
        ax[1].plot(x,y2,label='redshifted data')
        ax[1].plot(x, gauss3(x, *popt1),'r', label='guassianfit')
        ax[1].set_xlabel('Wavelength')
        ax[1].set_ylabel('Flux')       
        ax[1].legend()
        
        
    
    elif ((z>=2.092) and (z<=3.664)):
        p0=[med,0,10,1215,1,10,1549,1,10,1908,1]
        popt1,pcov1= curve_fit(gauss3,x,y2,p0=p0)
        
        FWHM_LyA=round(2*np.sqrt(2*np.log(2))*popt1[4])
        FWHM_CIV=round(2*np.sqrt(2*np.log(2))*popt1[7])
        FWHM_CIII=round(2*np.sqrt(2*np.log(2))*popt1[10])


        LyA_low=(popt1[0]+popt1[1]*popt1[3])
        LyA_high=(gauss3(popt1[3],*popt1))    
        LyA_dif=LyA_high-LyA_low

        CIV_low=(popt1[0]+popt1[1]*popt1[6])
        CIV_high=(gauss3(popt1[6],*popt1))    
        CIV_dif=CIV_high-CIV_low
        
        CIII_low=(popt1[0]+popt1[1]*popt1[9])
        CIII_high=(gauss3(popt1[9],*popt1))    
        CIII_dif=CIII_high-CIII_low
    
        print('name[i],LyA_dif,CIV_dif,CIII_dif,MGII_dif,FWHM_LyA,FWHM_CIV,FWHM_CIII,FWHM_MGII')
        value1=[name[i],LyA_dif,CIV_dif,CIII_dif,0,FWHM_LyA,FWHM_CIV,FWHM_CIII,0]
   # print(value1)
        ax[1].plot(x,y2,label='redshifted data')
        ax[1].plot(x, gauss3(x, *popt1),'r', label='guassianfit')
        ax[1].set_xlabel('Wavelength')
        ax[1].set_ylabel('Flux')       
        ax[1].legend()
        
    elif ((z>=0.9) and (z<=2.178)):
        p0=[med,0,10,1908,1,10,2800,1]
        popt1,pcov1= curve_fit(gauss2,x,y2,p0=p0)
                
        
        FWHM_CIII=round(2*np.sqrt(2*np.log(2))*popt1[4])
        FWHM_MgII=round(2*np.sqrt(2*np.log(2))*popt1[7])


        CIII_low=(popt1[0]+popt1[1]*popt1[3])
        CIII_high=(gauss2(popt1[3],*popt1))    
        CIII_dif=CIII_high-CIII_low

        MgII_low=(popt1[0]+popt1[1]*popt1[6])
        MgII_high=(gauss2(popt1[6],*popt1))    
        MgII_dif=MgII_high-MgII_low
    
        
    
        print('name[i],LyA_dif,CIV_dif,CIII_dif,MGII_dif,FWHM_LyA,FWHM_CIV,FWHM_CIII,FWHM_MGII')
        value1=[name[i],0,0,CIII_dif,MgII_dif,0,0,FWHM_CIII,FWHM_MgII]
   # print(value1)
        ax[1].plot(x,y2,label='redshifted data')
        ax[1].plot(x, gauss2(x, *popt1),'r', label='guassianfit')
        ax[1].set_xlabel('Wavelength')
        ax[1].set_ylabel('Flux')       
        ax[1].legend()
            
        
        
        
    elif ((z>=1.42) and (z<=3.6)):
        p0=[med,0,10,1549,1,10,1908,1]
        popt1,pcov1= curve_fit(gauss2,x,y2,p0=p0)
        
        FWHM_CIV=round(2*np.sqrt(2*np.log(2))*popt1[4])
        FWHM_CIII=round(2*np.sqrt(2*np.log(2))*popt1[7])


        CIV_low=(popt1[0]+popt1[1]*popt1[3])
        CIV_high=(gauss2(popt1[3],*popt1))    
        CIV_dif=CIV_high-CIV_low

        CIII_low=(popt1[0]+popt1[1]*popt1[6])
        CIII_high=(gauss2(popt1[6],*popt1))    
        CIII_dif=CIII_high-CIII_low
    
        
    
        print('name[i],LyA_dif,CIV_dif,CIII_dif,MGII_dif,FWHM_LyA,FWHM_CIV,FWHM_CIII,FWHM_MGII')
        value1=[name[i],0,CIV_dif,CIII_dif,0,0,FWHM_CIV,FWHM_CIII,0]
        #print(value1)
        ax[1].plot(x,y2,label='redshifted data')
        ax[1].plot(x, gauss2(x, *popt1),'r', label='guassianfit')
        ax[1].set_xlabel('Wavelength')
        ax[1].set_ylabel('Flux')       
        ax[1].legend()
    
    elif ((z>=0.34) and (z<=2.178)):
        p0=[med,0,10,2800,1]
        popt1,pcov1= curve_fit(gauss1,x,y2,p0=p0)
            
          
        FWHM_MgII=round(2*np.sqrt(2*np.log(2))*popt1[4])


        MgII_low=(popt1[0]+popt1[1]*popt1[3])
        MgII_high=(gauss1(popt1[3],*popt1))    
        MgII_dif=MgII_high-MgII_low
    
    
        print('name[i],LyA_dif,CIV_dif,CIII_dif,MGII_dif,FWHM_LyA,FWHM_CIV,FWHM_CIII,FWHM_MGII')
        value1=[name[i],0,0,0,MgII_dif,0,0,0,FWHM_MgII]
    #print(value1)
        ax[1].plot(x,y2,label='redshifted data')
        ax[1].plot(x, gauss1(x, *popt1),'r', label='guassianfit')
        ax[1].set_xlabel('Wavelength')
        ax[1].set_ylabel('Flux')       
    ax[1].legend()
    print(value)
    print(value1)
    with open('CSVFILE_bright.csv', 'a', newline='',) as f_object:
        writer_object = csv.writer(f_object)
        writer_object.writerow(value)  
        f_object.close()
    with open('CSVFILE_faint.csv', 'a', newline='',) as f_object:
        writer_object = csv.writer(f_object)
        writer_object.writerow(value1)  
        f_object.close()
    

    plt.savefig(path+name[i]+'.png')
    #plt.show()
    plt.clf()
