# -*- coding: utf-8 -*-
"""
Created on Fri Aug 28 12:27:56 2020

@author: akeating
"""

# multilayer thermal model
#!/usr/bin/python3
import numpy as np
import cmath
import yaml
import sys,os
import csv
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle  
import ctypes  # An included library with Python install.   
import time
from scipy.integrate import quad
from scipy import  LowLevelCallable
import scipy.integrate as si_integrate
from scipy import interpolate
import quadpy # install using conda install conda-forge::quadpy

# import argparse

# parser=argparse.ArgumentParser(
#     description='''A thermal model for plotting measured and simulated data as well as extracting thermal parameters from measured data obtained through a wide-band 3-omega measurement.  The model implements published models for bulk materials as well as thin films ''',
#     epilog="""Provide feedback to adrian.keating@uwa.edu.au.""")
# parser.add_argument('filename.yaml', type=int, default=42, help='Requires a yaml configuration file to setup the model')
# args=parser.parse_args()
__author__='Adrian Keating (UWA)'
__version__='0.3'

LAYERS='layers'
TEMPERATURE='measured'
AMPLIFIER='amplifier'
HEATER='heater'
MATERIALS='materials'

# next  TRY to reduce uncertainry by
# EXTEND heater capacitance optimization to borca and Kim models
#1) include effect of blackboddy radiatiojn
#2) add effect of top layer AIR

def is_number_tryexcept(s):
    """ Returns True is string is a number. """
    try:
        float(s)
        return True
    except (ValueError,TypeError):
        return False

    
class ThermalModels:
    def __init__(self,filedata):
        self.fileio = filedata 
        self.dB=filedata.dB
        self.glb_index=0  # keeps track of all the optimization parameeters
        self.LAYERS='layers'
        self.TEMPERATURE='measured'
        self.AMPLIFIER='amplifier'
        self.HEATER='heater'
        self.MATERIALS='materials'
        self.keywords=[self.MATERIALS,self.HEATER,self.AMPLIFIER,self.TEMPERATURE,self.LAYERS] #,'model','measured'] 
        self.optimseflag=False
        self.Prms=1
        self.l=1
        self.ThermalConductivity=1
        self.solveLoop=0
        self.b=1
        self.Layers=[]
        self.ThemralDiffusivity=1
        self.HeaterFlag=False
        self.opt={}
        self.k=[]
        self.D=[]
        self.cross_k=[]
        self.thick=[]
        self.yref=1
        self.font=15
        self.interation=0
        self.heaterlayer=None
        #self.optimize=False
        self.boundarycondition=[ 'isothermal','semi-infinite', 'adiabatic']
        self.boundarystate=0
        self.modelindex=0  # can have any number of models used to fit the data starting from 0 (1 layer model)
        self.modelnames=['Borca','Kim', 'Cahill'] # 'Bypass-Cahill']
        
    def plotlayers(self,ax,message,id_str):
        xw=.3
        yh=1
        x0=-0
        n=len(self.fileio.dB[self.LAYERS])
        scale=1/(n)
        space=scale*.2
        imgsize=scale-space
        y0=self.yref#-self.modelindex*space
        thick=[]
        k=[]
        D=[]
        righthandfont=self.font*.7
        Layers=self.Layers[-1]
        for ind,eachlayer in enumerate(Layers):
            if ('Prms' in eachlayer or 'prms' in eachlayer):
                # this layer contains the heater
                thick.append(float(eachlayer['thickness'])) 
            else:
                thick.append(float(eachlayer['thickness']))
                k.append(float(eachlayer['thermal_conductivity']))
                if ('thermal_diffusivity' in eachlayer):
                    D.append(float(eachlayer['thermal_diffusivity']))
                else:
                    D.append(diffusivity(float(eachlayer['thermal_conductivity']),float(eachlayer['cp'])*float(eachlayer['density'])))
        Dmax=max(D)
        kmax=max(k)
        ratio=.8  # 1-ratio is the amount  of the substrate shown
        maxt=sum(thick[:-1])
        deltay=0
        width=2*xw*.8
        y0=y0-space # leave a spac efor the Data and Layer stack text
        y00=0
        ax.text(x0,y0+space/4, id_str,size=righthandfont)
        ax.text(x0+width+space/3,y0-(righthandfont/(5*72))*message.count('\n') ,message[:-1],size=righthandfont) 
        for ind,eacht in enumerate(thick):
            height=imgsize*(eacht/maxt)*ratio           
            if(ind==len(thick)-1):
                pass
                height=(1-ratio)*imgsize
            if(ind==0):
                b=(width/5)
                deltax=width/2-b/2
                color='orange'
            else:
                k0=min(k[ind-1]/kmax,1)
                D0=min(abs(np.log(D[ind-1]/Dmax)),1)
                deltax=0
                b=width
                color=(0,k0,D0)
            y00+=height
            y0=y0-height 
            rect =Rectangle((x0+deltax, y0),  b, height, fc =color,  ec ='b', lw = 1)
            ax.add_patch( rect ) 
        self.yref=y0
        return 
        
    def buildopt_dict(self,eachlayer):
        if('optimize' in  eachlayer):
            self.opt=eachlayer['optimize'].split(',')
        else:
            self.opt={}
        return

    def resetkeyparameters(self):    
        self.thick=[]
        self.k=[]
        self.D=[]
        self.cross_k=[]
        self.HeaterFlag=False
        self.glb_index=0
 

    def checkoptimise(self,keyword,ind,eachlayer,default):
        if(keyword in self.opt):
            returnval=((float(self.args_opt[self.glb_index+self.opt.index(keyword)])))
        else:
            if (keyword in eachlayer):
                returnval=float(eachlayer[keyword])  # if not optimization is specified use the value given in the file
            else:
                returnval=default      # if no value given in the file, use an internall defined fault value
        return returnval

    def checklayeroptimization(self,Layers,*args_opt):
        self.args_opt=args_opt
        for ind,eachlayer in enumerate(Layers):
            self.buildopt_dict(eachlayer) # build the optimisation check dictionary
            if ('prms' in eachlayer or 'p' in eachlayer):
                # this layer contains the heater
                keyword='thermal_boundary_resistance'
                self.heaterlayer=ind
                self.dB[self.HEATER][keyword]=self.checkoptimise(keyword,ind,eachlayer,0)
                if (keyword in self.opt):
                    self.HeaterFlag=True 
                    self.glb_index+=1
            else:
                # TO DO Thsi section need to add the k and D send from the routine as these are used in te iteration FITTING
                self.thick.append(float(eachlayer['thickness']))
                layerD=self.checkoptimise('thermal_diffusivity',ind,eachlayer,1e-8)
                self.D.append(layerD)
                layerk=self.checkoptimise('thermal_conductivity',ind,eachlayer,1)
                self.k.append(layerk)
                layercross=self.checkoptimise('cross_thermal_conductivity',ind,eachlayer,1)
                self.cross_k.append(layercross)
                self.glb_index=self.glb_index+len(self.opt) # keeps track of optimization parameters for sdifferent structures
        return  


                
    def Bulk_fit2(self,Layers,freq, *args_opt,verbose=True):
        borca={}
        self.resetkeyparameters()        
        self.args_opt=args_opt
        self.checklayeroptimization(Layers, *args_opt)            
        print('Cahill Model forces Boundry COndition Set to ISOTHERMAL on last layer')  
        out=[]
        heater= self.dB[self.HEATER]
        P=float(heater['prms'])
        l=float(heater['l'])
        b=float(heater['b'])
        lowerlimit= self.dB[self.TEMPERATURE]['integration']['min']/b
        upperlimit= self.dB[self.TEMPERATURE]['integration']['max']/b
        ktoplayer=k[0]
        Dtoplayer=D[0]
        arglist=[P,l,b,0,ktoplayer,Dtoplayer] # the arg list must be passed as a tuple
        arglist=np.array([0,P,l,b,1/7,ktoplayer,Dtoplayer])
        print('Solving Loop',self.solveLoop)
        for ff in freq:
            omega = 2*np.pi*ff
            arglist[4]=omega #np.array([0,P,l,b,omega,ktoplayer,Dtoplayer]) #[3]=omega
            out1 = self.complex_integration(self.CahillModel , lowerlimit,upperlimit,args=arglist)
            out.append(out1)
        self.solveLoop+=1
        outlist=np.array(out)
        if(self.HeaterFlag): outlist=self.heater_capacitance_effect(fileio.freq,outlist)
        real_imag_data=np.append(np.real(outlist),np.imag(outlist))
        return real_imag_data

# see here on how to speed up with JIT and arbitary args
# https://stackoverflow.com/questions/60600672/how-to-reduce-integration-time-for-integration-over-2d-connected-domains
    def complex_integration(self,integrand, a, b,  args):
        # varg1,arg2,arg3,arg4,arg5,arg6
        def real_func(x,args):
            return np.real(integrand(x,args))
        def imag_func(x,args):
            return np.imag(integrand(x,args))
        real_integral =  si_integrate.quad(real_func, a, b, args=List(args))
        imag_integral = si_integrate.quad(imag_func, a, b, args=List(args))   
        return real_integral[0] + 1j*imag_integral[0]
    
    # see here on how to call jit from a class
    # https://stackoverflow.com/questions/41769100/how-do-i-use-numba-on-a-member-function-of-a-class
    # allowing arrays to be send and compelx returned     
    @staticmethod    
    #@jit(nopython=True, cache=True)
    def CahillModel(s,args):
        #s=args[0]
        Prms=args[1]
        l=args[2]
        b=args[3]
        omega=args[4]
        ThermalConductivity=args[5]
        ThemralDiffusivity=args[6]
        f = (Prms/l)
        f/=(np.pi * ThermalConductivity)
        f*= ((np.sin(s*b))**2)/(((s*b)**2) * np.sqrt(s**2 + (2j*omega/ThemralDiffusivity)))
        return f
        
# high speed integration
#  https://scicomp.stackexchange.com/questions/15907/why-does-matlabs-integral-outperform-integrate-quad-in-scipy
    def speedtest(self):
        t0 = time.time()
        omega = 10
        #vintegral = np.vectorize(integral)
        for i in range(300): 
            result = self.integral(self.f_integrand, 0, 30, omega)
        print ('Time Difference 300 runs :',time.time()-t0)
        print ('COMPLX RESULT',result)
    

    
    def integral(self,integrand, a, b,  args):
        def real_func(arg):
            return np.real(integrand(args))
        def imag_func(arg):
            return np.imag(integrand(args))
        real_integral =  si.quad(integrand, a, b, args=args)
        imag_integral = si.quad(integrand, a, b, args=args)   
        return real_integral[0] + 1j*imag_integral[0]
    
    # see here need for static for jit dfor use in a class
    # https://stackoverflow.com/questions/41769100/how-do-i-use-numba-on-a-member-function-of-a-class
    @staticmethod
    #@jit(complex128(float64, float64), nopython=True, cache=True)
    def f_integrand(s, omega):
        sigma = np.pi/(np.pi+2)
        xs = np.exp(-np.pi*s/(2*sigma))
        x1 = -2*sigma/np.pi*(np.log(xs/(1+np.sqrt(1-xs**2)))+np.sqrt(1-xs**2))
        x2 = 1-2*sigma/np.pi*(1-xs)
        zeta = x2+x1*1j
        Vc = 1/(2*sigma)
        theta =  -1*np.arcsin(np.exp(-np.pi/(2.0*sigma)*s))
        t1 = 1/np.sqrt(1+np.tan(theta)**2)
        t2 = -1/np.sqrt(1+1/np.tan(theta)**2)
        return np.real((t1-1j*t2)/np.sqrt(zeta**2-1))*np.exp(1j*omega*s/Vc);
    
    def heater_capacitance_effect(self, Frequency, Temperature):
        heater=self.dB[self.HEATER]
        omega = np.array(2*np.pi*Frequency);
        Thermal_resistance=heater['prms']/(2*heater['b']*heater['l'])  #W/ohm
        Rth=float(heater['thermal_boundary_resistance']) # K*m^2/W
        HeatCapacty= heater['density']*heater['specific_heat']  # (J/kg.K)*(kg/m^3)=(J/K*m^3)
        numerator= (Temperature + Rth*Thermal_resistance)
        denominator=(1 +  HeatCapacty*heater['thickness']*complex(0,2)*omega*(Rth + Temperature/Thermal_resistance))
        y=numerator/denominator
        return y

    def borca_thin_film_fit3_1D(self,Layers,freq,*args_opt):
        print('Solve Loop:'+f'{self.interation:3d} optimizing:'+' and '.join([f"{x:.4g}" for x in args_opt]), end="\r")
        out=self.borca_thin_film_fit3(Layers,freq, *args_opt,verbose=False) 
        outlist=np.append(np.real(out),np.imag(out))
        self.interation+=1
        return outlist 

    def checkboundaryconditions(self,Layers,defaultstate=0):
        #Substrate boundary condition,        
                      #  0 for semi-infinite, 
                      # +1 for an adiabatic condition 
                      # -1 for an isothermal boundary when the substrate is finite %  
        if ('boundary' in Layers[-1]):
            if (Layers[-1]['boundary'] in self.boundarycondition):
                state=int(self.boundarycondition.index(Layers[-1]['boundary']))-1
                bcmsg='Boundry condition found and set to:'+Layers[-1]['boundary']+' (State='+str(state)+')'
            else:
                raise KeyError('Within the database file, a boundary condition was defined in the last layer ['+self.dB[LAYERS][self.modelindex][-1]+'] but the keyword used was not one of isothermal, semi-infinite or adiabatic.  Please edit the file and correct')
                state =  defaultstate
                bcmsg='No boundry condition Found. Set to '+self.boundarycondition[defaultstate+1].upper()+' last layer'
        else:
            # default
            state =defaultstate
            bcmsg='Noboundry condition found. Set to '+self.boundarycondition[defaultstate+1].upper()+' last layer'
        if self.interation==0:
            print(bcmsg)
        return state


                            
    def borca_thin_film_fit3(self,Layers,freq, *args_opt,verbose=False):
        borca={}
        self.resetkeyparameters()
        self.args_opt=args_opt
        self.checklayeroptimization(Layers, *args_opt)      
        self.boundarystate=self.checkboundaryconditions(Layers,defaultstate=0)

        borca['n'] = len(Layers)-1  
        borca['d']= np.array(self.thick) #Thicknesses
        borca['Kx']       = np.array(self.k)*self.cross_k   
        borca['Ky']        = np.array(self.k)
        borca['Kxy']       = np.array(self.cross_k) #= borca['Kx'] /borca['Ky'] 
        borca['Dy']          = np.array(self.D)
        out=[]
        heater= self.dB[HEATER]
        lowerlimit= self.dB[TEMPERATURE]['integration']['min']/heater['b']
        upperlimit= self.dB[TEMPERATURE]['integration']['max']/heater['b']
        TempR=0
        Pdensity=heater['prms']/(heater['l']*heater['width'])
        if( 'thermal_boundary_resistance' in heater):
            #'K*m^2/W'
            TempR=heater['thermal_boundary_resistance']*heater['prms']/(heater['l']*heater['width'])
        for ff in freq:
            borca['omega'] = 2*np.pi*ff
            out1, err = quadpy.quad(lambda x: (self.borca_equation(x,borca,heater)), lowerlimit,upperlimit,epsabs=1e-10, epsrel=1e-10, limit=100000)
            out.append(out1)
        outlist=np.array(out) 
        if(heater['capacitance_flag'] ==True or self.HeaterFlag):   
            if(self.interation==0):               
                if(not(self.HeaterFlag)):
                    print('Including the effect of heater thermal heat capacity BUT NOT thermal boundary resistance')
                else:
                    print('Including Tthe effect of heater thermal heat capacity through the thermal_boundary_resistance parameter')
            outlist=self.HeaterEffect(self.dB, freq, outlist)
        return outlist

    def complex_borca(self,integrand, a, b,  args):
        def real_func(x,args):
            return np.real(integrand(x,args))
        def imag_func(x,args):
            return np.imag(integrand(x,args))
        real_integral =  si_integrate.quad(lambda x:(real_func(x,args)), a, b)
        imag_integral = si_integrate.quad(imag_func, a, b, args=(args))   
        return real_integral[0] + 1j*imag_integral[0]
    


    def borca_equation(self,ta,borca,heater):
        length=heater['l']
        T0=(-heater['prms']/(length*np.pi*borca['Ky'][0])) #(-P/(l*np.pi * Ky1)) 
        Aj0=self.Aj(0, ta, borca)
        invAj0=1/Aj0
        Bj0=self.Bj(0, ta, borca)
        invBj0=1/Bj0
        fun= T0*invAj0*invBj0*((np.sin(ta*heater['b']))**2)/((ta*heater['b'])**2)    
        return fun
    
    def Aj(self,j, ta, fun_structure):    
        n= fun_structure['n']-1
        if (j == n): 
            phase=(self.Bj(n, ta, fun_structure)*fun_structure['d'][n])
            y=np.array([])
            for phi in phase:
                y=np.append(y,-1*(cmath.tanh(phi))**self.boundarystate)
        else:
            # adding a np.finfo(float).eps to all numbers prevents an underfloat calculation result
            Ajplus1=self.Aj(j+1, ta, fun_structure)+np.finfo(float).eps*(1+1j)
            Bjplus1=self.Bj(j+1, ta, fun_structure)+np.finfo(float).eps*(1+1j)
            Bj0=self.Bj(j, ta, fun_structure)+np.finfo(float).eps*(1+1j)
            phase=fun_structure['d'][j]*Bj0
            gamma=fun_structure['Ky'][j]*Bj0
            beta=Ajplus1*fun_structure['Ky'][j+1]*Bjplus1  # generates FloatingPointError: underflow encountered in multiply from many zero values in Ajplau1 
            num = (beta/(gamma) - np.tanh(phase))
            den=1-(beta*np.tanh(phase)/(gamma))
            y=num/den
        return np.array(y)

    def Bj(self,j, ta, fun_structure):
        y= np.sqrt(fun_structure['Kxy'][j] *(ta**2) + 2j*fun_structure['omega']/fun_structure['Dy'][j])
        return y

    def HeaterEffect(self,dB,freq,out):   
        out=self.heater_capacitance_effect( freq, out)
        return out
            
    def kim_thin_film_fit_1D(self, Layers,freq,*args_opt):
        print('Solve Loop:'+f'{self.interation:3d} optimizing:'+' and '.join([f"{x:.4g}" for x in args_opt]), end="\r")
        out=self.kim_thin_film_fit(Layers,freq, *args_opt,verbose=False)
        outlist=np.append(np.real(out),np.imag(out))
        self.interation+=1
        return outlist 
       
    def kim_equation(self,tta,kim,heater):
        b=heater['b']
        T0=(heater['prms']/(2*heater['l']*np.pi*kim['Ky'][0])) #(-P/(l*np.pi * Ky1)) 
        places=10
        A=np.array((1,1))/2
        b=heater['b']
        omega=kim['omega']
        ta=tta
        fun=[]
        if(1):
            
            for mm in ta:
                m=mm
                if(self.heaterlayer>0):
                # CALC A
                    strM=''
                    uj=np.sqrt(m**2 -complex(0,2)*kim['omega']/kim['Dy'][0])
                    realphase=np.real(uj*kim['d'][0])
                    imagphase=np.imag(uj*kim['d'][0])
                    phase_j=complex(np.min([abs(realphase),100])*np.sign(realphase),imagphase)
                    exp_n_l=cmath.exp(phase_j)
                    gamma_j=uj*kim['Kx'][0]
                    M=np.array( ((1,0), (0, 1)) )
                    prop=np.array( ((1/exp_n_l,0), (0, exp_n_l)) )
                    for j,eachD in enumerate(kim['d'][0:self.heaterlayer]):  # go to n-1
#
                        ujplus=np.sqrt(m**2 -complex(0,2)*kim['omega']/kim['Dy'][j+1])
                        gamma_jplus=ujplus*kim['Kx'][j+1]
        
                        aa=np.round((gamma_jplus+gamma_j)/(2*gamma_jplus),places)
                        bb=np.round((gamma_jplus-gamma_j)/(2*gamma_jplus),places)
                        AB= np.array( ((aa,bb), (bb,aa)) )
                        Mj=np.dot(prop,AB)
                        #Mj=prop
                        M=np.dot(Mj,M)
                        # copy over old values
                        if(j<len(kim['d'][0:self.heaterlayer])):
                            realphase=np.real(ujplus*kim['d'][j+1])
                            imagphase=np.imag(ujplus*kim['d'][j+1])
                            phase_jplus=complex(np.min([abs(realphase),100])*np.sign(realphase),imagphase)
                            exp_nplus_l=cmath.exp(phase_jplus) 
                            gamma_j=gamma_jplus
                            prop=np.array( ((1/exp_nplus_l,0), (0, exp_nplus_l)) )
                    A= np.dot(M,np.array((1,0)))

                strM=''
                uj=np.sqrt(m**2 -complex(0,2)*kim['omega']/kim['Dy'][self.heaterlayer])
                realphase=np.real(uj*kim['d'][self.heaterlayer])
                imagphase=np.imag(uj*kim['d'][self.heaterlayer])
                phase_j=complex(np.min([abs(realphase),100])*np.sign(realphase),imagphase)
                exp_n_l=cmath.exp(phase_j)
                gamma_j=uj*kim['Kx'][self.heaterlayer]
                M=np.array( ((1,0), (0, 1)) )
                prop=np.array( ((1/exp_n_l,0), (0, exp_n_l)) )
                for j,eachD in enumerate(kim['d'][self.heaterlayer:-1]):  # go to n-1
        
                    #uj=np.sqrt(m**2 -complex(0,2)*kim['omega']/kim['Dy'][j])#
                    ujplus=np.sqrt(m**2 -complex(0,2)*kim['omega']/kim['Dy'][j+1])
                    gamma_jplus=ujplus*kim['Kx'][j+1]
    
                    aa=np.round((gamma_jplus+gamma_j)/(2*gamma_j),places)
                    bb=np.round((gamma_j-gamma_jplus)/(2*gamma_j),places)
                    AB= np.array( ((aa,bb), (bb,aa)) )
                    Mj=np.dot(prop,AB)
                    #Mj=prop
                    M=np.dot(M,Mj)
                    # copy over old values
                    if(j<len(kim['d'][self.heaterlayer:-1])):
                        realphase=np.real(ujplus*kim['d'][j+1])
                        imagphase=np.imag(ujplus*kim['d'][j+1])
                        phase_jplus=complex(np.min([abs(realphase),100])*np.sign(realphase),imagphase)
                        exp_nplus_l=cmath.exp(phase_jplus) 
                        gamma_j=gamma_jplus
                        prop=np.array( ((1/exp_nplus_l,0), (0, exp_nplus_l)) )
                B= np.dot(M,np.array((0,1)))
                scale= (B[0]+B[1])/(A[0]*B[1]-A[1]*B[0])
                if(m!=0):
                    fun.append(T0 * scale*((np.sin(m*b))**2)/(((m*b)**2) * np.sqrt(m**2 + (complex(0,2)*omega/kim['Dy'][0]))))
                else:
                    fun.append(T0 * scale*1/np.sqrt(m**2 +(complex(0,2)*omega/kim['Dy'][0])))
        return fun
    
    def kim_thin_film_fit(self, Layers,freq, *args_opt,verbose=True):
        kim={}
        thick=[]
        k=[]
        D=[]

        self.resetkeyparameters()
        self.args_opt=args_opt
        self.checklayeroptimization(Layers, *args_opt)  
        self.boundarystate=self.checkboundaryconditions(Layers,defaultstate=-1)
        #index=0

        kim['n'] = len(Layers)-1  
    
        reverse=False
        if(reverse==True):
            kim['Kx']       = np.array(self.k[::-1])  #reverse the list   
            kim['Ky']        = np.array(self.k[::-1])  #reverse the list  
            kim['Dy']          = np.array(self.D[::-1]) #reverse the list  
            kim['d']= np.array(self.thick[::-1]) #reverse the list  
        else:
            kim['Kx']       = np.array(self.k)  #DON'T reverse the list   
            kim['Ky']        = np.array(self.k)  #DON'T reverse the list  
            kim['Dy']          = np.array(self.D) #DON'T reverse the list  
            kim['d']= np.array(self.thick) #DON'T reverse the list  
        kim['Kxy']         = kim['Kx'] /kim['Ky'] 
        out=[]

        heater= self.dB[HEATER]
        lowerlimit= self.dB[TEMPERATURE]['integration']['min']/heater['b']
        upperlimit= self.dB[TEMPERATURE]['integration']['max']/heater['b']
        for ff in freq:
            kim['omega'] = 2*np.pi*ff
            tol=1/float(self.dB[TEMPERATURE]['integration']['max'])**2
            limit0=int(0.5/np.sqrt(tol))
            out1, err = quadpy.quad(lambda x: self.kim_equation(x,kim,heater), lowerlimit,upperlimit)
            out.append(out1)
        outlist=np.array(out)
        if(heater['capacitance_flag'] ==True or self.HeaterFlag): 
            if(self.interation==0):               
                if(not(self.HeaterFlag)):
                    print('Including the effect of heater thermal heat capacity BUT NOT thermal boundary resistance')
                else:
                    print('Including Tthe effect of heater thermal heat capacity through the thermal_boundary_resistance parameter')
            outlist=self.HeaterEffect(self.dB, freq, outlist)
        return outlist
    
    
    def Bulk_function(self,ta,Prms, l, b, omega, ThermalConductivity, ThemralDiffusivity):
        f = (Prms/(l * np.pi * ThermalConductivity)) * ((np.sin(ta*b))**2)/(((ta*b)**2) * np.sqrt(ta**2 + (2*complex(0,1)*omega/ThemralDiffusivity)))
        return f
    
    def Bulk_fit2_1D(self,Layers,freq,*args_opt):
        params=' and '.join([f"{x:1.4g}" for x in args_opt])+'      '
        print('Solve Loop:'+f'{self.interation:3d} optimizing:'+params, end="\r")
        out=self.Bulk_fit2(Layers,freq, *args_opt,verbose=False)
        outlist=np.append(np.real(out),np.imag(out))
        self.interation+=1
        return outlist 
    
    def Bulk_fit2(self,Layers,freq, *args_opt,verbose=True):
        # borca={}
        # thick=[]
        # k=[]
        # D=[]

        # glb_index=0 # a counter that cyces through the opt list
        # HeaterFlag=False
        
        self.resetkeyparameters()
        self.args_opt=args_opt
        self.checklayeroptimization(Layers, *args_opt) 
        
        # for ind,eachlayer in enumerate(Layers):
        #     if ('prms' in eachlayer or 'p' in eachlayer):
        #         # this layer contains the heater
        #         heatind=ind
        #         if('optimize' in  eachlayer):
        #             opt=eachlayer['optimize'].split(',')
        #         else:
        #             opt={}
        #         checkparam='thermal_boundary_resistance'
        #         if(checkparam in opt):
        #             self.dB[HEATER]['thermal_boundary_resistance']=((float(args_opt[glb_index+opt.index(checkparam)])))
        #             HeaterFlag=True
        #             glb_index=glb_index+1
        #     else:
        #         # TO DO Thsi section need to add the k and D send from the routine as these are used in te iteration FITTING
        #         thick.append(float(eachlayer['thickness']))
        #         if('optimize' in  eachlayer):
        #             opt=eachlayer['optimize'].split(',')
        #         else:
        #             opt={}
        #         checkparam='thermal_diffusivity'
        #         if(checkparam in opt):
        #             D.append((float(args_opt[glb_index+opt.index(checkparam)])))
        #         else:
        #             # if not on the optimization list just copy the existign value
        #             if ('thermal_diffusivity' in eachlayer):
        #                 D.append(float(eachlayer['thermal_diffusivity']))
        #             else:
        #                 D.append(diffusivity(float(eachlayer['thermal_conductivity']),float(eachlayer['cp'])*float(eachlayer['density'])))
        #         checkparam='thermal_conductivity'
        #         if(checkparam in opt):
        #              # check the order
        #             k.append((float(args_opt[glb_index+opt.index(checkparam)])))
        #         else:
        #             # if not on the optimization list just copy the existign value
        #             k.append(float(eachlayer['thermal_conductivity']))
        #         glb_index=glb_index+len(opt)
                
            
        if self.interation==0:
            print('WARNING: Cahill model forces the boundry condition to be set to ISOTHERMAL on last layer, independent of the adiabatic setting in the file')  
            # if(not(self.HeaterFlag)):
            #     print('Including the effect of heater thermal heat capacity BUT NOT thermal boundary resistance')
            # else:
            #     print('Including Tthe effect of heater thermal heat capacity through the thermal_boundary_resistance parameter')
        out=[]
        heater= self.dB[HEATER]
        P=float(heater['prms'])
        l=float(heater['l'])
        b=float(heater['b'])
        lowerlimit= self.dB[TEMPERATURE]['integration']['min']/b
        upperlimit= self.dB[TEMPERATURE]['integration']['max']/b
        ktoplayer=self.k[0]
        Dtoplayer=self.D[0]
        for ff in freq:
            omega = 2*np.pi*ff
            out1, err = quadpy.quad(lambda x: self.Bulk_function(x,P,l,b, omega, ktoplayer,Dtoplayer), lowerlimit,upperlimit)
            out.append(out1)
        outlist=np.array(out)
        # if(self.HeaterFlag): outlist=self.HeaterEffect(self.dB,freq,outlist)
        if(heater['capacitance_flag'] ==True or self.HeaterFlag): 
            if(self.interation==0):               
                if(not(self.HeaterFlag)):
                    print('Including the effect of heater thermal heat capacity BUT NOT thermal boundary resistance')
                else:
                    print('Including Tthe effect of heater thermal heat capacity through the thermal_boundary_resistance parameter')
            outlist=self.HeaterEffect(self.dB, freq, outlist)
    
        return outlist
    
     

        
class ThermalFiles:
    def __init__(self,yamlfile):
        self.yamlfile = yamlfile 
        self.dB={}
        self.LAYERS='layers'
        self.TEMPERATURE='measured'
        self.AMPLIFIER='amplifier'
        self.HEATER='heater'
        self.MATERIALS='materials'
        self.keywords=[self.MATERIALS,self.HEATER,self.AMPLIFIER,self.TEMPERATURE,self.LAYERS] #,'model','measured'] 
        self.getdata()
        self.allparameters={'Materials':[],'measured':[],'model':[],'Heater':['materixal','R','l','width','thickness','Prms','capacitance_flag'],'Amplifier':['gain_nominal',{'use_calibration':['Calibration_in','Calibration_out']}],'Measured':['downsample','output','range','ambientTinC','model','plot'],'Layers':[]}
        #self.checkparametersinfile(self.keywords,self.dB)
        self.downsample=int(self.temperature('downsample'))
        #print('self.downsample',self.downsample,self.dB[TEMPERATURE].keys())
        self.fmin_sim=self.temperature('fmin',1)
        self.fmax_sim=self.temperature('fmax',1e6)
        self.Nsim=int(self.temperature('samples',50))  
        self.Temperature=self.temperature('ambientTinC',25)
        self.internalf=np.logspace(np.log10(self.fmin_sim),np.log10(self.fmax_sim),self.Nsim)
        self.Nmax=int(self.temperature('range',1)  *len(self.internalf[0::self.downsample]))
        self.freq=np.array(self.internalf[0::self.downsample])[0:self.Nmax] # select vevery nth element
        self.scaleT=1 
        self.flist=np.append(self.freq,self.freq)
        self.setGain()


    def setGain(self):
        if('gain_nominal' in self.dB[AMPLIFIER]):
            # overwriteGainUSed laer if 'measured' is defined
            self.GainUsed=self.dB[AMPLIFIER]['gain_nominal']*np.ones(self.Nmax)
        else:
            self.GainUsed=np.ones(self.Nmax)
        return
    
    def checkparametersinfile(self,defineddict,dicttocheck):
        mainkeys=list(dicttocheck.keys())
        mainkeyslower=[key.lower() for key in mainkeys]
        for each in defineddict:
            if each.lower() in mainkeyslower:
                ind=mainkeyslower.index(each.lower())
                subkeys=dicttocheck[mainkeys[ind]]
                if type(defineddict[each])==list:
                    newdict=dict.fromkeys(defineddict[each],"") #{item:'' for item in defineddict[each]} # create a blank dict
                    self.checkparametersinfile(newdict,subkeys)
                elif type(defineddict[each])==dict:
                    newdict=defineddict[each]
                    self.checkparametersinfile(newdict,subkeys)
                else:
                    pass # do nothing and return
            else:
                print('File: ['+self.yamlfile+'[ is missing a main keyword ['+each+']')
        return
                
    def checkparametersinfilex(self):
        mainkeys=list(self.dB.keys())
        mainkeyslower=[key.lower() for key in mainkeys]
        for each in self.allparameters:
            if each.lower() in mainkeyslower:
                ind=mainkeyslower.index(each.lower())
                for word in self.allparameters[each]:
                    subkeys=list(self.dB[mainkeys[ind]])
                    subkeyslower=[key.lower() for key in subkeys]
                    if type(word)==str :
                        if word.lower() in subkeyslower:
                            pass
                        else:
                            print('File: ['+self.yamlfile+'[ under MAIN KEY ['+mainkeys[ind]+'] is missing a main keyword ['+word+']')
                    else:
                        print('TYPE',type(word),word)
            else:
                print('File: ['+self.yamlfile+'[ is missing a main keyword ['+each+']')
    def findmeasuredfiles(self):
        if ('data' in self.dB[self.TEMPERATURE]):
            if (type(self.dB[self.TEMPERATURE]['data'])!=list):
                filelist=self.dB[self.TEMPERATURE]['data'].split(',')  # allows file format of 'Glass.csv,PS77.csv' or just 'Glass.csv'
            else:
                filelist=self.dB[self.TEMPERATURE]['data']      
        else:
            filelist=[None]
        self.filelist=filelist
            
    def calibratemeasured(self):
        if ('data' in self.dB[self.TEMPERATURE]):
            CalibrationIn=self.three_omega_file_reading(self.dB[self.AMPLIFIER]['calibration_in'])
            CalibrationOut=self.three_omega_file_reading(self.dB[self.AMPLIFIER]['calibration_out'])
            self.V3wLockin_rms_I= (np.array(CalibrationIn['x']) +complex(0,1)* np.array(CalibrationIn['y']))
            self.V3wLockin_rms_O=(np.array(CalibrationOut['x']) +complex(0,1)* np.array(CalibrationOut['y']))
            self.CalibrationCoef= self.V3wLockin_rms_O/self.V3wLockin_rms_I

    def temperature(self,keyparam,default=1):
        if (keyparam in self.dB[TEMPERATURE]):
            variable=int(self.dB[TEMPERATURE][keyparam])
        else:
            variable=default # use default otherwise
        return variable

    def interpolatecp(self,tempRange,cpRange):
        if type(cpRange)==list:
            tck = interpolate.splrep(tempRange,cpRange, s=0)
            cp = interpolate.splev(self.Temperature+278, tck, der=0)
        else:
            cp=cpRange
        return cp
    
    def diffusivity(self,k,cv):
        if(cv==0):
            D=np.inf
        else:
            D=k/cv
        return (D)
    
    def material(self, each,keyword):
        # this keyword is optional and should return a default value
        if keyword[-1]=='*':
            return float(self.default(keyword))
        else:
            if keyword in self.dB[MATERIALS][each.lower()]:
                return self.dB[MATERIALS][each.lower()][keyword]
                # if type(self.dB[MATERIALS][each.lower()][keyword])==list:
                #     return self.dB[MATERIALS][each.lower()][keyword]
                # else:
                #     return float(self.dB[MATERIALS][each.lower()][keyword])
            else:
                raise KeyError('Error KEYWORD MISSING')   
                return
    def getdata(self):
        print('Reading in database['+self.yamlfile+']')
        with open(self.yamlfile,'r', encoding='utf8') as f:
        # utf-8 required to load accented text such as André-Marie Ampère
        #try:
            dMaterial={}
            indexlayers=0
            GradingConfigData = yaml.load(f, Loader=yaml.FullLoader)
            for key in self.keywords:                
                if not(key in GradingConfigData.keys()):
                    self.dB[key.lower()]=[]
                    # set list to empty if parameter no given
            for each in  (GradingConfigData.keys()):
                if (each.lower() in self.keywords):
                    self.dB[each.lower()]={}
                    if not(each in GradingConfigData):
                        raise Exception('The required keyword: ['+each+ '] is missing from the datafile: '+self.yamlfile+'\nPlease edit and ensure this parameter is specified')
                    elif(GradingConfigData[each]==None):
                        raise Exception('The required keyword: ['+each+ '] is specified in the datafile: '+self.yamlfile+' but contains no parameters.\nPlease edit and ensure this parameter is specified')
                    for elem in GradingConfigData[each]:
                        if(type(elem)==list):
                            if(is_number_tryexcept(elem)):
                                self.dB[each.lower()][indexlayers]=float(elem)
                            else:
                                self.dB[each.lower()][indexlayers]=elem
                            indexlayers=indexlayers+1
                        else:
                            allkeys=[]
                            self.dB[each.lower()][elem.lower()]={}
                            if (type(GradingConfigData[each][elem])==dict):
                                for keys in GradingConfigData[each][elem]:
                                    allkeys.append([*keys][0])
                                    if(is_number_tryexcept(GradingConfigData[each][elem][keys])):
                                        self.dB[each.lower()][elem.lower()][keys.lower()]=float(GradingConfigData[each][elem][keys])
                                    else:
                                        self.dB[each.lower()][elem.lower()][keys.lower()]=GradingConfigData[each][elem][keys]
                            else:
                                # just copy elements
                                if(is_number_tryexcept(GradingConfigData[each][elem])):
                                    self.dB[each.lower()][elem.lower()]= float(GradingConfigData[each][elem])
                                else:
                                    self.dB[each.lower()][elem.lower()]= (GradingConfigData[each][elem])
                            dMaterial[each]={}
                else:
                    
                    mystr=''
                    for ind,thestr in enumerate(self.keywords):
                        if(ind==len(self.keywords)-1):
                            mystr=mystr+' and ' +str(thestr)
                        else:
                            mystr=mystr+str(thestr)+', '
                    print('Keyword:'+each+' is unknown in the file: '+self.yamlfile+'\n Expected values are: '+mystr)    
                    sys.exit(0)
        # update the heater with the selected material
        self.dB[HEATER].update(self.dB[MATERIALS][self.dB[HEATER]['material'].lower()])
        if ('width' in self.dB[HEATER] or 'w' in self.dB[HEATER] ):
            if (not('b' in self.dB[HEATER])):
                self.dB[HEATER]['b']=float(self.dB[HEATER]['width'])/2
        print('Checking database parameters in ['+self.yamlfile+']')        
        for each in self.dB[MATERIALS]:
            if  not('thermal_diffusivity' in self.dB[MATERIALS][each]):
                cp=None
                if 'cp' in self.dB[MATERIALS][each]:
                    cp=self.dB[MATERIALS][each]['cp']
                if 'specific_heat' in self.dB[MATERIALS][each]:
                    cp=self.dB[MATERIALS][each]['specific_heat']
                    self.dB[MATERIALS][each]['cp']=cp
                if (cp!=None):
                    self.Temperature=self.temperature('ambientTinC',25)
                    if (type(cp)==list):
                        if(type(self.dB[MATERIALS][each]['temp'])==list):
                            cp=self.interpolatecp(self.dB[MATERIALS][each]['temp'],self.dB[MATERIALS][each]['cp'])
                        else:
                            raise KeyError('Material ['+each+'] requires both a [cp] and [temp] list but both lists were not found in the file: '+self.yamlfile)
                    else:
                        pass
                    if ('density' in self.dB[MATERIALS][each]):
                        cv=cp*self.dB[MATERIALS][each]['density']
                        if  not('thermal_conductivity' in self.dB[MATERIALS][each]):
                            raise KeyError('Material ['+each+'] is missing a defintion of thermal_conductivity in the file: '+self.yamlfile)
                        if  not('thermal_conductivity' in self.dB[MATERIALS][each]):
                            raise KeyError('Material ['+each+'] is missing a defintion of thermal_conductivity in the file: '+self.yamlfile)
                        self.dB[MATERIALS][each]['thermal_diffusivity'] =self.diffusivity(self.dB[MATERIALS][each]['thermal_conductivity'], cv)
                    else:
                        raise KeyError('Material ['+each+'] is missing a defintion of density required to calculate thermal_diffusivity - provide density and specific heat (cp) or thermal_diffusivity in the file: '+self.yamlfile)
        return self.dB

    def three_omega_file_reading(self,filename,validate=False):
        heater=self.dB[HEATER]
        data={}
        configfile=filename
        try:
            with open(filename, mode='r') as csv_file:
                csv_reader = csv.DictReader(csv_file)
                line_count = 0
                f=[]
                x=[]
                y=[]
                for row in csv_reader:
                    if line_count == 0:
                        line_count += 1
                    x.append(check_and_get_keyword(row,['InPhase X (V)']))
                    y.append(check_and_get_keyword(row,['OutPhase Y (V)']))
                    f.append(check_and_get_keyword(row,['Frequency (Hz)']))
                    data['b']=check_and_get_keyword(row,['Z2 Heater Width','ZHeaterWidth'])/2         
                    data['l']=check_and_get_keyword(row,['Z1 Heater Length','ZHeaterLength'])
                    data['P']=check_and_get_keyword(row,['Power'])
                    data['t']=check_and_get_keyword(row,['Thin Film Thickness'])
                    data['R']=check_and_get_keyword(row,['RHeater'])
                    line_count += 1
                data['x']=x
                data['y']=y
                data['f']=f
        
                if(validate==True):
                    if not('prms' in heater):
                        heater['prms']=data['P']
                        print('Using power of '+str(heater['prms'])+' W, defined in the File:'+filename+' ' )
                    else:
                        print('Using power of '+str(heater['prms'])+' W, defined in the config file' )
                    
                    if(data['R']==0):
                        print('No Power was defined in the File:'+filename+'\n Usiung Default power of '+str(heater['Prms'])+' W')
                    else:
                        heater['r']=data['R']
                    if(data['t']==0):
                        print('No thickness was defined in the File:'+filename+'\n Usiung Default thickness of '+str(heater['thickness'])+' W')
                    else:
                        pass
        except FileNotFoundError:
            print('ERROR - File not found ---------\nThe entry in the file of measure:'+filename+' does not refer to a measured file in the current directory.  Press enter, recheck the spelling and re-run')
            input()
            sys.exit()
        return data





def check_and_get_keyword(row,keyw):
    out=None
    for eachkeyw in keyw:
        if (eachkeyw in row):
            out=float(row[eachkeyw])
    if out==None:
        raise Exception('Missing column: ['+eachkeyw+'] in the input file: '+configfile+'\n.....Exiting')
    return out




 


configfile=''    
def RC_Correction(f,R1,R2,C):
    omega=2*np.pi*np.array(f)
    correction=(float(R2)+float(R1))/(float(R2)+float(R1)*(1+complex(0,1)*omega*float(R2)*float(C)))
    return correction




class FormatThermalData:
    def __init__(self,Layers,popt=[],pcov=[],fit=False):
        self.opt = {}
        self.popt=popt
        self.pcov=pcov
        self.fit=fit
        self.Layers=Layers
        self.ThermalVar={'thermal_conductivity':'k','thermal_boundary_resistance':'Rth','thermal_diffusivity':'D'}
        #self.ThermalUnit={'thermal_conductivity':'W/mK','thermal_boundary_resistance':'K*m^2/W','thermal_diffusivity':'$m^2/s$'}
        
        self.ThermalUnit={'thermal_conductivity':['W/mK',1],'thermal_boundary_resistance':['K*mm^2/W',1e6],'thermal_diffusivity':['$mm^2/s$',1e6]}
        self.msg=''
        self.logdata=[]
        self.gbl_index=0
        self.index=0
        self.eachopt=None
        self.BuildData()


    def BuildData(self):

        eachlayer=self.Layers[0:][0]  
        self.BuildThermalData('thermal_boundary_resistance',eachlayer)
        self.gbl_index=self.gbl_index+len(self.opt)
        self.Layers[0]['thermal_boundary_resistance']=self.eachopt
        for j,eachlayer in enumerate(self.Layers[1:],start=1):
            self.logdata.append(['Layer', str(j),''])
            if (self.fit):
                self.msg==self.msg+'Data#1 fit using '
            self.msg=self.msg+'Layer'+str(j)+'\n'
               
            self.BuildThermalData('thermal_conductivity',eachlayer)
            self.Layers[j]['thermal_conductivity']=self.logdata[-1][1]
            self.BuildThermalData('thermal_diffusivity',eachlayer)
            self.Layers[j]['thermal_diffusivity']=self.logdata[-1][1]
            self.logdata.append(['thickness', str(float(eachlayer['thickness'])*1e6),'microns',str(0)])                   
            #self.BuildStructuralData('porosity',eachlayer)
            self.BuildStructuralData('material',eachlayer)
            self.BuildStructuralData('boundary',eachlayer)
            self.gbl_index=self.gbl_index+len(self.opt)
                    
    def BuildStructuralData(self,keyword,eachlayer): 
        if('optimize' in  eachlayer):
            self.opt=eachlayer['optimize'].split(',')
        else:
            self.opt={} 
        if(keyword in self.opt): 
            self.logdata.append([keyword, str(eachlayer[keyword]),'',str(0)])
                       
    def BuildThermalData(self,keyword,eachlayer):   
        if('optimize' in  eachlayer):
            self.opt=eachlayer['optimize'].split(',')
        else:
            self.opt={} 
   
        if(keyword in self.opt and self.fit==True):  
            popt_ind=self.gbl_index+self.opt.index(keyword)
            self.eachopt=self.popt[popt_ind]
            Gain=self.ThermalUnit[keyword][1]
            val=Gain*self.eachopt 
            if self.pcov.size>0:
                uncertainty=(self.pcov[popt_ind,popt_ind]**.5)*2  # 95% coverage
                uncertainty_pc=100*uncertainty/self.eachopt
                self.logdata.append([keyword+'(opt)', val ,self.ThermalUnit[keyword][0],str(uncertainty_pc)])
                self.msg=self.msg+self.ThermalVar[keyword]+'(opt)='+str(np.round(val,2))+' '+self.ThermalUnit[keyword][0]+'+/-'+str(np.round(uncertainty_pc,2))+'%\n'
            else:
                self.logdata.append([keyword+'(opt)', self.eachopt*Gain ,self.ThermalUnit[keyword][0],str(0)])
                self.msg=self.msg+self.ThermalVar[keyword]+'(opt)='+str(np.round(val,2))+' '+self.ThermalUnit[keyword][0]+'\n'
            self.index+=1
        else:
            if(keyword in eachlayer):
                self.eachopt=eachlayer[keyword]
                Gain=self.ThermalUnit[keyword][1]
                val=Gain*self.eachopt 
                self.logdata.append([keyword,val,self.ThermalUnit[keyword][0],str(0)])
                self.msg=self.msg+self.ThermalVar[keyword]+'='+str(np.round(val,2))+' '+self.ThermalUnit[keyword][0]+'\n'
                self.index+=1
            else:
                self.eachopt=0
        return 
                        
def main(filename): 
    fileio=ThermalFiles(filename)
    fileio.findmeasuredfiles()
    fileio.calibratemeasured()
    model=ThermalModels(fileio)  # pss the file data along to the models

    numLayers=len(fileio.dB[LAYERS])
    Basefontsize=14
    model.font=Basefontsize/np.max([1,np.sqrt(np.sqrt(numLayers))])
    identifyDatafiles=''
    for fileid,eachfile in enumerate(fileio.filelist):
        if( eachfile!=None):
            identifyDatafiles+='Data'+str(fileid+1)+':'+os.path.basename(eachfile)+'\n'
            
    if(fileio.dB[TEMPERATURE]['plot'] ==True):
        if(1):
            # do this only once for all the models
            w=12;h=12
            splashmsg="Thermal thin film model, Version:"+str(__version__)+" by "+str(__author__)
            print(splashmsg)
            fig = plt.figure(splashmsg,figsize=(w,h))
            ax = plt.subplot2grid((1,2),(0, 0))
            ax2 = plt.subplot2grid((1,2),(0, 1))
            ax2.axis('off')  
            ax2.relim()  # Recompute the data limits based on current artists
            ax2.autoscale_view()
            ax.set_ylabel('$\Delta $Temperature(C)',style='oblique', weight='bold', fontsize=Basefontsize)
            ax.set_xlabel('Frequency(Hz), $10^x$',style='oblique', weight='bold', fontsize=Basefontsize)
            ax.set_title('Model/Fit -\n'+identifyDatafiles,  loc='left', weight='bold', fontsize=Basefontsize*.6,transform=ax.transAxes)
            ax2.set_title('Layers',style='oblique', weight='bold', fontsize=Basefontsize)
            ax2.set_ylim((0,1))
            ax2.set_xlim((0,1))
        logfplotflag=True
        lineartoxpercent=0
                
    print('Data files (if any) :',fileio.filelist)
    
    Datalist=[]
    Frequencylist=[]
    for fileid,eachfile in enumerate(fileio.filelist):
        if( eachfile!=None):
            Data=fileio.three_omega_file_reading(eachfile,validate=True) # valiate=True check essential parameters are in the file
            fileio.Nsim=len(Data['x'])  
            fileio.internalf=Data['f']
            if ('range' in fileio.dB[TEMPERATURE]):
                fileio.Nmax=int(float(fileio.dB[TEMPERATURE]['range'])*len(fileio.internalf[0::fileio.downsample]))
            else:
                fileio.Nmax=len(fileio.internalf[0::fileio.downsample])
            fileio.freq=np.array(fileio.internalf[0::fileio.downsample])[0:fileio.Nmax] # select vevery nth element
            scaleT=1#np.sqrt(freq) #1#(freq**.33)
            V3wLockin_rms0=(np.array(Data['x'][0::fileio.downsample][0:fileio.Nmax]) +complex(0,1)* np.array(Data['y'][0::fileio.downsample][0:fileio.Nmax]))           
            if (fileio.dB[AMPLIFIER]['use_calibration'] ==True):
                GainUsed=-1*fileio.CalibrationCoef[0::fileio.downsample][0:fileio.Nmax] 
                offset=0.1
                y0= 1*(V3wLockin_rms0/GainUsed) 
                c=1#c=RC_Correction(freq,50,29.4,0)
                y=y0/c  # this is done on the collected data already
                print('Using Calibration in files: \nInput file:'+fileio.dB[AMPLIFIER]['calibration_in']+'\nOutput file : '+fileio.dB[AMPLIFIER]['calibration_out'])
            else:
                GainUsed=fileio.dB[AMPLIFIER]['gain_nominal']*np.ones(len(V3wLockin_rms0))
                y= V3wLockin_rms0/fileio.dB[AMPLIFIER]['gain_nominal']
                print('Using Nomincal Calibration Gain of : '+str(fileio.dB[AMPLIFIER]['gain_nominal']))
                
            V_1w_heater_rms = np.sqrt(Data['P'] *Data['R'])  # A
            V_3w_heater_rms =y #((R_heater + R_heater + R_source)/ (R_heater));  % V
            VtoTempScale=(2 / (V_1w_heater_rms * fileio.dB[HEATER]['tcr']))
            
            Temperature_experimental = np.array(V_3w_heater_rms*VtoTempScale) # Equation In reference next to file.
            if(type(scaleT)==np.ndarray):
                scale1DT=np.append(scaleT,scaleT)
            else:
                scale1DT=scaleT
            if (1):
                #only fit the data associated with the first file specified
                #Temperature_experimental=(Temperature_experimental) # select vevery nth element
                Tlist=np.append(np.real(Temperature_experimental),np.imag(Temperature_experimental))
                flist=np.append(fileio.freq,fileio.freq)
                Datalist.append(Temperature_experimental)
                Frequencylist.append(fileio.freq)
        else:
            fileio.Nmax=100
            internalf=np.linspace(1, 1e4,fileio.Nmax)
            fileio.freq=np.array(internalf[0::fileio.downsample])[0:fileio.Nmax]
            fileio.setGain()

        if(logfplotflag):
            #if(1):
            dataindex=fileid
            ax.plot(np.log10(Frequencylist[dataindex]),scaleT*np.real(Temperature_experimental),'x',label='Data'+str(dataindex+1)+':real')
            ax.plot(np.log10(Frequencylist[dataindex]),scaleT*np.imag(Temperature_experimental),'x',label='Data'+str(dataindex+1)+':imag' )
        else:
            ax.plot(np.log10(Frequencylist[dataindex]),scaleT*np.abs(Temperature_experimental),'x',label='Data'+str(dataindex+1)+':mag')
            phi=np.array([cmath.phase(x) for x in Temperature_experimental])
            ax.plot(np.log10(Frequencylist[dataindex]),scaleT*phi,'x',label='Data'+str(dataindex+1)+':phase' )

    out=np.array([])
    for modelindex in fileio.dB[LAYERS]:
        model.modelindex=modelindex
        Layers=[]
        guess=[]
        for each in fileio.dB[LAYERS][modelindex]:
            if each.lower()==HEATER:
                Layers.append(fileio.dB[HEATER].copy())
                if('optimize' in  fileio.dB[HEATER]):
                    listopt=fileio.dB[HEATER]['optimize'].split(',')
                    for eachopt in listopt:
                            if(eachopt =='thermal_boundary_resistance'):
                                if('thermal_boundary_resistance' in fileio.dB[HEATER]):
                                    guess.append((fileio.dB[HEATER][eachopt]))
                                else:
                                    raise KeyError("Guess for thermal_boundary_resistance parameter in the heater is missing")
            else:
                if (each.lower() in fileio.dB[MATERIALS]):
                    if (not('thermal_diffusivity' in fileio.dB[MATERIALS][each.lower()])):
                        fileio.dB[MATERIALS][each.lower()]['thermal_diffusivity']=diffusivity(fileio.material(each,'thermal_conductivity'),interpolatecp(fileio.material(each,'temp'),fileio.material(each,'cp')*fileio.material(each,'density')))
                    Layers.append(fileio.dB[MATERIALS][each.lower()].copy())
                    
                    if('optimize' in fileio.dB[MATERIALS][each.lower()]):
                        listopt=fileio.dB[MATERIALS][each.lower()]['optimize'].split(',')
                        for eachopt in listopt:
                            if(eachopt =='thermal_diffusivity'):
                                guess.append((fileio.dB[MATERIALS][each.lower()][eachopt]))
                            else:
                                guess.append(fileio.dB[MATERIALS][each.lower()][eachopt])

                else:
                    raise KeyError("The Layer material name ["+each+ "] could not be found in the database file ["+fileio.yamlfile+"].  The materials found are ["+','.join(fileio.dB[MATERIALS].keys())+"].  Please either edit the Layer defintion of the materials defined i the database")
        if (len(guess)>0):   
            fit=True
            print('Importing Solvers due to presence of optimize parameter in file (this can take a while the first time).....')
            from scipy.optimize import leastsq
            print('Importing leastsq from scipy.optimize')
            from scipy.optimize import curve_fit
            print('Importing curve_fit from scipy.optimize')
            from scipy.optimize import minimize
            print('Importing minimize from scipy.optimize')
        else:
            fit=False
        pcov=np.array([[0,0],[0,0]])
        if(1):

            popt=guess
            FitModel=''
            if fit==True:
                minguess=[]
                maxguess=[]
                for each in guess:
                    minguess.append(each/10) #=[fileio.dB[TEMPERATURE]['guess']['k']['min'],fileio.dB[TEMPERATURE]['guess']['d']['min']]
                    maxguess.append(each*10) #=[k/100,fileio.dB[TEMPERATURE]['guess']['d']['max']]
                if(minguess==[] or maxguess==[]):
                    raise Exception('Exiting......No optimization for any Layer has been specified.  Please edit the input file and ensure the keyword:[optimize] is defined on at least one (1) layer with a specified parameter to optimize (such as thermal_conductivity, thermal_diffusivity, thermal_boundary_resistance)')
            else:
                print('Measured data present but no fitting specified...processing to plot based on current parameters')
            
            if (fileio.dB[TEMPERATURE]['model'].lower().find('borca')>-1):  #Options 'Kim or 'Cahill'
                FitModel='borca'
                
                if fit==True:
                    p_init = guess #[1.,-1.,.5,.5,2.]
                    Tabs=np.abs(Temperature_experimental)
                    model.optimseflag=True
                    model.interation=0
                    popt,pcov = curve_fit(lambda xx,*aa: scale1DT*model.borca_thin_film_fit3_1D(Layers,fileio.freq,*aa) ,flist,  scale1DT*Tlist,guess,bounds=(minguess, maxguess)  )
                    print('\n') # force new line after optimization
                else:
                    popt=[] 
                    pcov=[]    
                model.interation=0
                out=model.borca_thin_film_fit3(Layers,fileio.freq, *popt)

            elif (fileio.dB[TEMPERATURE]['model'].lower().find('kim')>-1):  #Options 'Kim or 'Cahill'
                FitModel='Kim'
                if fit==True:
                    minguess=[]
                    maxguess=[]
                    for each in guess:
                        minguess.append(each/10) #=[fileio.dB[TEMPERATURE]['guess']['k']['min'],fileio.dB[TEMPERATURE]['guess']['d']['min']]
                        maxguess.append(each*10) #=[k/100,fileio.dB[TEMPERATURE]['guess']['d']['max']]
                    #popt,pcov = curve_fit(lambda xx,*aa:  np.log(flist)*kim_thin_film_fit_1D(dB,Layers,freq,*aa) ,flist,  np.log(flist)*Tlist,guess,bounds=(minguess, maxguess)  )
                    model.interation=0
                    popt,pcov = curve_fit(lambda xx,*aa:  model.kim_thin_film_fit_1D(Layers,fileio.freq,*aa) ,flist,  Tlist,guess,bounds=(minguess, maxguess)  )
                    print('\n') # force new line after optimization
                model.interation=0
                out=model.kim_thin_film_fit(Layers,fileio.freq, *popt)

            elif (fileio.dB[TEMPERATURE]['model'].lower().find('cahill')>-1 or fileio.dB[TEMPERATURE]['model'].lower().find('cahil')>-1):  #Options 'Kim 
                FitModel='cahill'
                if fit==True:
                    minguess=[]
                    maxguess=[]
                    for each in guess:
                        minguess.append(each/10) #=[fileio.dB[TEMPERATURE]['guess']['k']['min'],fileio.dB[TEMPERATURE]['guess']['d']['min']]
                        maxguess.append(each*10) #=[k/100,fileio.dB[TEMPERATURE]['guess']['d']['max']]
                    model.interation=0
                    popt,pcov = curve_fit(lambda xx, *aa : scale1DT*model.Bulk_fit2_1D(Layers,fileio.freq,*aa),flist, scale1DT*Tlist,guess,bounds=(minguess, maxguess))
                    print('\n') # force new line after optimization
                else:
                    popt=[] 
                    pcov=[]
                model.interation=0
                out=model.Bulk_fit2(Layers,fileio.freq, *popt)
            if(len(out)==0):
                raise KeyError('No  known model was identified - check the databas file, valid models are '+'or'.model.modelnames)

            model.Layers.append(Layers)
            

            if(logfplotflag):
                ax.plot(np.log10(fileio.freq),scaleT*np.real(out),'-',label='Stack#'+str(modelindex+1)+':fit real')
                ax.plot(np.log10(fileio.freq),scaleT*np.imag(out),'-',label='Stack#'+str(modelindex+1)+':fit imag')
                ax.plot(np.log10(fileio.freq),np.zeros(len(fileio.freq)),'-.')
            else:

                ax.plot(np.log10(fileio.freq),scaleT*np.abs(out),'-',label='Model'+str(modelindex+1)+':fit mag')
                theta= np.array([cmath.phase(x) for x in out])
                ax.plot(np.log10(fileio.freq),scaleT*theta,'-',label='Stack#'+str(modelindex+1)+':fit phase')
                rediff=np.abs(out)-np.abs(Temperature_experimental)
                imdiff=theta-phi
                ax.plot(np.log10(fileio.freq),10*scaleT*(rediff),'-.',label='Stack#'+str(modelindex+1)+':DIFF magx10')
                ax.plot(np.log10(fileio.freq),10*scaleT*(imdiff),'-.',label='Stack#'+str(modelindex+1)+':DIFF phasex10')
                #ax.plot(np.log10(fileio.freq),10*scaleT*(imdiff+rediff),'-.',label='Model'+str(modelindex+1)+':DIFF+DIFF imagx10')
                ax.plot(np.log10(fileio.freq),np.zeros(len(fileio.freq)),'-.')

                # else:
                #     smax=max(np.abs(out))
                #     outx=[out[i] for i,each in enumerate(abs(out)) if each>smax*lineartoxpercent]
                #     newfreq=freq[0:len(outx)]
                #     if 'data' in fileio.dB[TEMPERATURE] :
                #         Temperature_experimentalx=Temperature_experimental[0:len(outx)]
                #         ax.plot( newfreq,scaleT*np.real(Temperature_experimentalx),'x',label='Data'+str(modelindex+1)+':real')
                #         ax.plot( newfreq ,scaleT*np.imag(Temperature_experimentalx),'x',label='Data'+str(modelindex+1)+':imag' )
                #     if fit==True:
                #         ax.plot( newfreq,scaleT*np.real(outx),'-',label='Stack#'+str(modelindex+1)+':fit real')
                #         ax.plot( newfreq,scaleT*np.imag(outx),'-',label='Stack#'+str(modelindex+1)+':fit imag')
                #         ax.plot( newfreq,np.zeros(len(newfreq)),'-.')

   
                #Call class to build all the parameters found or used for each layer
            ThermalMsg=FormatThermalData(Layers,popt,pcov,fit)  
            feature='output'

            vFlag=False
            if(feature in fileio.dB[TEMPERATURE]):
                if(fileio.dB[TEMPERATURE][feature].lower().find('v')>-1):vFlag=True
            outfile=configfile.split('.')[0]+'_'+FitModel+'-'+'Model'+str(modelindex+1)+'.txt'
            with open(outfile, mode='w') as csv_file:
                fieldnames = ['Parameters', 'Value',  'Unit','+/-Uncertainty%']
                writer = csv.DictWriter(csv_file, fieldnames=fieldnames, lineterminator="\n")
                writer.writeheader()

                for paramstr in ThermalMsg.logdata:
                    outdict={}
                    for jj,param in enumerate(paramstr):
                        outdict.update({fieldnames[jj]: str(param)})
                    writer.writerow(outdict)
                        #, fieldnames[1]:np.real(Temperature_experimental)[ii], fieldnames[2]:np.real(out)[ii],fieldnames[3]: np.imag(Temperature_experimental)[ii], fieldnames[4]:np.imag(out)[ii]})

                if vFlag:
                    if 'data' in fileio.dB[TEMPERATURE] :
                        fieldnames = ['frequency (Hz)', 'measured_real_v',  'model_real_v','measured_imag_v', 'model_imag_v']
                    else:
                        fieldnames = ['frequency (Hz)', 'model_real_v','model_imag_v']
                        VtoTempScale=1
                    Gain=GainUsed
                else:                    
                    if 'data' in fileio.dB[TEMPERATURE] :
                        fieldnames = ['frequency (Hz)', 'measured_real_T',  'model_real_T','measured_imag_T', 'model_imag_T']
                    else:
                        fieldnames = ['frequency (Hz)',   'model_real_T','model_imag_T']
                    Gain=np.ones(len(fileio.freq))
                    VtoTempScale=1
                writer = csv.DictWriter(csv_file, fieldnames=fieldnames, lineterminator="\n")
            
                writer.writeheader()
                ModelOut=out*Gain/VtoTempScale
                if 'data' in fileio.dB[TEMPERATURE] :
                    DataOut=Temperature_experimental #*Gain/VtoTempScale
                    for ii,freqi in enumerate(fileio.freq):
                        writer.writerow({fieldnames[0]: freqi, fieldnames[1]:np.real(DataOut)[ii], fieldnames[2]:np.real(ModelOut)[ii],fieldnames[3]: np.imag(DataOut)[ii], fieldnames[4]:np.imag(ModelOut)[ii]})
                else:
                    for ii,freqi in enumerate(fileio.freq):
                        writer.writerow({fieldnames[0]: freqi, fieldnames[1]:np.real(ModelOut)[ii], fieldnames[2]:np.imag(ModelOut)[ii]})

            print('Thermal Data:\n', ThermalMsg.msg)
            if(1):
                model.plotlayers(ax2,ThermalMsg.msg,'Layer stack#'+str(modelindex+1)+':'+'/'.join(fileio.dB[LAYERS][modelindex][1:])+' fit using model:'+FitModel)
                fig.canvas.draw()
                fig.canvas.flush_events()
    ax.legend( loc='best', fontsize=model.font)
    plt.show()
 
    print('Press Enter to Exit')            
    #input()


if __name__ == '__main__':
   yref=.9 # defines the position of the plotted stack starting from the top 
   print('Starting.....')
   np.seterr(all='raise')
   n = len(sys.argv) 
   #print('ARGS',n,sys.argv)
   if n>1:
       filename=sys.argv[1]
   else:
       filename='thin.yaml'
   main(filename)
   try:
       pass
       
   except yaml.YAMLError as e:
        msg="Error while parsing YAML file:"
        if hasattr(e, 'problem_mark'):
            if e.context != None:
                print ('  parser says\n' + str(e.problem_mark) + '\n  ' +
                    str(e.problem) + ' ' + str(e.context) +
                    '\nPlease correct data and retry.')
                msg=msg+'  parser says\n' + str(e.problem_mark) + '\n  ' +str(e.problem) + ' ' + str(e.context) +'\n\nMost likely error is that you did not leave a space after the colon following the keyword on the line indicated above\n\nPlease correct data and retry.'
            else:
                print ('  parser says\n' + str(e.problem_mark) + '\n  ' +
                    str(e.problem) + '\nPlease correct data and retry.')
                msg=msg+'  parser says\n' + str(e.problem_mark) + '\n  ' +str(e.problem) + '\n\nMost likely error is that you did not include 2-spaces before the keyword on the line indicated above\n\nPlease correct data and retry.'
        else:
            print ("Something went wrong while parsing yaml file")
            msg=msg+"\nSomething went wrong while parsing yaml file\n"+str(e)
        ctypes.windll.user32.MessageBoxW(0, msg, 'Error', 0)
   except KeyError as e:
        ctypes.windll.user32.MessageBoxW(0, 'File Key is incorrect.  The keyword: '+str(e)+' is missing', 'File Key is incorrect ', 0)
   except FileNotFoundError as e:
        ctypes.windll.user32.MessageBoxW(0, 'The Configuration File :\n'+str(configfile)+'\nis missing.  Please check the name and path', 'File missing', 0)
   except Warning as e:
        ctypes.windll.user32.MessageBoxW(0, 'A Warning has occured:\n'+str(e)+'\nExiting', 'Warning', 0)
   except Exception as e:
        ctypes.windll.user32.MessageBoxW(0, 'An Error has occured:\n'+str(e)+'\nExiting', 'Warning', 0)