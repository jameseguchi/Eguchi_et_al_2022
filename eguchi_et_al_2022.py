#Script to run C cycle model from Eguchi, Diamond, Lyons 2022, PNAS Nexus.
#Model was designed to investigate how differential mantle residence times of 
#organic C vs carbonates affect d13C of volcanic C outgassing and 
#marine carbonates


import numpy as np
import matplotlib.pyplot as plt


# time domain
# model time domain
t0 = 0       # Model start time [Myr]   
tf = 5000    # Model end time [Myr]
t  = np.linspace(t0,tf,(tf-t0+1))
tnum=len(t)
ticker=5000 #this just so years can be changed more easily
tchange1=ticker-2380#timing for perturbation for GOE and Lomagundi[Myr]
tchange2=ticker-1400#timing of second tectonic change - Breakup of Nuna [Myr]
tchange3=ticker-800#timing of third change - Break up of Rodinia for NOE [Myr]
crb_tau=30 #delay time for release of carbonates at arcs [Myr]
org_tau=310#delay time for release of organic C at OIBs at 1st perturbation[Myr]
org_tau2=400#delay time for release of organic C at OIBs at 1st perturbation[Myr]
org_tau3=230#delay time for release of organic C at OIBs at 1st perturbation[Myr]

# Constants
forg=0.2 #fraction of C from atmosphere-ocean as deposited as organics
xcrb=1-forg #fraction of C from atmosphere-ocean as deposited as carbonate
chi_org_i=0.2# fraction of organics that are subducted from the surface
chi_carb_i=chi_org_i # fraction of carbonates that are subducted from the surface
alpha_org=0.0 # fraction of subducted organics that come out at arcs
alpha_crb=1.0 # fraction of subducted carbonates that come out at arcs
orgc4=0.0 #fraction of subducted organics that remain in mantle
crbc4=0.0 # fraction that remains in the mantle


# initial conditions
c_atmi=0 # mass of atm-ocean C reservoir [g]
c_crbi=0 # mass of crustal carbonate C reservoir [g]
c_orgi=0 # mass of crustal organic C reservoir [g]
c_mcrbi=0 # mass of mantle carbonate C reservoir [g]
c_morgi=0 # mass of mantle organic C reservoir [g]
c_mntli=1e23 # mass of mantle primordial C reservoir [g]
d13C_atmi=-5 # d13c of C in atm-ocean [permill]
d13C_crbi=0  # d13c of C in carbonate [permill]
d13C_orgi=-25 # d13c of C in organic C [permill]
d13C_prim=-5 # d13c of C in primitive mantle C [permill]


F_mori=1e18#intial MORB C flux [g/Myr]
F_oibi=1e18#initial OIB C flux [g/Myr]
F_arci=1e18#initial arc C flux [g/Myr]

#Parameter changes to drive model perturbations

morb_change1=1e18 # MORB flux after Initial perturbation [g/Myr]
morb_change2=1e18 # MORB flux after 2nd perturbation [g/Myr]
morb_change3=1e18 # MORB flux after 3rd perturbation [g/Myr]


ki=1e-9 # weathering constant initial
k_change1=2e-3 #increased weathering constant at 1st perturbation
k_change2=1e-5 #decayed weathering constant after 1st perturbation
k_change3=2e-4 #increased weathering constant after 2nd pertubation
k_change4=2e-5 #decayed weathering constant after 2nd perturbation
k_change5=9e-4 #increased weathering constant after 3rd perturbation
k_change6=7e-5 #decayed weathering constant after 3rd perturbation
k_change7=1.2e-3 #final increase in k

decay_len1=350 #k decay length at 1st perturbation
decay_len2=300 #k decay length at 2nd perturbation
decay_len3=200 #k decay length at 3rd perturbation
decay_len4=100 #k decay length at for final increas in k


chi_org_i=0.1 #fraction of organics that are subducted from the surface prior to first perturbation
chi_carb_i=chi_org_i #fraction of carbonates that are subducted from the surface first pertubation
chi_carb_2=.999 #fraction of carbonates that are subducted from the surface during increased k periods
chi_org_2=.999 #fraction of organics that are subducted from the surface during increased k periods
chi_org_3=.1 #fraction of organics that are subducted from the surface during final k increase
chi_carb_3=.1 #fraction of organics that are subducted from the surface during final k increase


#------oxidation of org C at surface-------------------------------------------
k17=7.75e12 #weathering constant molC/yr
k17=k17*12*1e6 # g C/Myr
G=1.25e21 #initial org C [mol] from bergman COPSE
G=G*12#G [g]
U=1 #uplift rate
pres_o2=5.15e18*1000*.21 #mass of  present day O2


#-----------No Need to change anything after this line-------------------------
#initiate arrays for parameters

#initial fluxes of org c and carbonate deposition
F_orgi=forg*ki*c_atmi
F_crbi=xcrb*ki*c_atmi
#initial fluxes of org c and carbonate subduction
F_sorgi=0
F_scrbi=0

# reservoirs
c_atm=np.zeros(tnum) #atmosphere-ocean c reservoir [g]
c_atm[0]=c_atmi
c_crb=np.zeros(tnum) #continental carbonate reservoir [g]
c_crb[0]=c_crbi
c_org=np.zeros(tnum) #continental organic C reservoir [g]
c_org[0]=c_orgi
c_mcrb=np.zeros(tnum) #mantle subducted carbonate reservoir [g]
c_morg=np.zeros(tnum) #mantle subducted organic C reservoir [g]
c_mntl=np.zeros(tnum) #mantle primitive organic C reservoir [g]
c_mntl[0]=c_mntli

mol_orgc=np.zeros(tnum) #moles of organic C in all reservoirs [moles]
g_o2=np.zeros(tnum) #mass of o2 in ocean-atmosphere [g]
pal=np.zeros(tnum) #ocean-atmosphere o2 [PAL]
calc_o2=np.zeros(tnum) #mass of O2 calculated from inputs and sinks [g]

# isotopes
d13C_atm=np.zeros(tnum) #d13C of atm-ocean [per mil]
d13C_atm[0]=d13C_atmi
d13C_crb=np.zeros(tnum) #d13C of marine carbonate [per mil]
d13C_crb[0]=d13C_crbi
d13C_org=np.zeros(tnum) #d13C of organic C [per mil]
d13C_org[0]=d13C_orgi
d13C_prm=np.ones(tnum) #d13C primitive mantle C[per mil]
d13C_prm=d13C_prm*d13C_prim
d13C_arc=np.zeros(tnum) #d13C arc C [per mil]
d13C_arc[0]=d13C_prim 
d13C_oib=np.zeros(tnum) #d13C oib c[per mil]
d13C_oib[0]=d13C_prim
d13C_mor=np.ones(tnum) #d13C of MOR c [per mil]
d13C_mor=d13C_mor*d13C_prim


# fluxes
F_org=np.zeros(tnum) #flux of organic C production at surface [g/Myr]
F_org[0]=F_orgi
F_crb=np.zeros(tnum) #flux of carbonate production at surface [g/Myr]
F_crb[0]=F_crbi
F_sorg=np.zeros(tnum) #flux of organic C subduction [g/Myr]
F_sorg[0]=F_sorgi
F_scrb=np.zeros(tnum) #flux of carbonate subduction [g/Myr]
F_scrb[0]=F_scrbi
F_mor=np.ones(tnum) #flux of MOR subduction [g/Myr]
F_mor=F_mor*F_mori 
F_orgw=np.zeros(tnum) #flux of org C weathering [g/Myr]
F_oib=np.zeros(tnum) # flux of oib C [g/Myr]
F_oib[0]=F_oibi
F_arc=np.zeros(tnum) # flux of arc C [g/Myr]
F_arc[0]=F_arci
F_tot=np.zeros(tnum) # total flux of org C [g/Myr]
Foiborg=np.zeros(tnum) # flux of org C at oibs [g/Myr]
F_tot[0]=F_oib[0]+F_arc[0]+F_mor[0]

#strenghth of weathering feedback
k=np.zeros(tnum)

#fraction of C subducted
chi_carb=np.ones(tnum)
chi_carb=chi_carb*chi_carb_i
chi_org=np.ones(tnum)
chi_org=chi_org*chi_org_i



# iterate through model
for time in range(1,tnum):
           
    #this block of code changes values of model drivers during perturbations
    if (t0+time) <tchange1:#time before first perturbation
        k[time] = ki
        F_mor[time] = F_mori
        chi_org[time]=chi_org_i
        chi_carb[time]=chi_carb_i
    elif (t0+time) == (tchange1):#1st perturbation
        k[time]=k_change1
        F_mor[time]=morb_change1
        chi_org[time]=chi_org_2
        chi_carb[time]=chi_carb_2        
    elif (t0+time)>(tchange1) and (t0+time) < (tchange1+decay_len1):#decay time of 1st perturbation
        k[time]=k[time-1]-(k_change1-k_change2)/(decay_len1)
        F_mor[time]=morb_change1
        chi_carb[time]=chi_carb_2
        chi_org[time]=chi_org_2
    elif (t0+time)>=tchange1+decay_len1 and ((t0+time)<tchange2):#time between decay of 1st perturbation and 2nd perturbation
        k[time]=k_change2        
        F_mor[time]=morb_change1
        chi_org[time]=chi_org_i
        chi_carb[time]=chi_carb_i
    elif (t0+time)==tchange2: #time of 2nd perturbation
        k[time]=k_change3
        F_mor[time]=morb_change2
        org_tau=org_tau2
        chi_org[time]=chi_org_2
        chi_carb[time]=chi_carb_2
    elif (t0+time)>tchange2 and  (t0+time) < tchange2+decay_len2: #decay time of 2nd perturbation
        k[time]=k[time-1]-(k_change3-k_change4)/(decay_len2)        
        F_mor[time]=morb_change2
        org_tau=org_tau2
        chi_org[time]=chi_org_2
        chi_carb[time]=chi_carb_2
    elif(t0+time)>=tchange2+decay_len2 and (t0+time) < tchange3: #time between decay of 2nd perturbation and 3rd perturbation
        k[time]=k_change4       
        F_mor[time]=morb_change2
        org_tau=org_tau2
        chi_org[time]=chi_org_i
        chi_carb[time]=chi_carb_i
    elif(t0+time) == tchange3: #time of 3rd perturbation
        k[time]=k_change5        
        F_mor[time]=morb_change3
        org_tau=org_tau3
        chi_org[time]=chi_org_2
        chi_carb[time]=chi_carb_2
    elif (t0+time) > tchange3 and  (t0+time) <= tchange3 + decay_len3: #decay time of 3rd perturbation
        k[time]=k[time-1]-(k_change5-k_change6)/(decay_len3)      
        F_mor[time]=morb_change3
        org_tau=org_tau3
        chi_carb[time]=chi_carb_2
        chi_org[time]=chi_org_2
    elif (t0+time) > tchange3 + decay_len3 and (t0+time) <= tchange3 + decay_len3 + decay_len4: #increase time for last k increase
        k[time]=k[time-1]-(k_change6-k_change7)/(decay_len4) 
        F_mor[time]=morb_change3
        org_tau=org_tau3
        chi_carb[time]=chi_carb_3
        chi_org[time]=chi_org_3
    elif (t0+time) > tchange3 + decay_len4: #time after last k increase
        k[time]=k_change7     
        F_mor[time]=morb_change3
        org_tau=org_tau3
        chi_carb[time]=chi_carb_3
        chi_org[time]=chi_org_3
   
    
    
    atm_Fout=F_org[time-1]+F_crb[time-1]
    c_atm[time]=c_atm[time-1]+(F_tot[time-1]-(atm_Fout))+F_orgw[time-1]
    c_crb[time]=c_crb[time-1]+(F_crb[time-1]-F_scrb[time-1])
    c_org[time]=c_org[time-1]+(F_org[time-1]-F_sorg[time-1])-F_orgw[time-1]
    c_mcrb[time]=c_mcrb[time-1]+F_scrb[time-1]-(1-crbc4)*alpha_crb*F_scrb[time-crb_tau]-(1-crbc4)*(1-alpha_crb)*F_scrb[time-org_tau]
    c_morg[time]=c_morg[time-1]+F_sorg[time-1]-(1-orgc4)*(1-alpha_org)*F_sorg[time-org_tau]
    c_mntl[time]=c_mntl[time-1]-F_oibi-F_mor[time-1]-F_arci
    
     
    F_org[time]=forg*k[time]*c_atm[time]
    F_crb[time]=xcrb*k[time]*c_atm[time]
    F_scrb[time]=chi_carb[time]*F_crb[time]
    F_sorg[time]=chi_org[time]*F_org[time]    
    F_orgw[time]=(k17*U)*(c_org[time-1]/G)    
    if (t0 + time) <= (crb_tau):
        Farccrb=0
        Farcorg=0
    else:
        Farccrb=alpha_crb*F_scrb[time-crb_tau]
        Farcorg=alpha_org*F_sorg[time-crb_tau]   
    if (t0 + time) <= (org_tau):
        Foibcrb=0
        Foiborg[time]=0
    else :
        Foibcrb=(1-alpha_crb)*F_scrb[time-org_tau]
        Foiborg[time]=(1-alpha_org)*F_sorg[time-org_tau]    
    F_oib[time]=(Foibcrb+Foiborg[time]+F_oibi)
    F_arc[time]=F_arci+Farccrb+Farcorg    
    F_tot[time]=F_oib[time]+F_arc[time]+F_mor[time]

    d13C_oiborg=Foiborg[time]/F_oib[time]*d13C_org[time-org_tau]*(1-orgc4)
    d13C_oibcrb=Foibcrb/F_oib[time]*d13C_crb[time-org_tau]
    d13C_oibmntl=F_oibi/F_oib[time]*d13C_prm[time]
    d13C_oib[time]=d13C_oiborg+d13C_oibcrb+d13C_oibmntl
    d13C_arcorg=Farcorg/F_arc[time]*d13C_org[time-crb_tau]
    d13C_arccrb=Farccrb/F_arc[time]*d13C_crb[time-crb_tau]
    d13C_arcmntl=F_arci/F_arc[time]*d13C_prm[time]
    d13C_arc[time]=d13C_arcorg+d13C_arccrb+d13C_arcmntl
    d13C_atm[time]=F_oib[time]/F_tot[time]*d13C_oib[time]+F_arc[time]/F_tot[time]*d13C_arc[time]+F_mor[time]/F_tot[time]*d13C_mor[time]
    d13C_crb[time]=d13C_atm[time]+5
    d13C_org[time]=d13C_atm[time]-20
    
    mol_orgc[time]=(c_org[time]+c_morg[time])/12.011
    g_o2[time]=mol_orgc[time]*15.999*2
    pal[time]=g_o2[time]/pres_o2


   



tmin=0
tmax=ticker
tnew=ticker-t
tmax_plot=3500



# visualize data
#carbon isotopes
fig=plt.figure(1,[4,12])
plt.subplot(6,1,1)
plt.plot(tnew,d13C_crb,'orange',label='Marine Carb')
plt.plot(tnew,d13C_arc,label='Arc')
plt.plot(tnew,d13C_oib,label='OIB')
plt.xlim([tmin,tmax_plot])
ax = plt.gca()
ax.invert_xaxis()
plt.ylim([-20,20])
plt.legend(frameon=False)
plt.ylabel(r'$\delta^{13}C$')

#oxygen
plt.subplot(6,1,2)
plt.semilogy(tnew, pal)
plt.semilogy(tnew, calc_o2/pres_o2)
plt.xlim([tmin,tmax_plot])
plt.ylim([1e-7,1e1])
ax = plt.gca()
ax.invert_xaxis()
plt.yticks([1e-6, 1e-4, 1e-2, 1e0])
plt.ylabel(r'O$_2$ [PAL]')

#C fluxes
plt.subplot(6,1,3)
plt.semilogy(tnew,F_oib,'k',label='oib')
plt.semilogy(tnew,F_arc,'b',label='arc')
plt.semilogy(tnew,F_mor,'y',label='mor')
plt.semilogy(tnew,F_scrb+F_sorg, label='Subduction')
plt.semilogy(tnew,F_org+F_crb-(F_scrb+F_sorg), label='Continental')
plt.xlim([tmin,tmax_plot])
plt.ylim([1e11, 1e20])
ax = plt.gca()
ax.invert_xaxis()
plt.legend(frameon=False)
plt.ylabel('C flux [g/Myr]')

#C reservoirs
plt.subplot(6,1,4)
plt.semilogy(tnew,c_org,'r',label='c crustal org')
plt.semilogy(tnew,c_morg,'g',label='c mantle org')
plt.semilogy(tnew,c_crb,'m',label='crustal carb')
plt.semilogy(tnew,c_mcrb,'c',label='mantle carb')
plt.legend(frameon=False)
plt.xlim([tmin,tmax_plot])
plt.ylim([1e12, 1e23])
plt.ylabel('C reservoir [g]')
ax = plt.gca()
ax.invert_xaxis()

#k and chi (model perturbation drivers)
plt.subplot(6,1,5)
plt.semilogy(tnew,k)
plt.xlim([tmin,tmax_plot])
plt.ylim([1e-9, 1])
plt.xlabel('Million Years Ago')
ax = plt.gca()
ax2=ax.twinx()
ax.set_ylabel('k')
ax2.plot(tnew,chi_carb,'orange')
ax2.set_ylabel(r'$\chi$ Fraction of C subducted')
ax = plt.gca()
ax.invert_xaxis()




plt.show()




