import numpy as np


TYPE1='vespaBC_1_gals_SFH6HSsig_dr7_sp2'
TYPE2='LINksmhigh22_sSFRzcut_parent'

raGe,decGe,zGe,jackGe=np.loadtxt('{}_jack.dat'.format(TYPE1),unpack=True)
raGl,decGl,zGl,jackGl=np.loadtxt('{}_jack.dat'.format(TYPE2),unpack=True)
raRe,decRe,zRe,jackRe=np.loadtxt('/uufs/astro.utah.edu/common/uuastro/astro_data/zhengzheng/mccarthy/GravBinary_Data/random_{}_jack.dat'.format(TYPE1),unpack=True)
raRl,decRl,zRl,jackRl=np.loadtxt('random_{}_jack.dat'.format(TYPE2),unpack=True)

#lam=0.73
#mat=0.27

lam=0.692885
mat=0.307115

dGe=2000*(np.sqrt(lam+mat*(1+3*zGe))-1)/mat
dGl=2000*(np.sqrt(lam+mat*(1+3*zGl))-1)/mat
dRe=2000*(np.sqrt(lam+mat*(1+3*zRe))-1)/mat
dRl=2000*(np.sqrt(lam+mat*(1+3*zRl))-1)/mat

datae=np.vstack((raGe,decGe,dGe,jackGe)).T
datal=np.vstack((raGl,decGl,dGl,jackGl)).T

ra_Ge=raGe*np.pi/180.
dec_Ge=np.pi/2.-decGe*np.pi/180.
ra_Gl=raGl*np.pi/180.
dec_Gl=np.pi/2.-decGl*np.pi/180.
ra_Re=raRe*np.pi/180.
dec_Re=np.pi/2.-decRe*np.pi/180.
ra_Rl=raRl*np.pi/180.
dec_Rl=np.pi/2.-decRl*np.pi/180.

zmaxcut=0
if zmaxcut==1:
	cutg_early=dGe<=270.
	cutr_early=dRe<=270.
        cutg_late=dGl<=270.
        cutr_late=dRl<=270.
	dGe=dGe[cutg_early]
	dRe=dRe[cutr_early]
	ra_Ge=ra_Ge[cutg_early]
	ra_Re=ra_Re[cutr_early]
        dec_Ge=dec_Ge[cutg_early]
        dec_Re=dec_Re[cutr_early]
        dGl=dGl[cutg_late]
        dRl=dRl[cutr_late]
        ra_Gl=ra_Gl[cutg_late]
        ra_Rl=ra_Rl[cutr_late]
        dec_Gl=dec_Gl[cutg_late]
        dec_Rl=dec_Rl[cutr_late]
	jackGe=jackGe[cutg_early]
	jackGl=jackGl[cutg_late]
	jackRe=jackRe[cutr_early]
	jackRl=jackRl[cutr_late]

X_Ge=dGe*np.sin(dec_Ge)*np.cos(ra_Ge)
Y_Ge=dGe*np.sin(dec_Ge)*np.sin(ra_Ge)
Z_Ge=dGe*np.cos(dec_Ge)

X_Gl=dGl*np.sin(dec_Gl)*np.cos(ra_Gl)
Y_Gl=dGl*np.sin(dec_Gl)*np.sin(ra_Gl)
Z_Gl=dGl*np.cos(dec_Gl)

X_Re=dRe*np.sin(dec_Re)*np.cos(ra_Re)
Y_Re=dRe*np.sin(dec_Re)*np.sin(ra_Re)
Z_Re=dRe*np.cos(dec_Re)

X_Rl=dRl*np.sin(dec_Rl)*np.cos(ra_Rl)
Y_Rl=dRl*np.sin(dec_Rl)*np.sin(ra_Rl)
Z_Rl=dRl*np.cos(dec_Rl)

dataGe_jack=np.vstack((X_Ge,Y_Ge,Z_Ge,jackGe)).T
dataGl_jack=np.vstack((X_Gl,Y_Gl,Z_Gl,jackGl)).T
dataRe_jack=np.vstack((X_Re,Y_Re,Z_Re,jackRe)).T
dataRl_jack=np.vstack((X_Rl,Y_Rl,Z_Rl,jackRl)).T

dataGe=np.vstack((X_Ge,Y_Ge,Z_Ge)).T
dataGl=np.vstack((X_Gl,Y_Gl,Z_Gl)).T
dataRe=np.vstack((X_Re,Y_Re,Z_Re)).T
dataRl=np.vstack((X_Rl,Y_Rl,Z_Rl)).T


#np.savetxt('earlyLIN_SFH_CART.dat',dataGe, fmt=['%5e','%5e','%5e'])
np.savetxt('{}_CART.dat'.format(TYPE1),dataGe, fmt=['%5e','%5e','%5e'])
np.savetxt('{}_CART_jack.dat'.format(TYPE1),dataGe_jack, fmt=['%5e','%5e','%5e','%5d'])

#np.savetxt('lateLIN_SFH_CART.dat',dataGl, fmt=['%5e','%5e','%5e'])
#np.savetxt('{}_CART.dat'.format(TYPE2),dataGl, fmt=['%5e','%5e','%5e'])
#np.savetxt('{}_CART_jack.dat'.format(TYPE2),dataGl_jack, fmt=['%5e','%5e','%5e','%5d'])

#np.savetxt('random_SFH_earlyLIN_CART.dat',dataRe, fmt=['%5e','%5e','%5e'])
np.savetxt('/uufs/astro.utah.edu/common/uuastro/astro_data/zhengzheng/mccarthy/GravBinary_Data/random_{}_CART.dat'.format(TYPE1),dataRe, fmt=['%5e','%5e','%5e'])
np.savetxt('/uufs/astro.utah.edu/common/uuastro/astro_data/zhengzheng/mccarthy/GravBinary_Data/random_{}_CART_jack.dat'.format(TYPE1),dataRe_jack, fmt=['%5e','%5e','%5e','%5d'])

#np.savetxt('random_SFH_lateLIN_CART.dat',dataRl, fmt=['%5e','%5e','%5e'])
#np.savetxt('random_{}_CART.dat'.format(TYPE2),dataRl, fmt=['%5e','%5e','%5e'])
#np.savetxt('random_{}_CART_jack.dat'.format(TYPE2),dataRl_jack, fmt=['%5e','%5e','%5e','%5d'])
