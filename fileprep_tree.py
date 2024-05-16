import numpy as np


choice=int(input('For xi(s,u) or xi(rp,pi) (0 or 1) : '))
choice1=int(input('Random or data catalog (0 or 1): '))
choice2=int(input('Jackknife sample? (0=No 1=Yes): '))
data=raw_input('Catalog selection (early or late): ')
Bin=int(input('Bin number for s and rp (12 or 19 or 20): '))
eta=input('Eta: ')

if data=='early':
	TYPE='vespaBC_1_gals_SFH6HSsig_dr7_sp2'
else:
	TYPE='LINksmhigh22_sSFRzcut_parent'

TYPE1='dr7'

from datetime import datetime
print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\nProgram started at {}.\n---------------------------'.format(datetime.now().time()))


print('Started data load at {}'.format(datetime.now().time()))
if choice2==0:
	if choice1==0:
		dx,dy,dz=np.loadtxt('/uufs/astro.utah.edu/common/uuastro/astro_data/zhengzheng/mccarthy/GravBinary_Data/random_{}_CART.dat'.format(TYPE),unpack=True)
	if choice1==1:
        	dx,dy,dz=np.loadtxt('{}_CART.dat'.format(TYPE),unpack=True)
if choice2==1:
        if choice1==0:
                dx,dy,dz,jack=np.loadtxt('/uufs/astro.utah.edu/common/uuastro/astro_data/zhengzheng/mccarthy/GravBinary_Data/random_{}_CART_jack.dat'.format(TYPE),unpack=True)
        if choice1==1:
                dx,dy,dz,jack=np.loadtxt('{}_CART_jack.dat'.format(TYPE),unpack=True)

if choice==0:
	binNum=np.loadtxt('suCELLnum_{2}_{0}tree_{1}S.dat'.format(Bin,eta,TYPE1),unpack=False).astype(int)
	binNEI=np.loadtxt('suCELLnei_{2}_{0}tree_{1}S.dat'.format(Bin,eta,TYPE1),unpack=False)
if choice==1:
	binNum=np.loadtxt('rppiCELLnum_{2}_{0}tree_{1}Rp.dat'.format(Bin,eta,TYPE1),unpack=False).astype(int)
        binNEI=np.loadtxt('rppiCELLnei_{2}_{0}tree_{1}Rp.dat'.format(Bin,eta,TYPE1),unpack=False)
N=len(dx)
Nc=len(binNEI)

print('Finished data load at {}.\n---------------------------'.format(datetime.now().time()))

print('<< Digitizing... >>')
if choice==0:
# For su calculations
	if Bin==12:	
		lgSmax=1.4
	if Bin==19:
		lgSmax=1.8
	Lbox=np.power(10,lgSmax)
if choice==1:
# For rppi calculation
        if Bin==12:
                pi_max=40
                Lbox=pi_max
        if Bin==19:
                lgSmax=1.8
                Lbox=np.power(10,lgSmax)
	if Bin==20:
		lgSmax=20
		Lbox=lgSmax
Nbins=np.ceil(np.power(Nc,(1/3.))).astype(int)+1
rangeMin=-Nbins/2.*Lbox
rangeMax=Nbins/2.*Lbox
bins=np.linspace(rangeMin,rangeMax,(Nbins+1))
digx=np.digitize(dx,bins)
digy=np.digitize(dy,bins)
digz=np.digitize(dz,bins)
	
print('Started cell identification and sorting at {}.'.format(datetime.now().time()))
IND=np.zeros(N)
NgalID=np.linspace(0,N-1,N).astype(int)
for i in range(0,Nc):
	if i < Nc:
		digytmp=digy[digx==binNum[i,0]]
 		digztmp=digz[digx==binNum[i,0]]
		dig=NgalID[digx==binNum[i,0]]
		digztmp2=digztmp[digytmp==binNum[i,1]]
		dig2=dig[digytmp==binNum[i,1]]
		dig3=dig2[digztmp2==binNum[i,2]]
		IND[dig3]=(i+1)

if choice2==0:		
	coordsort=np.vstack((dx,dy,dz,IND)).T
if choice2==1:
        coordsort=np.vstack((dx,dy,dz,jack,IND)).T

coordsort=coordsort[coordsort[:,-1].argsort()]
dx=coordsort[:,0]
dy=coordsort[:,1]
dz=coordsort[:,2]
if choice2==0:
	IND=coordsort[:,3].astype(int)
	coords=np.vstack((dx,dy,dz)).T
if choice2==1:
	jack=coordsort[:,3].astype(int)
        IND=coordsort[:,4].astype(int)
        coords=np.vstack((dx,dy,dz,jack)).T

if choice2==0:
	if choice1==1:
		if choice==0:
        		np.savetxt('{0}_suCALC_{1}tree_{2}.dat'.format(TYPE,Bin,TYPE1),coords)
		if choice==1:
        		np.savetxt('{0}_rppiCALC_{1}tree_{2}.dat'.format(TYPE,Bin,TYPE1),coords)
	if choice1==0:
        	if choice==0:
                	np.savetxt('/uufs/astro.utah.edu/common/uuastro/astro_data/zhengzheng/mccarthy/GravBinary_Data/random_{0}_suCALC_{1}tree_{2}.dat'.format(TYPE,Bin,TYPE1),coords)
        	if choice==1:
                	np.savetxt('/uufs/astro.utah.edu/common/uuastro/astro_data/zhengzheng/mccarthy/GravBinary_Data/random_{0}_rppiCALC_{1}tree_{2}.dat'.format(TYPE,Bin,TYPE1),coords)
if choice2==1:
        if choice1==1:
                if choice==0:
                        np.savetxt('{0}_jack_suCALC_{1}tree_{2}.dat'.format(TYPE,Bin,TYPE1),coords,fmt=['%7f','%7f','%7f','%7d'])
                if choice==1:
                        np.savetxt('{0}_jack_rppiCALC_{1}tree_{2}.dat'.format(TYPE,Bin,TYPE1),coords,fmt=['%7f','%7f','%7f','%7d'])
        if choice1==0:
                if choice==0:
                        np.savetxt('/uufs/astro.utah.edu/common/uuastro/astro_data/zhengzheng/mccarthy/GravBinary_Data/random_{0}_jack_suCALC_{1}tree_{2}.dat'.format(TYPE,Bin,TYPE1),coords,fmt=['%7f','%7f','%7f','%7d'])
                if choice==1:
                        np.savetxt('/uufs/astro.utah.edu/common/uuastro/astro_data/zhengzheng/mccarthy/GravBinary_Data/random_{0}_jack_rppiCALC_{1}tree_{2}.dat'.format(TYPE,Bin,TYPE1),coords,fmt=['%7f','%7f','%7f','%7d'])

print('Finished cell identification and sorting at {}.\nPrinted file.\n---------------------------'.format(datetime.now().time()))

print('Started cell indexing at {}.'.format(datetime.now().time()))
#Ngalc=len(np.unique(IND))
galIND=np.zeros((Nc,2))
i=0
k=0
while(True):
	galIND[(IND[k]-1),0]=k+1
	l=len(IND[IND==IND[k]])-1
	galIND[(IND[k]-1),1]=k+1+l
	s=1
	while((k+s)!=N and IND[k]==IND[k+s]):
		s=s+1
	k=k+s
	i=i+1
	if (k==N): break

galIND=galIND.astype(int)

if choice1==1:
	if choice==0:
		np.savetxt('{0}_suGALind_{1}tree_{2}.dat'.format(TYPE,Bin,TYPE1),galIND, fmt='%7d')
	if choice==1:
       		np.savetxt('{0}_rppiGALind_{1}tree_{2}.dat'.format(TYPE,Bin,TYPE1),galIND, fmt='%7d')
if choice1==0:
      	if choice==0:
    		np.savetxt('random_{0}_suGALind_{1}tree_{2}.dat'.format(TYPE,Bin,TYPE1),galIND, fmt='%7d')
        if choice==1:
                np.savetxt('random_{0}_rppiGALind_{1}tree_{2}.dat'.format(TYPE,Bin,TYPE1),galIND, fmt='%7d')

print('Finished cell indexing at {}.\nPrinted file.\n---------------------------'.format(datetime.now().time()))
	
print('Program finished at {}.\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'.format(datetime.now().time()))


