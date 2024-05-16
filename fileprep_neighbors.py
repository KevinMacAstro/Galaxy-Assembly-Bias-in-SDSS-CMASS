import numpy as np

choice=int(input('For xi(s,u) or xi(rp,pi) (0 or 1) : '))

TYPE1='LINearly_sSFR_parent'


from datetime import datetime
print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\nProgram started at {}.\n---------------------------'.format(datetime.now().time()))

print('Load data and find bounds started at {}.\n'.format(datetime.now().time()))
coords1=np.loadtxt('/uufs/astro.utah.edu/common/uuastro/astro_data/zhengzheng/mccarthy/ObsAB_Data/random_{}_CART.dat'.format(TYPE1),unpack=False)
coords3=np.loadtxt('{}_CART.dat'.format(TYPE1),unpack=False)


Max=np.amax([np.amax(coords1),np.amax(coords3)])
Min=np.amin([np.amin(coords1),np.amin(coords3)])
MAX=np.amax([abs(Max),abs(Min)])
print('Load data and find bounds finished at {}.\n---------------------------'.format(datetime.now().time()))

print('Started bin and bin neighbor identification at {}.\n'.format(datetime.now().time()))
if choice==0:
# For su calculations
        lgSmax=1.4 
        S=np.power(10,lgSmax)
if choice==1:
# For rppi calculation
	pi_max=80
        S=pi_max

###Layers of neighboring cells###
##eta==1 has one layer of neighbors##
eta=1
Nbins=2*np.ceil(np.ceil(MAX*2/(eta*S))/2).astype(int)
rangeMin=-Nbins/2.*eta*S
rangeMax=Nbins/2.*eta*S

bins=np.linspace(rangeMin,rangeMax,(Nbins+1))
Lhalf=abs(bins[1]-bins[0])/2.

Nc=(Nbins)**3
binNUM=np.zeros((Nc,6))
p=0
for i in range(1,(Nbins+1)):
	for j in range(1,(Nbins+1)):
		for k in range(1,(Nbins+1)):
			binNUM[p,0]=i
			binNUM[p,1]=j
			binNUM[p,2]=k
			binNUM[p,3]=bins[(i-1)]+Lhalf
			binNUM[p,4]=bins[(j-1)]+Lhalf
			binNUM[p,5]=bins[(k-1)]+Lhalf
			p=p+1
binNum=binNUM[:,0:3].astype(int)
binNum=binNum.astype(int)
if choice==0:
        np.savetxt('suCELLnum_{}.dat'.format(TYPE),binNum,fmt='%7d')
if choice==1:
        np.savetxt('rppiCELLnum_{}.dat'.format(TYPE),binNum,fmt='%7d')
print('Finished bin identification at {}.\n'.format(datetime.now().time()))
j=2
i=0
while(True):
	while((binNum[i,0]!=j or binNum[i,1]!=j or binNum[i,2]!=j)):
		i=i+1
	if np.sqrt((binNUM[0,3]-binNUM[i,3])**2+(binNUM[0,4]-binNUM[i,4])**2+(binNUM[0,5]-binNUM[i,5])**2)>S:
		l=np.sqrt((binNUM[0,3]-binNUM[i,3])**2+(binNUM[0,4]-binNUM[i,4])**2+(binNUM[0,5]-binNUM[i,5])**2)
		break
	j=j+1

binIND=np.linspace(1,Nc,Nc).astype(int)
binTMP=binNUM[:,4:6][binNUM[:,3]>0]
indTMP=binIND[binNUM[:,3]>0]
binTMP2=binTMP[:,1][binTMP[:,0]>0]
indTMP2=indTMP[binTMP[:,0]>0]
indTMP3=indTMP2[binTMP2>0]
ind=indTMP3[0]-1

Neighbors=len(binNUM[:,0][np.sqrt((binNUM[:,3]-binNUM[ind,3])**2+(binNUM[:,4]-binNUM[ind,4])**2+(binNUM[:,5]-binNUM[ind,5])**2)<=l])

binNEI=np.zeros((Nc,Neighbors))
#perm=np.array([[-1,-1,-1], [-1,-1,0], [-1,-1,1], [-1,0,-1], [-1,0,0], [-1,0,1], [-1,1,-1], [-1,1,0], [-1,1,1], [0,-1,-1], [0,-1,0], [0,-1,1], [0,0,-1], [0,0,1], [0,1,-1], [0,1,0], [0,1,1], [1,-1,-1], [1,-1,0], [1,-1,1], [1,0,-1], [1,0,0], [1,0,1], [1,1,-1] ,[1,1,0], [1,1,1]])
#f=0

for i in range(0,Nc):
        if i < Nc:
		print('{}'.format(i))
		binNUMtmp=binNUM[:,0:3][np.sqrt((binNUM[:,3]-binNUM[i,3])**2+(binNUM[:,4]-binNUM[i,4])**2+(binNUM[:,5]-binNUM[i,5])**2)<=l]
		binINDtmp=binIND[np.sqrt((binNUM[:,3]-binNUM[i,3])**2+(binNUM[:,4]-binNUM[i,4])**2+(binNUM[:,5]-binNUM[i,5])**2)<=l]
		s=len(binINDtmp)
		binNEI[i,:s]=binINDtmp
print('Finished bin neighbor identification at {}.\n---------------------------'.format(datetime.now().time()))

binNEI=binNEI.astype(int)
if choice==0:
        np.savetxt('suCELLnei_{}.dat'.format(TYPE),binNEI,fmt='%7d')
if choice==1:
        np.savetxt('rppiCELLnei_{}.dat'.format(TYPE),binNEI,fmt='%7d')

print('Program finished at {}.\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'.format(datetime.now().time()))
