import numpy as np
import matplotlib.pyplot as plt


####### From https://gist.github.com/eteq/4599814
try:
	from scipy.spatial import cKDTree as KDT
except ImportError:
	from scipy.spatial import KDTree as KDT

def spherematch(ra1, dec1, ra2, dec2, tol=None, nnearest=1):
    """
    Finds matches in one catalog to another.
    Parameters
    ra1 : array-like
        Right Ascension in degrees of the first catalog
    dec1 : array-like
        Declination in degrees of the first catalog (shape of array must match `ra1`)
    ra2 : array-like
        Right Ascension in degrees of the second catalog
    dec2 : array-like
        Declination in degrees of the second catalog (shape of array must match `ra2`)
    tol : float or None, optional
        How close (in degrees) a match has to be to count as a match.  If None,
        all nearest neighbors for the first catalog will be returned.
    nnearest : int, optional
        The nth neighbor to find.  E.g., 1 for the nearest nearby, 2 for the
        second nearest neighbor, etc.  Particularly useful if you want to get
        the nearest *non-self* neighbor of a catalog.  To do this, use:
        ``spherematch(ra, dec, ra, dec, nnearest=2)``
    Returns
    -------
    idx1 : int array
        Indecies into the first catalog of the matches. Will never be
        larger than `ra1`/`dec1`.
    idx2 : int array
        Indecies into the second catalog of the matches. Will never be
        larger than `ra1`/`dec1`.
    ds : float array
        Distance (in degrees) between the matches
    """

    ra1 = np.array(ra1, copy=False)
    dec1 = np.array(dec1, copy=False)
    ra2 = np.array(ra2, copy=False)
    dec2 = np.array(dec2, copy=False)

    if ra1.shape != dec1.shape:
        raise ValueError('ra1 and dec1 do not match!')
    if ra2.shape != dec2.shape:
        raise ValueError('ra2 and dec2 do not match!')

    x1, y1, z1 = _spherical_to_cartesian(ra1.ravel(), dec1.ravel())

    # this is equivalent to, but faster than just doing np.array([x1, y1, z1])
    coords1 = np.empty((x1.size, 3))
    coords1[:, 0] = x1
    coords1[:, 1] = y1
    coords1[:, 2] = z1

    x2, y2, z2 = _spherical_to_cartesian(ra2.ravel(), dec2.ravel())

    # this is equivalent to, but faster than just doing np.array([x1, y1, z1])
    coords2 = np.empty((x2.size, 3))
    coords2[:, 0] = x2
    coords2[:, 1] = y2
    coords2[:, 2] = z2

    kdt = KDT(coords2)
    if nnearest == 1:
        idxs2 = kdt.query(coords1)[1]
    elif nnearest > 1:
        idxs2 = kdt.query(coords1, nnearest)[1][:, -1]
    else:
        raise ValueError('invalid nnearest ' + str(nnearest))

    ds = _great_circle_distance(ra1, dec1, ra2[idxs2], dec2[idxs2])

    idxs1 = np.arange(ra1.size)

    if tol is not None:
        msk = ds < tol
        idxs1 = idxs1[msk]
        idxs2 = idxs2[msk]
        ds = ds[msk]

    return idxs1, idxs2, ds


def _spherical_to_cartesian(ra, dec):
    """
    (Private internal function)
    Inputs in degrees.  Outputs x,y,z
    """
    rar = np.radians(ra)
    decr = np.radians(dec)

    x = np.cos(rar) * np.cos(decr)
    y = np.sin(rar) * np.cos(decr)
    z = np.sin(decr)

    return x, y, z


def _great_circle_distance(ra1, dec1, ra2, dec2):
    """
    (Private internal function)
    Returns great circle distance.  Inputs in degrees.
    Uses vicenty distance formula - a bit slower than others, but
    numerically stable.
    """
    from numpy import radians, degrees, sin, cos, arctan2, hypot

    # terminology from the Vicenty formula - lambda and phi and
    # "standpoint" and "forepoint"
    lambs = radians(ra1)
    phis = radians(dec1)
    lambf = radians(ra2)
    phif = radians(dec2)

    dlamb = lambf - lambs

    numera = cos(phif) * sin(dlamb)
    numerb = cos(phis) * sin(phif) - sin(phis) * cos(phif) * cos(dlamb)
    numer = hypot(numera, numerb)
    denom = sin(phis) * sin(phif) + cos(phis) * cos(phif) * cos(dlamb)
    return degrees(arctan2(numer, denom))
########################################

################# From Z.Z.
def binarySearch(alist, item):
	first = 0
	last = len(alist)-1
	found = False
	while first<=last and not found:
	        midpoint = (first + last)//2
	        if alist[midpoint] == item:
	            found = True
	        else:
	            if item < alist[midpoint]:
	                last = midpoint-1
	            else:
	                first = midpoint+1
	if found:
		return found,midpoint
	else:
		return found, 00
##########################################


TYPE1='vespaBC_1_gals_SFH6HSsig_dr7_sp2'
TYPE2='LINksmhigh22_sSFRzcut_parent'

# Randoms from Hong, early/late created by Kevin
raR,decR,jackR=np.loadtxt('random.dat',unpack=True)
data_early=np.loadtxt('{}.dat'.format(TYPE1),unpack=False)
data_late=np.loadtxt('{}.dat'.format(TYPE2),unpack=False)

#TYPE1='LINearly_sSFR'
#TYPE2='LINlate_sSFR'

raGe=data_early[:,0]
decGe=data_early[:,1]
zGe=data_early[:,2]

raGl=data_late[:,0]
decGl=data_late[:,1]
zGl=data_late[:,2]

# Check footprint 
plt.scatter(raR,decR,color='y')
plt.scatter(raGl,decGl,color='b',s=1,alpha=0.5)
plt.show()

plt.scatter(raR,decR,color='y')
plt.scatter(raGe,decGe,color='r',s=1,alpha=0.5)
plt.show()

# Assign jack id to early/late galaxies
indGe,indRe,ds=spherematch(raGe,decGe,raR,decR,tol=1,nnearest=1)
indGl,indRl,ds=spherematch(raGl,decGl,raR,decR,tol=1,nnearest=1)

jackGe=jackR[indRe]
jackGl=jackR[indRl]

# Assign redshift to random catalogs
Nr=len(raR)

znum=100
histe,edgee=np.histogram(zGe, znum)

randmulte=Nr/len(zGe)
hist1e=np.around(histe*randmulte)
binwidthe=edgee[1]-edgee[0]

histl,edgel=np.histogram(zGl, znum)

randmultl=Nr/len(zGl)
hist1l=np.around(histl*randmultl)
binwidthl=edgel[1]-edgel[0]


zrande=[]
for i in range(0,znum):
        if i <znum:
                Ne=int(hist1e[i])
                for j in range(0,Ne):
                        if j < Ne:
                                ztmpe=np.random.rand()*binwidthe+edgee[i]
                                zrande.append(ztmpe)

se=len(zrande)
raRe=raR[0:se]
decRe=decR[0:se]
jackRe=jackR[0:se]
catRe=np.vstack((raRe,decRe,zrande,jackRe)).T
#np.savetxt('/uufs/astro.utah.edu/common/uuastro/astro_data/zhengzheng/mccarthy/GravBinary_Data/random_{}_jack.dat'.format(TYPE1),catRe,fmt=['%5e','%5e','%5e','%5d'])


zrandl=[]
for i in range(0,znum):
	if i <znum:
		Nl=int(hist1l[i])
		for j in range(0,Nl):
			if j < Nl:
				ztmpl=np.random.rand()*binwidthl+edgel[i]
				zrandl.append(ztmpl)

sl=len(zrandl)
raRl=raR[0:sl]
decRl=decR[0:sl]
jackRl=jackR[0:sl]
catRl=np.vstack((raRl,decRl,zrandl,jackRl)).T
#np.savetxt('/uufs/astro.utah.edu/common/uuastro/astro_data/zhengzheng/mccarthy/GravBinary_Data/random_{}_jack.dat'.format(TYPE2),catRl,fmt=['%5e','%5e','%5e','%5d'])

plt.hist(zrandl,bins=100,normed=True,label='randoms')
plt.hist(zGl,bins=100,normed=True,alpha=0.5,label='late')
plt.legend(loc='upper right')
plt.show()

plt.hist(zrande,bins=100,normed=True,label='randoms')
plt.hist(zGe,bins=100,normed=True,alpha=0.5,label='early')
plt.legend(loc='upper right')
plt.show()



earlyjack=np.vstack((raGe,decGe,zGe,jackGe)).T
latejack=np.vstack((raGl,decGl,zGl,jackGl)).T
#np.savetxt('{}_jack.dat'.format(TYPE1),earlyjack)
#np.savetxt('{}_jack.dat'.format(TYPE2),latejack)
