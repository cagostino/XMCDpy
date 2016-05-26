#lets gooooo. Author: Christopher Agostino
import numpy as np
import os
import osgeo
import scipy
import gdal
import matplotlib.pyplot as plt
from lmfit.models import LorentzianModel

#import pyfits
import gc
plt.rc('text',usetex=True)
#make these plots look dank
plt.rc('font',family='serif')
l = os.listdir("left/")
r =os.listdir("right/")
l.sort()
r.sort()
def read_dat(x,y,px, area,l):
	fil = np.loadtxt(str(x)+'_'+str(y)+'/'+str(x)+'_'+str(y)+'_'+str(area)+'_area_'+str(px)+'dir_'+l+'.csv', skiprows=1)
	return fil.transpose()
plt.ion()
