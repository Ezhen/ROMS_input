
import numpy as np
from numpy import shape
import pyroms
import pyroms_toolbox
from datetime import datetime
import subprocess
import os
import commands
import _iso

from nco import Nco
import _remapping

def flood(varz, Bgrd, Bpos='t', irange=None, jrange=None, \
          spval=-32768, dmax=0, cdepth=0, kk=0):
    #print varz
    #print type(varz)
    varz = varz.copy()
    #print type(varz)
    assert len(varz.shape) == 3, 'var must be 3D'

    # replace spval by nan
    idx = np.where(varz.mask==True)
    #print np.where(varz==-32767)
    varz[idx] = np.nan
    #print varz

    x = Bgrd.lon_t
    y = Bgrd.lat_t
    h = Bgrd.h
    if Bpos is 't':
        mask = Bgrd.mask_t[0,:,:]
    elif Bpos is 'uv':
        mask = Bgrd.mask_uv[0,:,:]

    nlev, Mm, Lm = varz.shape
    #print varz.shape

    if irange is None:
        irange = (0,Lm)
    else:
        assert varz.shape[2] == irange[1]-irange[0], \
               'var shape and irange must agreed'

    if jrange is None:
        jrange = (0,Mm)
    else:
        assert varz.shape[1] == jrange[1]-jrange[0], \
               'var shape and jrange must agreed'

    x = x[jrange[0]:jrange[1], irange[0]:irange[1]]
    y = y[jrange[0]:jrange[1], irange[0]:irange[1]]
    h = h[jrange[0]:jrange[1], irange[0]:irange[1]]
    mask = mask[jrange[0]:jrange[1], irange[0]:irange[1]]
    #print np.shape(x), np.shape(y), np.shape(h), np.shape(mask)

    # Finding nearest values in horizontal
    # critical deph => no change if depth is less than specified value
    cdepth = abs(cdepth)
    if cdepth != 0:
        idx = np.where(h >= cdepth)
        msk = np.zeros(mask.shape)
        msk[idx] = 1
    else:
        msk = mask.copy()
    for k in range(nlev-1,0,-1):
        c1 = np.array(msk, dtype=bool)
	#print 	c1
        c2 = np.isnan(varz[k,:,:]) == 1
	#print c2
        if kk == 0:
            c3 = np.ones(mask.shape).astype(bool)
        else:
            c3 = np.isnan(varz[min(k-kk,0),:,:]) == 0
	#print np.shape(c1)
	#print np.shape(c2)
	#print np.shape(c3)
        c = c1 & c2 & c3
        idxnan = np.where(c == True)
        idx = np.where(c2 == False)
        if list(idx[0]):
            wet = np.zeros((len(idx[0]),2))
            dry = np.zeros((len(idxnan[0]),2))
            wet[:,0] = idx[0]+1
            wet[:,1] = idx[1]+1
            dry[:,0] = idxnan[0]+1
            dry[:,1] = idxnan[1]+1
	    #print np.round(varz[:,30,30],2)
            varz[k,:] = _remapping.flood(varz[k,:], wet, dry, x, y, dmax)

    # drop the deepest values down
    idx = np.where(np.isnan(varz) == 1)
    varz[idx] = spval
    bottom = pyroms.utility.get_bottom(varz[::-1,:,:], mask, spval=spval)
    bottom = (nlev-1) - bottom
    for i in range(Lm):
        for j in range(Mm):
            if mask[j,i] == 1:
                varz[bottom[j,i]:,j,i] = varz[bottom[j,i],j,i]
    return varz
