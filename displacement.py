#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
 Check of displacement between 
           DEM and AVNIR-2 with orthorectification
   
    2016.7.7 by Y.Iikura
"""

# cd /Volumes/Transcend/AlosAVNIR2
# displacement.py ALAV2A017062800_20060521 dem.tif

import sys
import os
import numpy as np
import cv2
import convert_util as ut
import proj_util as pr

'''
os.chdir('..')
reload(ut)
os.chdir(fscene)
'''

param=sys.argv
if len(param)!=4:
    print 'Usage: displacement.py scene_name dem.tif '

fscene=param[1]
fname=param[2]

'''
fscene='ALAV2A017062800_20060521'
fname='dem.tif'
'''

gt,wkt,dem=pr.read_tif(fname)
ymax,xmax=dem.shape
pr.xs=gt[0]
pr.dx=gt[1]
pr.xe=gt[0]+xmax*gt[1]
pr.ye=gt[3]
pr.dy=-gt[5]
pr.ys=gt[3]-ymax*gt[5]

print gt
print dem.shape,dem.dtype
print xmax,ymax


#demx=cv2.resize(dem,(1000,200))
#cv2.imshow('dem',demx/np.max(demx))



os.chdir(fscene)
sname=fscene.split('_')[0]

flag=1
sat=ut.original(sname)
if flag:
  print sat.jmax,sat.imax
  print sat.xs,sat.dx
  print sat.ye,sat.dy
  print sat.offset
  print sat.gain

sat.chokka()
if flag:
  print sat.pv
  print sat.qv
  print sat.xn,sat.yn

dem[dem<0.0]=0.0
inc=ut.incident(dem,sat)
#incx=cv2.resize(inc,(1200,300))
#cv2.imshow('incx',incx/np.max(incx))
#cv2.destroyWindow('incx')


#sat.display('sat',600,400)
#cv2.destroyWindow('sat')

sat.read_band(4)
conv=ut.convert(sat,xmax,ymax)

ss=ut.s_file(pr.xs,pr.ye,ymax,xmax,sat)
difx=dem*np.tan(ss*np.pi/180.0)*sat.pv[0]
dify=dem*np.tan(ss*np.pi/180.0)*sat.pv[1]
difx=difx.astype(np.float32)
dify=dify.astype(np.float32)

'''
new=conv.convert(0.0,0.0)
newx=cv2.resize(new,(1200,300))
cv2.imshow('newx',newx)
cv2.destroyWindow('newx')
np.corrcoef(inc.flatten(),new.flatten())


new=conv.convert2(0.0,0.0,0.0,0.0)
new32=new.astype(np.float32)
dd,dt=cv2.phaseCorrelate(inc,new32); print dd,dt
new=conv.convert2(-pr.dx*dd[0],pr.dy*dd[1],0.0,0.0)
new32=new.astype(np.float32)
dd,dt=cv2.phaseCorrelate(inc,new32); print dd,dt

new=conv.convert2(0.0,0.0,-difx,dify)
new32=new.astype(np.float32)
dd,dt=cv2.phaseCorrelate(inc,new32); print dd,dt
new=conv.convert2(-pr.dx*dd[0],pr.dy*dd[1],-difx,dify)
new32=new.astype(np.float32)
dd,dt=cv2.phaseCorrelate(inc,new32); print dd,dt
'''


print("*** Phase Only Correlation ***")
for band in [1,2,3,4]:
  print ' Band '+str(band)+':'
  sat.read_band(band)
  new=conv.convert2(0.0,0.0,0.0,0.0)
  new32=new.astype(np.float32)
  dd,dt=cv2.phaseCorrelate(inc,new32); print dd,dt
  new=conv.convert2(-pr.dx*dd[0],pr.dy*dd[1],0.0,0.0)
  new32=new.astype(np.float32)
  dd,dt=cv2.phaseCorrelate(inc,new32); print dd,dt
  new=conv.convert2(0.0,0.0,-difx,dify)
  new32=new.astype(np.float32)
  dd,dt=cv2.phaseCorrelate(inc,new32); print dd,dt
  new=conv.convert2(-pr.dx*dd[0],pr.dy*dd[1],-difx,dify)
  new32=new.astype(np.float32)
  dd,dt=cv2.phaseCorrelate(inc,new32); print dd,dt
  #pr.write_tif('../'+fnew+'/band'+str(band)+'.tif',new,1)

#ut.gwrite(sat,fnew)
exit()
