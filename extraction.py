#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
 Extraction of ALOS AVNIR-2 Image with orthorectification
	Original data : CEOS format
		Specified Region from DEM Geotif Image region 

        2016.7.6 by Y. Iikura	
"""
# cd /Volumes/Transcend/LandsatETM+
# extraction.py ALAV2A017062800_20060521 dem.tif NEW

import sys
import os
import numpy as np
import cv2
import convert_util as ut
import proj_util as pr

'''
os.chdir('..')
reload(ut)
'''

param=sys.argv
if len(param)!=4 and len(param)!=6:
    print 'Usage: extraction.py scene_name dem.tif new_folder (dx=0.0 dy=0.0)'
    exit()

fscene=param[1]
fname=param[2]
fnew=param[3]
if len(param)==6:
  dx=float(param[4])
  dy=float(param[5])
else:
  dx=0.0
  dy=0.0


if os.path.isdir(fnew) == False: os.mkdir(fnew)  

'''
fscene='ALAV2A017062800_20060521'
fname='dem.tif'
fnew='NEW'
dx=0.0
dy=0.0
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

ut.gwrite(sat,fnew)

exit()

for band in [1,2,3,4]:
  sat.read_band(band)
  conv=ut.convert(sat,xmax,ymax)
  new=conv.convert(dx,dy)
  pr.write_tif('../'+fnew+'/band'+str(band)+'.tif',new,1)

exit()
