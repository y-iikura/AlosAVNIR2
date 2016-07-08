#!/usr/bin/env python
# -*- coding: utf-8 -*-
# cd /Volumes/Transcend/JAXA土地被覆

import numpy as np
from osgeo import gdal
from osgeo import osr

a=6378137.0
b=6356752.314
ee=(a*a-b*b)/a/a
aee=a*(1.0-ee)
eee=ee/(1.0-ee)
zone=54
a0=1.005052501813147
a1=0.005063108622327
a2=0.000010627590328
a3=0.000000020820407
a4=0.000000000039332
a5=0.000000000000073
n=(1.0-np.sqrt(1-ee))/(1.0+np.sqrt(1-ee))
b0=(1+n)/(1+n**2/4+n**4/4)
b1=3*n/2-27*n**3/32+269*n**5/512
b2=21*n**2/16-55*n**4/32
b3=151*n**3/96-417*n**5/128
b4=1097*n**4
b5=8011*n**5/2560
dtor=np.pi/180.0

lat0=0.0 ; lon0=0.0
dlat=0.0 ; dlon=0.0

xs=0.0 ; xe=0.0
ye=0.0 ; ys=0.0
dx=0.0 ; dy=0.0

# for utm zone54
wkt='PROJCS["WGS 84 / UTM zone 54N",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.0174532925199433,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",141],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",0],UNIT["metre",1,AUTHORITY["EPSG","9001"]],AXIS["Easting",EAST],AXIS["Northing",NORTH],AUTHORITY["EPSG","32654"]]'

# for lat-lon
wkt2='GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0],UNIT["degree",0.0174532925199433],AUTHORITY["EPSG","4326"]]'

def slat(lat):
  rlat=lat*dtor
  s2=np.sin(2*rlat)
  s4=np.sin(4*rlat)
  s6=np.sin(6*rlat)
  s8=np.sin(8*rlat)
  s10=np.sin(10*rlat)
  temp=a0*rlat-a1*s2/2+a2*s4/4-a3*s6/6+a4*s8/8-a5*s10/10
  return aee*temp

def tlat(y):
  rlat=b0*y/a
  s2=np.sin(2*rlat)
  s4=np.sin(4*rlat)
  s6=np.sin(6*rlat)
  s8=np.sin(8*rlat)
  s10=np.sin(10*rlat)
  temp=rlat+b1*s2+b2*s4+b3*s6+b4*s8+b5*s10
  return temp

def bl2utm(ll):
  lat=ll[0]
  lon=ll[1]
  lon=lon-((zone-30.0)*6-3.0)
  y0=slat(lat)
  lat=lat*dtor
  lon=lon*dtor
  #print,lat,y0
  cf=np.cos(lat)
  sf=np.sin(lat)
  tf=sf/cf
  tt=tf*tf
  ttt=tt*tt
  uu=eee*cf*cf
  uuu=uu*uu
  N=a/np.sqrt(1-ee*sf*sf)
  x1=lon*cf
  x2=x1*x1
  x=x1*(1.0+x2*((1-tt+uu)+x2*(5-18*tt+ttt+14*uu-58*tt*uu)/20)/6)
  y=y0+N*x2*tf*(1+x2*((5-tt+9*uu+4*uuu)+x2*(61-58*tt+ttt+270*uu-330*tt*uu)/30)/12)/2
  x=0.9996*N*x+500000
  y=0.9996*y
  return [x,y]

def utm2bl(xy):
  x0=(xy[0]-500000.0)/0.9996
  y0=xy[1]/0.9996
  f0=tlat(y0)
  cf=np.cos(f0)
  sf=np.sin(f0)
  tf=np.tan(f0)
  tt=tf*tf
  ttt=tt*tt
  uu=eee*cf*cf
  N=a/np.sqrt(1-ee*sf*sf)
  RR=aee/(1-ee*sf*sf)**1.5
  r1=x0/N
  r2=r1*r1
  r=r1*(1.0-r2*((1.0+2*tt+uu)-r2*(5+28*tt+24*ttt)/20)/6)
  #f=r2*(1.0-r2*((5.0+3*tt+uu-9*tt*uu-4*uu**2)-r2*(61+90*tt+45*ttt)/30)/12)/2
  f=r2*(1.0-r2*((5.0+3*tt+uu-9*tt*uu)-r2*(61+90*tt+45*ttt)/30)/12)/2
  lon=((zone-30)*6-3.0)+r/cf/dtor	
  lat=(f0-N*tf*f/RR)/dtor
  return[lat,lon]

def read_tif(fname):
  src = gdal.Open(fname, gdal.GA_Update)
  pdem = src.GetRasterBand(1)
  gt = src.GetGeoTransform()
  image = pdem.ReadAsArray()
  proj = osr.SpatialReference()
  proj.ImportFromWkt(src.GetProjectionRef())
  wkt=proj.ExportToWkt()
  return [gt,wkt,image]

def write_tif(dname,data,select):
  driver = gdal.GetDriverByName('GTiff')
  #wkt_projection=proj.ExportToWkt()
  y_pixels,x_pixels=data.shape
  dataset = driver.Create(
    dname,
    x_pixels,
    y_pixels,
    1,
    gdal.GDT_Byte, )
  dataset.SetGeoTransform((
    xs,
    dx,
    0, 
    ye,
    0,
    -dy))  
  #dataset.SetProjection(wkt_projection)
  if select==1: dataset.SetProjection(wkt)
  else: dataset.SetProjection(wkt2)
  dataset.GetRasterBand(1).WriteArray(data)
  dataset.FlushCache()

if __name__=='__main__':
    import sys
    param=sys.argv
    if len(param) !=3:
      print 'Usaage: jaxa_util.py latitude longitude'
      exit()
    #print param[1:]
    lat=float(sys.argv[1])
    lon=float(sys.argv[2])
    x,y=bl2utm([lat,lon])
    #print('{0:6.2f} {1:6.2f}'.format(x,y))
    print('x={0:6.2f},  y={1:6.2f}'.format(x,y))
    lat2,lon2=utm2bl([x,y])
    print('lat2={0:9.7f},  lon2={1:9.7f}'.format(lat2,lon2))

    

