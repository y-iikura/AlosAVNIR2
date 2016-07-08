#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
 Satellie Image Class for Geometric Correction
	Original data : AVNIR-2 : CEOS format
	New data : converted image with specified region 	
 7/5/2016 copied from ETM+
"""
import numpy as np
from numpy.linalg import inv
import cv2
#from osgeo import gdal
#from osgeo import osr
import proj_util as pr
import array

class original:
  def __init__(self,fname):
    self.name=fname
    f=open('LED-'+fname+'-O1B2G_U','r')
    disc=f.read(4680)
    header=f.read(4680)
    proj1=f.read(1916)
    data=f.read(48)
    proj2=f.read(2716)
    radio=f.read(4680)
    f.close()
    #gain = np.zeros(11,dtype='float32')
    #offset = np.zeros(11,dtype='float32')
    gain=[]
    offset=[]
    lat=float(str(header[212:228])) # scene center
    lon=float(str(header[228:244]))
    line=float(str(header[244:260]))
    pixl=float(str(header[260:276]))
    jmax=int(2*line-1)
    imax=int(2*pixl-1)
    pangle=float(str(header[372:388]))
    self.sun_el=float(str(header[458:461]))
    self.sun_az=float(str(header[463:466]))
    str(header[1732:1748]) # upper left
    str(header[1748:1763])
    #dx0=str(proj1[44:60])
    #dy0=str(proj1[60:76])
    dx0=10.0
    dy0=10.0
    y0=float(str(proj1[140:156]))*1000.0
    x0=float(str(proj1[156:172]))*1000.0
    t0=float(str(proj1[204:220]))*180.0/np.pi
    tt=float(str(proj1[700:716]))*180.0/np.pi
    str(radio[2694:2766])
    gain=[]
    offset=[]
    for i in range(4):
      gain.append(float(str(radio[(2702+16*i):(2710+16*i)])))
      offset.append(float(str(radio[(2710+16*i):(2718+16*i)])))
    print imax,jmax
    self.gain=gain
    self.offset=offset
    #self.image=image
    self.imax=imax
    self.jmax=jmax
    self.dx=dx0
    self.dy=dy0
    self.pc=self.imax/2.0
    self.lc=self.jmax/2.0
    self.xc=x0
    self.yc=y0
    self.xs=x0-self.pc*dx0
    self.ye=y0+self.lc*dx0
    self.t0=t0
    self.pangle=pangle
    self.hsat=691650.00
    self.earth=6378137.0
    self.read_band(1)
  def read_band(self,band):
    sname2='IMG-{:02d}'.format(band)+'-'+self.name+'-O1B2G_U'
    temp=array.array('B')
    f=open(sname2)
    nbytes=(self.imax+100)*(self.jmax+1)
    temp.fromfile(f,nbytes)
    f.close()
    temp2=np.array(temp,dtype=np.uint8).reshape(self.jmax+1,self.imax+100)
    self.image=temp2[1:,34:(self.imax+34)]
  def display(self,name,xmax,ymax):
    old=self.image
    min=np.min(old) ; max=np.max(old)
    oldx=np.uint8(255.0*(old-min)/(max-min))
    oldy=cv2.resize(oldx,(xmax,ymax))
    cv2.imshow(name,oldy)
  def xy2pl(self,x,y):
    aa=np.cos(self.t0)/self.dx
    bb=-np.sin(self.t0)/self.dy
    p=aa*(x-self.xs)+bb*(y-self.ye)
    l=bb*(x-self.xs)-aa*(y-self.ye)
    return [p,l]
  def pl2xy(self,p,l):
    aa=np.cos(self.t0)*self.dx
    bb=-np.sin(self.t0)*self.dy
    x=aa*p+bb*l+self.xs
    y=bb*p-aa*l+self.ye
    return [x,y]
  def coverage(self):
    ul=self.pl2xy(0.0,0.0) 
    ur=self.pl2xy(float(self.imax),0.0)
    ll=self.pl2xy(0.0,float(self.jmax))
    lr=self.pl2xy(float(self.imax),float(self.jmax)) 
    left=int(ll[0]/10000.0)
    right=int(np.ceil(ur[0]/10000.0)) 
    lower=int(lr[1]/10000.0) 
    upper=int(np.ceil(ul[1]/10000.0))
    return [left,right,lower,upper]
  def chokka(self):
    mask=np.zeros((self.jmax,self.imax),dtype=np.uint8)
    temp=np.where(self.image > 0) 
    mask[temp]=1
    xtemp = temp[1]
    ytemp = temp[0]
    set1=np.min(ytemp)
    set2=np.max(xtemp)
    set3=np.max(ytemp)
    set4=np.min(xtemp)
    state=0
    position=[0]
    for i in range(self.imax):
      if mask[set1,i] != state:
        position.append(i)
        state=1-state
    state=0
    for j in range(self.jmax):
      if mask[j,set2] != state:
        position.append(j)
        state=1-state
    state=0
    for i in range(self.imax):
      if mask[set3,i] != state:
        position.append(i)
        state=1-state
    state=0
    for j in range(self.jmax):
      if mask[j,set4] != state:
        position.append(j)
        state=1-state
    tq1=float(position[1]-set4)/(position[7]-set1)
    tp1=float(position[3]-set1)/(set2-position[2])
    tq2=float(set2-position[6])/(set3-position[4])
    tp2=float(set3-position[8])/(position[5]-set4)
    tp=np.arctan((tp1+tp2)/2.0)
    tq=np.arctan((tq1+tq2)/2.0)
    self.pv=[np.cos(tp),-np.sin(tp)]
    self.qv=[-np.sin(tq),-np.cos(tq)]
    #pq=transpose([[pv],[qv]])
    #pinv=invert(pq)
    angle=self.pangle*np.pi/180.0
    angle2=np.arcsin((self.earth+self.hsat)*np.sin(angle)/self.earth)
    dn=self.earth*(angle2-angle)
    self.xn=self.xc-dn*self.pv[0]
    self.yn=self.yc-dn*self.pv[1]
  #6/14/2016 copied from apm_util.py
  def xdist(self,xxx,yyy,pqinv):
    xy=np.array([[xxx-self.xn,yyy-self.yn]])
    #xy=np.array([xxx-self.xc,yyy-self.yc])
    temp=pqinv.dot(xy)
    return temp[0]
  def angle(self,px):
    a=px/self.earth
    s=(self.earth+self.hsat)*np.sin(a)
    c=(self.earth+self.hsat)*np.cos(a)
    f=np.arctan(s/(c-self.earth))
    return -np.degrees(f)
  def sangle(self,xxx,yyy):
    pq=np.array([self.pv,self.qv])
    pqt = pq.T
    pqinv = inv(pqt)
    pxy=self.xdist(xxx,yyy,pqinv)
    return self.angle(pxy)

class convert:
  def __init__(self,old,xmax,ymax):
    self.old=old
    self.xs=pr.xs
    self.ye=pr.ye
    self.dx=pr.dx
    self.dy=pr.dy
    self.xmax=xmax
    self.ymax=ymax
  def convert(self,tx,ty):
    x=self.xs+np.arange(self.xmax)*self.dx
    y=self.ye-np.arange(self.ymax)*self.dy
    xso=self.old.xs+tx
    yeo=self.old.ye+ty
    xc=(x-xso)/self.old.dx
    yc=(yeo-y)/self.old.dy
    xx,yy=np.meshgrid(xc,yc)
    xx=np.float32(xx)
    yy=np.float32(yy)
    return cv2.remap(self.old.image,xx,yy,cv2.INTER_LINEAR)
  def convert2(self,tx,ty,difx,dify):
    x=self.xs+np.arange(self.xmax)*self.dx
    y=self.ye-np.arange(self.ymax)*self.dy
    xso=self.old.xs+tx
    yeo=self.old.ye+ty
    xc=(x-xso)/self.old.dx
    yc=(yeo-y)/self.old.dy
    xx,yy=np.meshgrid(xc,yc)
    xx=np.float32(xx)+difx/self.old.dx
    yy=np.float32(yy)+dify/self.old.dy
    return cv2.remap(self.old.image,xx,yy,cv2.INTER_LINEAR)
  def hyouka(self,parm):
    new=self.convert(parm[0],parm[1])
    res=np.corrcoef(inc[100:500,100:500].flat,new[100:500,100:500].flat)
    return -res[0,1]

def gwrite(sat,fnew):
  f=open('../'+fnew+'/gparm.txt','w')
  f.write(" * datum : WGS84\n")
  f.write(" image size:\n")  f.write("  pixel= "+str(sat.imax)+"\n")  f.write("  line= "+str(sat.jmax)+"\n")  f.write("image scene center:\n")  f.write("  p0= "+str(sat.pc)+"\n")  f.write("  l0= "+str(sat.lc)+"\n")  f.write("utm scene center:\n")  #f.write("  x0= "+str(sat.xc)+"\n")  #f.write("  y0= "+str(sat.yc)+"\n")  f.write('x0= {:12.2f} \n'.format(sat.xc))  f.write('y0= {:12.2f} \n'.format(sat.yc))  f.write("spatial resolution:\n")  f.write("  dx0= "+str(sat.dx)+"\n")  f.write("  dy0= "+str(sat.dy)+"\n")  f.write("orientation angle:\n")  f.write("  t0= "+str(sat.t0)+"\n")
  f.write("pointing angle:\n")
  f.write("  angle= "+str(sat.pangle)+"\n")  f.write("sun position:\n")  f.write("  el= "+str(sat.sun_el)+"\n")  f.write("  az= "+str(sat.sun_az)+"\n")
  pvx=['{:10.6f}'.format(x) for x in sat.pv]
  qvx=['{:10.6f}'.format(x) for x in sat.qv]
  f.write("scanning direction:\n")
  f.write("  pv= "+" ".join(pvx)+"\n")
  f.write("orbit direction:\n")
  f.write("  qv= "+" ".join(qvx)+"\n")
  f.write("nadir utm point:\n")
  #f.write("  xn=  "+str(sat.xn)+"\n")
  #f.write("  yn=  "+str(sat.yn)+"\n")
  f.write('xn= {:12.2f} \n'.format(sat.xn))
  f.write('yn= {:12.2f} \n'.format(sat.yn))
  f.write("  hsat= "+str(sat.hsat)+"\n")
  left,right,lower,upper=sat.coverage()
  f.write("coverage:"+"\n")
  f.write("  left = "+str(left)+"\n")
  f.write("  right= "+str(right)+"\n")
  f.write("  lower= "+str(lower)+"\n")
  f.write("  upper= "+str(upper)+"\n")
  offsetx=['{:7.2f}'.format(x) for x in sat.offset]
  gainx=['{:7.5f}'.format(x) for x in sat.gain]
  f.write("calibration:"+"\n")
  f.write("  offset= "+" ".join(offsetx)+"\n")
  f.write("  gain=   "+" ".join(gainx)+"\n")
  f.close()

def incident(dem,sat) :
  el=np.pi*sat.sun_el/180 ; az=np.pi*sat.sun_az/180
  jmax,imax=dem.shape
  a=(np.roll(dem,-1,1)-np.roll(dem,1,1))/60.0
  a[:,0]=a[:,1] ; a[:,imax-1]=a[:,imax-2] 
  b=(np.roll(dem,1,0)-np.roll(dem,-1,0))/60.0
  b[0,:]=b[1,:] ; b[jmax-1,:]=b[jmax-2,:]
  temp=-a*np.cos(el)*np.sin(az)-b*np.cos(el)*np.cos(az)+np.sin(el)
  return temp/np.sqrt(1+a**2+b**2)

#6/14/2016 copied from apm_util.py
def s_file(xs,ye,jmax,imax,sat):
  res=[]
  for y in range(jmax):
    xx=xs+np.arange(imax)*30
    yy=ye-np.ones(imax)*30*y
    s_ang=sat.sangle(xx,yy)
    res.append(s_ang)
  res = np.array(res)
  return res.reshape(jmax,imax)


