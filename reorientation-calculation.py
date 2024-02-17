import numpy as np
import cv2
import copy

AANUMAA="1"
clat=20.    #crater center latitude
clon=90.    #crater center east longitude

title="T"+AANUMAA+"E_vol.txt"
with open(title) as f:
    im_name = [s.strip() for s in f.readlines()]

xx=[]
yy=[]
zz=[]
vol=[]
for l in im_name :
    lo,la,v,x,y,z = l.split()
    xx.append(float(x))
    yy.append(float(y))
    zz.append(float(z))
    vol.append(float(v))

del lo,la,im_name,x,y,z,v


num=len(vol)
Ixx=0.0
Iyy=0.0
Izz=0.0
Ixy=0.0
Ixz=0.0
Iyz=0.0

for i in range(num):
    x=-xx[i] /1188000.
    y=-yy[i] /1188000.
    z=zz[i] /1188000.
    v=vol[i]/10000.
    Ixx=Ixx+(y*y+z*z)*v
    Iyy=Iyy+(x*x+z*z)*v
    Izz=Izz+(y*y+x*x)*v
    Ixy=Ixy -x*y*v
    Ixz=Ixz -x*z*v
    Iyz=Iyz -y*z*v
    


del xx,yy,zz,vol

moi=np.array([ [Ixx,Ixy,Ixz],[Ixy,Iyy,Iyz],[Ixz,Iyz,Izz] ])
w,v = np.linalg.eig(moi)

t0=0
t1=1
t2=2
if w[t0] > w[t1] :
   tem=t0
   t0 = t1
   t1 = tem

if w[t1] > w[t2] :
   tem=t1
   t1 = t2
   t2 = tem

if w[t0] > w[t1] :
   tem=t0
   t0 = t1
   t1 = tem

if w[t1] > w[t2] :
   tem=t1
   t1 = t2
   t2 = tem

#t0,t1,t2

output=[]
output.append([AANUMAA])
output.append([w[t0],w[t1],w[t2],w[t1]/w[t0],w[t2]/w[t0]])

res=[]
t=t0
x=v[0][t]
y=v[1][t]
z=v[2][t]
lat=180.*np.arcsin(z/(x*x+y*y+z*z)**0.5)/np.pi
lon=180.*np.arctan2(y,x)/np.pi
if lon < 0 :
   lon=lon+360.

if abs(clon-lon) > 90. :
   lat=-lat
   lon=lon+180.

if lon > 360. :
   lon=lon-360.

res.append([lat,lon])

t=t1
x=v[0][t]
y=v[1][t]
z=v[2][t]
lat=180.*np.arcsin(z/(x*x+y*y+z*z)**0.5)/np.pi
lon=180.*np.arctan2(y,x)/np.pi
if lon < 0 :
   lon=lon+360.

res.append([lat,lon])

t=t2
x=v[0][t]
y=v[1][t]
z=v[2][t]
lat=180.*np.arcsin(z/(x*x+y*y+z*z)**0.5)/np.pi
lon=180.*np.arctan2(y,x)/np.pi
if lon < 0 :
   lon=lon+360.

res.append([lat,lon])

output.append(res)

###################

img = []
title="T"+AANUMAA+"E.bmp"
img.append(cv2.imread(title))
ylength = img[0].shape[0]
xlength = img[0].shape[1]
outimg=copy.copy(img[0])
# left top corner = (y=0,x=0) 

col=[]
col.append( [ 254, 254,  254]  )
col.append( [ 128, 128,  128]  )
col.append( [ 1, 1, 1 ]  )

for j in range(3):
    yy=int(round(90-res[j][0]))
    xx=int(round(res[j][1]))
    for i in range(-2,3) :
     if yy+i > 180 or yy+i < 0 or xx+i < 0 or xx+i > 359 or xx < 0 or xx >359 or yy > 180 or yy < 0:
       continue
     outimg[yy+i][xx][0]=col[j][0]
     outimg[yy+i][xx][1]=col[j][1]
     outimg[yy+i][xx][2]=col[j][2]
     outimg[yy][xx+i][0]=col[j][0]
     outimg[yy][xx+i][1]=col[j][1]
     outimg[yy][xx+i][2]=col[j][2]

for j in range(3):
    yy=int(round(90+res[j][0]))
    tem=180+res[j][1]
    if tem > 360. :
       tem=tem-360.
    xx=int(round(tem))
    for i in range(-2,3) :
     if yy+i > 180 or yy+i < 0 or xx+i < 0 or xx+i > 359 or xx < 0 or xx >359 or yy > 180 or yy < 0:
       continue
     outimg[yy+i][xx][0]=col[j][0]
     outimg[yy+i][xx][1]=col[j][1]
     outimg[yy+i][xx][2]=col[j][2]
     outimg[yy][xx+i][0]=col[j][0]
     outimg[yy][xx+i][1]=col[j][1]
     outimg[yy][xx+i][2]=col[j][2]

ti="T"+AANUMAA+"F.bmp"
cv2.imwrite(ti, outimg)

##################

lat= res[0][0]
lon= res[0][1]
shi=(180.-lon)*np.pi/180.0
fai=-lat*np.pi/180.0

z_rot=np.matrix([ [np.cos(shi),-np.sin(shi),0],[np.sin(shi),np.cos(shi),0],[0,0,1] ])
y_rot=np.matrix([ [np.cos(fai),0,np.sin(fai)],[0,1,0],[-np.sin(fai),0,np.cos(fai)] ])
x=v[0][t2]
y=v[1][t2]
z=v[2][t2]
tem=np.matrix([[x],[y],[z]])
tem2=y_rot*(z_rot*tem)
x=float(tem2[0])
y=float(tem2[1])
z=float(tem2[2])
lat=np.arcsin(z/(x*x+y*y+z*z)**0.5)
ram=(np.pi*0.5-lat)

cx=np.cos(clon*np.pi/180.)*np.cos(clat*np.pi/180.)
cy=np.sin(clon*np.pi/180.)*np.cos(clat*np.pi/180.)
cz=np.sin(clat*np.pi/180.)
ctem=np.matrix([[cx],[cy],[cz]])
ctem2=y_rot*(z_rot*ctem)
if float(ctem2[2]) > 0.01 :
    ram=ram+np.pi

x_rot=np.matrix([ [1,0,0],[0,np.cos(ram),-np.sin(ram)],[0,np.sin(ram),np.cos(ram)] ])


ctem3=x_rot*(y_rot*(z_rot*ctem))
clat2=180.*np.arcsin(float(ctem3[2]))/np.pi
clon2=180.*np.arctan2(float(ctem3[1]),float(ctem3[0]))/np.pi
if clon2 < 0. :
   clon2=clon2+360.
output.append([clat,clon,clat2,clon2])
print(output)

invx=x_rot**-1
invy=y_rot**-1
invz=z_rot**-1

outimg2=copy.copy(img[0])
for la in range(ylength):
  for lo in range(xlength):
    lat=np.pi*(90-la)/180.
    lon=np.pi*lo/180.
    x=np.cos(lon)*np.cos(lat)
    y=np.sin(lon)*np.cos(lat)
    z=np.sin(lat)
    tem3=np.matrix([[x],[y],[z]])
    tem4=invz*(invy*(invx*tem3))
    x=float(tem4[0])
    y=float(tem4[1])
    z=float(tem4[2])
    lat=180.*np.arcsin(z/(x*x+y*y+z*z)**0.5)/np.pi
    lon=180.*np.arctan2(y,x)/np.pi
    if lon < 0 :
       lon=lon+360.
    yy=int(round(90-lat))
    xx=int(round(lon))
    if xx == 360:
       xx = 0
    outimg2[la][lo][0]=img[0][yy][xx][0]
    outimg2[la][lo][1]=img[0][yy][xx][1]
    outimg2[la][lo][2]=img[0][yy][xx][2]

yy=90
xx=180
for i in range(-2,3):
    outimg2[yy+i][xx][0]=col[0][0]
    outimg2[yy+i][xx][1]=col[0][1]
    outimg2[yy+i][xx][2]=col[0][2]
    outimg2[yy][xx+i][0]=col[0][0]
    outimg2[yy][xx+i][1]=col[0][1]
    outimg2[yy][xx+i][2]=col[0][2]

yy=90
xx=0
for i in range(-2,3):
    if xx +i < 0 or xx+i >359:
       continue
    outimg2[yy+i][xx][0]=col[0][0]
    outimg2[yy+i][xx][1]=col[0][1]
    outimg2[yy+i][xx][2]=col[0][2]
    outimg2[yy][xx+i][0]=col[0][0]
    outimg2[yy][xx+i][1]=col[0][1]
    outimg2[yy][xx+i][2]=col[0][2]

yy=90
xx=359
for i in range(-2,3):
    if xx +i < 0 or xx+i >359:
       continue
    outimg2[yy+i][xx][0]=col[0][0]
    outimg2[yy+i][xx][1]=col[0][1]
    outimg2[yy+i][xx][2]=col[0][2]
    outimg2[yy][xx+i][0]=col[0][0]
    outimg2[yy][xx+i][1]=col[0][1]
    outimg2[yy][xx+i][2]=col[0][2]

yy=90
xx=90
for i in range(-2,3):
    outimg2[yy+i][xx][0]=col[1][0]
    outimg2[yy+i][xx][1]=col[1][1]
    outimg2[yy+i][xx][2]=col[1][2]
    outimg2[yy][xx+i][0]=col[1][0]
    outimg2[yy][xx+i][1]=col[1][1]
    outimg2[yy][xx+i][2]=col[1][2]

yy=90
xx=270
for i in range(-2,3):
    outimg2[yy+i][xx][0]=col[1][0]
    outimg2[yy+i][xx][1]=col[1][1]
    outimg2[yy+i][xx][2]=col[1][2]
    outimg2[yy][xx+i][0]=col[1][0]
    outimg2[yy][xx+i][1]=col[1][1]
    outimg2[yy][xx+i][2]=col[1][2]

ti="T"+AANUMAA+"G.bmp"
cv2.imwrite(ti, outimg2)

exit()
