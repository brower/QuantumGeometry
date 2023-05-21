reset session
plot ARG1

### Function
r(x,y,z) = sqrt(2.0*x*x + ayy*y*y + azz*z*z + axy*x*y  + ayz*y*z  + azx*z*x )
f(x,y,z) = a*exp( - b *r(x,y,z)) /(r(x,y,z)**c) + d

### Initial values
a=4.0
b=0.05
c=1.0
d = 40.0
## axx = 2.0
ayy = 2.0
azz = 2.0
axy = 0.2
ayz = 0.2
azx = 0.2

### Fitting
set dummy x, y, z
fit f(x,y,z)  ARG1 using  3:4:5:2 via a,b,c,d, ayy, azz, axy, ayz, azx
replot f(x,0,0), f(0,x,0), f(0,0,x)

### Prepare output graph.
#set term postscript color
# set output  "affine16.eps"
# replot
 
# set term x11
 
