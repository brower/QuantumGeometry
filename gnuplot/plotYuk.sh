reset session
plot ARG1

### Function
f(x) = a*exp( - b *x)/(x**c) + d

### Initial values
a=2.0
b=0.5
c=1.0
d = 0.5

### Fitting
fit f(x)  ARG1 via a,b,c,d
replot f(x)
pause -1

### Prepare output graph.
#set term postscript color
# set output  "yuk16a.eps"
# replot
 
# set term x11
 
