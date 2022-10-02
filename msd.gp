reset
set encoding utf8
set log
set key Left
set key left

#set tics scale 3.0
set border linewidth 2.5
set format y "10^{%T}"
set format x "10^{%T}"
set key font ",18"
set key reverse

LW = 5
PS = 2
LT = 7
EV = 1

dataFile = "MeanSqd_1_0.04405_1_0_1.dat"

set term post enh col font "Latin_Modern_Roman" 22 size 20cm,16cm
set out "conta_msd.eps"

phi=0
aux=0.04405
#aux=1.7425

g=10.0
#k=aux/phi
k=0.04405
q=0.1
y=1.0
b=1.0

S = 2.*q*k*(y+2.*k)*(y+k) / (g+2.*q*k*(y+2.*k)*(y+k))
P = 1./(y+2.*k)

#set xrange [S/10:1000*P]
#set xrange [0.0001:1000]

set xlabel "{/:Italic {/Symbol D} T}"
set ylabel "{/:Italic MSD({/Symbol D} T)}"


f(x) = (g/((y+2.*k)*(y+k)))*(x - (1./(y+2.*k))*(1.-exp(-(y+2.*k)*x)) ) + 2.*q*k*x #MSD Normal


msd(x) = (2*(b*b)/((y+k)*(y+k))) * (-1/(k*k)-1/(2*(y+k)*(y+2*k)) - exp(-(y+2*k)*x)/(y*k) + exp(-k*x)/(k*k) + exp(-(y+k)*x)/(k*(y+k)) + exp(-2*(y+k)*x)/(2*y*(y+k)) + x/k + exp(-(y+2*k)*x)/(y*(y+2*k))) + f(x)

plot dataFile every EV u 1:2 w lp lt LT ps PS lw LW lc rgb "blue" title sprintf("k=%g b=%g",k,b),\
     msd(x) lw LW lc rgb "black" title "Analytical Biased MSD",\
     f(x) lw LW lc rgb "red" title "Analytical Unbiased MSD"
