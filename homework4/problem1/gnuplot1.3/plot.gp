# generate latex formatted images with the following
# 
# produce a '.tex' and '.eps' file for latex to interpret
#
#   gnuplot -p plot.gp
#
# latex will combine the two files into a single dvi file and
#
#   latex image.tex
#
# dvips generates a postscript file and open for analysis
#
#   dvips -o image.ps image.dvi
#
# or just run the jobs file
#
#   /bin/bash jobs
#
set terminal epslatex size 4.5, 3.2 standalone color colortext 10
set out 'image.tex'

t(x,k) = 4/(pi*(2*k-1))*sin((2*k-1)*x)
p(x,k) = ( (1/(k+1)) * ((-1)**(k+1) - 1) + (1/(1-k))*((-1)**(k+1) - 1) ) * cos(k*x)

#an(x,k) = ( (cos(pi*(1+k)) - 1) / (1+k) + (cos(pi*(1-k)) - 1) / (1-k) ) * cos(k*x)
an(x,k) = cos(2*k*x) / ((2*k)**2 - 1)
f(x) = (x>0 ? sin(x) : 0)
test(k) = k

#series(x,n) = (n>0 ? t(x,n) + series(x,n-1) : 0)
#all(x, n) = (1 / pi) + 0.5*sin(x) - 2*((cos(2*x)/(2**2 - 1)) + (cos(4*x)/(4**2 - 1)) + (cos(6*x)/(6**2 - 1))) / pi # - series(x, n)

series(x,n) = (n>0 ? an(x,n) + series(x,n-1) : 0)
all(x, n) = (1 / pi) + 0.5*sin(x) - 2*series(x, n) / pi

set key right top
set style line 1 lw 1 ps 2 pt 6 lc 0
set style line 2 lt 1 pi 5 pt 11 lc 0
set style line 3 lt 1 pi 1 pt 6 ps 1 lc 0
set style line 4 lw 4 lt 2 lc 0

set title "1.3 Fourier Expansion" font ", 24"
set grid
set xlabel '$x$'
set ylabel '$f(x)$'
set mytics 5
#set xtics 0 0.01
set mxtics 5

set xrange [-pi:pi]

set samples 100
#set table "output.dat"
#plot f(x) with line linestyle 4, series(x,10) with line linestyle 1
plot all(x,2) with line linestyle 1, f(x) with line linestyle 4
#unset table

#plot "../data/output.dat" using 1:2 with line linestyle 4 title '$I = \int f dx$', \
#     "../data/output.dat" using 1:3 with linespoints linestyle 1 title '$I = w h \braket{f}$'

set out
