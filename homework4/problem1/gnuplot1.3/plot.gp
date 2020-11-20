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

# ANALYTICAL FUNCTION
f(x) = (x>0 ? sin(x) : 0)

#series(x,n) = (n>0 ? t(x,n) + series(x,n-1) : 0)
#all(x, n) = (1 / pi) + 0.5*sin(x) - 2*((cos(2*x)/(2**2 - 1)) + (cos(4*x)/(4**2 - 1)) + (cos(6*x)/(6**2 - 1))) / pi # - series(x, n)

series(x,n) = (n>0 ? an(x,n) + series(x,n-1) : 0)
all(x, n) = (1 / pi) + 0.5*sin(x) - 2*series(x, n) / pi

# LEGEND OPTIONS
set key right top
set key spacing 4

# LINE STYLES
set style line 1 lw 1 ps 2 pt 6 lc 0
set style line 2 lw 2 lt 1 pi 5 pt 11 lc 1
set style line 3 lt 1 pi 1 pt 6 ps 1 lc 0
set style line 4 lw 4 lt 1 lc 0

# PLOT STYLE
set title "1.3 Fourier Expansion" font ", 24"
#set grid
set xlabel '$t$'
set ylabel '$f(t)$'
set mytics 5
#set xtics 0 0.01
set mxtics 5
set xrange [-pi:pi]

plot all(x,2) with line linestyle 1 title '$\frac{1}{\pi} + \frac{1}{2}\sin(t) + \sum_{n=1}^{\infty} \frac{-2}{\pi((2n)^2 - 1)} \cos(2n t)$', f(x) with line linestyle 2

#plot "../data/output.dat" using 1:2 with line linestyle 4 title '$I = \int f dx$', \
#     "../data/output.dat" using 1:3 with linespoints linestyle 1 title '$I = w h \braket{f}$'

set out
