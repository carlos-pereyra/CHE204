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

# FOURIER EXPANSION SERIES
p(x,k) = (2*sin(k*pi*x)*(-1)**(k+1)) / (pi*k)
series(x,n) = (n>0 ? p(x,n) + series(x,n-1) : 0)

# ANALYTICAL FUNCTION
f(x) = x # conditional function = (x>0 ? 1 : -1)

# LEGEND OPTIONS
set key right top
set key spacing 4

# LINE STYLES
set style line 1 lw 1 ps 2 pt 6 lc 0
set style line 2 lw 2 lt 1 pi 5 pt 11 lc 1
set style line 3 lt 1 pi 1 pt 6 ps 1 lc 0
set style line 4 lw 4 lt 1 lc 0

set title "1.2 Fourier Expansion" font ", 24"
#set grid
set xlabel '$x$'
set ylabel '$f(x)$'
set mytics 5
#set xtics 0 0.01
set mxtics 5

set xrange [-1.1:1.1]
plot series(x,5) with line linestyle 1 title '$\frac{2}{n\pi} \sum_{n=1}^{5} \left(-\cos(n\pi) + \frac{2\sin(n\pi)}{n\pi}\right) \sin(n\pi x)$', series(x,100) with line linestyle 4 title '$\frac{2}{n\pi} \sum_{n=1}^{100} \left(-\cos(n\pi) + \frac{2\sin(n\pi)}{n\pi}\right) \sin(n\pi x)$', f(x) with line linestyle 2

#set yrange [4:5.5]
#rotate by 90 offset 0, graph -0.2
#set label 1 at  1, 2 '\hl{\small force comparison}'

#plot "../data/output.dat" using 1:2 with line linestyle 4 title '$I = \int f dx$', \
#     "../data/output.dat" using 1:3 with linespoints linestyle 1 title '$I = w h \braket{f}$'

set out
