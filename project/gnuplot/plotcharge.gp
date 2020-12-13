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
set terminal epslatex size 5.5, 5.5 standalone color colortext 10 header "\\newcommand{\\ft}[0]{\\huge}"
set out 'chargeplot.tex'

# LEGEND OPTIONS
#set key right top

# LINE STYLES 
set style line 1 lw 1 ps 2 pt 6 lc 0
set style line 2 lt 1 pi 5 pt 11 lc 0
set style line 3 lt 1 pi 1 pt 6 ps 1 lc 0
set style line 4 lw 4 lt 2 lc 0

# MULTIPLOT
#set multiplot layout 1,2

# PLOT STYLES
set size square
#set title '$f(x,y) = \sum_{n,m=-4}^{4} e^{\left(-(x-nr_x^{(1)} -mr_x^{(2)})^2/\sigma -\right)}$'
set title 'Charge Density $\rho(x,y)$'
#set grid
set xlabel '\large $x$'
set ylabel '\large $y$'
set mytics 5
#set xtics 0 0.01
set mxtics 5

# HEAT MAP
set pm3d map

# CONTOUR
#set contour base
#set cntrparam levels 21
#set cntrparam levels incremental -100,10
#set view map
#unset surface
#set contour

# SHOW
splot "../data/charge0.dat" u 1:2:3 title ''
#unset view
#unset pm3d

