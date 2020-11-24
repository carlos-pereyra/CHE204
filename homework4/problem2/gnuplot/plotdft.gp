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
set out 'imagedft.tex'

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
set title '$\log\left(|| \mathcal{F}[f(x,y)] ||\right)$'
set xlabel '\large $j$'
set ylabel '\large $k$'

# HEAT MAP
#set pm3d map
#set pal gray;

# CONTOUR
#set auto
#set surface
#set view map
#unset key
#unset surface
#set contour base
set contour
set cntrparam levels 20
#set cntrparam levels incremental -100,10
# 2.4 SHOW
set hidden3d
#set dgrid3d 50,50 qnorm 2
set dgrid3d 50,50
set log z
#set view map
#unset surface
#set contour
splot "../data/dft2d_abs.dat" u 1:2:3 with lines title ''

