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
set terminal epslatex size 4.5, 3.2 standalone color colortext 10
set out 'image.tex'

set key right top
set style line 1 lw 1 ps 2 pt 6 lc 0
set style line 2 lt 1 pi 5 pt 11 lc 0
set style line 3 lt 1 pi 1 pt 6 ps 1 lc 0
set style line 4 lw 4 lt 2 lc 0

set title "Integration Result" font ", 24"
set grid
set logscale x 10
set xlabel '$n$'
set ylabel '$I$' 
set mytics 5
#set xtics 0 0.01
set mxtics 5

set yrange [4:5.5]
#rotate by 90 offset 0, graph -0.2
#set label 1 at  1, 2 '\hl{\small force comparison}'

plot "../data/output.dat" using 1:2 with line linestyle 4 title '$I = \int f dx$', \
     "../data/output.dat" using 1:3 with linespoints linestyle 1 title '$I = w h \braket{f}$'

set out
