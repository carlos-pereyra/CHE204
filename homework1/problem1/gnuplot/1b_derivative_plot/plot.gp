# generate latex formatted images with the following
# 
# produce a '.tex' and '.eps' file for latex to interpret
#
#   gnuplot -p plot.gp
#
# latex will combine the two files into a single dvi file and 
#
#   latex doc.tex
#
# dvips generates a postscript file and open for analysis
#
#   dvips -o image.ps doc.dvi
#
set terminal epslatex size 4.5, 3.2 standalone color colortext 10
set out 'image.tex'

set key right top
set style line 1 lw 1      pt 8 lc 0
set style line 2 lt 1 pi 5 pt 11 lc 0
set style line 3 lt 1 pi 1 pt 14 lc 0
set style line 4 lw 3 lt 5 lc 0

set title 'Negative First Derivative Potential $-\nabla_iV(r)$' font ", 24"
set xlabel '$r$' offset 0
set ylabel '$F$' 
set mytics 5
#set xtics 0 0.01
set mxtics 5

set yrange [-2:10]
#rotate by 90 offset 0, graph -0.2
e(x) = 0.45

#set label 1 at  .028, 0.445 '\hl{\small kinetic energy at T = 0.3}'

plot "../../data/output.dat" using 1:3 with points linestyle 2 title '$-\nabla_iV(r)$'

set out
