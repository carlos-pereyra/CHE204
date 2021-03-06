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
#set terminal epslatex size 5.5, 5.5 standalone color colortext 10 header "\\newcommand{\\ft}[0]{\\huge}"

# ANIMATION
set terminal gif animate delay 2
set output 'tmp.gif'

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
set title 'potential, gamma=2'
#set grid
#set xlabel '\large $x$'
#set ylabel '\large $y$'
set xlabel 'X'
set ylabel 'Y'
set mytics 5
#set xtics 0 0.01
set mxtics 5

# HEAT MAP
#set pm3d map
set pm3d at b

# CONTOUR
#set contour base
#set cntrparam levels 21
#set cntrparam levels incremental -100,10
#set view map
#unset surface
#set contour

# SHOW
#set out 'fig0.tex'
#splot "../data/potential0.dat" u 1:2:3 title ''

#set out 'fig1.tex'
#splot "../data/potential1.dat" u 1:2:3 title ''
#do for [i=0:200] { set out sprintf('fig%d.tex', i); splot sprintf('../data/potential%d.dat', i) using 1:2:3 title sprintf('time = %d', i) with lines}
#do for [i=0:200] { set out sprintf('fig%d.tex', i); splot sprintf('../data/potential%d.dat', i) using 1:2:3 title sprintf('time = %d', i) with lines}

#set zrange [0:0.016];
#set cbrange [0:0.016]
do for [i=1:327] {
    splot sprintf('../data/pottmp%d.dat', i) using 1:2:3 title sprintf('iteration = %d', i) with dots
}

