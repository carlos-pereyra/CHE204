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
set style line 3 lt 1 pi 1 pt 6 ps 1 lc 0
set style line 4 lw 5 lt 5 lc 1

#set title " Morse Potential V(r), \n and negative first derivative, F(r)" font ", 24"
set xlabel '$r$' offset 0
#set ylabel '$\braket{E}$' 
set mytics 5
#set xtics 0 0.01
set mxtics 5

set yrange [-2:10]
#rotate by 90 offset 0, graph -0.2
e(x) = 0.45

#set label 1 at  1, 2 '\hl{\small force comparison}'

set multiplot
plot "../../data/output.dat" using 1:4 with points linestyle 1 title '$-\nabla_iV(r)$ Forward Diff.', \
     "../../data/output.dat" using 1:5 with points linestyle 2 title '$-\nabla_iV(r)$ Backward Diff.', \
     "../../data/output.dat" using 1:6 with points linestyle 3 title '$-\nabla_iV(r)$ Central Diff.', \
     "../../data/output.dat" using 1:3 with line linestyle 4 title '$-\nabla_iV(r)$ Analytical' 

# smaller overlayed plot
set size 0.45, 0.45       # set size of inset
set origin 0.52, 0.22     # move bottom left corner of inset
set xrange [0.5:0.8]
set yrange [2:8]
set xtics 0 0.1
unset key
unset xlabel
plot "../../data/output.dat" using 1:4 with points linestyle 1 title '$-\nabla_iV(r)$ Forward Diff.', \
     "../../data/output.dat" using 1:5 with points linestyle 2 title '$-\nabla_iV(r)$ Backward Diff.', \
     "../../data/output.dat" using 1:6 with points linestyle 3 title '$-\nabla_iV(r)$ Central Diff.', \
     "../../data/output.dat" using 1:3 with line linestyle 4 title '$-\nabla_iV(r)$ Analytical'
unset multiplot
set out
