#set term postscript landscape enhanced color 
#set terminal png  enhanced truecolor
set terminal pdfcairo size 6.75, 2 enhanced truecolor
set output 'fig6.pdf' # nom du fichier de sortie
set ticslevel 0.0
set border lw 1.5
set multiplot # plusieurs panel
unset key

#panel A
set origin -0.02,0.04
set size 0.29,0.9
set xtics 60 in nomirror font "Helvetica,17" offset -0.,0,0
set ytics 0.04 in nomirror font "Helvetica,16 offset 0,0,0
set xlabel 'Dissipation rate, T{/Symbol s} (kJ/l/h)' offset 0,0,0  font "Helvetica,16"
set ylabel 'J_{bex}' offset 2,0,0  font "Helvetica,20" rotate by 90
set border 3 front lc rgb "black"
set label "A" at -27,-0.055 font "Helvetica,24"
set label "B" at 73,-0.055 font "Helvetica,24"
set label "C" at 180,-0.055 font "Helvetica,24"
set label "D" at 262,-0.055 font "Helvetica,24"
set label "T{/symbol s}^{opt}{/Symbol \273}17" at 11,0.005 font "Helvetica,14"
set label "{/Symbol b}^{opt}{/Symbol \273}0.5" at 17,-0.034 font "Helvetica,14"
set label "{/Symbol b}^{max}{/Symbol \273}0.66" at 1,-0.044 font "Helvetica,14"
set palette defined (0. 'red', 0.5 'yellow',  1. 'blue' )
plot [-1:60][0.:-0.05] 'fig56_samp.dat' u 1:4 w p pt 7 ps 0.3 lc rgb 'orange','../MANIFOLD/jbexglu_fig5.dat' u 7:3 w l lw 6 lc rgb 'dark-gray','circ.dat' w p ps 1 pt 7 lc rgb 'black','line_fig6a.dat'  w l lw 5 dt 2 lc 'black'

#panel D
set border 31
unset label
set size 0.25,1.08
set origin 0.76,-0.08
set xtics 10 in nomirror font "Helvetica,14" offset 0.,0.,0 
set ytics 20 in nomirror font "Helvetica,14" offset 0,0,0 
set ztics 1 in nomirror font "Helvetica,14" offset 0.5,0,0
set xlabel 'S_{ATP,gly}' offset  2.,0.,0  font "Helvetica,17" rotate parallel
set ylabel 'S_{ATP,ares}' offset -2,0.,0   font "Helvetica,17" rotate parallel
set zlabel '{/Symbol b}^{opt}' offset 4,0.,0  font "Helvetica,19" rotate parallel
set palette defined (0. 'blue', 0.5 'yellow',  1. 'red' )
set cbrange [0:30]
set cbtics 30 font 'helvetica,15' offset 0,0.2
set colorbox horizontal user origin 0.81,0.89 size 0.15,0.025
set cblabel 'T{/symbol s}^{opt} (kJ/l/h)' font 'helvetica, 15' offset 0,1
set view 65,45
set pm3d interpolate 0,0
splot [22:0] [0:52] [0:1][0:30]  'fig6d.dat' u 1:2:3:8 w pm3d palette,'c3_pand.dat' w p pt 7 ps 0.7 lc rgb 'dark-blue','c2_pand.dat' w p pt 7 ps 1 lc rgb 'orange','c1_pand.dat' w p pt 7 ps 1 lc rgb 'dark-red','circ_panB.dat' w p pt 6 ps 1 lc rgb 'black'

#panel B
unset tics
unset label
unset xlabel
unset ylabel
set size 0.16,0.395
set origin 0.285,0.62
set ylabel "J_{bex}"  font 'helvetica, 15'
set object circle at 6,-0.036 front size 2.5 fillcolor rgb "black" fillstyle solid
plot [-1:60][0:-0.05] 'fig6b1.dat' u 1:4 w p ps 0.1 pt 7 lc rgb 'dark-blue','../MANIFOLD/jbexglu_fig5.dat' u 7:3 w p pt 7 ps 0.3 lc rgb 'dark-gray'
unset object

set origin 0.285,0.33
set object circle at 16,-0.03 front size 2.5 fillcolor rgb "black" fillstyle solid
plot [-1:60][0:-0.05] 'fig56_samp.dat' u 1:4 w p ps 0.1 pt 7 lc rgb 'orange','../MANIFOLD/jbexglu_fig5.dat' u 7:3 w p pt 7 ps 0.3 lc rgb 'dark-gray
set xlabel "T{/symbol s}" font 'helvetica, 15'
unset object

set origin 0.285,0.0
set size 0.16,0.435
set object circle at 43,-0.0043 front size 2.5 fillcolor rgb "black" fillstyle solid
plot [-1:60][0:-0.05] 'fig6b3.dat' u 1:4 w p ps 0.1 pt 7 lc rgb 'dark-orange','../MANIFOLD/jbexglu_fig5.dat' u 7:3 w p pt 7 ps 0.3 lc rgb 'dark-gray'

set xlabel '{/symbol  h}_{gly}' offset 0,1 font 'helvetica,20'
set ylabel '{/symbol  h}_{ares}' offset 1,0 font 'helvetica,20'
set xtics 1
set ytics 1
set origin 0.53,0.0
set size 0.23,0.71
plot [0:1][0:1] 'fig6b3.dat' u 14:15 w p ps 0.1 pt 7 lc rgb 'dark-orange','fig6b1.dat' u 14:15 w p ps 0.1 pt 7 lc rgb 'dark-blue','fig56_samp.dat' u 14:15 w p ps 0.1 pt 7 lc rgb 'orange','fig6c.dat' u 4:5 w p ps 0.5 pt 7 lc 'green'

unset label
unset xlabel
unset ylabel
unset xtics
unset ytics
unset key
unset border
set origin 0.39,0.02
set size 0.16,0.99
plot "fig6_panelb.png" binary filetype=png with rgbalpha
set origin 0.57,0.62
set size 0.21,0.42
plot "fig6_paneld.png" binary filetype=png with rgbalpha
