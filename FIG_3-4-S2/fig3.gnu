#set term postscript landscape enhanced color 
#set terminal png  enhanced truecolor
set terminal pdfcairo size 6.5, 3.75 enhanced truecolor
set output 'fig3.pdf' # nom du fichier de sortie
set ticslevel 0.0
set border lw 1.5
set multiplot # plusieurs panel

unset key
unset tics
unset cblabel
unset xlabel
unset ylabel

set size 0.135,0.17
unset colorbox
set palette defined (-1 "slateblue1", 0 "white", 1 "salmon")
set origin 0.03,0.45
plot 'mat3_0c.data' matrix with image pixels
unset label
set origin 0.115,0.45
plot 'mat3_02a.data' matrix with image pixels
set origin 0.2,0.45
plot 'mat3_036a.data' matrix with image pixels

unset border
unset label
set size 0.15, 0.2
set origin 0.41, 0.84
set label at 1.8,0 "n_{m}=1" font "Helvetica,14"
plot [-2:2][-1:1] 'circle1.dat' w lp pt 7 ps 1 lc rgb "black"
unset label
set origin 0.75, 0.84
set label at 1.8,0 "n_{m}=2" font "Helvetica,14"
plot [-2:2][-1.3:1.3] 'circle2.dat' w p pt 7 ps 1 lc rgb "black",'circle3.dat' w l smooth csplines lc rgb "black",'circle4.dat' w l smooth csplines lc rgb "black"
unset label

set origin 0.01,0.6
set size 0.3,0.42
plot "fig3_panel_a.png" binary filetype=png with rgbalpha
set origin 0.34,0.0
set size 0.35,0.5
plot "fig3_panel_f.png" binary filetype=png with rgbalpha
set border

#PANEL BCG
set colorbox
set size 0.37,0.46   # taille du plot
set cbtics 1   font 'Helvetica,13' offset -1,0,0
set xlabel 'Entropy production rate, {/Symbol s}/{/Symbol s}^*' offset 0,0.3,0  font "Helvetica,14"
set ylabel 'Exchange flux, J_3' offset 1.5,0,0  font "Helvetica,15"
set palette defined ( 0 'black', 0.2 'blue', 0.5 'green', 0.7 'yellow', 0.85 'orange', 1. 'red' )
set cbrange [0:1]
set xtics 1 in nomirror font "Helvetica,13" offset 0,0,0
set ytics 0.5 in nomirror font "Helvetica,13" offset 0,0,0

#PANEL Bs
set origin 0.33,0.46
set label "{/Symbol s}={/Symbol s}*" at 0.8,0.61 font "Helvetica,14"
set label '{/Symbol s}^i/{/Symbol s}' font "Helvetica,14" at 1.08,0.82
set label "-J_3=J^*" at 0.3,-0.41 font "Helvetica,15"
set label "A" at -2.7,0.97 font "Helvetica,23"
set label "D" at -2.7,-0.6 font "Helvetica,23"
set label "E" at -2.7,-1.7 font "Helvetica,23"
set label "B" at -0.5,0.97 font "Helvetica,23"
plot [-0.1:1.1][-0.5:0.7][-1.8:0.2] 'fig3b_nm1.dat' u 1:4:5 w p pt 7 ps 0.2 palette,'manifold_fig3.dat' u 5:3 w p pt 5 ps 0.2 lc 0 lt 1 ,'fig3_circ.dat' w p pt 7 ps 1 lc rgb 'violet'
unset label

#PANEL C
set origin 0.66,0.47
set label "{/Symbol s}={/Symbol s}*" at 0.8,0.61 font "Helvetica,14"
set label '{/Symbol s}^i/{/Symbol s}' font "Helvetica,14" at 1.08,0.82
#set label "-J_3=J^*" at 0.3,-0.41 font "Helvetica,15"
set label "C" at -0.5,0.96 font "Helvetica,22"
plot [-0.1:1.1][-0.5:0.7][-1.8:0.2] 'fig3c_nm2_50K.dat' u 1:4:5 w p pt 7 ps 0.2 palette,'manifold_fig3.dat' u 5:3 w p pt 5 ps 0.2 lc 0 lt 1 ,'fig3_circ.dat' w p pt 7 ps 1 lc rgb 'violet'

#PANEL G
set origin 0.66,-0.01
set size 0.32,0.47
unset label
set label "-J_3=J^*" at 0.3,-0.41 font "Helvetica,15"
set label "G" at -0.5,0.84 font "Helvetica,23"
set label "F" at -2.3,0.84 font "Helvetica,23"

plot [-0.1:1.1][-0.5:0.7][-1.8:0.2] 'fig3g_nb2_100K.dat' u 1:4:5 w p pt 7 ps 0.1 lc rgb 'gray','manifold_fig3.dat' u 5:3 w p pt 5 ps 0.2 lc 0 lt 1 notitle,'fig3g.dat' u 1:4 w p ps 0.2 lc rgb "red",'fig3_circ.dat' w p pt 7 ps 1 lc rgb 'violet' notitle
unset key


#panel E (distribution)
unset label
set size 0.33,0.35
set origin 0.01,-0.02
set xlabel 'Maximum anabolic flux, (-J_3)^{max}' offset 0,0.4,0  font "Helvetica,15"
set ylabel 'Distribution' offset 0.4,0,0  font "Helvetica,15"
set xtics 0.2 in nomirror font "Helvetica,13" offset 0,0.3,0
unset ytics # 50 in nomirror font "Helvetica,13" offset 0,0,0
set style fill transparent solid 0.5 noborder
set key samplen 0.5
set label "J^*" at 0.38,0.54 font "Helvetica,15"
plot [-0.:0.41][0:4]  'fig3e_nec2_nm2.dat' u 1:2 smooth kdensity bandwidth 0.02 with filledcurves y=0 lt 7 title "n_{m}=2,n_{bc}=1",'fig3e_nec1_nm2.dat' u 1:($2*1.4) smooth kdensity bandwidth 0.02 with filledcurves y=0 lt 9 title "n_{m}=2,n_{bc}=2",'fig3e_nec1_nm1.dat' u 1:2 smooth kdensity bandwidth 0.02 with filledcurves y=0 lt 4 title "n_{m}=1,n_{bc}=1",'fig3e_c.res' w p pt 6 lc rgb 'black' notitle
#title "(-J_3)^{max}=J^*
unset key

#PANEL D (examples)
unset label
unset xlabel
unset ylabel
unset tics
set size 0.135,0.17
set origin 0.03,0.34
set label "(-J_3)^{max}{/Symbol \273}0" at  -0.15,-0.7 font "Helvetica,11"
plot [-0.1:1.1][-0.43:0.6][-1.8:0.2]  'pol_0c.data' u ($1*0.15):2 w p ps 0.2 lc rgb "grey",'manifold_fig3.dat' u 5:3 w p pt 5 ps 0.1 lc 0 lt 1
set origin 0.115,0.34
unset label
set label "(-J_3)^{max}>0" at  -0.1,-0.7 font "Helvetica,11"
plot [-0.1:1.1][-0.43:0.6][-1.8:0.2]'pol_02a.data' u ($1*0.15):2 w p ps 0.2 lc rgb "grey",'manifold_fig3.dat' u 5:3 w p pt 5 ps 0.1 lc 0 lt 1
set origin 0.2,0.34
unset label
set label "(-J_3)^{max}{/Symbol \273}J^*" at -0.,-0.7 font "Helvetica,11"
plot [-0.1:1.1][-0.43:0.6][-1.8:0.2] 'pol_036a.data' u ($1*0.15):2 w p ps 0.2 lc rgb "grey",'manifold_fig3.dat' u 5:3 w p pt 5 ps 0.1 lc 0 lt 1
