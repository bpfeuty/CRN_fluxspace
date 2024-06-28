#set term postscript landscape enhanced color 
#set terminal png  enhanced truecolor
set terminal pdfcairo size 6.5, 3.75 enhanced truecolor
set output 'fig4.pdf' # nom du fichier de sortie
set ticslevel 0.0
set border lw 1.5
set multiplot # plusieurs panel

unset key
unset tics
unset border
unset label

set size 0.15, 0.2
set origin 0.4, 0.84
set label at 1.8,0 "n_{m}=1" font "Helvetica,15"
plot [-2:2][-1:1] 'circle1.dat' w lp pt 7 ps 1 lc rgb "black"
unset label
set origin 0.75, 0.84
set label at 1.8,0 "n_{m}=2" font "Helvetica,15"
plot [-2:2][-1.3:1.3] 'circle2.dat' w p pt 7 ps 1 lc rgb "black",'circle3.dat' w l smooth csplines lc rgb "black",'circle4.dat' w l smooth csplines lc rgb "black"
unset label

# PANEL A & G
unset label
set origin 0.0,0.56
set size 0.35,0.44
plot "fig4_panel_a.png" binary filetype=png with rgbalpha
set origin 0.33,0.0
set size 0.35,0.46
plot "fig4_panel_e.png" binary filetype=png with rgbalpha

set border
set size 0.4,0.475   # taille du plot
set cbtics 1   font 'Helvetica,13' offset -1,0,0
set xlabel 'Exchange flux, J_1' offset 0,0,0  font "Helvetica,15"
set ylabel 'Exchange flux, J_2' offset 1.5,0,0  font "Helvetica,15"
set palette defined ( 0 'black', 0.2 'blue', 0.5 'green', 0.7 'yellow', 0.85 'orange', 1. 'red' )
set cbrange [0:1]
set xtics 0.5 in nomirror font "Helvetica,13" offset 0,0,0
set ytics 0.5 in nomirror font "Helvetica,13" offset 0,0,0

set origin 0.3,0.46
set label "{/Symbol s}={/Symbol s}*" at 0.65,0.81 font "Helvetica,15"
set label '{/Symbol s}^i/{/Symbol s}' font "Helvetica,15" at 0.98,1.15 

set label "A" at -3.1,1.27 font "Helvetica,22"
set label "B" at -0.9,1.27 font "Helvetica,22"
plot [-0.5:1.][-0.5:1.][-1.8:0.2] 'fig4b_nm1.dat' u 2:3:5 w p pt 7 ps 0.2 palette,'manifold_fig4.dat' u 1:2 w p pt 5 ps 0.2 lc 0 lt 1 ,'z.dat' w p pt 7 ps 1 lc rgb 'violet'
unset label

set origin 0.65,0.46
set label "C" at -0.9,1.27 font "Helvetica,22"
set label '{/Symbol s}^i/{/Symbol s}' font "Helvetica,14" at 0.98,1.15 	
set label "{/Symbol s}={/Symbol s}*" at 0.65,0.81 font "Helvetica,15"
plot [-0.5:1.][-0.5:1.][-1.8:0.2] 'fig4c_nm2.dat' u 2:3:5 w p pt 7 ps 0.2 palette,'manifold_fig4.dat' u 1:2 w p pt 5 ps 0.2 lc 0 lt 1 ,'z.dat' w p pt 7 ps 1 lc rgb 'violet'

# panel G
unset label
unset cblabel
set cbrange [-2:2]
set cbtics 2
set origin 0.65,-0.03
set label "{/Symbol m}_4-{/Symbol m}_5" at 0.98,1.17 font "Helvetica,13"
set label "G" at -0.9,1.22 font "Helvetica,22"
set label "D" at -5.6,1.86 font "Helvetica,22"
set label "E" at -5.6,0.44 font "Helvetica,22"
set label "F" at -3.3,1.22 font "Helvetica,22"
set label "{/ZapfChancery V} {/Symbol \273} 0.82" at  0.,1.12 font "Helvetica,15"
plot [-0.5:1.][-0.5:1.][-1.5:1.5] 'fig4g.dat' u 2:3:6 w p pt 7 ps 0.4 palette,'manifold_fig4.dat' u 1:2 w p pt 5 ps 0.2 lc 0 lt 1 

unset label
unset cblabel
unset xlabel
unset ylabel
unset tics

# panel C
set size 0.14,0.175
unset colorbox
set palette defined (-1 "slateblue1", 0 "white", 1 "salmon")
set origin 0.02,0.45
plot 'matrixtt02b.data' matrix with image pixels
set origin 0.11,0.45
plot 'matrixtt05a.data' matrix with image pixels
set origin 0.2,0.45
plot 'matrixtt08a.data' matrix with image pixels

set size 0.14,0.175
set origin 0.02,0.33
set label "{/ZapfChancery V} {/Symbol \273} 0.1" at  -0.4,-0.9 font "Helvetica,13"
plot [-0.5:1.][-0.5:1.][-1.8:0.2] 'polyto02b.data' w p pt 7 ps 0.1 lc rgb "grey", 'manifold_fig4.dat' u 1:2 w p pt 5 ps 0.1 lc 0 lt 1
set origin 0.11,0.33
unset label
set label "{/ZapfChancery V} {/Symbol \273} 0.5" at  -0.4,-0.9 font "Helvetica,13"
plot [-0.5:1.][-0.5:1.][-1.8:0.2] 'polyto05a.data' w p pt 7 ps 0.1 lc rgb "grey",'manifold_fig4.dat' u 1:2 w p pt 5 ps 0.1 lc 0 lt 1 
set origin 0.20,0.33
unset label
set label "{/ZapfChancery V} {/Symbol \273} 0.8" at -0.4,-0.9 font "Helvetica,13"
plot [-0.5:1.][-0.5:1.][-1.8:0.2] 'polyto08a.data' w p pt 7 ps 0.1 lc rgb "grey",'manifold_fig4.dat' u 1:2 w p pt 5 ps 0.1 lc 0 lt 1

#panel D
unset label
set size 0.33,0.3
set origin 0.01,0.01
set xlabel 'Solution space size, {/ZapfChancery V}' offset 0,0.7  font "Helvetica,15"
set ylabel 'Distribution' offset 1.,0,0  font "Helvetica,15"
set xtics 0.2 in nomirror font "Helvetica,13" offset 0,0.3,0
unset ytics # 50 in nomirror font "Helvetica,13" offset 0,0,0
set style fill transparent solid 0.5 noborder
set style function filledcurves y1=0
set key samplen 0.5 left
plot [0.05:0.85][0.:8] 'fig4e_ncc2_nm2.dat' u 1:2 smooth kdensity bandwidth 0.02 with filledcurves lt 2 title 'n_m=2,n_{cc}=2','fig4e_ncc2_nm1.dat' u 1:2 smooth kdensity bandwidth 0.02 with filledcurves lt 9 title 'n_m=1,n_{cc}=2','fig4e_ncc0_nm2.dat' u 1:2 smooth kdensity bandwidth 0.02 with filledcurves lt 5 title 'n_m=2,n_{cc}=0','c1.dat' w p ps 0.7 pt 6 lc rgb 'black' notitle 
