#set term postscript landscape enhanced color 
#set terminal png  enhanced truecolor
set terminal pdfcairo size 6.75, 3.5 enhanced truecolor
set output 'fig2.pdf' # nom du fihier de sortie
set ticslevel 0.0
set multiplot # plusieurs panel

unset key
unset border
unset tics
unset label

set origin 0.1,0.57
set size 0.24,0.43
plot "panel_a.png" binary filetype=png with rgbalpha
set origin 0.44,-0.02
set size 0.54,0.97
plot "panel_cd.png" binary filetype=png with rgbalpha
set border

set border lw 1.5
set key center top
unset ytics
set style fill transparent solid 0.5 noborder
set xtics 5. out nomirror font "Helvetica,14" offset 0,0.3,0
set xlabel 'Dissipation rate, T{/Symbol s}' offset 0,0.8,0  font "Helvetica,15"

set origin 0.07,0.01
set size 0.34,0.285
set key samplen 0.5
plot [-0.:10.][0:8.8]  'dis_n8_10.dat'  smooth kdensity bandwidth 0.2 with filledcurves y=0 lt 4 title "k_r=10",'dis_n8_1.dat' u 1:2 smooth kdensity bandwidth 0.2 with filledcurves y=0 lt 7 title "k_r=1",'dis_n8_01.dat'  smooth kdensity bandwidth 0.2 with filledcurves y=0 lt 9 title "k_r=0.1" ,'dis_n8_001.dat'  smooth kdensity bandwidth 0.2 with filledcurves y=0 lt 2 title "k_r=0.01"
unset label

set size 0.361,0.22
set origin 0.05,0.23
unset  xlabel
unset xtics
set ylabel 'Sampled distribution' offset 0.5,0,0  font "Helvetica,16"
plot [-0.1:10][0:1.5] 'dis_n16_nm1.dat' u 1:2 smooth kdensity bandwidth 0.2 with filledcurves y=0 lt 7 title "n_s=16,n_r=20,n_m=1",'dis_n8_nm2.dat' u 1:2 smooth kdensity bandwidth 0.2 with filledcurves y=0 lt 9 title "n_s=8,n_r=12,n_m=2"

set size 0.34,0.22
set origin 0.07,0.388
unset ylabel
set label "B" at -3,1.6 font "Helvetica, 24"
set label "A" at -3,5.6 font "Helvetica, 24"
set label "C" at 12.,5.6 font "Helvetica, 24"
set arrow from 9.74,-3.3 to 9.74,1.7 nohead  lw 2 dt 2 lc 'black' 
set label "T{/Symbol s}*" at 9.4,1.85 font "Helvetica, 15"
unset xtics
plot [-0.1:10][0:1.5] 'dis_n8_nm1.dat' u 1:2 smooth kdensity bandwidth 0.2 with filledcurves y=0 lt 4 title "n_s=8,n_r=12,n_m=1"

unset arrow
unset key
unset border
unset tics
unset label
unset xlabel
unset ylabel

set size 0.15, 0.2
set origin 0.5, 0.83
set label at 2.,0 "n_{m}=1" font "Helvetica,14"
plot [-2:2][-1:1] 'circle1.dat' w lp pt 7 ps 1 lc rgb "black"
unset label
set origin 0.73, 0.83
set label at 2.,0 "n_{m}=2" font "Helvetica,14"
plot [-2:2][-1.3:1.3] 'circle2.dat' w p pt 7 ps 1 lc rgb "black",'circle3.dat' w l smooth csplines lc rgb "black",'circle4.dat' w l smooth csplines lc rgb "black"
unset label