reset
#
set style line 1 lw 4 pt 7 ps 1.8 lc "black"
set style line 2 lw 4 pt 7 ps 1.8 lc "red"
set style line 3 lw 4 pt 7 ps 1.8 lc "blue"
set style line 4 lw 4 pt 3 ps 1.5 lc "red" dt (2,4)
set style line 5 lw 4 pt 3 ps 1.5 lc "blue" dt (2,4)
set style line 6 lw 2 pt 5 ps 1.5 lc "red"
set style line 7 lw 2 pt 5 ps 1.5 lc "blue"
#
set size 1.5,1.8
#
set ter pos eps color 32 enhanced
#
set yr [10:0]
set xl "s_{xx} (kPa)"
set yl "depth (m)"
#
set mytics 2
set ticscale 3
#
set key spacing 1.5
#
set yzeroaxis lw 4
#
set out "eps/comp.eps"
p "earth_pressure_depth_05w.dat" u ($2*0.001):1 t "T_{ini}" w lp ls 1, \
    "" u ($3*0.001):1 t "M50% T_{act}" w lp ls 2, \
    "" u ($4*0.001):1 t "T_{pas}" w lp ls 3, \
"earth_pressure_depth_025w.dat" u ($3*0.001):1 t "M25% T_{act}" w lp ls 6, \
        "" u ($4*0.001):1 t "T_{pas}" w lp ls 7, \
"igarashi_result.dat" u 4:(-$1) t "(ref.) T_{act}" w lp ls 4, \
    "" u 3:(-$1) t "T_{pas}" w lp ls 5
