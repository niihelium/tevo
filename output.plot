set term png
#enhanced color 14

set output 'output.png'

set tics in

set style function lines
set size 1.0, 1.0
set origin 0.0, 0.0

set autoscale

#set xrange [1e3:1e6]
#set yrange [1e-20:1e-10]

set log y
#set format x "10^{%2.T}"
set format y "10^{%2.T}"

set xlabel 'z'
#Спектральный поток
set ylabel 'k, cm^{3} s^{-1} '

#set xtics 0,50,1000
#set ytics 0,50,1000

 #[-5:5] f(x) title "f(x)=cos(x)", integral_f(x)

set key bottom right

plot "./output.out" u 1:10 w l lt 10 lc 1 title "total", \
    "./output.out" u 1:2 w l lt 10 lc 2 title "HI", \
    "./output.out" u 1:3 w l lt 10 lc 3 title "HII", \
    "./output.out" u 1:($9/1e8) w l lt 10 lc 4 title "T", \
    "./output.out" u 1:5 w l lt 10 lc 7 title "HeI", \
    "./output.out" u 1:6 w l lt 10 lc 8 title "HeII", \
    "./output.out" u 1:7 w l lt 10 lc 9 title "HeIII", \
    "./output.out" u 1:4 w l lt 10 lc 5 title "H", \
    "./output.out" u 1:8 w l lt 10 lc 6 title "He"
