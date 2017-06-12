set term png
#enhanced color 14

set output 'concentration.png'

set tics in

set style function lines
set size 1.0, 1.0
set origin 0.0, 0.0

set autoscale

#set xrange [1e3:1e7]
#set yrange [1e-35:1e-18]

set log y
#set format x "10^{%2.T}"
set format y "10^{%2.T}"

set xlabel 'Z'
#Спектральный поток
set ylabel 'n, cm^{-3} '

#set xtics 0,50,1000
#set ytics 0,50,1000

 #[-5:5] f(x) title "f(x)=cos(x)", integral_f(x)

#set key bottom right

plot "./concentration.out" u 1:2 w l lt 10 lc 1  title "Full concentration"
