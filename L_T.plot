set term png
#enhanced color 14

set output 'l-t.png'

set tics in

set style function lines
set size 1.0, 1.0
set origin 0.0, 0.0

set autoscale

set xrange [1e4:1e8]
set yrange [1e-25:1e-20]

set log xy
set format x "10^{%2.T}"
set format y "10^{%2.T}"

set xlabel 'T, K'
#Спектральный поток
set ylabel 'k, cm^{3} s^{-1} '

#set xtics 0,50,1000
#set ytics 0,50,1000

 #[-5:5] f(x) title "f(x)=cos(x)", integral_f(x)

set key bottom right

plot "./output.out" u 9:10 w l lt 10 lc 1 title "L"
