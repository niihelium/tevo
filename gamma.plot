set term png
#enhanced color 14

set output 'gamma.png'

set tics in

set style function lines
set size 1.0, 1.0
set origin 0.0, 0.0

set autoscale

#set xrange [12:1e4]
#set yrange [1e-29:1e-20]

set log y
#set format x "10^{%2.T}"
set format y "10^{%2.T}"

set xlabel 'z'
#Спектральный поток
set ylabel 'Г, s^{-1}'

#set xtics 0,50,1000
#set ytics 0,50,1000

 #[-5:5] f(x) title "f(x)=cos(x)", integral_f(x)

plot "./gamma_H.out" u 1:2 w l lt 10 lc 1 title "H", \
    "./gamma_HeI.out" u 1:2 w l lt 10 lc 2 title "HeI", \
    "./gamma_HeII.out" u 1:2 w l lt 10 lc 3 title "HeII", \
    #"./data/hm_gal.dat" u 1:31 w l lt 10 lc 4 title "z = 3", \
    #"./data/hm_gal.dat" u 1:36 w l lt 10 lc 5 title "z = 5", \
    # "./data/hm_gal.dat" u 1:43 w l lt 10 lc 6 title "z = 6"
