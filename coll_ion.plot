set term png
#enhanced color 14

set output 'Collisional.png'

set tics in

set style function lines
set size 1.0, 1.0
set origin 0.0, 0.0

set autoscale

set xrange [1e3:1e6]
set yrange [1e-20:1e-7]

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

plot "./coll_H.out" u 1:2 w l lt 10 lc 1 title "H", \
    "./coll_HeI.out" u 1:2 w l lt 10 lc 2 title "HeI", \
    "./coll_HeII.out" u 1:2 w l lt 10 lc 3 title "HeII"
    #"./data/hm_gal.dat" u 1:31 w l lt 10 lc 4 title "z = 3", \
    #"./data/hm_gal.dat" u 1:36 w l lt 10 lc 5 title "z = 5", \
    # "./data/hm_gal.dat" u 1:43 w l lt 10 lc 6 title "z = 6"
