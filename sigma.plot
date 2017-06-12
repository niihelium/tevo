set term png
#enhanced color 14

set output 'sigma.png'

set tics in

set style function lines
set size 1.0, 1.0
set origin 0.0, 0.0

set autoscale

#set xrange [12:1e4]
#set yrange [1e-29:1e-20]

set log xy
set format x "10^{%2.T}"
set format y "10^{%2.T}"

set xlabel 'E, eV'
#Спектральный поток
set ylabel 'J, erg s^{-1} cm^{-2} Hz^{-1}'

#set xtics 0,50,1000
#set ytics 0,50,1000

 #[-5:5] f(x) title "f(x)=cos(x)", integral_f(x)

plot "./sigma.out" u 1:2 w l lt 10 lc 1
    #"./data/hm_gal.dat" u 1:16 w l lt 10 lc 2 title "z = 1", \
    #"./data/hm_gal.dat" u 1:25 w l lt 10 lc 3 title "z = 2", \
    #"./data/hm_gal.dat" u 1:31 w l lt 10 lc 4 title "z = 3", \
    #"./data/hm_gal.dat" u 1:36 w l lt 10 lc 5 title "z = 5", \
    # "./data/hm_gal.dat" u 1:43 w l lt 10 lc 6 title "z = 6"
