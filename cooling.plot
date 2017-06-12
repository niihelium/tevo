set term png
#enhanced color 14

set output 'cooling.png'

set tics in

set style function lines
set size 1.0, 1.0
set origin 0.0, 0.0

set autoscale

set xrange [1e3:1e7]
set yrange [1e-35:1e-18]

set log xy
set format x "10^{%2.T}"
set format y "10^{%2.T}"

set xlabel 'T, K'
#Спектральный поток
set ylabel 'k, cm^{3} s^{-1} '

#set xtics 0,50,1000
#set ytics 0,50,1000

 #[-5:5] f(x) title "f(x)=cos(x)", integral_f(x)

#set key bottom right

plot "./cooling.out" u 1:2 w l lt 10 lc 1  title "H_excitation", \
  "./cooling.out" u 1:3 w l lt 10 lc 2 title "He_excitation_11s", \
  "./cooling.out" u 1:4 w l lt 10 lc 3 title "He_excitation_23s", \
  "./cooling.out" u 1:5 w l lt 10 lc 4 title "HeII_excitation", \
  "./cooling.out" u 1:6 w l lt 10 lc 5 title "H_collisional_ionization", \
  "./cooling.out" u 1:7 w l lt 10 lc 6 title "He_collisional_ionization", \
  "./cooling.out" u 1:8 w l lt 10 lc 7 title "HeII_collisional_ionization", \
  "./cooling.out" u 1:9 w l lt 10 lc 8 title "HII_recombination", \
  "./cooling.out" u 1:10 w l lt 10 lc 9 title "HeII_recombination", \
  "./cooling.out" u 1:11 w l lt 10 lc 10 title "HeII_recombination_dielectronic", \
  "./cooling.out" u 1:12 w l lt 10 lc 11 title "HeIII_recombination", \
  "./cooling.out" u 1:13 w l lt 10 lc 12 title "compton_cooling", \
  "./cooling.out" u 1:14 w l lt 10 lc 13  lw 5 title "bremsstrahlung"
