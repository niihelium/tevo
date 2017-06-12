set term png
#enhanced color 14

set output 'lambda.png'

set tics in

set style function lines
set size 1.0, 1.0
set origin 0.0, 0.0

set autoscale

set xrange [1e4:1e8]
set yrange [1e-25:1e-22]

set log xy
set format x "10^{%2.T}"
set format y "10^{%2.T}"

set xlabel 'T, K'
#Спектральный поток
set ylabel '\Lambda, erg cm^{3} s^{-1}  '

#set xtics 0,50,1000
#set ytics 0,50,1000

 #[-5:5] f(x) title "f(x)=cos(x)", integral_f(x)

#set key bottom right

plot "./lambda.out" u 1:2 w l lt 10 lc 1  title "H excitation", \
  "./lambda.out" u 1:3 w l lt 10 lc 2 title "He excitation 11s", \
  "./lambda.out" u 1:4 w l lt 10 lc 3 title "He excitation 23s", \
  "./lambda.out" u 1:5 w l lt 10 lc 4 title "HeII excitation", \
  "./lambda.out" u 1:6 w l lt 10 lc 5 title "H collisional ionization", \
  "./lambda.out" u 1:7 w l lt 10 lc 6 title "He collisional ionization", \
  "./lambda.out" u 1:8 w l lt 10 lc 7 title "HeII collisional ionization", \
  "./lambda.out" u 1:9 w l lt 10 lc 8 title "HII recombination", \
  "./lambda.out" u 1:10 w l lt 10 lc 9 title "HeII recombination", \
  "./lambda.out" u 1:11 w l lt 10 lc 10 title "HeII recombination dielectronic", \
  "./lambda.out" u 1:12 w l lt 10 lc 11 title "HeIII recombination", \
  "./lambda.out" u 1:13 w l lt 10 lc 12 title "compton cooling", \
  "./lambda.out" u 1:14 w l lt 10 lc 13 lw 5 title "bremsstrahlung", \
  "./lambda.out" u 1:15 w l lt 10 lc 14 title "total"
