set terminal pngcairo size 400,600
set output 'bulk_modulus_vs_temperature.png'
set title "Bulk Modulus vs Temperature"
set xlabel "Temperature"
set ylabel "Bulk Modulus"
set xrange [0:300]

plot 'bulk_modulus-temperature.dat' using 1:2 with lines title 'Bulk Modulus'

