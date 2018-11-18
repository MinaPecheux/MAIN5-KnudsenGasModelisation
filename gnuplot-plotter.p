reset
P = 100
T = 20

set title 'Heat approximation - Euler implicit scheme'
set xlabel 'P'
set ylabel 'Values'
set xrange [0:P-1]
set yrange [0:4]
set size 1.0, 1.0
set term gif animate delay 25
set output "my-plot.gif"
do for [i=0:T-1] {
  plot 'output.dat' every ::(i*P)::((i+1)*P-1) using 1:3 with lines title sprintf("Exact (t = %d)", i), \
  'output.dat' every ::(i*P)::((i+1)*P-1) using 1:2 with lines title sprintf("Scheme")
}
set size 1,1
