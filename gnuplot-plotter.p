reset
set title 'Heat approximation - Euler implicit scheme'
set xlabel 'P'
set ylabel 'Values'
set xrange [1:20]
set yrange [0:4]
set size 1.0, 1.0
set term gif animate delay 25
set output "my-plot.gif"
do for [i=0:19] {
  plot 'output.dat' every ::(i*20)::((i+1)*20-1) using 1:3 with lines title sprintf("Exact (t = %d)", i), \
  'output.dat' every ::(i*20)::((i+1)*20-1) using 1:2 with lines title sprintf("Scheme")
}
set size 1,1
