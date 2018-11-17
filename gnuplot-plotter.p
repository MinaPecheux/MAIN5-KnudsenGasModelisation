reset
set title 'Heat approximation - Euler implicit scheme'
set xlabel 'P'
set ylabel 'Values'
set xrange [0:20]
set yrange [0:4]
set size 1.0, 1.0
set term gif animate delay 25
set output "my-plot.gif"
block_len = 20
do for [i=0:19] {
  plot 'output.dat' every ::(i*block_len)::((i+1)*block_len-1) using 1:3 with lines title sprintf("Exact (t = %d)", i), \
  'output.dat' every ::(i*block_len)::((i+1)*block_len-1) using 1:2 with lines title sprintf("Scheme")
}
set size 1,1
