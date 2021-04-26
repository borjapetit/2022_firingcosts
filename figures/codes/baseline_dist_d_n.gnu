
reset

set style line 1 lc rgb '#006400' lt 1 lw 5
set style line 2 lc rgb '#4ca34c' lt 1 lw 5
set style line 3 lc rgb '#71c471' lt 1 lw 5
set style line 4 lc rgb '#97d897' lt 1 lw 5
set style line 5 lc rgb '#cceacc' lt 1 lw 5

md = "results/baseline_dist_d.txt"
mn = "results/baseline_dist_n.txt"

set terminal pdfcairo size 10,3.5 enhanced color font 'arial,28' lw 2
set output 'figures/baseline_dist_d_n.pdf'
set multiplot layout 1,2

set key top right
set title 'Productivity distribution'
set xlabel 'Log productivity'
set format x "%4.1f"
set format y "%3.2f"
set yrange[-0.01:0.08]
plot md using 'd':'Dd' smooth csplines ls 1 lw 3 title 'Population' , \
     md using 'd':'D0' smooth csplines ls 3 lw 3 title 'Entrants'
unset key

set title 'Size distribution'
set xlabel 'Num. employees'
set format x "%3.0f"
set format y "%3.2f"
set yrange[-0.01:0.15]
plot mn using  'n':'Dn' smooth acsplines ls 1 lw 3 title ' '

unset multiplot
unset output
