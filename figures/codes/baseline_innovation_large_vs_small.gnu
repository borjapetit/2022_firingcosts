
reset

set style line 1 lc rgb '#006400' lt 1 lw 5
set style line 2 lc rgb '#4ca34c' lt 1 lw 5
set style line 3 lc rgb '#71c471' lt 1 lw 5
set style line 4 lc rgb '#97d897' lt 1 lw 5
set style line 5 lc rgb '#cceacc' lt 1 lw 5

m = "results/graphs/baseline_innovation_large_vs_small.txt"

set xzeroaxis
set yzeroaxis

set format x "%3.1f"
set format y "%4.2f"
set xlabel 'Productivity growth'

################################################################################

set terminal pdfcairo size 10,3.5 enhanced color font 'arial,21' lw 2
set output 'figures/baseline_innovation_large_vs_small.pdf'
set multiplot layout 1,2
set title 'Low productivity firm (p25)'
set xrange [-2:2]
plot m using 'dlow':'td0low' w l ls 4 lw 3 title ' Default ' , \
     m using 'dlow':'td1low' w l ls 1 lw 3 title ' Chosen  '
unset key
set title 'High productivity firm (p75)'
set xrange [-2:2]
plot m using 'dhigh':'td0high' w l ls 4 lw 3 title ' ' , \
     m using 'dhigh':'td1high' w l ls 1 lw 3 title ' '
unset multiplot
unset output

################################################################################
