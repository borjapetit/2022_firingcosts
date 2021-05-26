
reset

set style line 1 lc rgb '#006400' lt 1 lw 5
set style line 2 lc rgb '#4ca34c' lt 1 lw 5
set style line 3 lc rgb '#71c471' lt 1 lw 5
set style line 4 lc rgb '#97d897' lt 1 lw 5
set style line 5 lc rgb '#cceacc' lt 1 lw 5

m  = "results/graphs/baseline_optimalpi.txt"

set format x  "%3.1f"
set format y  "%4.2f"
set format y2 "%4.2f"
set xlabel '(log) Productivity'

################################################################################
# FIGURE 1

set terminal pdfcairo size 10,3.5 enhanced color font 'arial,21' lw 2
set output 'figures/baseline_optimalpi.pdf'

set size 0.9,1
set lmargin 18
set xrange [-3:2]
set yrange [0:0.25]
set ytics (0.00,0.05,0.10,0.15,0.20,0.25)
set y2range [0:0.10]
set y2tics (0.00,0.02,0.04,0.06,0.08,0.10)
#set title "Choice of next period's productivity distribution"
plot m using 'prod':'eta'   w filledcurves ls 5 lw 3 title 'Benchamark dist.' , \
     m using 'prod':'value' w l ls 3 lw 3 title 'Value function' axes x1y2 , \
     m using 'prod':'pi'    w l ls 1 lw 3 title 'Chosen dist.'
unset output
