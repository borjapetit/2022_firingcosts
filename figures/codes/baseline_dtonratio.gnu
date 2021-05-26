
reset

set style line 1 lc rgb '#006400' lt 1 lw 5
set style line 2 lc rgb '#4ca34c' lt 1 lw 5
set style line 3 lc rgb '#71c471' lt 1 lw 5
set style line 4 lc rgb '#97d897' lt 1 lw 5
set style line 5 lc rgb '#cceacc' lt 1 lw 5

m  = "results/graphs/baseline_dtonratio.txt"

set format x "%3.0f"
set format y "%4.0f"
set xlabel '(log) Productivity'

################################################################################
# FIGURE 1

set terminal pdfcairo size 10,3.5 enhanced color font 'arial,21' lw 2
set output 'figures/baseline_dtonratio.pdf'

set size 0.7,1
set lmargin 30
set yrange [0:6]
set xrange [-3.5:6]
unset key
set ytics (0,1,2,3,4,5,6)
set title "Productivity-to-labor ratio"
plot m using 'prod':'ratio' w l ls 1 lw 3

unset output
