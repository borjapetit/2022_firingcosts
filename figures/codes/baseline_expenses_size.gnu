
reset

set style line 1 lc rgb '#006400' lt 1 lw 5
set style line 2 lc rgb '#4ca34c' lt 1 lw 5
set style line 3 lc rgb '#71c471' lt 1 lw 5
set style line 4 lc rgb '#97d897' lt 1 lw 5
set style line 5 lc rgb '#cceacc' lt 1 lw 5

ma  = "results/baseline_expenses_size.txt"

set format x "%3.0f"
set format y "%4.0f"
set xlabel 'Productivity growth'

################################################################################
# FIGURE 1

set terminal pdfcairo size 10,3.5 enhanced color font 'arial,28' lw 2
set output 'figures/baseline_expenses_size.pdf'
set multiplot layout 1,2

set xzeroaxis
unset key
set xlabel 'firm size'
set xrange [0:60]
set ytics (0,10,20,30,40,50,60)

set yrange [0:5]
set ytics (0,1,2,3,4,5)
set title "Innovation expenses"
plot ma using 'Lab':'Exp' w l ls 1 lw 3

set yrange [0:65]
set ytics (0,10,20,30,40,50,60)
set title "Innovation expenses as a share of firm's output"
plot ma using 'Lab':'Share' w l ls 1 lw 3

unset multiplot
unset output

################################################################################
