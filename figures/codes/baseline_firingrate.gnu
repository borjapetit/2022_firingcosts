
reset

set style line 1 lc rgb '#006400' lt 1 lw 5
set style line 2 lc rgb '#4ca34c' lt 1 lw 5
set style line 3 lc rgb '#71c471' lt 1 lw 5
set style line 4 lc rgb '#97d897' lt 1 lw 5
set style line 5 lc rgb '#cceacc' lt 1 lw 5

datam = "results/graphs/baseline_firingchoice.txt"

set format x "%3.0f"
set format y "%4.0f"
set xlabel '(log) Productivity'

################################################################################
# FIGURE 1

set terminal pdfcairo size 10,3.5 enhanced color font 'arial,21' lw 2
set output 'figures/baseline_firingchoice.pdf'

set size 0.7,1
set lmargin 25
set yrange [1:30]
set xrange [-3:3.5]
set ytics (0,5,10,15,20,25,30)
set title "Firing rate"
set xlabel '(log) Productivity'
set ylabel 'Initial size'
set pm3d map
set palette rgbformulae 22,13,-31
splot datam

unset output
