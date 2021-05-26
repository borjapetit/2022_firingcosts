
reset

set style line 1 lc rgb '#006400' lt 1 lw 5
set style line 2 lc rgb '#4ca34c' lt 1 lw 5
set style line 3 lc rgb '#71c471' lt 1 lw 5
set style line 4 lc rgb '#97d897' lt 1 lw 5
set style line 5 lc rgb '#cceacc' lt 1 lw 5

md = "results/graphs/baseline_distortion.txt"

set terminal pdfcairo size 10,3.5 enhanced color font 'arial,21' lw 2
set output 'figures/baseline_distortion.pdf'
set multiplot layout 1,2

set format x "%2.1f"
set format y "%2.1f"
set format y2 "%2.1f"
set xlabel 'Log productivity'


set key top right
set xrange [-1:2]
set xtics (1.0,-0.5,0.0,0.5,1.0,1.5,2.0)
set yrange [0:0.8]
set ytics nomirror
set ytics (0.0,0.2,0.4,0.6,0.8)
set y2range [0:0.4]
set y2tics (0.0,0.1,0.2,0.3,0.4)
set title 'Firm with (d,n) = (0.5,10)'
plot md using ($1):($4) w l ls 1 lw 3 title 'Firing rate (left)' , \
     md using ($1):($5) w l ls 3 lw 3 title 'Prob. (right)' axes x1y2
unset key

set key top right
set title 'Firm with (d,n) = (-0.5,3)'
set xrange [-2:1]
set xtics (-2.0,-1.5,-1.0,-0.5,0.0,0.5,1.0)
set yrange [0:0.8]
set ytics nomirror
set ytics (0.0,0.2,0.4,0.6,0.8)
set y2range [0:0.4]
set y2tics (0.0,0.1,0.2,0.3,0.4)
plot md using ($6):($9)  w l ls 1 lw 3 title 'Firing rate' , \
     md using ($6):($10) w l ls 3 lw 3 title 'Prob.' axes x1y2
unset key

unset multiplot
unset output
