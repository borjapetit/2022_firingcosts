
reset

set style line 1 lc rgb '#006400' lt 1 lw 5
set style line 2 lc rgb '#4ca34c' lt 1 lw 5
set style line 3 lc rgb '#71c471' lt 1 lw 5
set style line 4 lc rgb '#97d897' lt 1 lw 5
set style line 5 lc rgb '#cceacc' lt 1 lw 5

ma  = "results/graphs/experiment_innovation_d_fc0.txt"
mb  = "results/graphs/experiment_innovation_d_2fc0.txt"
mc  = "results/graphs/experiment_innovation_d_one.txt"
na  = "results/graphs/experiment_innovation_n_fc0.txt"
nb  = "results/graphs/experiment_innovation_n_2fc0.txt"
nc  = "results/graphs/experiment_innovation_n_one.txt"

set format x "%3.1f"
set format y "%4.1f"
set xlabel 'Productivity growth'

################################################################################
# FIGURE 1

set terminal pdfcairo size 10,3.5 enhanced color font 'arial,21' lw 2
set output 'figures/experiment_innovation_d_02.pdf'
set multiplot layout 1,2

set xzeroaxis
unset key
set xlabel 'log productivity'
set xrange [-2:2]

set yrange [-4:0]
set ytics (0.0,-0.5,-1,-1.5,-2.0,-2.5,-3.0)
set title 'Differential expected growth rate (p.p.)'
plot ma using 'd':(100*$4) w l ls 1 lw 3

set yrange [-4:1]
set ytics (0.0,-0.5,-1,-1.5,-2.0,-2.5,-3.0)
set title 'Change in sd of innovation (%)'
plot ma using 'd':(100*$7/$8) w l ls 1 lw 3

#($14-$13)

unset multiplot
unset output

################################################################################
# FIGURE 2

set terminal pdfcairo size 5,3.5 enhanced color font 'arial,21' lw 2
set output 'figures/experiment_innovation_d_02_ratio.pdf'

set xzeroaxis
unset key
set xlabel 'log productivity'
set xrange [-2:2]

set yrange [-1:2]
set ytics (2.0,1.5,1.0,0.5,0.0,-0.5,-1)
set title 'Change in growth relative to change in volatility'
#plot ma using 'd':(100*$4)/(100*$7/$8) w l ls 1 lw 3
plot ma using 'd':($5/$7)/($6/$8) w l ls 1 lw 3

unset output

################################################################################
# FIGURE 3

set terminal pdfcairo size 10,3.5 enhanced color font 'arial,28' lw 2
set output 'figures/experiment_innovation_d_04.pdf'
set multiplot layout 1,2

set xzeroaxis
unset key
set xlabel 'log productivity'
set xrange [-2:2]

set yrange [-6:0]
set ytics (0,-1,-2,-3,-4,-5,-6)
set title 'Differential expected growth rate (p.p.)'
plot mb using 'd':(100*$4) w l ls 1 lw 3

set yrange [-4:0.0]
set ytics (0.0,-0.5,-1,-1.5,-2.0,-2.5,-3.0,-3.5,-4.0)
set title 'Change in sd of innovation (%)'
plot mb using 'd':(100*$7/$8) w l ls 1 lw 3

unset multiplot
unset output

################################################################################
# FIGURE 3

set terminal pdfcairo size 10,3.5 enhanced color font 'arial,28' lw 2
set output 'figures/experiment_innovation_d_10.pdf'
set multiplot layout 1,2

set format x "%3.1f"
set format y "%3.0f"

set xzeroaxis
unset key
set xlabel 'log productivity'
set xrange [-2:2]

set yrange [-10:0]
set ytics (0,-2,-4,-6,-8,-10)
set title 'Differential expected growth rate (p.p.)'
plot mc using 'd':(100*$4)  w l ls 1 lw 3

set yrange [-7:0.0]
set ytics (0,-1,-2,-3,-4,-5,-6,-7)
set title 'Change in sd of innovation (%)'
plot mc using 'd':(100*$7/$8) w l ls 1 lw 3

unset multiplot
unset output

################################################################################
# FIGURE 4

set terminal pdfcairo size 10,3.5 enhanced color font 'arial,28' lw 2
set output 'figures/experiment_innovation_d_all.pdf'
set multiplot layout 1,2

set format x "%3.1f"
set format y "%3.0f"

set xzeroaxis
set key bottom right
set xlabel 'log productivity'
set xrange [-2:2]

set yrange [-10:0]
set ytics (0,-2,-4,-6,-8,-10)
set title 'Differential expected growth rate (p.p.)'
plot ma using 'd':(100*$4) w l ls 1 lw 3 title ' Firing cost = 0.2' , \
     mb using 'd':(100*$4) w l ls 3 lw 3 title ' Firing cost = 0.4' , \
     mc using 'd':(100*$4) w l ls 5 lw 3 title ' Firing cost = 1.0'

unset key
set yrange [-7:0.0]
set ytics (0,-1,-2,-3,-4,-5,-6,-7)
set title 'Change in sd of innovation (%)'
plot ma using 'd':(100*$7/$8) w l ls 1 lw 3 title ' Firing cost = 0.2', \
     mb using 'd':(100*$7/$8) w l ls 3 lw 3 title ' Firing cost = 0.4', \
     mc using 'd':(100*$7/$8) w l ls 5 lw 3 title ' Firing cost = 1.0'

unset multiplot
unset output

################################################################################
# FIGURE 5
#
# set terminal pdfcairo size 10,7 enhanced color font 'arial,28' lw 2
# set output 'figures/experiment_innovation_n_all.pdf'
# set multiplot layout 2,2
#
# set format x "%3.1f"
# set format y "%3.0f"
#
# set xzeroaxis
# set key bottom right
# set xlabel 'firm size'
#
# set boxwidth 10
# set style data histograms
# set style histogram cluster
#
# set yrange [-3:1]
# set xrange [1:6]
# set ytics (1,0,-1,-2,-3,-4)
# set title 'Differential expected growth rate'
# set style fill solid
# plot na using 'n':(100*$4):2 w boxes ls 1 lw 3 title ' Firing cost = 0.2'
#
# set yrange [-3:1]
# set xrange [1:2]
# set ytics (1,0,-1,-2,-3,-4)
# set title 'Differential expected growth rate'
# set style fill solid
# plot nb using 'n':(100*$4):2 w boxes ls 3 lw 3 title ' Firing cost = 0.4'
#
# set yrange [-3:1]
# set xrange [1:2]
# set ytics (1,0,-1,-2,-3,-4)
# set title 'Differential expected growth rate'
# set style fill solid
# plot nc using 'n':(100*$4):2 w boxes ls 5 lw 3 title ' Firing cost = 1.0'
#
# set key bottom left
# set yrange [-3:0]
# set xrange [0:7]
# set ytics (1,0,-1,-2,-3,-4)
# set title 'Differential expected growth rate'
# set style fill solid
# plot na using 1:(100*$4) w boxes ls 1 lw 3 title ' Firing cost = 0.2' , \
#      nb using 1:(100*$4) w boxes ls 3 lw 3 title ' Firing cost = 0.4' , \
#      nc using 1:(100*$4) w boxes ls 5 lw 3 title ' Firing cost = 1.0'
#
# unset key
# unset multiplot
# unset output
#
################################################################################
