
reset

set style line 1 lc rgb '#006400' lt 1 lw 5
set style line 2 lc rgb '#4ca34c' lt 1 lw 5
set style line 3 lc rgb '#71c471' lt 1 lw 5
set style line 4 lc rgb '#97d897' lt 1 lw 5
set style line 5 lc rgb '#cceacc' lt 1 lw 5

pe    = "results/experiment_3_pe.txt"
ge    = "results/experiment_3_ge.txt"
geng  = "results/experiment_3_ge_ni.txt"

set xtics (0.0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4003)

set format x "%4.2f" ; set xzeroaxis ; set xrange [0.0:0.4003]
set format y "%3.1f" ; set yzeroaxis ;

################################################################################
# FIGURE 1, % CHANGE IN AGGREGATE VARIABLES

set terminal pdfcairo size 10,3.5 enhanced color font 'arial,28' lw 2
set output 'figures/experiment_change_1.pdf'
set multiplot layout 1,2

set key left bottom
set xlabel 'Firing cost parameter'
set title 'Aggregate productivity (% w.r.t. frictionless)'
plot pe using 'fc':'tfp' w l ls 4 lw 3 title 'Partial Equil.' , \
     ge using 'fc':'tfp' w l ls 1 lw 3 title 'General Equil.' ; unset key
set title 'Average firm productivity (% w.r.t. frictionless)'
plot pe using 'fc':'d' w l ls 4 lw 3 title 'Partial Equil.' , \
     ge using 'fc':'d' w l ls 1 lw 3 title 'General Equil.'

unset multiplot
unset output

################################################################################
# FIGURE 3, % CHANGE IN AGGREGATE VARIABLES

set terminal pdfcairo size 10,3.5 enhanced color font 'arial,28' lw 2
set output 'figures/experiment_change_2.pdf'
set multiplot layout 1,2

set key left bottom
set xlabel "Firing cost parameter"
set title 'Aggregate productivity (% w.r.t. frictionless)'
plot geng using 'fc':'tfp' w l ls 4 lw 3 title 'Without innovation' , \
     ge   using 'fc':'tfp' w l ls 1 lw 3 title 'With innovation' ; unset key
set title 'Average firm productivity  (% w.r.t. frictionless)'
plot geng using 'fc':'d' w l ls 4 lw 3 title 'Without innovation' , \
     ge   using 'fc':'d' w l ls 1 lw 3 title ' innovation'

unset multiplot
unset output
