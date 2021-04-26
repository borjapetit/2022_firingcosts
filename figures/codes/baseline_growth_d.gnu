
################################################################################

reset

set style line 1 lc rgb '#006400' lt 1 lw 5
set style line 2 lc rgb '#4ca34c' lt 1 lw 5
set style line 3 lc rgb '#71c471' lt 1 lw 5
set style line 4 lc rgb '#97d897' lt 1 lw 5
set style line 5 lc rgb '#cceacc' lt 1 lw 5

md = "results/baseline_growth_d_mean.txt"
sd = "results/baseline_growth_d_sd.txt"
rd = "results/baseline_growth_d_rho.txt"

unset key
set format x "%3.1f"
set format y "%3.2f"
set xrange [-2.2:2.2]
set xlabel 'Log productivity'
set xzeroaxis
set yzeroaxis

################################################################################

set terminal pdfcairo size 10,3.5 enhanced color font 'arial,28' lw 2
set output 'figures/baseline_growth_d_mean_sd.pdf'
set multiplot layout 1,2
set title 'Expected growth rate'
plot md using 'd':'dist' smooth csplines ls 3 lw 3 axes x1y2 title 'Distribution of firms', \
     md using 'd':'agd'  smooth csplines ls 1 lw 3 title 'Ave. growth rate'
set title 'Std growth'
plot sd using 'd':'dist' smooth csplines ls 3 lw 3 axes x1y2 , \
     sd using 'd':'agd'  smooth csplines ls 1 lw 3
unset multiplot
unset output

################################################################################

set terminal pdfcairo size 10,3.5 enhanced color font 'arial,28' lw 2
set output 'figures/baseline_growth_d_diff_innov_noninnov_1.pdf'
set multiplot layout 1,2
set title 'Expected growth rate'
set key top right
plot md using 'd':'ngd' smooth csplines ls 3 lw 3 title 'Non-innovators', \
     md using 'd':'igd' smooth csplines ls 1 lw 3 title 'Innovators'
set title 'Std growth'
unset key
plot sd using 'd':'ngd' smooth csplines ls 3 lw 3, \
     sd using 'd':'igd' smooth csplines ls 1 lw 3
unset multiplot
unset output

################################################################################

set terminal pdfcairo size 10,7 enhanced color font 'arial,28' lw 2
set output 'figures/baseline_growth_d_diff_innov_noninnov_2.pdf'
set multiplot layout 2,2
set title 'Expected growth rate'
set key top right
plot md using 'd':'ngd' smooth csplines ls 3 lw 3 title 'Default', \
     md using 'd':'igd' smooth csplines ls 1 lw 3 title 'Chosen'
set title 'Std. deviation growth'
unset key
plot sd using 'd':'ngd' smooth csplines ls 3 lw 2 , \
     sd using 'd':'igd' smooth csplines ls 1 lw 3
set title 'Prob. of innovation'
plot rd using 'd':'rho0' smooth csplines ls 3 lw 3 , \
     rd using 'd':'rho'  smooth csplines ls 1 lw 3
unset multiplot
unset output

################################################################################

set terminal pdfcairo size 10,3.5 enhanced color font 'arial,28' lw 2
set output 'figures/baseline_growth_d_diff_innov_average.pdf'
set multiplot layout 1,2
set title 'Expected growth rate'
set key top right
plot md using 'd':'agd' smooth csplines ls 3 lw 3 title 'Average', \
     md using 'd':'igd' smooth csplines ls 1 lw 3 title 'Innovators'
set title 'Std growth rate'
unset key
plot sd using 'd':'agd' smooth csplines ls 3 lw 3, \
     sd using 'd':'igd' smooth csplines ls 1 lw 3
unset multiplot
unset output

################################################################################
