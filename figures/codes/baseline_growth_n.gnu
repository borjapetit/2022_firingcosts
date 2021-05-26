
################################################################################

reset

set style line 1 lc rgb '#006400' lt 1 lw 5
set style line 2 lc rgb '#4ca34c' lt 1 lw 5
set style line 3 lc rgb '#71c471' lt 1 lw 5
set style line 4 lc rgb '#97d897' lt 1 lw 5
set style line 5 lc rgb '#cceacc' lt 1 lw 5

md = "results/graphs/baseline_growth_n_mean.txt"
sd = "results/graphs/baseline_growth_n_sd.txt"
rd = "results/graphs/baseline_growth_n_rho.txt"

unset key
set format x "%3.0f"
set format y "%3.2f"
set xlabel 'Emplyoment'
set xzeroaxis

################################################################################

set terminal pdfcairo size 10,3.5 enhanced color font 'arial,21' lw 2
set output 'figures/baseline_growth_n_mean_sd.pdf'
set multiplot layout 1,2
set title 'Expected growth rate'
plot md using 'n':'dist' smooth sbezier ls 3 lw 3 axes x1y2 title 'Distribution of firms', \
     md using 'n':'agd'  smooth sbezier ls 1 lw 3 title 'Ave. growth rate'
set title 'Std growth'
plot sd using 'n':'dist' smooth sbezier ls 3 lw 3 axes x1y2 , \
     sd using 'n':'agd'  smooth sbezier ls 1 lw 3
unset multiplot
unset output

################################################################################

set terminal pdfcairo size 10,3.5 enhanced color font 'arial,21' lw 2
set output 'figures/baseline_growth_n_diff_innov_noninnov_1.pdf'
set multiplot layout 1,2
set title 'Expected growth rate'
set key top right
plot md using 'n':'ngd' smooth sbezier ls 3 lw 3 title 'Non-innovators', \
     md using 'n':'igd' smooth sbezier ls 1 lw 3 title 'Innovators'
set title 'Std growth'
plot sd using 'n':'ngd' smooth sbezier ls 3 lw 3, \
     sd using 'n':'igd' smooth sbezier ls 1 lw 3
unset multiplot
unset output

################################################################################

set terminal pdfcairo size 10,7 enhanced color font 'arial,21' lw 2
set output 'figures/baseline_growth_n_diff_innov_noninnov_2.pdf'
set multiplot layout 2,2
set title 'Expected growth rate'
set xzeroaxis
set key top right
plot md using 'n':'ngd' smooth sbezier ls 3 lw 3 title 'Default', \
     md using 'n':'igd' smooth sbezier ls 1 lw 3 title 'Chosen'
set title 'Std. deviation growth'
unset key
plot sd using 'n':'ngd' smooth sbezier ls 3 lw 2 , \
     sd using 'n':'igd' smooth sbezier ls 1 lw 3
set title 'Prob. of innovation'
plot rd using 'n':'rho0' smooth sbezier ls 3 lw 3 , \
     rd using 'n':'rho'  smooth sbezier ls 1 lw 3
unset multiplot
unset output

################################################################################

set terminal pdfcairo size 10,3.5 enhanced color font 'arial,21' lw 2
set output 'figures/baseline_growth_n_diff_innov_average.pdf'
set multiplot layout 1,2
set title 'Expected growth rate'
set xzeroaxis
set key top right
plot md using 'n':'agd' smooth sbezier ls 3 lw 3 title 'Average', \
     md using 'n':'igd' smooth sbezier ls 1 lw 3 title 'Innovators'
set title 'Std growth rate'
unset xzeroaxis
unset key
plot sd using 'n':'agd' smooth sbezier ls 3 lw 3, \
     sd using 'n':'igd' smooth sbezier ls 1 lw 3
unset multiplot
unset output

################################################################################
