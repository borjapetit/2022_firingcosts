

################################################################################
# MAIN SETTING
################################################################################

reset

# Load the dataset
m  = "txtfiles/graphs/sector.txt"

# Set terminal: pdf file, with arial font of size 21
set terminal pdfcairo size 10,3.5 enhanced color font 'arial, 21' lw 2

# Name resulting graph
set output "figures/evidence_sector.pdf"

# Set multiplot environment: 1 row, 2 columns
set multiplot layout 1,2

# Set line styles
set style line 1 lc rgb '#006400' lt 1 lw 3
set style line 2 lc rgb '#4ca34c' lt 1 lw 3
set style line 3 lc rgb '#71c471' lt 1 lw 3
set style line 4 lc rgb '#97d897' lt 1 lw 3
set style line 5 lc rgb '#cceacc' lt 1 lw 3

# Set some general options
set format x "%5.2f"    # Format of y axis ticks: #.#
set format y "%3.1f"    # Format of x axis ticks: #.#

unset key

set fit quiet
f(x) = p*x + a
fit f(x) m using 3:4 via p,a
g(x) = n*x + b
fit g(x) m using 5:6 via n,b

################################################################################
# GRAPH 1
################################################################################

set title 'Employment growth'
set xlabel 'Average'
set ylabel 'Standard deviation'
plot m using 3:4 w p ls 4 lw 3 , f(x) w l ls 1 lw 3

################################################################################
# GRAPH 2
################################################################################

set format x "%5.1f"    # Format of y axis ticks: #.#

set fit quiet
set xzeroaxis
set yzeroaxis

set title 'Revenue growth'
set xlabel 'Average'
set ylabel 'Standard deviation'
plot m using 5:6 w p ls 4 lw 3 , g(x) w l ls 1 lw 3

################################################################################

unset multiplot
unset output
