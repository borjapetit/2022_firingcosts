

################################################################################
# MAIN SETTING
################################################################################

reset

# Load the dataset
m  = "txtfiles/graphs/firmsize.txt"

# Set terminal: pdf file, with arial font of size 21
set terminal pdfcairo size 10,3.5 enhanced color font 'arial, 21' lw 2

# Name resulting graph
set output "figures/evidence_firmsize_1.pdf"

# Set multiplot environment: 1 row, 2 columns
set multiplot layout 1,2

# Set line styles
set style line 1 lc rgb '#006400' lt 1 lw 3
set style line 2 lc rgb '#4ca34c' lt 1 lw 3
set style line 3 lc rgb '#71c471' lt 1 lw 3
set style line 4 lc rgb '#97d897' lt 1 lw 3
set style line 5 lc rgb '#cceacc' lt 1 lw 3

# Set some general options
set format x "%5.0f"    # Format of y axis ticks: #.#
set format y "%4.2f"    # Format of x axis ticks: #.#

unset key

set fit quiet
f(x) =  q*x*x + p*x + a
fit f(x) m using 1:2 via q,p,a
g(x) = s*x*x + n*x + b
fit g(x) m using 1:3 via s,n,b

set xrange [0:30]

################################################################################
# GRAPH 1
################################################################################

set title 'Average growth rates'
set xlabel 'Number of workers'
unset ylabel
plot m using 1:2 w p ls 4 lw 3 , f(x) w l ls 1 lw 3

################################################################################
# GRAPH 2
################################################################################

set fit quiet

set title 'Standard deviation of growth rates'
set xlabel 'Number of workers'
plot m using 1:3 w p ls 4 lw 3 , g(x) w l ls 1 lw 3

################################################################################

unset multiplot
unset output


################################################################################
# MAIN SETTING
################################################################################

reset

# Load the dataset
m  = "txtfiles/firmsize.txt"

# Set terminal: pdf file, with arial font of size 21
set terminal pdfcairo size 10,3.5 enhanced color font 'arial, 24' lw 2

# Name resulting graph
set output "figures/evidence_firmsize_2.pdf"

# Set multiplot environment: 1 row, 2 columns
set multiplot layout 1,2

# Set line styles
set style line 1 lc rgb '#006400' lt 1 lw 3
set style line 2 lc rgb '#4ca34c' lt 1 lw 3
set style line 3 lc rgb '#71c471' lt 1 lw 3
set style line 4 lc rgb '#97d897' lt 1 lw 3
set style line 5 lc rgb '#cceacc' lt 1 lw 3

# Set some general options
set format x "%5.0f"    # Format of y axis ticks: #.#
set format y "%4.2f"    # Format of x axis ticks: #.#

unset key

set xrange [0:60]

set fit quiet
f(x) =  q*x*x + p*x + a
fit f(x) m using 1:2 via q,p,a
g(x) = s*x*x + n*x + b
fit g(x) m using 1:3 via s,n,b

set xrange [0:60]

################################################################################
# GRAPH 1
################################################################################

set title 'Average growth rates'
set xlabel 'Number of workers'
unset ylabel
plot m using 1:2 w p ls 4 lw 3 , f(x) w l ls 1 lw 3 , m using 1:2 smooth sbezier ls 3 lw 3

################################################################################
# GRAPH 2
################################################################################

set fit quiet

set title 'Standard deviation of growth rates'
set xlabel 'Number of workers'
plot m using 1:3 w p ls 4 lw 3 , g(x) w l ls 1 lw 3 , m using 1:3 smooth sbezier ls 3 lw 3

################################################################################

unset multiplot
unset output
