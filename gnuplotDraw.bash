#!/bin/bash

#############################
#   FUNCTIONS
#############################

# Check if a directory is available. If it's a file then we have to delete it before and if it doesn't exist we create it
checkDirAvailability(){
    if ! [ -d "$1" ] 
    then
        if [ -e "$1" ]
        then
            echo "ERROR: you have to delete the file '$1' before running this script"
            exit 1
        else
            mkdir "$1"
        fi
    fi
}

#############################
#   MAIN PART OF THIS SCRIPT
#############################

checkDirAvailability "graphs"

# We initialize the needed variables
gnuplot_instructionsPhi=""
gnuplot_instructionsPhi2=""
potentialType="$(tail -n 2 data/phi.dat | head -n 1)"
energy="$(tail -n 1 data/phi.dat)"
line_style="with lines"

# We get the data needed to create the graph
gnuplot_instructionsPhi="'data/phi.dat' u 1:2 title '{/Symbol f}(x)' $line_style lc rgb 'blue', 'data/potential.dat' u 1:2 title 'V(x)' $line_style lw 4 lc rgb 'black'"
gnuplot_instructionsPhi2="'data/phi_square.dat' u 1:2 title '|{/Symbol f}(x)|²' $line_style lc rgb 'blue', 'data/potential.dat' u 1:2 title 'V(x)' $line_style lw 4 lc rgb 'black'"

# We create the graph
gnuplot -e "reset; set terminal jpeg size 1600, 900; set grid; set ylabel ''; set xlabel 'x (nm)'; set title '{/Symbol f}(x) for a $potentialType and E = $energy eV'; set output 'graphs/schrodingerGraph.jpg'; plot $gnuplot_instructionsPhi"
gnuplot -e "reset; set terminal jpeg size 1600, 900; set grid; set ylabel ''; set xlabel 'x (nm)'; set title '|{/Symbol f}(x)|² for a $potentialType and E = $energy eV'; set output 'graphs/schrodingerGraphProbability.jpg'; plot $gnuplot_instructionsPhi2"
