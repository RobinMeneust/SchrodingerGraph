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
gnuplot_instructions=""
gnuplot_instructionsALL=""
potentialType="$(tail -n 2 data/phi.dat | head -n 1)"
energy="$(tail -n 1 data/phi.dat)"
line_style="with lines"
border_style="set xtics nomirror; set ytics nomirror; set border 3 lw 1; set arrow from graph 0,0 to graph 0,1 filled lw 4; set arrow from graph 1,0 to graph 1,1 filled lw 4"

# We get the data needed to create the graph
gnuplot_instructions="'data/phi.dat' u 1:2 title 'phi(x)' $line_style lc rgb 'blue', 'data/potential.dat' u 1:2 title 'V(x) (volts)' $line_style lw 4 lc rgb 'black'"

# We create the graph
gnuplot -e "reset; set terminal jpeg size 1600, 900; $border_style; $label_style; set ylabel ''; set xlabel 'x (nm)'; set title 'phi(x) for a $potentialType and E = $energy '; set output 'graphs/schrodingerGraph.jpg'; plot $gnuplot_instructions"