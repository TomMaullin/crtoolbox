#!/bin/bash
fslpython -c "from MSmp.join_and_plot_Mmu import *; joinAndPlot('$1',$2)"