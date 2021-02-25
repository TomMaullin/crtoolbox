#!/bin/bash
fslpython -c "from TwoSmp.join_and_plot_2mu import *; joinAndPlot('$1',$2)"