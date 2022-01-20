#!/usr/bin/env python
# coding: utf-8


#!/usr/bin/env python
import os, sys, shutil
from numpy import *
from pylab import *
from scipy.interpolate import interp1d
import ROOT
from plotlibrary import *
from plotting_helpers import *

if len(sys.argv) < 2:
    print("Usage: python /path/plotsinglemodel.py /path/analysis_output.root /path/destination/")
    print("Optional arguments: histogram_number(integer) groomed_plots(boolean) scale_au(boolean) var_legend(boolean)")
    sys.exit(1)
    
infile = sys.argv[1]
destination = sys.argv[2]
if len(sys.argv) > 3:
    histnumber = int(sys.argv[3])
else:
    histnumber = -1

# Check whether the specified path exists or not, create it if it does not exist
isExist = os.path.exists(destination)
if not isExist:
    os.makedirs(destination)

f0 = ROOT.TFile(infile,"READ")
keys = f0.GetListOfKeys()
    
if histnumber >= 0 and histnumber <= len(keys):
    histoname = GetAllHistoNames(infile)[0][histnumber]
    if 'pid' in histoname:
        print('Will do special PID plot')
        MakePIDPlot(infile, destination, histnumber)
    else:
        MakeSinglePlot(infile, destination, histnumber, groomed_plots = False, scale_au = True, var_legend = False)
        
else:
    for i in range(0,len(keys)):
        histoname = GetAllHistoNames(infile)[0][i]
        if 'pid' in histoname:
            print('Will do special PID plot')
            MakePIDPlot(infile, destination, histnumber = i)
        else:
            MakeSinglePlot(infile, destination, histnumber = i, groomed_plots = False, scale_au = True, var_legend = False)