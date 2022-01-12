#!/usr/bin/env python
# coding: utf-8


#!/usr/bin/env python
import os, sys, shutil
from numpy import *
from pylab import *
from scipy.interpolate import interp1d
import ROOT
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


def MakeSinglePlot(infile = '', destination = './', histnumber = int, groomed_plots = False, scale_au = True, var_legend = True):
    rootfile1 = infile
    histoname = GetAllHistoNames(infile)[0][histnumber]
    histotitle = GetAllHistoNames(infile)[1][histnumber]
    plot_title = "%s%s" %(destination,histoname)
    legend_positions = ['r','r','r','r','r','r','l','c','r','r','r','r','r','r','r','r','r','r']

    ROOT.gROOT.SetBatch(True)
    c=ROOT.TCanvas("c","c",800,600)
    f0 = ROOT.TFile(rootfile1,"READ")
    
    # Get the histogram
    h0=f0.Get(histoname)
    n_bins = h0.GetNbinsX()
    divider = 20 #larger divider gives more bins
    while n_bins%divider != 0 or divider%2 != 0:
        divider += 1
    if divider == n_bins:
        divider = divider/2
    h0.Rebin(int(n_bins/divider))
    #print(n_bins, int(n_bins/divider), h0.GetNbinsX())
    # The following lines will normalise the histogram if scale_au=True and groomed_plots=False
    if scale_au == True and groomed_plots == False:
        norm = h0.GetEntries()
        scale = 1/h0.Integral();
        h0.Scale(scale);
    
    c.cd()
    #Set axis title, you can use latex as demonstrated
    h0.GetXaxis()
    # Thie following lines will set the Y-axis title depending on whether it was normalised
    if scale_au == True and groomed_plots == False:
        h0.GetYaxis().SetTitle('A.U.')
    else:
        h0.GetYaxis()
    #switch off the ugly stats box on the upper right corner
    h0.SetStats(0)
    
    #set label sizes
    h0.GetYaxis().SetLabelSize(0.04)
    h0.GetYaxis().SetTitleSize(0.04)
    h0.GetXaxis().SetTitleSize(0.04)
    h0.GetXaxis().SetLabelSize(0.04)
    h0.SetTitle("")
    #set linewidth and color of the lines
    h0.SetLineWidth(2)
    h0.SetLineColor(ROOT.kBlue-2)
    #actually draw the histogram
    h0.Draw('hist')
    h0.Draw('E1 same')
    
    #set axis offsets so that labels don't get messed up
    h0.GetXaxis().SetTitleOffset(1.4)
    h0.GetYaxis().SetTitleOffset(1.4)
    #set axis ranges
    max_x = h0.GetXaxis().GetXmax()
    min_x = min(h0.GetXaxis().GetXmin(),0)
    h0.GetXaxis().SetRangeUser(min_x, max_x)
    max_y = h0.GetBinContent(h0.GetMaximumBin())*1.1
    h0.GetYaxis().SetRangeUser(0, max_y)
    #set margins on canvas
    c.SetRightMargin(0.09)
    c.SetLeftMargin(0.15)
    c.SetBottomMargin(0.15)

    # Add a legend and write sample name, information etc. Position dependent on variable legend choice
    latex = ROOT.TLatex ()
    latex.SetNDC(ROOT.kTRUE)
    if var_legend == True:
        legendp = GetLegendPlacement(legend_positions[histnumber])
    else: 
        legendp = GetLegendPlacement('r')
    a,b = legendp
    latex.SetTextSize(0.04)
    legend = ROOT.TLegend (a,0.64 ,b,0.74)
    latex.DrawText(a+0.01, 0.85 , "Dark Shower Sample")
    latex.SetTextSize(0.03)
    latex.DrawLatex(a+0.01, 0.8 , "#bf{pp events at #sqrt{s}=13 TeV}")
    latex.DrawLatex(a+0.01, 0.76 , "#bf{m_{Z'}=800 GeV, N_{c}=2, N_{f}=2}")
    latex.DrawLatex(a+0.01, 0.6 , "#bf{There are %i events}"%(h0.GetEntries()))
    latex.DrawLatex(a+0.01, 0.56 , "#bf{mean = %0.2f std = %0.2f}"%(h0.GetMean(),h0.GetStdDev()))
    legend.AddEntry (histoname, histotitle)
    legend.SetLineWidth (0)
    legend.Draw("same")
    
    # Print statistics
    stat_string = GetStatistics(h0,histotitle)
    print(stat_string)
           
    # Draws groomed histograms in same plot if groomed_plots=True and the groomed histograms exist
    if groomed_plots == True:
        try:
            GetGroomedPlots(f0,c,histoname,h0.GetNbinsX(),legend)
            plot_title = plot_title + "_groomed"
        except ValueError:
            print(sys.exc_info()[0], "occured: No groomed histograms exist of this kind. Will continue without groomed histograms.")
            
    # Save the canvas
    plot_title = plot_title + ".png"
    c.Print(plot_title)
    

if histnumber >= 0 and histnumber <= 29:
    MakeSinglePlot(infile, destination, histnumber, groomed_plots = False, scale_au = True, var_legend = True)
else:
    for hist in range(0,30):
        #if hist == 23:
        #    continue
        #else:
        #print(hist)
        MakeSinglePlot(infile, destination, histnumber = hist, groomed_plots = False, scale_au = True, var_legend = False)
