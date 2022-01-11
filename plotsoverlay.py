#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import os, sys, shutil
from numpy import *
from pylab import *
from scipy.interpolate import interp1d
import ROOT
from plotting_helpers import *


if len(sys.argv) < 5:
    print("Usage: python /path/plotsinglemodel.py /path/output.root")
    print("Optional arguments: histogram_number(integer) groomed_plots(boolean) scale_au(boolean) var_legend(boolean)")
    sys.exit(1)
    
infile1 = sys.argv[1]
infile2 = sys.argv[2]
id1 = sys.argv[3]
id2 = sys.argv[4]
if len(sys.argv) > 5:
    histnumber = int(sys.argv[5])
else:
    histnumber = -1


def MakeSinglePlot(infile1 = '', infile2 = '', id1 = '', id2 = '', histnumber = int, groomed_plots = False, var_legend = True):
    rootfile1 = infile1
    rootfile2 = infile2
    histoname1 = GetAllHistoNames(infile1)[0][histnumber]
    histotitle1 = GetAllHistoNames(infile1)[1][histnumber] + ' %s'%(id1)
    histoname2 = GetAllHistoNames(infile2)[0][histnumber]
    histotitle2 = GetAllHistoNames(infile2)[1][histnumber] + ' %s'%(id2)
    plot_title = "Plots/%s" %(histoname1)
    #legend_positions = ['r','r','r','r','r','r','l','c','r','r','r','r','r','r','r','r','r','r']
   
    ROOT.gROOT.SetBatch(True)
    c=ROOT.TCanvas("c","c",800,600)
    f1 = ROOT.TFile(rootfile1,"READ")
    f2 = ROOT.TFile(rootfile2,"READ")
    
    # Get the histogram, you can rebin the histogram with the commented out rebin option
    h1=f1.Get(histoname1)
    n_bins = h1.GetNbinsX()
    h2=f2.Get(histoname2)
    n_bins = h2.GetNbinsX()
    divider = 25 #larger divider gives more bins
    while n_bins%divider != 0 :
        divider += 1
    if divider == n_bins:
        divider = divider/2
    h1.Rebin(int(n_bins/divider))
    h2.Rebin(int(n_bins/divider))
    #print(n_bins, int(n_bins/divider), h1.GetNbinsX())

    
    c.cd()
    #Set axis title, you can use latex as demonstrated
    h1.GetXaxis()
    h1.GetYaxis()
    #switch off the ugly stats box on the upper right corner
    h1.SetStats(0)
    
    #set label sizes
    h1.GetYaxis().SetLabelSize(0.04)
    h1.GetYaxis().SetTitleSize(0.04)
    h1.GetXaxis().SetTitleSize(0.04)
    h1.GetXaxis().SetLabelSize(0.04)
    h1.SetTitle("")
    #set linewidth and color of the lines
    h1.SetLineWidth(2)
    h1.SetLineColor(ROOT.kBlue-2)
    #actually draw the histogram
    h1.Draw('hist')
    h1.Draw('E1 same')
    #set linewidth and color of the lines
    h2.SetLineWidth(2)
    h2.SetLineColor(ROOT.kMagenta-3)
    #actually draw the histogram
    h2.Draw('hist same')
    h2.Draw('E1 same')
    
    #set axis offsets so that labels don't get messed up
    h1.GetXaxis().SetTitleOffset(1.4)
    h1.GetYaxis().SetTitleOffset(1.4)
    #set axis ranges
    max_x = max(h1.GetXaxis().GetXmax(),h2.GetXaxis().GetXmax())
    min_x = min(h1.GetXaxis().GetXmin(),h2.GetXaxis().GetXmin())
    h1.GetXaxis().SetRangeUser(min_x, max_x)
    max_y = max(h1.GetBinContent(h1.GetMaximumBin()),h2.GetBinContent(h2.GetMaximumBin()))*1.1
    h1.GetYaxis().SetRangeUser(0, max_y)
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
    legend.AddEntry (h1, histotitle1)
    legend.AddEntry (h2, histotitle2)
    legend.SetLineWidth (0)
    legend.Draw("same")
    
    # Print statistics
    #stat_string = GetStatistics(h1,histotitle1)
    #print(stat_string)
           
    # Draws groomed histograms in same plot if groomed_plots=True and the groomed histograms exist
    if groomed_plots == True:
        try:
            GetGroomedPlots(f1,c,histoname1,h1.GetNbinsX(),legend)
            plot_title = plot_title + "_groomed"
        except ValueError:
            print(sys.exc_info()[0], "occured: No groomed histograms exist of this kind. Will continue without groomed histograms.")
            
    # Save the canvas
    plot_title = plot_title + "_overlay_%s_%s.png"%(id1,id2)
    c.Print(plot_title)
    

if histnumber >= 0 and histnumber <= 28:
    MakeSinglePlot(infile1, infile2, id1, id2, histnumber, groomed_plots = False, scale_au = True, var_legend = True)
else:
    for hist in range(0,29):
        #if hist == 23:
        #    continue
        #else:
        print(hist)
        MakeSinglePlot(infile1, infile2, id1, id2, histnumber = hist, groomed_plots = False, var_legend = False)
