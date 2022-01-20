#!/usr/bin/env python
# coding: utf-8


#!/usr/bin/env python
import os, sys, shutil
from numpy import *
from pylab import *
from scipy.interpolate import interp1d
import ROOT
from plotting_helpers import *


def MakeSinglePlot(infile = '', destination = './', histnumber = int, groomed_plots = False, scale_au = True, var_legend = True):
    rootfile1 = infile
    histoname = GetAllHistoNames(infile)[0][histnumber]
    histotitle = GetAllHistoNames(infile)[1][histnumber]
    plot_title = "%s%s" %(destination,histoname)
    #legend_positions = ['r','r','r','r','r','r','l','c','r','r','r','r','r','r','r','r','r','r']

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
    if 'diff' not in histoname and 'rinv' not in histoname:
        h0.Rebin(int(n_bins/divider))
    #print(n_bins, int(n_bins/divider), h0.GetNbinsX())
    # The following lines will normalise the histogram if scale_au=True and groomed_plots=False
    if scale_au == True and groomed_plots == False and h0.Integral()!=0:
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
    h0.SetLineColorAlpha(ROOT.kBlue, 0.35)
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
        legendp = GetLegendPlacement('l')
    a,b = legendp
    latex.SetTextSize(0.04)
    legend = ROOT.TLegend (a,0.64 ,b,0.74)
    latex.DrawText(a+0.01, 0.85 , "Dark Shower Sample")
    latex.SetTextSize(0.03)
    latex.DrawLatex(a+0.01, 0.8 , "#bf{pp events at #sqrt{s}=13 TeV}")
    latex.DrawLatex(a+0.01, 0.76 , "#bf{m_{Z'}=800 GeV, N_{c}=3, N_{f}=3}")
    latex.DrawLatex(a+0.01, 0.6 , "#bf{There are %i events}"%(h0.GetEntries()))
    latex.DrawLatex(a+0.01, 0.56 , "#bf{mean = %0.4f std = %0.4f}"%(h0.GetMean(),h0.GetStdDev()))
    legend.AddEntry (histoname, histotitle)
    legend.SetLineWidth (0)
    legend.Draw("same")
    
    # Print statistics
    stat_string = GetStatistics(h0,histotitle)
    #print(stat_string)
           
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
    c.Close()

    
    

def MakePIDPlot(infile = '', destination = './', histnumber = int):
    rootfile1 = infile
    histoname = GetAllHistoNames(infile)[0][histnumber]
    histotitle = GetAllHistoNames(infile)[1][histnumber]
    plot_title = "%s%s" %(destination,histoname)
    #legend_positions = ['r','r','r','r','r','r','l','c','r','r','r','r','r','r','r','r','r','r']

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
    #print(n_bins, int(n_bins/divider), h0.GetNbinsX())
    # The following lines will normalise the histogram if scale_au=True and groomed_plots=False
    #if h0.Integral()!=0:
    #    norm = h0.GetEntries()
    #    scale = 1/h0.Integral();
    #    h0.Scale(scale);
    
    c.cd()
    #Set axis title, you can use latex as demonstrated
    h0.GetXaxis()

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
    h0.SetBarWidth(2)
    h0.SetFillColorAlpha(ROOT.kBlue, 0.35)
    h0.Draw('bar')
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
    content_list = []
    for b in range(n_bins):
        content = h0.GetBinContent(b)
        if content != 0:
            content_list.append(content)
            b_value = min_x + b*((max_x-min_x)/n_bins)
            print(content, h0.GetXaxis().FindBin(content), b, b_value)
    print(content_list)
    #set margins on canvas
    c.SetRightMargin(0.09)
    c.SetLeftMargin(0.15)
    c.SetBottomMargin(0.15)

    # Add a legend and write sample name, information etc. Position dependent on variable legend choice
    latex = ROOT.TLatex ()
    latex.SetNDC(ROOT.kTRUE)
    legendp = GetLegendPlacement('c')
    a,b = legendp
    latex.SetTextSize(0.04)
    legend = ROOT.TLegend (a-0.1,0.64 ,b,0.70)
    latex.DrawText(a-0.09, 0.85 , "Dark Shower Sample")
    latex.SetTextSize(0.03)
    latex.DrawLatex(a-0.09, 0.8 , "#bf{pp events at #sqrt{s}=13 TeV}")
    latex.DrawLatex(a-0.09, 0.76 , "#bf{m_{Z'}=800 GeV, N_{c}=3, N_{f}=3}")
    if len(content_list)!=0:
        latex.DrawLatex(0.24, 0.84 , "#bf{-213}")
        latex.DrawLatex(0.27, 0.4 , "#bf{-211}")
        latex.DrawLatex(0.66, 0.84 , "#bf{113}")
        latex.DrawLatex(0.63, 0.4, "#bf{111}")
        latex.DrawLatex(0.785, 0.84 , "#bf{213}")
        latex.DrawLatex(0.76, 0.4, "#bf{211}")
    legend.AddEntry (histoname, histotitle)
    legend.SetLineWidth (0)
    legend.Draw("same")
    
    # Save the canvas
    plot_title = plot_title + ".png"
    c.Print(plot_title)
    c.Close()