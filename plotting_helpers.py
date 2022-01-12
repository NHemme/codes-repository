#!/usr/bin/env python
# coding: utf-8

import numpy as np
import ROOT

def GetAllHistoNames(infile = ''):
    rootfile1 = infile
    f0 = ROOT.TFile(rootfile1,"READ")
    
    listofnames = []
    listoftitles = []
    
    for obj in f0.GetListOfKeys():
        listofnames.append(obj.GetName())
        listoftitles.append(obj.GetTitle())
    
    return listofnames, listoftitles


def GetLegendPlacement(position = ''):
        if position == 'r':
            placement = [0.62,0.9]
        if position == 'c':
            placement = [0.4,0.67]
        if position == 'l':
            placement = [0.19,0.48]
        return placement

def GetStatistics(histogram,histotitle):
    h0 = histogram
    q = np.array([0.0,0.0,0.0])
    p = np.array([0.25,0.5,0.75])
    h0.GetQuantiles(3,q,p)
    x_unit = h0.GetXaxis().GetTitle()
    if x_unit.find("trk") != -1:
        x_unit = x_unit
    else:
        index = x_unit.find("[")
        x_unit = x_unit[index+1:-1]
    stat_string = """
    Some statistics for %s: 
    There are %i events in the histogram. They are binned into %i bins. 
    The mean is %0.2f %s and the standard deviation is %0.2f %s. 
    Half of the events are contained between %0.2f and %0.2f %s.
    """%(histotitle,h0.GetEntries(), h0.GetNbinsX(),h0.GetMean(),x_unit,h0.GetStdDev(),x_unit,q[0],q[2],x_unit)
    return stat_string


    
def GetGroomedPlots(rootfile, canvas, histoname, n_bins,legend):
    h1_name = histoname + "_softdrop"
    h2_name = histoname + "_trimmed"
    # The following line checks if the groomed histograms exist and will raise an error and leave the function if not
    if rootfile.GetListOfKeys().Contains(h1_name) == False or rootfile.GetListOfKeys().Contains(h2_name) == False:
        raise ValueError
        return
    # Enables drawing on the same canvas
    canvas.cd()
    h1=rootfile.Get(h1_name)
    h2=rootfile.Get(h2_name)
    # Binning depends on the n_bins of the un-groomed histogram so that they have equal number of bins
    h1_bins = h1.GetNbinsX()
    binning = h1_bins/n_bins
    h1.Rebin(int(binning))
    h2_bins = h2.GetNbinsX()
    binning = h2_bins/n_bins
    h2.Rebin(int(binning))
    # Sets options for drawing and draws
    h1.SetLineWidth(2)
    h1.SetLineColor(ROOT.kMagenta-3)
    h1.DrawCopy('hist same')
    h1.DrawCopy('E1 same')
    h2.SetLineWidth(2)
    h2.SetLineColor(ROOT.kCyan-3)
    h2.DrawCopy('hist same')
    h2.DrawCopy('E1 same')
    # Defined titles and adds legends
    h1_title = h1.GetTitle()
    index = h1_title.find('with')
    h1_title = h1_title[:index] + "soft dropped"
    h2_title = h2.GetTitle()
    index = h2_title.find('with')
    h2_title = h2_title[:index] + "trimmed"
    legend.AddEntry(h1,h1_title)
    legend.AddEntry(h2,h2_title)
    print(GetStatistics(h1,h1_name))
    print(GetStatistics(h2,h2_name))
    
