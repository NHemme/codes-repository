#!/usr/bin/env python
# coding: utf-8
# ================================================================================================================
# Authors and user guide
# ================================================================================================================
###################################################################################################################
# Authors and contacts:
# Guillaume Albouy: guillaume.albouy@etu.univ-grenoble-alpes.fr
# Akanksha Singh: akki153209@gmail.com
# Harikrishnan Nair: hunair1996@gmail.com
# Nicoline Hemme: nicolinehemme@gmail.com

# This code needs helpers.py which are stored in the same folder as this file
# We assume that different branches corresponding to different jet radii are defined in root files
# In this code jet clustering radius is fixed to 1.4, can be changed at line 59
# This code will analyse input sample and create following reconstructed level normalized distributions:
# 1) pt of leading/subleading jet
# 2) dijet invarint mass
# 3) missing energy
# 4) transverse mass of dijet and met system
# 5) transverse ratio
# 6) delta phi between missin energy and leading/subleading jet
# 7) 2D histo of track pT of leading jet
# 8) delta eta between leading and subleading jet
# 9) pt and invariant mass for trimmed and SoftDropped leading/subleading jets

# command: python /path_of_code/analysis.py /path_of_rootfile/name_of_rootfile.root /path_of_rootfile/output_name.root
# Takes the rootfile as input and computes the defined variables and fills the respective histograms.
###################################################################################################################


import sys
import numpy as np
import ROOT
from array import array
from helpers import *


try:
  input = raw_input
except:
  pass

if len(sys.argv) < 3:
  print(" Usage: python analysis_gen.py /path/delphes_file.root /path/output.root")
  sys.exit(1)

ROOT.gSystem.Load("libDelphes")

try:
        ROOT.gInterpreter.Declare('#include "classes/DelphesClasses.h"')
        ROOT.gInterpreter.Declare('#include "external/ExRootAnalysis/ExRootTreeReader.h"')
except:
        pass

    
# ================================================================================================================
# Define the parameters of selection cuts
# ================================================================================================================
# Radius of jets (0.4, 1.0, 1.4) :
R = 1.4
# Events selection (pT in GeV)
pT_min_jet1 = 250
pT_min_jet2 = 250
eta_max = 2.5
M_PI = 3.14


# ================================================================================================================
# Get files from input and create ROOT objects
# ================================================================================================================
inputFile = sys.argv[1]
print("Input file :")
print(inputFile)

outputFile = sys.argv[2]
print("output file :")
print(outputFile)

# Create chain of root trees
chain = ROOT.TChain("Delphes")
chain.Add(inputFile)

# Create object of class ExRootTreeReader
treeReader = ROOT.ExRootTreeReader(chain)
numberOfEntries = treeReader.GetEntries()

# Get pointer to branches used in this analysis
# R-jet branches : 04, 08, 10, 12, 14
R_jet = str(int(R*10))
if R<1.0: R_jet = '0' + R_jet

invpid  = [12, 14, 16, -12, -14, -16, 6000111,-6000111] #Neutrinoes (& anti) and??
piPID = [4900111, 4900211]
rhoPID = [4900113, 4900213]



# ================================================================================================================
# Getting the required branches from Delphes ROOT file.
# ================================================================================================================
###############################
# Reconstruction level branches
###############################
branchPtcl = treeReader.UseBranch("Particle")
branchJet = treeReader.UseBranch("ParticleFlowJet%s"%R_jet)
branchMET = treeReader.UseBranch("MissingET")
branchtrack = treeReader.UseBranch("Track")

###########################
# Generator level branches
###########################
branchGenMET = treeReader.UseBranch("GenMissingET")
branchGenJet = treeReader.UseBranch("GenJet%s"%R_jet)
branchGenTrack = treeReader.UseBranch("GenTrack")


# ================================================================================================================
# Book histograms at generator level
# ================================================================================================================
Nbins = 150
hist1GenJet1PT = ROOT.TH1F("gen_jet1_pt", "Gen Lead jet p_{T} with R=%.1f;p_{T} GeV;A.U."%(R), Nbins, pT_min_jet1, 2000.0)
hist1GenJet2PT = ROOT.TH1F("gen_jet2_pt", "Gen sub-Lead jet p_{T} with R=%.1f;p_{T} GeV;A.U."%(R), Nbins, pT_min_jet2, 2000.0)
histGenJetmJJ = ROOT.TH1F("gen jet_mJJ", "Gen Invariant mass m_{JJ} with R=%.1f"%(R), Nbins, 0.0, 3500.0)
histGenMET = ROOT.TH1F("gen_met", "Gen Missing transverse energy", Nbins, .0, 1500.0) # to be implemented
histGenmT = ROOT.TH1F("gen_jet_met_mt" , "Gen transverse mass of jet + MET" , Nbins, 0.0 , 3000.0)
histGenJetdphi = ROOT.TH1F("gen_jet_dphi", "dphi between Gen jets with R=%.1f"%(R), Nbins, 0.0, 6.0)

Nbins = 50
histGenrinv = ROOT.TH1F("gen_rinv", "Gen invisible ratio", Nbins, 0.0, 1.0)
histGenratT = ROOT.TH1F("gen_jet_met_rt" , "Gen transverse ratio of jet + MET" , Nbins, 0.0 , 1.0)

#Nbins = 40
#histGenJet1Ntrk
#histGenJet2Ntrk

Nbins = 50
histGenPions = ROOT.TH1F("gen_pions", "Gen number of pions", Nbins, 0.0, 50)
histGenRhos = ROOT.TH1F("gen_rhos", "Gen number of rhos", Nbins, 0.0, 50)

histGendphi1 = ROOT.TH1F("gen_dPhi1", "Gen dPhi_jet1_met", Nbins, -4.0, 4.0)
histGendphi2 = ROOT.TH1F("gen_dPhi2", "Gen dPhi_jet2_met", Nbins, -4.0, 4.0)
histGendelEta = ROOT.TH1F("gen_delta_eta", "Gen delta_eta_jet1_and_jet2", Nbins, -4.0, 4.0)

#Nbins = 150
#histGen2Jet1PT = ROOT.TH1F("jet1_pt_softdrop", "Lead jet p_{T} with R=%.1f;p_{T} GeV;A.U."%(R), Nbins, 0.0, 2000.0)
#histGen2Jet2PT = ROOT.TH1F("jet2_pt_softdrop", "Sub jet p_{T} with R=%.1f;p_{T} GeV;A.U."%(R), Nbins, 0.0, 2000.0)
#hist2JetmJJ = ROOT.TH1F("jet_mJJ_softdrop", " Invariant mass m_{JJ} with R=%.1f;m_{JJ} GeV;A.U."%(R), Nbins, 0.0, 3500.0)

#histGen3Jet1PT = ROOT.TH1F("jet1_pt_trimmed", "Lead jet p_{T} with R=%.1f;p_{T} GeV;A.U."%(R), Nbins, 0.0, 2000.0)
#histGen3Jet2PT = ROOT.TH1F("jet2_pt_trimmed", "Sub jet p_{T} with R=%.1f;p_{T} GeV;A.U."%(R), Nbins, 0.0, 2000.0)
#hist3JetmJJ = ROOT.TH1F("jet_mJJ_trimmed", " Invariant mass m_{JJ} with R=%.1f;m_{JJ} GeV;A.U."%(R), Nbins, 0.0, 3500.0)


# ================================================================================================================
# Book histograms at reconstruction level
# ================================================================================================================
Nbins = 150
hist1Jet1PT = ROOT.TH1F("jet1_pt", "Lead jet p_{T} with R=%.1f;p_{T} GeV;A.U."%(R), Nbins, pT_min_jet1, 2000.0)
hist1Jet2PT = ROOT.TH1F("jet2_pt", "Sub jet p_{T} with R=%.1f;p_{T} GeV;A.U."%(R), Nbins, pT_min_jet2, 2000.0)
histMET = ROOT.TH1F("jet_met", "Missing transverse energy;MET GeV;A.U.", Nbins, .0, 3000.0)
hist1JetmJJ = ROOT.TH1F("jet_mJJ", "Invariant mass m_{JJ} with R=%.1f;m_{JJ} GeV;A.U."%(R), Nbins, 0.0, 3500.0)
histmT = ROOT.TH1F("jet_met_mt" , "Transverse mass of jet + MET" , Nbins, 0.0 , 3000.0)

Nbins = 100
histratT = ROOT.TH1F("jet_met_rt" , "Transverse ratio of jet + MET" , Nbins, 0.0 , 1.0)

Nbins = 40
hist1Jet1Ntrk = ROOT.TH1F("jet1_ntrk", "Lead jet N_{trk} with R=%.1f;N_{trk};A.U."%(R), Nbins, .0, 170.0)
hist1Jet2Ntrk = ROOT.TH1F("jet2_ntrk", "Sub jet N_{trk} with R=%.1f;N_{trk};A.U."%(R), Nbins, .0, 170.0)

Nbins = 50
histdphi1 = ROOT.TH1F("dPhi1", "dPhi_jet1_met", Nbins, -4.0, 4.0)
histdphi2 = ROOT.TH1F("dPhi2", "dPhi_jet2_met", Nbins, -4.0, 4.0)
histdelEta = ROOT.TH1F("delta_eta", "delta_eta_jet1_and_jet2", Nbins, -4.0, 4.0)


#Softdrop and trimmed histograms
Nbins = 150
hist2Jet1PT = ROOT.TH1F("jet1_pt_softdrop", "Softdrop lead jet p_{T} with R=%.1f;p_{T} GeV;A.U."%(R), Nbins, 0.0, 2000.0)
hist2Jet2PT = ROOT.TH1F("jet2_pt_softdrop", "Softdrop sub jet p_{T} with R=%.1f;p_{T} GeV;A.U."%(R), Nbins, 0.0, 2000.0)
hist2JetmJJ = ROOT.TH1F("jet_mJJ_softdrop", "Softdrop invariant mass m_{JJ} with R=%.1f;m_{JJ} GeV;A.U."%(R), Nbins, 0.0, 3500.0)

hist3Jet1PT = ROOT.TH1F("jet1_pt_trimmed", "Trimmed lead jet p_{T} with R=%.1f;p_{T} GeV;A.U."%(R), Nbins, 0.0, 2000.0)
hist3Jet2PT = ROOT.TH1F("jet2_pt_trimmed", "Trimmed sub jet p_{T} with R=%.1f;p_{T} GeV;A.U."%(R), Nbins, 0.0, 2000.0)
hist3JetmJJ = ROOT.TH1F("jet_mJJ_trimmed", "Trimmed invariant mass m_{JJ} with R=%.1f;m_{JJ} GeV;A.U."%(R), Nbins, 0.0, 3500.0)

#Ignore 2D histogram for now
#histtrackpt = ROOT.TH2F("track_pt", "track PT vs Ntrk1", 20,0.0,100.0 , 40,0.0,160.0)


# ================================================================================================================
# Perform analysis - loop over all events
# ================================================================================================================

for entry in range(0, numberOfEntries):
    # Load selected branches with data from specified event
    treeReader.ReadEntry(entry)
    totpx, totpy = 0, 0
    npi, nrho = 0, 0 
    allmesons, stablemesons = 0,0
    npi_unstable, nrho_unstable = 0,0
    for iptcl in range(branchPtcl.GetEntries()):
        if abs(branchPtcl.At(iptcl).PID) not in piPID+rhoPID:
            continue
        if branchPtcl.At(iptcl).Status not in [1,83,84]:
            print('Unexpected status of dark particle. Status is %d.'%(branchPtcl.At(iptcl).Status))
        allmesons = allmesons + 1
        if branchPtcl.At(iptcl).Status in [83,84]:
            if branchPtcl.At(iptcl).PID in piPID:
                #print('Unstable pion')
                npi_unstable = npi_unstable + 1
            if branchPtcl.At(iptcl).PID in rhoPID:
                #print('Unstable rho')
                #print(branchPtcl.At(iptcl).PID)
                nrho_unstable = nrho_unstable + 1
        #if branchPtcl.At(iptcl).Status in [83,84]:
        #    daut = branchPtcl.At(iptcl).D1
        #    print(branchPtcl.At(daut).Status)
        if branchPtcl.At(iptcl).Status == 1:
            stablemesons = stablemesons + 1
            if abs(branchPtcl.At(iptcl).PID) in piPID: npi = npi + 1
            if abs(branchPtcl.At(iptcl).PID) in rhoPID: nrho = nrho + 1
        #totpx = branchPtcl.At(iptcl).Px + totpx
        #totpy = branchPtcl.At(iptcl).Py + totpy     
    if entry%1000 == 0:
        print('Number of unstable pions: %i. Number of unstable rhos: %i'%(npi_unstable,nrho_unstable))
    #    print(allmesons, stablemesons, stablemesons/allmesons)
    if allmesons != 0:
        rinv = stablemesons/allmesons
        histGenrinv.Fill(rinv)
    
    #histGenMET.Fill(np.sqrt(totpx**2 + totpy**2))
    histGenPions.Fill(npi)
    histGenRhos.Fill(nrho)
        
    #######################
    # Generator level plots
    #######################
    if branchGenJet.GetEntries() > 1: 
        # Take the two leading jets
        genjet1 = branchGenJet.At(0)
        genjet2 = branchGenJet.At(1)
        
        if genjet1.PT > pT_min_jet1 and np.abs(genjet1.Eta) < eta_max and genjet2.PT > pT_min_jet2 and np.abs(genjet2.Eta) < eta_max :
            genjet1vec = ROOT.TLorentzVector()
            genjet2vec = ROOT.TLorentzVector()
            genjet1vec = GetJetVector(genjet1.PT , genjet1.Eta , genjet1.Phi , genjet1.Mass)
            genjet2vec = GetJetVector(genjet2.PT , genjet2.Eta , genjet2.Phi , genjet2.Mass)
            genmet = branchMET.At(0)
            genm = genmet.MET
            genMETPhi = genmet.Phi
            genMETx = genm* np.cos(genMETPhi)
            genMETy = genm* np.sin(genMETPhi)
            genvecmet = ROOT.TLorentzVector()
            genvecmet.SetPxPyPzE(genMETx , genMETy , 0 , genm)

            genmt = (genjet1vec + genjet2vec + genvecmet).Mt()
            gendphi1 = GetdPhi(genjet1.Phi, genmet.Phi)
            gendphi2 = GetdPhi(genjet2.Phi , genmet.Phi)
            gendelEta = GetdEta(genjet1.Eta , genjet2.Eta)
            if (genmt!=0): genrt = genm/genmt


            hist1GenJet1PT.Fill(branchGenJet.At(0).PT)
            hist1GenJet2PT.Fill(branchGenJet.At(1).PT)
            histGenMET.Fill(branchGenMET.At(0).MET)
            histGenJetmJJ.Fill((genjet1vec+genjet2vec).M())
            histGenJetdphi.Fill(GetdPhi(branchGenJet.At(0).Phi,branchGenJet.At(1).Phi))
            histGenmT.Fill(genmt)
            histGenratT.Fill(genrt)
            histGendphi1.Fill(gendphi1)
            histGendphi2.Fill(gendphi2)
            histGendelEta.Fill(gendelEta)


    ############################
    # Reconstruction level plots
    ############################
    if branchJet.GetEntries() > 1:
        # Take the two leading jets
        jet1 = branchJet.At(0)
        jet2 = branchJet.At(1)

        #Selecting events
        if jet1.PT > pT_min_jet1 and np.abs(jet1.Eta) < eta_max and jet2.PT > pT_min_jet2 and np.abs(jet2.Eta) < eta_max :

            # Defining the TLorentz vector for leading jet
            vec1 = GetJetVector(jet1.PT , jet1.Eta , jet1.Phi , jet1.Mass)

            # Defining TLorentz vector for sub - leading jet
            vec2 = GetJetVector(jet2.PT , jet2.Eta , jet2.Phi , jet2.Mass)

            # Defining TLorentz vector for missing energy branch
            #{
            met = branchMET.At(0)
            m = met.MET
            METPhi = met.Phi
            METx = m* np.cos(METPhi)
            METy = m* np.sin(METPhi)
            vecmet = ROOT.TLorentzVector()
            vecmet.SetPxPyPzE(METx , METy , 0 , m)
            #}

            # Computing transverse mass for jet1+jet2+met
            mt = (vec1 + vec2 + vecmet).Mt()

            # Computing dPhi between jet1 and met
            dPhi1 = GetdPhi(jet1.Phi, met.Phi)

            # Computing dphi between jet2 and met
            dPhi2 = GetdPhi(jet2.Phi , met.Phi)

            # Computing transverse ratio
            if (mt!=0): rt = m/mt

            # SoftDropped jets
            #{
            # Getting the Soft dropped jet1 four momenta
            jet1_softdrop = jet1.SoftDroppedP4[0]
            jet1_pt_softdrop = jet1_softdrop.Pt()

            # Getting the Soft dropped Jet2 four momenta
            jet2_softdrop = jet2.SoftDroppedP4[0]
            jet2_pt_softdrop = jet2_softdrop.Pt()

            # Computing invariant mass using Soft dropped leading and sub-leading jets
            invmass_softdrop = Getinvmass(jet1_softdrop , jet2_softdrop)
            #}

            # Trimmed jets
            #{
            # Getting the Trimmed jet1 four momenta
            jet1_trimmed = jet1.TrimmedP4[0]
            jet1_pt_trimmed = jet1_trimmed.Pt()

            # Getting the Trimmed jet2 four momenta
            jet2_trimmed = jet2.TrimmedP4[0]
            jet2_pt_trimmed = jet2_trimmed.Pt()

            # Computing invariant mass using Trimmed leading and sub-leading jets
            invmass_trimmed = Getinvmass(jet1_trimmed , jet2_trimmed)
            #}

            #Calculating average track transverse momentum
            #{
            delEta = GetdEta(jet1.Eta , jet2.Eta)
            track = []
            selectedtrack = []
            trackpt = []
            #Store the tracks to array called "track"
            for i in range(0 , branchtrack.GetEntries()):
                 track.append(branchtrack.At(i))

            #Loop over array "track" and calculate delta phi and delta R
            for j in range(0 , len(track)):
                dPhi = GetdPhi(track[j].Phi , jet1.Phi)

                dEta = GetdEta(track[j].Eta , jet1.Eta)

                DeltaR = np.sqrt(dEta*dEta + dPhi*dPhi)

                #Store the tracks which satisfy the condition on delta R to array "selectedtrack"
                if DeltaR < 1.4 :
                    selectedtrack.append(track[j])

            #Loop over selected tracks to calculate average momentum of the tracks
            for k in range(0 , len(selectedtrack)):
                trackpt.append(selectedtrack[k].PT)

            averagept = np.average(trackpt)
            #}

            #Fill histograms
            hist1Jet1PT.Fill(jet1.PT)
            hist1Jet2PT.Fill(jet2.PT)
            hist1JetmJJ.Fill((jet1.P4() + jet2.P4()).M())
            hist1Jet1Ntrk.Fill(jet1.NCharged)
            hist1Jet2Ntrk.Fill(jet2.NCharged)
            histMET.Fill(branchMET.At(0).MET)
            histmT.Fill(mt)
            histdphi1.Fill(dPhi1)
            histdphi2.Fill(dPhi2)
            histratT.Fill(rt)
            histdelEta.Fill(delEta)
            #histtrackpt.Fill(averagept,jet1.NCharged)
            if jet1_pt_softdrop > pT_min_jet1 and jet2_pt_softdrop > pT_min_jet2:
                hist2Jet1PT.Fill(jet1_pt_softdrop)
                hist2Jet2PT.Fill(jet2_pt_softdrop)
                hist2JetmJJ.Fill(invmass_softdrop)
            if jet1_pt_trimmed > pT_min_jet1 and jet2_pt_trimmed > pT_min_jet2:
                hist3Jet1PT.Fill(jet1_pt_trimmed)
                hist3Jet2PT.Fill(jet2_pt_trimmed)
                hist3JetmJJ.Fill(invmass_trimmed)



# ================================================================================================================
# Normalising histograms
# ================================================================================================================
if histGenrinv.GetSumw2N()==0: histGenrinv.Sumw2(True)
if hist1Jet1PT.GetSumw2N()==0: hist1Jet1PT.Sumw2(True)
if hist1Jet2PT.GetSumw2N()==0: hist1Jet2PT.Sumw2(True)
if hist1JetmJJ.GetSumw2N()==0: hist1JetmJJ.Sumw2(True)
if hist1Jet1Ntrk.GetSumw2N()==0: hist1Jet1Ntrk.Sumw2(True)
if hist1Jet2Ntrk.GetSumw2N()==0: hist1Jet2Ntrk.Sumw2(True)
if histMET.GetSumw2N()==0: histMET.Sumw2(True)
if histmT.GetSumw2N()==0: histmT.Sumw2(True)
if histdphi1.GetSumw2N()==0: histdphi1.Sumw2(True)
if histdphi2.GetSumw2N()==0: histdphi2.Sumw2(True)
if histratT.GetSumw2N()==0: histratT.Sumw2(True)
if hist2Jet1PT.GetSumw2N()==0: hist2Jet1PT.Sumw2(True)
if hist2Jet2PT.GetSumw2N()==0: hist2Jet2PT.Sumw2(True)
if hist2JetmJJ.GetSumw2N()==0: hist2JetmJJ.Sumw2(True)
if hist3Jet1PT.GetSumw2N()==0: hist3Jet1PT.Sumw2(True)
if hist3Jet2PT.GetSumw2N()==0: hist3Jet2PT.Sumw2(True)
if hist3JetmJJ.GetSumw2N()==0: hist3JetmJJ.Sumw2(True)
if histdelEta.GetSumw2N()==0: histdelEta.Sumw2(True)


# ================================================================================================================
# Creating a list and saving the histograms to the list
# ================================================================================================================
histlist = ROOT.TList()
histlist.Add(hist1GenJet1PT)
histlist.Add(hist1GenJet2PT)
histlist.Add(histGenMET)
histlist.Add(histGenPions)
histlist.Add(histGenRhos)
histlist.Add(histGenrinv)
histlist.Add(histGenJetmJJ)
histlist.Add(histGenJetdphi)
histlist.Add(histGendphi1)
histlist.Add(histGendphi2)
histlist.Add(histGenmT)
histlist.Add(histGenratT)
histlist.Add(histGendelEta)
histlist.Add(hist1Jet1PT)
histlist.Add(hist1Jet2PT)
histlist.Add(hist1JetmJJ)
histlist.Add(hist1Jet1Ntrk)
histlist.Add(hist1Jet2Ntrk)
histlist.Add(histMET)
histlist.Add(histmT)
histlist.Add(histdphi1)
histlist.Add(histdphi2)
histlist.Add(histratT)
histlist.Add(hist2Jet1PT)
histlist.Add(hist2Jet2PT)
histlist.Add(hist2JetmJJ)
histlist.Add(hist3Jet1PT)
histlist.Add(hist3Jet2PT)
histlist.Add(hist3JetmJJ)
#histlist.Add(histtrackpt)
histlist.Add(histdelEta)


# ================================================================================================================
# Final steps
# ================================================================================================================
#outputFile = inputFile[:-5] + "_combined_R14.root"
rootFile = ROOT.TFile(outputFile, "RECREATE")
histlist.Write()
rootFile.Close()

