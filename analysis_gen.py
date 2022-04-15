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
# In this code jet clustering radius is fixed to 1.4, can be changed below 
# This code will analyse input sample and create several distributions

# command: python /path_of_code/analysis.py /path_of_rootfile/name_of_rootfile.root /path_of_output/output_name.root
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
        ROOT.gInterpreter.Declare('#include "modules/FastJetFinder.h"')
except:
        pass

    
# ================================================================================================================
# Define the parameters of selection cuts
# ================================================================================================================
# Radius of jets (0.4, 0.8, 1.0, 1.2 1.4) and get pointer to branches used in this analysis
R = 0.4
R_jet = str(int(R*10))
if R<1.0: R_jet = '0' + R_jet
    
# Events selection (pT in GeV)
pT_min_jet1 = 30
pT_min_jet2 = 30
eta_max = 3
dphi_max = 0.6 # For dphi bump investigation

# Define constants
M_PI = 3.14

# Model-dependent parameters that may be used in code - change as needed
mz_model = 1000 
jetFinder = ROOT.FastJetFinder()
jetFinder.fParameterR = R
jetFinder.JetAlgorithm = 6
jetFinder.JetPTMin = pT_min_jet1

# Relevant PID lists - add new PIDs as needed when flavours change or new mesons are introduced
invPID  = [12, 14, 16, -12, -14, -16, 6000111,-6000111] #Neutrinoes (& anti) and??
piPID = [4900111, 4900221, 4900331, 4900211, 4900311, 4900321, -490111, -4900221, -4900331, -4900211, -4900311, -4900321]
rhoPID = [4900113, 4900223, 4900333, 4900213, 4900313, 4900323, -4900113, -4900223, -4900333, -4900213, -4900313, -4900323]
deltaPID = [4901114, -4901114]
quarkPID = [4900101, 4900102, 4900103,-4900101, -4900102, -4900103]
darkPID =  piPID+rhoPID+deltaPID


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


# ================================================================================================================
# Getting the required branches from Delphes ROOT file.
# ================================================================================================================
branchPtcl = treeReader.UseBranch("Particle")

# Reconstruction level branches
branchJet = treeReader.UseBranch("ParticleFlowJet%s"%R_jet)
branchMET = treeReader.UseBranch("MissingET")
branchtrack = treeReader.UseBranch("Track")

# Generator level branches
branchGenMET = treeReader.UseBranch("GenMissingET")
branchGenJet = treeReader.UseBranch("GenJet%s"%R_jet)
branchGenTrack = treeReader.UseBranch("GenTrack")


# ================================================================================================================
# Book histograms for particle analysis 
# ================================================================================================================
histGenZmJJ = ROOT.TH1F("gen_Z_mJJ", "Invariant mass of Z' decay products;m_{JJ} [GeV];N_{events}", 150, 0.0, 1500.0)
histGenDarkMET = ROOT.TH1F("gen_dark_met", "Gen MET from dark particles ;Dark Particles MET [GeV];N_{events}", 250, 0.0, 1500.0)
histGenrinv = ROOT.TH1F("gen_rinv", "Gen r_inv;r_{inv};N_{events}", 50, 0.0, 1)
histGenPions = ROOT.TH1F("gen_pions", "Gen number of pions;Number of dark pions;N_{events}", 50, 0.0, 50)
histGenRhos = ROOT.TH1F("gen_rhos", "Gen number of rhos;Number of dark rhos;N_{events}", 50, 0.0, 50)
histGenInvPtcls1 = ROOT.TH1F("gen_ninv_jet1", "Number of invisible particles in leading jet (GenJet, R=%.1f);Percentage of particles in leading jet that are invisible (GenJet, R=%.1f) [%%];N_{events}"%(R,R), 90, 0, 30)
histGenInvPtcls2 = ROOT.TH1F("gen_ninv_jet2", "Number of invisible particles in sub-leading jet (GenJet, R=%.1f);Percentage of particles in sub-leading jet that are invisible (GenJet, R=%.1f) [%%];N_{events}"%(R,R), 90, 0, 30)
histGenInvPtcls3 = ROOT.TH1F("gen_ninv_jet3", "Number of invisible particles in sub-leading jet (GenJet, R=%.1f);Percentage of particles in third jet that are invisible (GenJet, R=%.1f) [%%];N_{events}"%(R,R), 90, 0, 30)
histGenDeltaCharge = ROOT.TH1F("gen_delta_charge_dark", "Gen total charge difference in dark particles pr event;#Delta charge;N_{events}", 80, -20, 20)
histGenPID = ROOT.TH1F("gen_pid", "Gen dark PID;PID;N_{events}", 400, -400, 400)

# Testing dphi bump
histGenZpT = ROOT.TH1F("gen_Z_pT", "Transverse momentum of Z' boson;p_{T} [GeV];N_{events}", 100, 0.0, 1000.0)
histGenqDarkdPhi = ROOT.TH1F("gen_qdark_dphi", "#Delta#phi of Z' decay dark quark pair;#Delta#phi [rad];N_{events}", 50, -3.2,3.2)
histGenUpperpT = ROOT.TH1F("gen_upper_pT", "Transverse momentum of hemisphere centered around leading jet;p_{T} [GeV];N_{events}", 200, 0, 500)
histGenLowerpT = ROOT.TH1F("gen_lower_pT", "Transverse momentum of hemisphere opposite leading jet;p_{T} [GeV];N_{events}", 200, 0.0, 500)
histGenInvPtclsUpper = ROOT.TH1F("gen_ninv_upper", "Percentage of particles that are invisible in hemisphere centered around leading jet (GenJet, R=%.1f);Percentage of particles that are invisible [%%];N_{events}"%(R), 200, 0, 5)
histGenInvPtclsLower = ROOT.TH1F("gen_ninv_lower", "Percentage of particles that are invisible in hemisphere opposite leading jet (GenJet, R=%.1f);Percentage of particles that are invisible [%%];N_{events}"%(R), 200, 0, 5)


# ================================================================================================================
# Book histograms at generator level
# ================================================================================================================
# All jet and MET variables
hist1GenJet1PT = ROOT.TH1F("gen_jet1_pt", "Gen leading jet p_{T} R=%.1f;p_{T} [GeV] (GenJet, R=%.1f);N_{events}"%(R,R), 150, pT_min_jet1, 1500.0)
hist1GenJet2PT = ROOT.TH1F("gen_jet2_pt", "Gen sub-leading jet p_{T} R=%.1f;p_{T} [GeV] (GenJet, R=%.1f);N_{events}"%(R,R), 150, pT_min_jet2, 1500.0)
hist1GenJet3PT = ROOT.TH1F("gen_jet3_pt", "Gen third jet p_{T} R=%.1f;p_{T} [GeV] (GenJet, R=%.1f);N_{events}"%(R,R), 150, pT_min_jet2, 1500.0)
histGenJetmJJ = ROOT.TH1F("gen_jet_mJJ", "Gen di-jet invariant mass R=%.1f;m_{JJ} [GeV] (GenJets, R=%.1f);N_{events}"%(R,R), 150, 0.0, 3000.0)
histGenMET = ROOT.TH1F("gen_met", "Gen MET R=%.1f;MET [GeV] (GenJets, R=%.1f);N_{events}"%(R,R), 150, .0, 1500.0)
histGenmT = ROOT.TH1F("gen_jet_met_mt" , "Gen transverse mass of jet + MET R=%.1f; M_{T} [GeV] (GenJets, R=%.1f);N_{events}"%(R,R) , 250, 0.0 , 2000.0)
histGenJetdphi = ROOT.TH1F("gen_jet_dphi", "Gen dphi between leading and sub-leading jets R=%.1f;#Delta#phi [rad] (GenJets, R=%.1f);N_{events}"%(R,R), 50, -3.2, 3.2)
histGenJetdphimin = ROOT.TH1F("gen_jet_dphi_min", "Gen minimum dphi between 2 leading jets and MET;#Delta#phi_{min} [rad] (GenJets, R=%.1f);N_{events}"%(R), 50, -3.2, 3.2)
histGenJet13dphi = ROOT.TH1F("gen_jet_dphi_13", "Gen dphi between leading and third jets R=%.1f;#Delta#phi [rad] (GenJets, R=%.1f);N_{events}"%(R,R), 50, -3.2, 3.2)
histGenJet23dphi = ROOT.TH1F("gen_jet_dphi_23", "Gen dphi between sub-leading and third jets R=%.1f;#Delta#phi [rad] (GenJets, R=%.1f);N_{events}"%(R,R), 50, -3.2, 3.2)
histGenAllJetsMT = ROOT.TH1F("gen_alljets_mt", "Gen all jets mass R=%.1f;Mass [GeV] (GenJets, R=%.1f);N_{events}"%(R,R), 250, 0, 1000.0)
histGenAllJetsPT = ROOT.TH1F("gen_alljets_pt", "Gen all jets p_{T} R=%.1f;p_{T} [GeV] (GenJets, R=%.1f);N_{events}"%(R,R), 200, pT_min_jet1, 500.0)
histGenNjets = ROOT.TH1F("gen_njets", "Gen number of jets in event with R=%.1f;Number of jets (GenJets, R=%.1f);N_{events}"%(R,R), 20, 0, 20.0)
histGenJetMET = ROOT.TH1F("gen_jet_MET", "Gen MET from negative jet p_{T} R=%.1f;MET [GeV] (GenJets, R=%.1f);N_{events}"%(R,R), 150, 0, 1500.0)
histGenratT = ROOT.TH1F("gen_jet_met_rt" , "Gen transverse ratio MET/M_T R=%.1f; R_{T} (GenJets, R=%.1f);N_{events}"%(R,R), 50, 0.0 , 1.0)
histGendphi1 = ROOT.TH1F("gen_dPhi1", "Gen dPhi_jet1_met;#Delta#phi (GenJet, R=%.1f);N_{events}", 50, -3.2, 3.2)
histGendphi2 = ROOT.TH1F("gen_dPhi2", "Gen dPhi_jet2_met;#Delta#phi (GenJet, R=%.1f);N_{events}", 50, -3.2, 3.2)
histGendelEta = ROOT.TH1F("gen_delta_eta", "Gen delta_eta_jet1_and_jet2;#Delta#eta;N_{events}", 50, -4, 4)


# MET differences 
histMETdiffGenDark = ROOT.TH1F("met_diff_gen_dark", "Difference in MET from gen and counted dark particles;MET Diff [GeV];N_{events}", 100, -100, 100)
histMETdiffGenRec = ROOT.TH1F("met_diff_gen_rec", "Difference in MET from gen and rec;MET Diff [GeV];N_{events}", 100, -100, 100)



# ================================================================================================================
# Book histograms at reconstruction level
# ================================================================================================================
Nbins = 150
hist1Jet1PT = ROOT.TH1F("jet1_pt", "Lead jet p_{T} with R=%.1f;p_{T} GeV;N_{events}"%(R), Nbins, pT_min_jet1, 2000.0)
hist1Jet2PT = ROOT.TH1F("jet2_pt", "Sub jet p_{T} with R=%.1f;p_{T} GeV;N_{events}"%(R), Nbins, pT_min_jet2, 2000.0)
histMET = ROOT.TH1F("jet_met", "Missing transverse energy;MET [GeV];N_{events}", Nbins, .0, 1500.0)
hist1JetmJJ = ROOT.TH1F("jet_mJJ", "Invariant mass m_{JJ} with R=%.1f;m_{JJ} GeV;N_{events}"%(R), Nbins, 0.0, 2500.0)
histmT = ROOT.TH1F("jet_met_mt" , "Transverse mass of jet + MET" , Nbins, 0.0 , 3000.0)

Nbins = 100
histratT = ROOT.TH1F("jet_met_rt" , "Transverse ratio of jet + MET" , Nbins, 0.0 , 1.0)

Nbins = 50
hist1Jet1Ntrk = ROOT.TH1F("jet1_ntrk", "Lead jet N_{trk} with R=%.1f;N_{trk};N_{events}"%(R), Nbins, .0, 100.0)
hist1Jet2Ntrk = ROOT.TH1F("jet2_ntrk", "Sub jet N_{trk} with R=%.1f;N_{trk};N_{events}"%(R), Nbins, .0, 100.0)

Nbins = 50
histdphi1 = ROOT.TH1F("dPhi1", "dPhi_jet1_met;#Delta #phi [rad]; N_{events}", Nbins, -4.0, 4.0)
histdphi2 = ROOT.TH1F("dPhi2", "dPhi_jet2_met;#Delta #phi [rad]; N_{events}", Nbins, -4.0, 4.0)
histdelEta = ROOT.TH1F("delta_eta", "delta_eta_jet1_and_jet2;#Delta#eta; N_{events}", Nbins, -4.0, 4.0)


#Softdrop and trimmed histograms
Nbins = 150
hist2Jet1PT = ROOT.TH1F("jet1_pt_softdrop", "Softdrop lead jet p_{T} with R=%.1f;p_{T} GeV;N_{events}"%(R), Nbins, 0.0, 2000.0)
hist2Jet2PT = ROOT.TH1F("jet2_pt_softdrop", "Softdrop sub jet p_{T} with R=%.1f;p_{T} GeV;N_{events}"%(R), Nbins, 0.0, 2000.0)
hist2JetmJJ = ROOT.TH1F("jet_mJJ_softdrop", "Softdrop invariant mass m_{JJ} with R=%.1f;m_{JJ} GeV;N_{events}"%(R), Nbins, 0.0, 2500.0)

hist3Jet1PT = ROOT.TH1F("jet1_pt_trimmed", "Trimmed lead jet p_{T} with R=%.1f;p_{T} GeV;N_{events}"%(R), Nbins, 0.0, 2000.0)
hist3Jet2PT = ROOT.TH1F("jet2_pt_trimmed", "Trimmed sub jet p_{T} with R=%.1f;p_{T} GeV;N_{events}"%(R), Nbins, 0.0, 2000.0)
hist3JetmJJ = ROOT.TH1F("jet_mJJ_trimmed", "Trimmed invariant mass m_{JJ} with R=%.1f;m_{JJ} GeV;N_{events}"%(R), Nbins, 0.0, 2500.0)

#Ignore 2D histogram for now
#histtrackpt = ROOT.TH2F("track_pt", "track PT vs Ntrk1", 20,0.0,100.0 , 40,0.0,160.0)


# ================================================================================================================
# Perform analysis - loop over all events
# ================================================================================================================
# Define variables across ALL events
N_EVENTS = treeReader.GetEntries()
rinv_sum, dark_events_counter = 0,0 
delta_charge_sum, delta_charge_events = 0,0
percent_pi_sum, percent_rho_sum, percent_delta_sum = 0,0,0 # Only for statistics
rho_diag, rho_off_diag = 0, 0 # Only for statistics
pi_diag, pi_off_diag = 0, 0 # Only for statistics
n_deltabaryon_list = [] 

# Start loop through every event
for event in range(0, N_EVENTS):
    # Load selected branches with data from specified event
    treeReader.ReadEntry(event)
    
    # Define variables for this event used across different analysis
    N_PTCLS = branchPtcl.GetEntries()
    dark_trigger, gen_trigger, rec_trigger = False, False, False # For storing MET differences between methods 
        
    ##############################################
    # Generator level analysis
    ##############################################
    genjet_px, genjet_py = 0,0
    N_JETS = branchGenJet.GetEntries()
    histGenNjets.Fill(N_JETS)
    
    if N_JETS > 1: 
        # Take the two leading jets
        genjet1 = branchGenJet.At(0)
        genjet2 = branchGenJet.At(1)
        dphi = GetdPhi(branchGenJet.At(0).Phi,branchGenJet.At(1).Phi)
            
        if genjet1.PT > pT_min_jet1 and np.abs(genjet1.Eta) < eta_max and genjet2.PT > pT_min_jet2 and np.abs(genjet2.Eta) < eta_max and abs(dphi) < dphi_max:
            gen_trigger = True
            # Define leading and sub-leading jets and MET as vectors
            genjet1_vec, genjet2_vec, genmet_vec = ROOT.TLorentzVector(), ROOT.TLorentzVector(), ROOT.TLorentzVector()
            genjet1_vec = GetJetVector(genjet1.PT , genjet1.Eta , genjet1.Phi , genjet1.Mass)
            genjet2_vec = GetJetVector(genjet2.PT , genjet2.Eta , genjet2.Phi , genjet2.Mass)
            genmet = branchGenMET.At(0)
            genmet_m = genmet.MET
            genmet_x, genmet_y = genmet_m* np.cos(genmet.Phi), genmet_m* np.sin(genmet.Phi)
            genmet_vec.SetPxPyPzE(genmet_x , genmet_y , 0 , genmet_m)

            # Calculate observables
            gen_mt = (genjet1_vec + genjet2_vec + genmet_vec).Mt()
            if (gen_mt!=0): gen_rt = genmet_m/gen_mt
            gendphi1 = GetdPhi(genjet1.Phi, genmet.Phi)
            gendphi2 = GetdPhi(genjet2.Phi , genmet.Phi)
            gendphimin = min(gendphi1, gendphi2)
            gendelEta = GetdEta(genjet1.Eta , genjet2.Eta)

            # Fill histograms
            hist1GenJet1PT.Fill(branchGenJet.At(0).PT)
            hist1GenJet2PT.Fill(branchGenJet.At(1).PT)
            #histGenMET.Fill(branchGenMET.At(0).MET)
            histGenJetmJJ.Fill((genjet1_vec+genjet2_vec).M())
            histGenJetdphi.Fill(GetdPhi(branchGenJet.At(0).Phi,branchGenJet.At(1).Phi))
            histGenmT.Fill(gen_mt)
            histGenratT.Fill(gen_rt)
            histGendphi1.Fill(gendphi1)
            histGendphi2.Fill(gendphi2)
            histGenJetdphimin.Fill(gendphimin)
            histGendelEta.Fill(gendelEta)
            
            histGenInvPtcls1.Fill(GetJetNinv(genjet1, R, darkPID, branchPtcl))
            histGenInvPtcls2.Fill(GetJetNinv(genjet2, R, darkPID, branchPtcl))
            
            # Look at third jet and fill histograms for jet3
            if N_JETS > 2 and branchGenJet.At(2).PT > pT_min_jet1:
                genjet3 = branchGenJet.At(2)
                genjet3_vec = ROOT.TLorentzVector()
                genjet3_vec = GetJetVector(genjet3.PT , genjet3.Eta , genjet3.Phi , genjet3.Mass)
                histGenInvPtcls3.Fill(GetJetNinv(genjet3, R, darkPID, branchPtcl)) 
                hist1GenJet3PT.Fill(genjet3.PT)
                histGenJet13dphi.Fill(GetdPhi(genjet1.Phi,genjet3.Phi))
                histGenJet23dphi.Fill(GetdPhi(genjet2.Phi,genjet3.Phi))
      
    
            upper_ninv_per, lower_ninv_per = GetHemisphereNinv(genjet1, darkPID, branchPtcl)
            histGenInvPtclsUpper.Fill(upper_ninv_per)
            histGenInvPtclsLower.Fill(lower_ninv_per)
            
            upper_pt, lower_pt = GetHemispherePt(genjet1, pT_min_jet1, branchGenJet)
            histGenUpperpT.Fill(upper_pt)
            histGenLowerpT.Fill(lower_pt)               
            
            # Get sum of pT and mass of all jets in event
            for i in range(branchGenJet.GetEntries()):
                current_genjet = branchGenJet.At(i)
                x, y = np.cos(current_genjet.Phi), np.sin(current_genjet.Phi)
                px, py = current_genjet.PT*x, current_genjet.PT*y
                genjet_px = genjet_px + px
                genjet_py = genjet_py + py
                if current_genjet.PT > pT_min_jet1:
                    current_genjet_vec = ROOT.TLorentzVector()
                    current_genjet_vec = GetJetVector(current_genjet.PT , current_genjet.Eta , current_genjet.Phi , current_genjet.Mass)
                    genjet_mass = current_genjet_vec.M()
                    genjet_pt = current_genjet.PT
                    histGenAllJetsMT.Fill(genjet_mass)
                    histGenAllJetsPT.Fill(genjet_pt)

    ###################################################
    # Reconstruction level analysis
    ###################################################
    if branchJet.GetEntries() > 1:
        # Take the two leading jets
        jet1 = branchJet.At(0)
        jet2 = branchJet.At(1)
        dphi = GetdPhi(jet1.Phi,jet2.Phi)

        #Selecting events
        
        #if jet1.PT > pT_min_jet1 and jet2.PT > pT_min_jet2 :    
        if jet1.PT > pT_min_jet1 and np.abs(jet1.Eta) < eta_max and jet2.PT > pT_min_jet2 and np.abs(jet2.Eta) < eta_max and abs(dphi) < dphi_max:
            rec_trigger = True
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
            invmass_softdrop = GetInvMass(jet1_softdrop , jet2_softdrop)
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
            invmass_trimmed = GetInvMass(jet1_trimmed , jet2_trimmed)
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
                if DeltaR < R :
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
            #histMET.Fill(branchMET.At(0).MET)
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

    
    histMET.Fill(branchMET.At(0).MET)
    histGenMET.Fill(branchGenMET.At(0).MET)
    histGenJetMET.Fill(np.sqrt(genjet_px**2+genjet_py**2))
        
        
    ###################################
    # Analysis of dark particles      #
    ###################################
    m_trigger = False                 # Triggered if Z' invariant mass is within 10% of Z' mass
    delta_charge = 0                  # Total charge difference for dark particles
    dark_px, dark_py = 0,0            # To calculate MET originating from dark particles
    allmesons, stablemesons = 0,0     # To calculate rinv
    npi, npi_unstable = 0,0           # Remaining definitions are for keeping track of number of particles of different kinds
    nrho, nrho_unstable = 0,0         # ------------
    ndelta = 0                        # ------------
    ninv = 0                          # ------------
    nrho_to_SM, nrho_to_DS = 0,0      # ------------
#    fastjetparticles = []
#    jets = []

#     for iptcl in range(N_PTCLS):
#         #jet = jet_def(ROOT.PseudoJet(branchPtcl.At(iptcl).Px, branchPtcl.At(iptcl).Py, branchPtcl.At(iptcl).Pz, branchPtcl.At(iptcl).E))
#         #jet = jetFinder(branchPtcl.At(iptcl))

#         pid = branchPtcl.At(iptcl).PID
#         status = branchPtcl.At(iptcl).Status
#         if status !=1:
#             continue 
#         if pid in invPID:
#             continue
#         px, py, pz, e = branchPtcl.At(iptcl).Px, branchPtcl.At(iptcl).Py, branchPtcl.At(iptcl).Pz, branchPtcl.At(iptcl).E
#         particle_vec = ROOT.TLorentzVector()
#         particle_vec.SetPxPyPzE(px, py, pz, e)
#         fastjetparticles.append(particle_vec)
        
        
#     jetFinder.InputArray = fastjetparticles
#     jetFinder.OutputArray = jets
#     jetFinder = ROOT.FastJetFinder()
    
#     for i in range(len(jets)):
#         print(i)
#         if abs(jets[i].eta()) < eta_max and jets[i].pt() > pT_min_jet1:
#             particle = part(0, 0, 0, jets[i].pt(), jets[i].eta(), jets[i].phi(), 0)
#             print('Yes')
    
    # Check for Z-prime boson in event, find daughters and get invariant mass - triggers m_trigger if within mass range
    # m_trigger needed for 2nd loop over N_PTCLS
    for iptcl in range(N_PTCLS):
        if branchPtcl.At(iptcl).PID == 4900023:
            daut1 = branchPtcl.At(iptcl).D1
            daut2 = branchPtcl.At(iptcl).D2
            if branchPtcl.At(daut1).PID in quarkPID or branchPtcl.At(daut2).PID in quarkPID:
                daut1_vec, daut2_vec = ROOT.TLorentzVector(), ROOT.TLorentzVector()
                px1, py1, pz1, e1 = branchPtcl.At(daut1).Px, branchPtcl.At(daut1).Py, branchPtcl.At(daut1).Pz, branchPtcl.At(daut1).E
                px2, py2, pz2, e2 = branchPtcl.At(daut2).Px, branchPtcl.At(daut2).Py, branchPtcl.At(daut2).Pz, branchPtcl.At(daut2).E
                daut1_vec.SetPxPyPzE(px1, py1, pz1, e1)
                daut2_vec.SetPxPyPzE(px2, py2, pz2, e2)     
                mz_inv = GetInvMass(daut1_vec, daut2_vec)

                
                
                # Get dphi of the two dark quarks - use GetLastCopy to use the quarks right before they fragment/hadronize
                daut1 = GetLastCopy(daut1,branchPtcl)
                daut2 = GetLastCopy(daut2,branchPtcl)
                qdark_dphi = GetdPhi(branchPtcl.At(daut1).Phi,branchPtcl.At(daut2).Phi)
                
                # Fill the following histograms 
                histGenZmJJ.Fill(mz_inv)
                # Only fill the histograms if the events contain 2 leading jets fulfilling pT, eta and dphi requirements
                if gen_trigger == True:
                    histGenZpT.Fill(branchPtcl.At(iptcl).PT)
                    histGenqDarkdPhi.Fill(qdark_dphi)
                
                dark_events_counter = dark_events_counter + 1
                # Trigger the m_trigger in this event only if mz_inv is within 10% of the Z' mass
                if 0.9*mz_model < mz_inv < 1.1*mz_model: m_trigger = True
    
    
    # Loop through every particle in event for dark particle analysis
    for iptcl in range(N_PTCLS):
        # Skip this particle if it is not a dark meson/baryon
        if branchPtcl.At(iptcl).PID not in darkPID:
            continue 
        
        pid = branchPtcl.At(iptcl).PID
        status = branchPtcl.At(iptcl).Status
        dark_trigger = True
        
        # Warn user if the dark particle status is not as expected - add statuses as needed
        if status not in [1,83,84,91]:
            print('Unexpected status of dark particle. Status is %d.'%(status))
       
        # Going through decayed dark mesons and counting allmesons, stablemesons and whether meson is diagonal or not
        elif status in [83,84,91]:
            if pid in piPID:
                allmesons = allmesons + 1
                pi_diag, pi_off_diag = CountDiagonality(pid, pi_diag, pi_off_diag)
                daut1, daut2 = branchPtcl.At(iptcl).D1, branchPtcl.At(iptcl).D2
                # For using branching ratios to determine number of stable/unstable pions PID 53 is decay product of stable pions
                if any(branchPtcl.At(daut1).PID == 53, branchPtcl.At(daut2).PID == 53):
                    stablemesons = stablemesons + 1
                else:
                    npi_unstable = npi_unstable + 1
            elif pid in rhoPID:
                rho_diag, rho_off_diag = CountDiagonality(pid, rho_diag, rho_off_diag)
                nrho_unstable = nrho_unstable + 1
                # Checks if rho has decayed only within the dark sector as it will then not count towards total number of dark hadrons
                daut_pids = GetDaughtersPID(branchPtcl.At(iptcl),branchPtcl)
                if all(item in darkPID for item in daut_pids):
                    nrho_to_DS = nrho_to_DS + 1
                elif any(item in darkPID for item in daut_pids):
                    allmesons = allmesons + 1
                    #nrho_to_DS = nrho_to_DS + 0.5
                    #nrho_to_SM = nrho_to_SM + 0.5
                else:
                    allmesons = allmesons + 1
                    nrho_to_SM = nrho_to_SM + 1
                
        # Going through final-state dark mesons and counting allmesons and whether meson is diagonal or not. Fill PID histogram
        elif branchPtcl.At(iptcl).Status == 1:
            allmesons = allmesons + 1
            stablemesons = stablemesons + 1
            dark_px = dark_px + branchPtcl.At(iptcl).Px 
            dark_py = dark_py + branchPtcl.At(iptcl).Py
            if branchPtcl.At(iptcl).PID in piPID:
                npi = npi + 1
                pi_diag, pi_off_diag = CountDiagonality(pid, pi_diag, pi_off_diag)
            elif branchPtcl.At(iptcl).PID in rhoPID:
                nrho = nrho + 1
                rho_diag, rho_off_diag = CountDiagonality(pid, rho_diag, rho_off_diag)
            elif branchPtcl.At(iptcl).PID in deltaPID:
                ndelta = ndelta + 1
            if m_trigger == True:
                if pid > 0:
                    short_pid = pid - 4900000
                    row, column = int(str(pid)[-3:-2]), int(str(pid)[-2:-1])
                    if row != column: delta_charge = delta_charge + 1
                else:
                    short_pid = pid + 4900000
                    row, column = int(str(pid)[-3:-2]), int(str(pid)[-2:-1])
                    if row != column: delta_charge = delta_charge - 1
                histGenPID.Fill(short_pid)
    
    # Calculate delta_charge only if m_trigger was triggered in event and sum delta_charge. Fill delta_charge histogram
    if m_trigger == True:
        delta_charge_sum = delta_charge_sum + delta_charge # only for statistics purposes
        delta_charge_events = delta_charge_events + 1 # only for statistics purposes
        histGenDeltaCharge.Fill(delta_charge)
    
    if (allmesons!=0): rinv = stablemesons/allmesons
    if rinv > 1:
        raise ValueError('r_inv > 1! r_inv=%0.3f'%(rinv))
    rinv_sum = rinv_sum + rinv
    
    histGenrinv.Fill(rinv)
    histGenPions.Fill(npi)
    histGenRhos.Fill(nrho)
        
        
        
    # Will store MET calculated directly from dark particles and difference between reconstructed and generated MET only if all triggers were triggered in event
    if dark_trigger == True and gen_trigger == True and rec_trigger == True:
        histGenDarkMET.Fill(np.sqrt(dark_px**2 + dark_py**2))
        MET_diff_gen_dark = branchGenMET.At(0).MET - np.sqrt(dark_px**2 + dark_py**2)
        MET_diff_gen_rec = branchGenMET.At(0).MET - branchMET.At(0).MET
        histMETdiffGenDark.Fill(MET_diff_gen_dark)
        histMETdiffGenRec.Fill(MET_diff_gen_rec)
        
###################################################
# Statistics and printing statements
###################################################
    # To do within events
    if stablemesons != 0:    
        percent_pi = npi/stablemesons
        percent_rho = nrho/stablemesons
        percent_delta = ndelta/stablemesons
        percent_pi_sum = percent_pi_sum + percent_pi
        percent_rho_sum = percent_rho_sum + percent_rho
        percent_delta_sum = percent_delta_sum + percent_delta
        n_deltabaryon_list.append(percent_delta)
        
    # Print some numbers at every 1000 events - comment out if not needed
    if event%1000 == 0:
        print('Event %i: Unstable pions: %i. Unstable rhos to SM: %i. Unstable rhos to DS: %i.'%(event, npi_unstable, nrho_to_SM, nrho_to_DS))
        print('             Stable mesons: %i. Total dark mesons: %i. r_inv = %0.3f'%(stablemesons,allmesons,rinv))
        print('             Final-state percentage of: Pions %i%% Rhos %i%% Deltas %i%%'%(percent_pi*100, percent_rho*100, percent_delta*100))
        if rho_diag !=0: print('             Ratio of off-diagonal to diagonal rhos is %0.3f'%(rho_off_diag/rho_diag))
        if pi_diag !=0: print('             Ratio of off-diagonal to diagonal pions is %0.3f'%(pi_off_diag/pi_diag))
        
        
# After all events have been looped through
# Print means for true mean
rinv_mean = rinv_sum/(dark_events_counter)
delta_charge_mean = delta_charge_sum/(delta_charge_events)
percent_pi_mean = percent_pi_sum/(dark_events_counter)
percent_rho_mean = percent_rho_sum/(dark_events_counter)
percent_delta_mean = percent_delta_sum/(dark_events_counter)
print('r_inv = %0.3f   delta_charge = %0.3f'%(rinv_mean,delta_charge_mean))
#print('Pions = %0.3f   Rhos = %0.3f   Deltas = %0.3f'%(percent_pi_mean, percent_rho_mean, percent_delta_mean))
if percent_delta_sum != 0:
    percent_delta_mean = percent_delta_sum/(dark_events_counter)
    deltabaryon_stddev_innersum = 0 
    deltabaryon_stddev = 0
    for perc in n_deltabaryon_list:
        #print(perc)
        deltabaryon_stddev_innersum = deltabaryon_stddev_innersum + (perc - percent_delta_mean)**2
    deltabaryon_stddev = np.sqrt(1/(len(n_deltabaryon_list))*deltabaryon_stddev_innersum)
    print('Pions = %0.5f   Rhos = %0.5f   Deltas = %0.5f \u00B1 %0.5f'%(percent_pi_mean, percent_rho_mean, percent_delta_mean, deltabaryon_stddev))
    pion_terror = (percent_pi_mean - 0.225)/(0.225)*100
    rho_terror = (percent_rho_mean - 0.675)/(0.675)*100
    delta_terror = (percent_delta_mean - 0.1)/(0.1)*100
    perc_sum = percent_pi_mean + percent_rho_mean + percent_delta_mean
    print('Pion error = %0.2f%%   Rho error = %0.2f%%   Delta error = %0.2f%%   Check sum = %0.2f'%(pion_terror, rho_terror, delta_terror, perc_sum))
    
print('Ratio of off-diagonal to diagonal rhos is %0.3f'%(rho_off_diag/rho_diag))
n_rhos = rho_diag + rho_off_diag
n_pis = pi_diag + pi_off_diag
print('Ratio of rhos produced is %0.3f'%(n_rhos/(n_rhos+n_pis)))


# ================================================================================================================
# Normalising histograms
# ================================================================================================================
if histGenrinv.GetSumw2N()==0: histGenrinv.Sumw2(True)
if histGenPID.GetSumw2N()==0: histGenPID.Sumw2(True)
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
histlist.Add(hist1GenJet3PT)
histlist.Add(histGenMET)
histlist.Add(histGenDarkMET)
histlist.Add(histGenAllJetsMT)
histlist.Add(histGenAllJetsPT)
histlist.Add(histGenUpperpT)
histlist.Add(histGenLowerpT)
histlist.Add(histGenNjets)
histlist.Add(histGenJetMET)
histlist.Add(histGenPions)
histlist.Add(histGenRhos)
histlist.Add(histGenInvPtcls1)
histlist.Add(histGenInvPtcls2)
histlist.Add(histGenInvPtcls3)
histlist.Add(histGenInvPtclsUpper)
histlist.Add(histGenInvPtclsLower)
histlist.Add(histGenrinv)
histlist.Add(histGenDeltaCharge)
histlist.Add(histGenPID)
histlist.Add(histGenZmJJ)
histlist.Add(histGenZpT)
histlist.Add(histGenqDarkdPhi)
histlist.Add(histGenJetmJJ)
histlist.Add(histGenJetdphi)
histlist.Add(histGenJetdphimin)
histlist.Add(histGenJet13dphi)
histlist.Add(histGenJet23dphi)
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
#histlist.Add(histtrackpt) # 2D histogram
histlist.Add(histdelEta)
histlist.Add(histMETdiffGenDark)
histlist.Add(histMETdiffGenRec)


# ================================================================================================================
# Final steps
# ================================================================================================================
#outputFile = inputFile[:-5] + "_combined_R14.root"
rootFile = ROOT.TFile(outputFile, "RECREATE")
histlist.Write()
rootFile.Close()

