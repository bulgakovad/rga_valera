#----------------------------------------Cell 1 -------------------------------------------------------
import sys
if len(sys.argv) > 1: method = sys.argv[1]
from ROOT import gRandom, TH1, TH1D, TCanvas, cout
import ROOT

import numpy as np
import pandas as pd
from array import array
import math 



intial_path = "/w/hallb-scshelf2102/clas12/bulgakov/projects/rga_valera/analysis/"

#old data without Valerii suggestion for continue
fData = ROOT.TFile.Open("/lustre24/expphy/volatile/clas12/bulgakov/rga_valera/analysis_output/Data_xsec_MCcorr_SMfix_Jul22.root","READ")
fSim = ROOT.TFile.Open("/lustre24/expphy/volatile/clas12/bulgakov/rga_valera/analysis_output/SIM_xsec_MCcorr_SMfix_Jul22.root","READ")

#fData = ROOT.TFile.Open("/w/hallb-scshelf2102/clas12/bulgakov/projects/rga_valera/analysis/Data_xsec_MCcorr_SMfix_Jul22.root","READ")
#fSim = ROOT.TFile.Open("/w/hallb-scshelf2102/clas12/bulgakov/projects/rga_valera/analysis/SIM_xsec_MCcorr_SMfix_Jul22.root","READ")


def lumi(charge):
    RD=57.1 
    qe    = 1.602177E-19     # Electron charge ---------- Coulomb 
    rho   = 0.0701          # Density of H2 at 20K ----- g/cm3
    A0    = 6.0221367E23     # Avogado number ----------- mol^-1
    MH    = 1.00794          # Atomic mass of hydrogen -- amu =g/mol
    LT    = 5.0              # Lenght of target --------- cm
    CMB   = 1e30             # cm^2 to microbarn --------mcbarn/cm^2
    chargeConvered = charge/1e9            # nanofarad to farad
    np           = LT*rho*A0/MH     # number of target nuclei per cm2 (Density of target )
    ne           = chargeConvered/qe       # Number of electron hitting the target

    factor = (ne*np)/CMB 
    return factor 

# q -bin number 0-20 
def normalize(histIN, histOUT, q, charge):
    
    deltaW = 0.05 
    e =  2.718 
    
    # bin Q2 width
    nQ2 = 50 
    q2Min = 1 
    q2Max = 2500 
    deltaQ2 = np.log(q2Max/q2Min)/nQ2 
    
    adjust = lumi(charge)/6.0 
    
    q2 =  (q2Min*math.exp( q*deltaQ2) + q2Min* math.exp((q + 1) *deltaQ2))/2 
    binQ2 = q2Min* math.exp( (q + 1) *deltaQ2) - q2Min* math.exp( q*deltaQ2) 
    
    for t in range(1, 1 + histIN.GetNbinsX(),1):
        histOUT.SetBinContent(t, histIN.GetBinContent(t)/deltaW/binQ2/adjust) 
        histOUT.SetBinError(t, histIN.GetBinError(t)/deltaW/binQ2/adjust) 
        
def ApplyCorr(histIN, histOUT, q, df):
    
    if (q < 6 or q > 14):
        print ('no corrections for the hist out of Q2 range')
        for iXbin in range(1, 1 + histIN.GetNbinsX(),1):
            histOUT.SetBinContent(iXbin, 0.)
            histOUT.SetBinError(iXbin, 0.)
        return 
    
    nQ2 = 50 
    q2Min = 1 
    q2Max = 2500 
    deltaQ2 = np.log(q2Max/q2Min)/nQ2 
    
    # center of the bin
    q2 =  (q2Min* math.exp( q*deltaQ2) + q2Min* math.exp( (q + 1) *deltaQ2))/2
    #nW bins:
    nW = df.shape[0]/9
    
    q2binsDF = pd.unique(df['Q2'])
    df = df[df['Q2'] == q2binsDF[q-6]]
    
    if (abs(q2 - q2binsDF[q-6]) >  0.2):
        print ('Q2 range Error')
        return
    
    for iXbin in range(1, 1 + histIN.GetNbinsX(),1):

        if (histIN.GetBinCenter(iXbin)<1.14  or histIN.GetBinCenter(iXbin) > 2.5):
            continue
        corrFound = False
        for index, row in df.iterrows():
            if (abs(row['W'] - histIN.GetBinCenter(iXbin)) < 0.01):
                corrFound = True
                corr = row['Combination']
                evYield = histIN.GetBinContent(iXbin)
                histOUT.SetBinContent(iXbin, evYield * corr)
                histOUT.SetBinError(iXbin, histIN.GetBinError(iXbin) * corr)
                break
                
        if (not corrFound):
            
            if (histIN.GetBinCenter(iXbin)>2.45):
                print(row['W'], histIN.GetBinCenter(iXbin))
            histOUT.SetBinContent(iXbin, 0.)
            histOUT.SetBinError(iXbin, 0.)

def GetCorrFactorFCC(curr):
    if (curr > 0):
        return (1 - (curr - 45.)* 0.0041)
    else:
        return 1.

fcc_55 = 4.67104e+06 * GetCorrFactorFCC(55.)
fcc_50 = 3.02021e+06 * GetCorrFactorFCC(50.)
fcc_52 = 365707 * GetCorrFactorFCC(52.)
fcc_53 = 198954 * GetCorrFactorFCC(53.)
fcc_54 = 363850 * GetCorrFactorFCC(54.)

# 45nA runs
fcc_data = 5.85024e+06 * 9.8088 / 10.61360 + 6.00024e+06  + 9.37112e+06
fcc_data +=(fcc_55 + fcc_50 + fcc_52 + fcc_53 + fcc_54)
fccConst = fcc_data
print(fccConst)
print("end cell 1")

#----------------------------------------End Cell 1 -------------------------------------------------------
#----------------------------------------Cell 2 -------------------------------------------------------
#2D unfolding rooUnfold
def getResponse2D(fileSim, MSL, TSML, RSL, ORSL):
    for sec in range(6):
        fold = 'Deconv2D/S' + str(sec) + '/'
        endR = '_my_secBySec_S' + str(sec)
        endR2 = '_secBySec_S' + str(sec)
        
        MSL.append(fileSim.Get(fold + 'Sep_measured2D'+ endR))
        nBinsX = MSL[-1].GetNbinsX()
        TSML.append(fileSim.Get(fold + 'Sep_truth2D'+ endR))
        RSL.append(fileSim.Get(fold + 'R2D'+'_S'+str(sec)))      
        ORSL.append(ROOT.RooUnfoldResponse(MSL[-1],TSML[-1], RSL[-1]))
        
#2D to 1D
def AntiLinear2D(histIN_1, listHistsOUT, nameEnd = ''):
    nQ2bins = 11
    
    wMin_fold = histIN_1.GetXaxis().GetXmin()
    wMax_fold = histIN_1.GetXaxis().GetXmax()
    nBinsX = histIN_1.GetNbinsX()
    
    for iQ2bin in range(nQ2bins):

        histN = 'wUnf_Q2_' + str(iQ2bin + 1)
        listHistsOUT.append( ROOT.TH1D(histN + nameEnd,histN + nameEnd, nBinsX, wMin_fold, wMax_fold))
        for iBin in range(1 , 1 + nBinsX,1):
            shift = 0
            listHistsOUT[-1].SetBinContent(iBin, histIN_1.GetBinContent(iBin, iQ2bin+1))
            listHistsOUT[-1].SetBinError(iBin, histIN_1.GetBinError(iBin, iQ2bin+1))
            
            
        
# GetXSEC using 2D rooUnfold methods        
def GetXSEC_roo2D(fileSim, fileData, outXSEC, addN, chargeFC_tmp, regParam, method = 'bayes'):
    # regParam is nIter for bayes or i_comp for SVD
    if (method!= 'bayes' and method!= 'bbb' and method!= 'svd'):
        print('unsupported deconvolution method is used, please put bayes or bbb or svd')
        return
    
    #sectors
    nSec = 6
    #number of Q2 bins
    nQ2bins = 11
    
    # initialization of internal rooUnfold matrices and struct.
    tmp_measuredL = []
    tmp_truthL = []
    tmp_respMatL = []
    tmp_respObjL = []

    # read root file Sim
    getResponse2D(fileSim, tmp_measuredL, tmp_truthL, tmp_respMatL, tmp_respObjL)

    # measured from Data
    tmp_hMeas = []

    # read root file Data
    for sec in range(nSec):
        fold = 'Deconv2D/S' + str(sec) + '_my/'
        tmp_hMeas.append(fileData.Get(fold + 'measured2D_DATA_S' + str(sec)))
        

    #RooUnfoldObjs
    tmp_unfold = []
    tmp_hUnfold = []

        
    for sec in range(nSec):
        
        # method selection
        if (method == 'bayes'):
            tmp_unfold.append(ROOT.RooUnfoldBayes(tmp_respObjL[sec],tmp_hMeas[sec], regParam))
            
        if (method == 'bbb'):            
            tmp_unfold.append(ROOT.RooUnfoldBinByBin(tmp_respObjL[sec],tmp_hMeas[sec]))
        if (method == 'svd'):
            tmp_unfold.append(ROOT.RooUnfoldSvd(tmp_respObjL[sec],tmp_hMeas[sec], regParam))
            
        #####################################
        # errors accounting
        tmp_unfold[-1].IncludeSystematics(0)
        tmp_unfold[-1].SetOverflow(0)
        #tmp_unfold[-1].SetVerbose(2)
        ######################################
        
        # Actual deconvolution
        tmp_hUnfold.append(tmp_unfold[-1].Hunfold())
                   
        
    # Convert 2D hist into list of 1D hists
    wYield1D = []
    
    for sec in range(nSec): 
        wYield1D.append([])
        AntiLinear2D(tmp_hUnfold[sec], wYield1D[-1], addN + 'AL' + str(sec)) 
    
    # read file with all the corrections factors
    #corrDF = pd.read_csv('/home/valerii/Clas12/Inclusive/xsecs/rc_bcc_et_corr.dat') Not need?
    
    #corrDF = pd.read_csv('/home/valerii/clas12/Inclusive/xsecs/rc_bcc_et_corr_targetFix.dat') COMMENTED OUT FOR NOW

    
    
    # W yield corrected for Lumin. 
    yieldNorm = []

    for iQ2 in range(nQ2bins):

        yieldNorm.append([])
        outXSEC.append([])

        for sec in range(nSec):
            # Q2 convert binning to match binning in ana12.
            q = iQ2 + 5

            WcurrentBin = wYield1D[sec][iQ2]
            yieldNorm[-1].append(ROOT.TH1D('nrY' + str(iQ2)+ ' ' + str(sec)+ addN ,
                                           'nrY' + str(iQ2)+ ' ' + str(sec)+ addN,
                                           WcurrentBin.GetNbinsX(),
                                           WcurrentBin.GetXaxis().GetXmin(), 
                                           WcurrentBin.GetXaxis().GetXmax()))

            outXSEC[-1].append(ROOT.TH1D('xsec' + str(iQ2)+ ' ' + str(sec) + addN,
                                         'xsec' + str(iQ2)+ ' ' + str(sec)+ addN,
                                         WcurrentBin.GetNbinsX(),
                                         WcurrentBin.GetXaxis().GetXmin(), 
                                         WcurrentBin.GetXaxis().GetXmax()))

            normalize(WcurrentBin, yieldNorm[-1][-1], q, chargeFC_tmp)
            
            #ApplyCorr(yieldNorm[-1][-1], outXSEC[-1][-1], q, corrDF)     COMMENTED OUT FOR NOW




    out_file = ROOT.TFile(addN + "_unfolded_output.root", "RECREATE")

    # Save unfolded 2D hists
    out_file.mkdir("Unfold2D")
    out_file.cd("Unfold2D")
    for sec in range(nSec):
        h2 = tmp_hUnfold[sec]
        h2.SetName(f"hUnfold2D_sec{sec}")
        h2.Write()

    # Save 1D W-yield hists
    out_file.mkdir("WYield1D")
    out_file.cd("WYield1D")
    for sec in range(nSec):
        out_file.mkdir(f"WYield1D/sec{sec}")
        out_file.cd(f"WYield1D/sec{sec}")
        for iQ2 in range(nQ2bins):
            h1 = wYield1D[sec][iQ2]
            h1.SetName(f"wYield1D_sec{sec}_Q2bin{iQ2}")
            h1.Write()

    out_file.Close()
    print(f"Saved output ROOT file: {addN}_unfolded_output.root")
    
def IntOverSec(listIn, listOut, endN, isTenBin = False):
    nQ2bins = 11
    if (isTenBin):
        nQ2bins = 10
    for iQ2 in range(nQ2bins):
        listOut.append(ROOT.TH1D('intCutFunc' + str(iQ2) + endN,
                                     'intCutFunc' + str(iQ2) + endN,
                                     listIn[0][0].GetNbinsX(),
                                     listIn[0][0].GetXaxis().GetXmin(),
                                     listIn[0][0].GetXaxis().GetXmax()))
        
        for iXbin in range(1, 1 + listIn[0][0].GetNbinsX(),1):
            value = 0
            for sec in range(6):
                value += listIn[iQ2][sec].GetBinContent(iXbin)
            listOut[-1].SetBinContent(iXbin, value/6)

def GetTitle(npad):
    nQ2 = 50;
    q2Min = 1;
    q2Max = 2500;
    deltaQ2 = np.log(q2Max/q2Min)/nQ2;
    Q2bin_min = q2Min*math.exp((npad + 5)*deltaQ2)
    Q2bin_max = q2Min*math.exp( (npad + 6)*deltaQ2)
    return  "             " + str(Q2bin_min)[0:4] + " < Q^{2} < " +  str(Q2bin_max)[0:4] + " GeV^{2}"
print("end cell 2")
#----------------------------------------End Cell 2 -------------------------------------------------------
#----------------------------------------Cell 3 -----------------------------------------------------------
epochs = [1,2,3]

xsec_allS_roo2D_bayes_1 = []
xsec_allS_roo2D_bayes_2 = []
xsec_allS_roo2D_bayes_3 = []

GetXSEC_roo2D(fSim, fData, xsec_allS_roo2D_bayes_2, 'bayes2_2D_2', fccConst, epochs[1], method = 'bayes')

xsec_integrated_2D_bayes_1 = []
xsec_integrated_2D_bayes_2 = []
xsec_integrated_2D_bayes_3 = []

IntOverSec(xsec_allS_roo2D_bayes_2, xsec_integrated_2D_bayes_2, 'xsec_integrated_2D_bayes_2')

"""
xsec_allS_roo2D_bayes_2_lastQ2 = []
GetXSEC_roo2D(fSim_lastQ2, fData_lastQ2, xsec_allS_roo2D_bayes_2_lastQ2, 'bayes2_2D_2_lastQ2', fccConst, epochs[1], method = 'bayes')
xsec_integrated_2D_bayes_2_lastQ2 = []
IntOverSec(xsec_allS_roo2D_bayes_2_lastQ2, xsec_integrated_2D_bayes_2_lastQ2, 'xsec_integrated_2D_bayes_2_lastQ2')
"""

print("end cell 3")
#----------------------------------------End Cell 3 -------------------------------------------------------
#----------------------------------------Cell 4 -------------------------------------------------------
#
##statDF = pd.read_csv('/home/valerii/Clas12/Inclusive/fromLaptop/Unfolding/XSEC_final/statErrF.dat')
#statDF = pd.read_csv('stat_combL.dat')
#
#def AddSysUNC(histIN_L, df):
#    
#    for iQ2 in range(1,10,1):
#        
#        histIN = histIN_L[iQ2]
#        
#        # convert binning
#        q = iQ2 + 5
#
#        nQ2 = 50 
#        q2Min = 1 
#        q2Max = 2500 
#        deltaQ2 = np.log(q2Max/q2Min)/nQ2 
#        # center of the bin
#        q2 =  (q2Min* math.exp( q*deltaQ2) + q2Min* math.exp( (q + 1) *deltaQ2))/2
#
#        q2binsDF = pd.unique(df['Q2'])
#        #print(q-6, q2binsDF)
#        oneQ2_df = df[df['Q2'] == q2binsDF[q-6]]
#
#        if (abs(q2 - q2binsDF[q-6]) >  0.2):
#            print ('Q2 range Error')
#            return
#
#        for iXbin in range(1, 1 + histIN.GetNbinsX(),1):
#            corrFound = False
#            for index, row in oneQ2_df.iterrows():
#                if (abs(row['W'] - histIN.GetBinCenter(iXbin)) < 0.01):
#                    corrFound = True
#                    #print(row['W'], histIN.GetBinCenter(iXbin))
#                    corr = row['stat_combL']
#                    evYield = histIN.GetBinContent(iXbin)
#                    histIN.SetBinError(iXbin, evYield * corr/100)
#                    break
#
#            if (not corrFound):
#                histIN.SetBinContent(iXbin, 0.)
#                histIN.SetBinError(iXbin, 0.)
#
#AddSysUNC(xsec_integrated_2D_bayes_2, statDF)
##AddSysUNC(xsec_integrated_2D_bayes_2_lastQ2, statDF)
#print("end cell 4")
#----------------------------------------End Cell 4 ---- NO file so far---------------------------------------------------