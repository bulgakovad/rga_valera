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
fData = ROOT.TFile.Open("Data_xsec_MCcorr_SMfix_Jul22.root","READ")
fSim = ROOT.TFile.Open("SIM_xsec_MCcorr_SMfix_Jul22.root","READ")

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
# Linearized 1D unfolding, but now reading the SAME 2D folders/objects as Valerii's original code:
#   SIM : Deconv2D/S{sec}/Sep_measured2D_my_secBySec_S{sec}
#         Deconv2D/S{sec}/Sep_truth2D_my_secBySec_S{sec}
#         Deconv2D/S{sec}/R2D_S{sec}
#   DATA: Deconv2D/S{sec}_my/measured2D_DATA_S{sec}
#
# We do NOT do TH2 RooUnfold directly.
# Instead:
#   2D hist  --> correctly linearized 1D hist
#   2D resp  --> correctly linearized 2D response matrix
#   then 1D RooUnfold on the linearized objects

LIN_NSEC = 6
LIN_NQ2 = 11
LIN_WMIN = 0.875
LIN_WMAX = 2.875

def _must_get(root_file, key):
    obj = root_file.Get(key)
    if not obj:
        raise RuntimeError(f"Cannot find object: {key}")
    return obj

def _clone_detached(h, new_name):
    hc = h.Clone(new_name)
    hc.SetDirectory(0)
    return hc

def _infer_2d_binning_from_hist(h2, expected_nq2=LIN_NQ2):
    nW = h2.GetNbinsX()
    nQ2 = h2.GetNbinsY()
    if nQ2 != expected_nq2:
        raise RuntimeError(
            f"Unexpected number of Q2 bins in 2D hist: got {nQ2}, expected {expected_nq2}"
        )
    return nW, nQ2

def _global_bin_1based(iW_0based, iQ2_0based, nW):
    # linearization convention:
    # global_bin_0based = iQ2 * nW + iW
    return iQ2_0based * nW + iW_0based + 1

def LinearizeTH2(hist2, new_name):
    """
    Convert TH2(W,Q2) into TH1(global_bin), where global bins run over all Q2 blocks.
    X axis of TH2 is assumed to be W, Y axis Q2.
    """
    nW, nQ2 = _infer_2d_binning_from_hist(hist2)
    nGlobal = nW * nQ2

    h1 = ROOT.TH1D(new_name, new_name, nGlobal, -0.5, nGlobal - 0.5)
    h1.SetDirectory(0)

    for iQ2 in range(nQ2):
        for iW in range(nW):
            gb = _global_bin_1based(iW, iQ2, nW)
            h1.SetBinContent(gb, hist2.GetBinContent(iW + 1, iQ2 + 1))
            h1.SetBinError(gb,   hist2.GetBinError(iW + 1, iQ2 + 1))

    return h1

def LinearizeResponseTH2(resp2, new_name):
    """
    Convert TH2 response matrix in flattened-bin coordinates into another detached TH2.
    This is mainly a safe detached clone, but we verify dimensions.
    In the original files, R2D_S{sec} is already in global-bin coordinates.
    """
    h = _clone_detached(resp2, new_name)
    return h

def UnlinearizeWQ2(histIN_lin, listHistsOUT, nameEnd='', nQ2bins=LIN_NQ2, wMin=LIN_WMIN, wMax=LIN_WMAX):
    nGlobal = histIN_lin.GetNbinsX()
    if nGlobal % nQ2bins != 0:
        raise RuntimeError(
            f"Linearized histogram has {nGlobal} bins, which is not divisible by nQ2bins={nQ2bins}"
        )
    nW = nGlobal // nQ2bins

    for iQ2 in range(nQ2bins):
        histN = f'wUnf_Q2_{iQ2+1}'
        h_out = ROOT.TH1D(histN + nameEnd, histN + nameEnd, nW, wMin, wMax)
        h_out.SetDirectory(0)

        for iW in range(nW):
            gb = _global_bin_1based(iW, iQ2, nW)
            h_out.SetBinContent(iW + 1, histIN_lin.GetBinContent(gb))
            h_out.SetBinError(iW + 1,   histIN_lin.GetBinError(gb))

        listHistsOUT.append(h_out)

def getResponseFromOriginal2D_asLinearized1D(fileSim, MSL, TSML, RSL, ORSL):
    """
    Read the SAME objects as Valerii's original 2D code,
    but convert measured/truth 2D histograms into linearized 1D histograms.
    The response matrix R2D_S{sec} is already in flattened global-bin coordinates.
    """
    for sec in range(LIN_NSEC):
        fold = f'Deconv2D/S{sec}/'
        endR = f'_my_secBySec_S{sec}'

        h_meas2 = _must_get(fileSim, fold + 'Sep_measured2D' + endR)
        h_truth2 = _must_get(fileSim, fold + 'Sep_truth2D' + endR)
        h_resp2 = _must_get(fileSim, fold + 'R2D_S' + str(sec))

        h_meas1 = LinearizeTH2(h_meas2,  f'lin_measured_sec{sec}')
        h_truth1 = LinearizeTH2(h_truth2, f'lin_truth_sec{sec}')
        h_resp_lin = LinearizeResponseTH2(h_resp2, f'lin_response_sec{sec}')

        MSL.append(h_meas1)
        TSML.append(h_truth1)
        RSL.append(h_resp_lin)
        ORSL.append(ROOT.RooUnfoldResponse(MSL[-1], TSML[-1], RSL[-1]))

def GetXSEC_rooLIN1D_fromOriginal2D(fileSim, fileData, outXSEC, addN, chargeFC_tmp, regParam, method='bayes', use_corr=True):
    """
    Same overall downstream logic as before,
    but input comes from Valerii's original 2D folders/objects.
    """
    if method not in ['bayes', 'bbb', 'svd']:
        print('unsupported deconvolution method is used, please put bayes or bbb or svd')
        return

    nSec = LIN_NSEC
    nQ2bins = LIN_NQ2

    # MC response from original 2D folders
    tmp_measuredL = []
    tmp_truthL = []
    tmp_respMatL = []
    tmp_respObjL = []
    getResponseFromOriginal2D_asLinearized1D(fileSim, tmp_measuredL, tmp_truthL, tmp_respMatL, tmp_respObjL)

    # measured DATA from original 2D folders, then linearize
    tmp_hMeas = []
    for sec in range(nSec):
        fold = f'Deconv2D/S{sec}_my/'
        h_data2 = _must_get(fileData, fold + f'measured2D_DATA_S{sec}')
        h_data1 = LinearizeTH2(h_data2, f'lin_meas_data_sec{sec}_{addN}')
        tmp_hMeas.append(h_data1)

    # actual unfolding
    tmp_unfold = []
    tmp_hUnfold_lin = []

    for sec in range(nSec):
        if method == 'bayes':
            unf = ROOT.RooUnfoldBayes(tmp_respObjL[sec], tmp_hMeas[sec], regParam)
        elif method == 'bbb':
            unf = ROOT.RooUnfoldBinByBin(tmp_respObjL[sec], tmp_hMeas[sec])
        elif method == 'svd':
            unf = ROOT.RooUnfoldSvd(tmp_respObjL[sec], tmp_hMeas[sec], regParam)

        unf.IncludeSystematics(0)
        unf.SetOverflow(0)

        h_unf_lin = unf.Hunfold()
        h_unf_lin = _clone_detached(h_unf_lin, f'hUnfoldLIN1D_sec{sec}_{addN}')

        tmp_unfold.append(unf)
        tmp_hUnfold_lin.append(h_unf_lin)

    # unlinearize back into W histograms per Q2 bin
    wYield1D = []
    for sec in range(nSec):
        wYield1D.append([])
        UnlinearizeWQ2(tmp_hUnfold_lin[sec], wYield1D[-1], nameEnd=f'_{addN}_sec{sec}')

    # correction table
    corrDF = pd.read_csv('rc_bcc_et_corr_targetFix.dat')

    # normalize and apply corrections
    yieldNorm = []

    for iQ2 in range(nQ2bins):
        yieldNorm.append([])
        outXSEC.append([])

        for sec in range(nSec):
            q = iQ2 + 5

            WcurrentBin = wYield1D[sec][iQ2]

            yieldNorm[-1].append(ROOT.TH1D(
                'nrY' + str(iQ2) + ' ' + str(sec) + addN,
                'nrY' + str(iQ2) + ' ' + str(sec) + addN,
                WcurrentBin.GetNbinsX(),
                WcurrentBin.GetXaxis().GetXmin(),
                WcurrentBin.GetXaxis().GetXmax()
            ))
            yieldNorm[-1][-1].SetDirectory(0)

            outXSEC[-1].append(ROOT.TH1D(
                'xsec' + str(iQ2) + ' ' + str(sec) + addN,
                'xsec' + str(iQ2) + ' ' + str(sec) + addN,
                WcurrentBin.GetNbinsX(),
                WcurrentBin.GetXaxis().GetXmin(),
                WcurrentBin.GetXaxis().GetXmax()
            ))
            outXSEC[-1][-1].SetDirectory(0)

            normalize(WcurrentBin, yieldNorm[-1][-1], q, chargeFC_tmp)

            if use_corr:
                ApplyCorr(yieldNorm[-1][-1], outXSEC[-1][-1], q, corrDF)
            else:
                for iBin in range(1, WcurrentBin.GetNbinsX() + 1):
                    outXSEC[-1][-1].SetBinContent(iBin, yieldNorm[-1][-1].GetBinContent(iBin))
                    outXSEC[-1][-1].SetBinError(iBin, yieldNorm[-1][-1].GetBinError(iBin))

    # save output
    out_file = ROOT.TFile(addN + "_unfolded_output.root", "RECREATE")

    out_file.mkdir("UnfoldLIN1D")
    out_file.cd("UnfoldLIN1D")
    for sec in range(nSec):
        h1 = tmp_hUnfold_lin[sec]
        h1.SetName(f"hUnfoldLIN1D_sec{sec}")
        h1.Write()

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

def IntOverSec(listIn, listOut, endN, isTenBin=False):
    nQ2bins = 11
    if isTenBin:
        nQ2bins = 10
    for iQ2 in range(nQ2bins):
        listOut.append(ROOT.TH1D('intCutFunc' + str(iQ2) + endN,
                                 'intCutFunc' + str(iQ2) + endN,
                                 listIn[0][0].GetNbinsX(),
                                 listIn[0][0].GetXaxis().GetXmin(),
                                 listIn[0][0].GetXaxis().GetXmax()))
        listOut[-1].SetDirectory(0)

        for iXbin in range(1, 1 + listIn[0][0].GetNbinsX(), 1):
            value = 0
            for sec in range(6):
                value += listIn[iQ2][sec].GetBinContent(iXbin)
            listOut[-1].SetBinContent(iXbin, value / 6)

def GetTitle(npad):
    nQ2 = 50
    q2Min = 1
    q2Max = 2500
    deltaQ2 = np.log(q2Max/q2Min) / nQ2
    Q2bin_min = q2Min * math.exp((npad + 5) * deltaQ2)
    Q2bin_max = q2Min * math.exp((npad + 6) * deltaQ2)
    return "             " + str(Q2bin_min)[0:4] + " < Q^{2} < " + str(Q2bin_max)[0:4] + " GeV^{2}"

print("end cell 2")
#----------------------------------------End Cell 2 -------------------------------------------------------
#----------------------------------------Cell 3 -----------------------------------------------------------
epochs = [1,2,3]

xsec_allS_rooLIN1D_bayes_1 = []
xsec_allS_rooLIN1D_bayes_2 = []
xsec_allS_rooLIN1D_bayes_3 = []

GetXSEC_rooLIN1D_fromOriginal2D(
    fSim,
    fData,
    xsec_allS_rooLIN1D_bayes_2,
    'bayes2_LIN1D_fromOriginal2D_2',
    fccConst,
    epochs[1],
    method='bayes',
    use_corr=False
)

xsec_integrated_LIN1D_bayes_1 = []
xsec_integrated_LIN1D_bayes_2 = []
xsec_integrated_LIN1D_bayes_3 = []

IntOverSec(
    xsec_allS_rooLIN1D_bayes_2,
    xsec_integrated_LIN1D_bayes_2,
    'xsec_integrated_LIN1D_bayes_2'
)

"""
xsec_allS_rooLIN1D_bayes_2_lastQ2 = []
GetXSEC_rooLIN1D(
    fSim_lastQ2,
    fData_lastQ2,
    xsec_allS_rooLIN1D_bayes_2_lastQ2,
    'bayes2_LIN1D_2_lastQ2',
    fccConst,
    epochs[1],
    method='bayes'
)
xsec_integrated_LIN1D_bayes_2_lastQ2 = []
IntOverSec(
    xsec_allS_rooLIN1D_bayes_2_lastQ2,
    xsec_integrated_LIN1D_bayes_2_lastQ2,
    'xsec_integrated_LIN1D_bayes_2_lastQ2'
)
"""

print("end cell 3")
#----------------------------------------End Cell 3 -------------------------------------------------------
#----------------------------------------Cell 4 -----------------------------------------------------------
from pathlib import Path

def _load_optional_table(candidates, required_cols=None):
    for path in candidates:
        p = Path(path)
        if not p.exists():
            continue

        # Try normal CSV first, then whitespace-separated
        for read_kwargs in ({}, {'sep': r'\s+'}):
            try:
                df = pd.read_csv(p, **read_kwargs)
                if required_cols is None or all(col in df.columns for col in required_cols):
                    return df, str(p)
            except Exception:
                pass

    return None, None


# Try local file first, then Valera-like absolute path
statDF, stat_path = _load_optional_table(
    [
        'stat_combL.dat',
        '/home/valerii/Clas12/Inclusive/fromLaptop/Unfolding/XSEC_final/statErrF.dat'
    ],
    required_cols=['Q2', 'W', 'stat_combL']
)

if statDF is None:
    print('stat file not found -> keeping central values unchanged and leaving current bin errors as they are')
else:
    print(f'Loaded stat uncertainties from: {stat_path}')


def AddSysUNC(histIN_L, df):

    # No stat file -> do nothing
    if df is None:
        return

    for iQ2 in range(1,10,1):

        histIN = histIN_L[iQ2]

        # convert binning
        q = iQ2 + 5

        nQ2 = 50
        q2Min = 1
        q2Max = 2500
        deltaQ2 = np.log(q2Max/q2Min)/nQ2
        # center of the bin
        q2 = (q2Min * math.exp(q * deltaQ2) + q2Min * math.exp((q + 1) * deltaQ2)) / 2

        q2binsDF = pd.unique(df['Q2'])
        oneQ2_df = df[df['Q2'] == q2binsDF[q-6]]

        if (abs(q2 - q2binsDF[q-6]) > 0.2):
            print('Q2 range Error')
            return

        for iXbin in range(1, 1 + histIN.GetNbinsX(), 1):
            corrFound = False
            for index, row in oneQ2_df.iterrows():
                if (abs(row['W'] - histIN.GetBinCenter(iXbin)) < 0.01):
                    corrFound = True
                    corr = row['stat_combL']
                    evYield = histIN.GetBinContent(iXbin)
                    histIN.SetBinError(iXbin, evYield * corr/100)
                    break

            # Keep Valera's behavior when file exists:
            if (not corrFound):
                histIN.SetBinContent(iXbin, 0.)
                histIN.SetBinError(iXbin, 0.)


AddSysUNC(xsec_integrated_LIN1D_bayes_2, statDF)


print("end cell 4")
#----------------------------------------End Cell 4 -------------------------------------------------------
#----------------------------------------Cell 5 -----------------------------------------------------------
def GetSysHists(xsecIN, sysOut, df, nameH):

    for iQ2 in range(10):
        histD = xsecIN[iQ2]
        sysOut.append(ROOT.TH1D('intSysF' + str(iQ2) + nameH,
                                'intSysF' + str(iQ2) + nameH,
                                histD.GetNbinsX(),
                                histD.GetXaxis().GetXmin(),
                                histD.GetXaxis().GetXmax()))
        sysOut[-1].SetDirectory(0)

        # No sys file -> leave zero histogram and continue
        if df is None:
            continue

        # convert binning
        q = iQ2 + 5

        nQ2 = 50
        q2Min = 1
        q2Max = 2500
        deltaQ2 = np.log(q2Max/q2Min)/nQ2
        # center of the bin
        q2 = (q2Min * math.exp(q * deltaQ2) + q2Min * math.exp((q + 1) * deltaQ2)) / 2

        q2binsDF = pd.unique(df['Q2'])
        oneQ2_df = df[df['Q2'] == q2binsDF[q-6]]

        if (abs(q2 - q2binsDF[q-6]) > 0.2):
            print('Q2 range Error')
            continue

        for iXbin in range(1, 1 + histD.GetNbinsX(), 1):
            corrFound = False
            for index, row in oneQ2_df.iterrows():
                if (abs(row['W'] - histD.GetBinCenter(iXbin)) < 0.01):
                    corrFound = True
                    corr = row['Total']
                    sysOut[-1].SetBinContent(iXbin, 1000 * histD.GetBinContent(iXbin) * corr/100)
                    break

            if (not corrFound):
                sysOut[-1].SetBinContent(iXbin, 0.)
                sysOut[-1].SetBinError(iXbin, 0.)


sysUncHistL = []
sysUncHistL_lastQ2 = []

sysDF, sys_path = _load_optional_table(
    [
        'totalSysReadable_Jan2024_scale.dat',
        '/home/valerii/Clas12/Inclusive/systematics/totalSysReadable.dat'
    ],
    required_cols=['Q2', 'W']
)

if sysDF is not None:
    if 'Total_withScale' in sysDF.columns:
        sysDF = sysDF[['Q2', 'W', 'Total_withScale']]
        sysDF = sysDF.rename(columns={"Total_withScale": "Total"})
    elif 'Total' in sysDF.columns:
        sysDF = sysDF[['Q2', 'W', 'Total']]
    else:
        print('systematics file found, but no Total/Total_withScale column -> using zero systematics for now')
        sysDF = None

if sysDF is None:
    print('systematics file not found -> creating zero systematic histograms')
else:
    print(f'Loaded systematics from: {sys_path}')

GetSysHists(xsec_integrated_LIN1D_bayes_2, sysUncHistL, sysDF, 'sys_Name')
#GetSysHists(xsec_integrated_1D_bayes_2_lastQ2, sysUncHistL_lastQ2, sysDF, 'sys_Name_lastQ2')

print("end cell 5")
#----------------------------------------End Cell 5 -------------------------------------------------------

#----------------------------------------Cell 6 -----------------------------------------------------------
import os
import re
import glob
import math
import numpy as np
import pandas as pd
import ROOT
from ROOT import TCanvas, TLegend, TGraph, TGraphErrors

ROOT.gStyle.SetOptStat(0)

def _extract_q2_from_filename(path):
    m = re.search(r'Q2=([0-9.]+)\.dat$', os.path.basename(path))
    if not m:
        raise ValueError(f"Cannot extract Q2 from filename: {path}")
    return float(m.group(1))

def _load_exp_graph(dat_file, xmin=1.15, xmax=2.5):
    # columns expected: W eps sigma error sys_error
    df = pd.read_csv(dat_file, sep=r"\s+", comment="#")
    df = df[(df["W"] >= xmin) & (df["W"] <= xmax)].copy()

    x = df["W"].to_numpy(dtype=float)
    y = df["sigma"].to_numpy(dtype=float) * 1e-3

    # total error = sqrt(error^2 + sys_error^2) * 1e-3
    ey = np.sqrt(
        df["error"].to_numpy(dtype=float)**2 +
        df["sys_error"].to_numpy(dtype=float)**2
    ) * 1e-3

    ex = np.zeros_like(x, dtype=float)
    
    #to get rid of error bars
    ey = np.zeros_like(x, dtype=float)

    gr = TGraphErrors(
        len(x),
        np.array(x, dtype='float64'),
        np.array(y, dtype='float64'),
        np.array(ex, dtype='float64'),
        np.array(ey, dtype='float64')
    )
    gr.SetMarkerStyle(20)              # red circles
    gr.SetMarkerSize(0.7)
    gr.SetMarkerColor(ROOT.kRed + 1)
    gr.SetLineColor(ROOT.kRed + 1)
    gr.SetLineWidth(1)
    return gr, df

def _q2_center_from_iQ2(iQ2):
    q = iQ2 + 5
    nQ2 = 50
    q2Min = 1.0
    q2Max = 2500.0
    deltaQ2 = np.log(q2Max / q2Min) / nQ2
    return 0.5 * (
        q2Min * math.exp(q * deltaQ2) +
        q2Min * math.exp((q + 1) * deltaQ2)
    )

def _hist_to_graph(hist, xmin=1.15, xmax=2.5):
    x_vals = []
    y_vals = []

    for ib in range(1, hist.GetNbinsX() + 1):
        x = hist.GetBinCenter(ib)
        if x < xmin or x > xmax:
            continue
        y = hist.GetBinContent(ib)
        if y == 0:
            continue
        x_vals.append(x)
        y_vals.append(y)

    gr = TGraph(len(x_vals), np.array(x_vals, dtype='float64'), np.array(y_vals, dtype='float64'))
    gr.SetMarkerStyle(20)              # black circles
    gr.SetMarkerSize(0.7)
    gr.SetMarkerColor(ROOT.kBlack)
    gr.SetLineColor(ROOT.kBlack)
    gr.SetLineWidth(2)
    return gr

def PlotIntegratedXsecVsW_withExp(hist_list, exp_dir, out_dir='xsec_vs_W_overlay_png', xmin=1.15, xmax=2.5):
    os.makedirs(out_dir, exist_ok=True)

    exp_files = sorted(glob.glob(os.path.join(exp_dir, "*.dat")))
    if not exp_files:
        raise RuntimeError(f"No .dat files found in: {exp_dir}")

    exp_map = {}
    for f in exp_files:
        q2 = _extract_q2_from_filename(f)
        exp_map[q2] = f

    exp_q2_vals = sorted(exp_map.keys())

    # only Q2 bins 1..9
    for iQ2 in range(1, 10):
        h_in = hist_list[iQ2]
        if not h_in:
            print(f"Skip Q2 bin {iQ2}: histogram is missing")
            continue

        q2_center = _q2_center_from_iQ2(iQ2)

        # match nearest experimental file in Q2
        q2_match = min(exp_q2_vals, key=lambda q: abs(q - q2_center))
        exp_file = exp_map[q2_match]

        gr_exp, df_exp = _load_exp_graph(exp_file, xmin=xmin, xmax=xmax)
        gr_me = _hist_to_graph(h_in, xmin=xmin, xmax=xmax)

        # determine plotting range
        y_vals = []

        for ib in range(1, h_in.GetNbinsX() + 1):
            x = h_in.GetBinCenter(ib)
            if xmin <= x <= xmax:
                y = h_in.GetBinContent(ib)
                if y > 0:
                    y_vals.append(y)

        if len(df_exp) > 0:
            y_exp = df_exp["sigma"].to_numpy(dtype=float) * 1e-3
            ey_exp = np.sqrt(
                df_exp["error"].to_numpy(dtype=float)**2 +
                df_exp["sys_error"].to_numpy(dtype=float)**2
            ) * 1e-3
            y_vals.extend(y_exp + ey_exp)

        ymax = max(y_vals) if y_vals else 1.0

        c = TCanvas(f'c_xsec_Q2_{iQ2}', f'c_xsec_Q2_{iQ2}', 900, 700)
        c.SetMargin(0.12, 0.05, 0.12, 0.08)
        c.SetGrid()

        frame = ROOT.TH1D(f"frame_Q2_{iQ2}", "", 100, xmin, xmax)
        frame.SetMinimum(0.0)
        frame.SetMaximum(1.20 * ymax)
        frame.SetTitle(GetTitle(iQ2))
        frame.GetXaxis().SetTitle("W (GeV)")
        frame.GetYaxis().SetTitle("Cross section")
        frame.GetXaxis().CenterTitle()
        frame.GetYaxis().CenterTitle()
        frame.Draw()

        gr_me.Draw("P SAME")
        gr_exp.Draw("P SAME")

        leg = TLegend(0.12, 0.72, 0.42, 0.88)
        leg.SetBorderSize(0)
        leg.SetFillStyle(0)
        leg.SetTextSize(0.05) 
        leg.AddEntry(gr_me, "My 1D unfolded w/ lin.bin WITH CORR", "p")
        leg.AddEntry(gr_exp, "Valerii published (2D unfolded)", "p")
        leg.Draw()

        out_name = os.path.join(out_dir, f"xsec_vs_W_Q2bin_{iQ2}.png")
        c.SaveAs(out_name)

        print(f"Q2 bin {iQ2}: theory Q2 center = {q2_center:.4f}, matched exp file Q2 = {q2_match:.4f}")
        print(f"Saved: {out_name}")

# Set this to the folder with Valera's .dat files
exp_data_dir = "../../HarryLeeDCC/paper_plots/exp_data"

PlotIntegratedXsecVsW_withExp(
    xsec_integrated_LIN1D_bayes_2,
    exp_dir=exp_data_dir,
    out_dir='NEW_ATTEMPT_xsec_vs_W_overlay_LIN1D_png_NO_CORRECTIONS',
    xmin=1.1,
    xmax=2.6
)
print("end cell 6")
#----------------------------------------End Cell 6 -------------------------------------------------------
#----------------------------------------Cell 7 -----------------------------------------------------------
#import os
#import re
#import glob
#import math
#import numpy as np
#import pandas as pd
#import ROOT
#from ROOT import TCanvas, TLegend, TGraph
#
#ROOT.gStyle.SetOptStat(0)
#
#def _extract_q2_from_filename(path):
#    m = re.search(r'Q2=([0-9.]+)\.dat$', os.path.basename(path))
#    if not m:
#        raise ValueError(f"Cannot extract Q2 from filename: {path}")
#    return float(m.group(1))
#
#def _load_exp_df(dat_file, xmin=1.15, xmax=2.5):
#    # expected columns: W eps sigma error sys_error
#    df = pd.read_csv(dat_file, sep=r"\s+", comment="#")
#    df = df[(df["W"] >= xmin) & (df["W"] <= xmax)].copy()
#    df["sigma_plot"] = df["sigma"].astype(float) * 1e-3
#    return df
#
#def _q2_center_from_iQ2(iQ2):
#    q = iQ2 + 5
#    nQ2 = 50
#    q2Min = 1.0
#    q2Max = 2500.0
#    deltaQ2 = np.log(q2Max / q2Min) / nQ2
#    return 0.5 * (
#        q2Min * math.exp(q * deltaQ2) +
#        q2Min * math.exp((q + 1) * deltaQ2)
#    )
#
#def _hist_value_at_x(hist, x):
#    ibin = hist.GetXaxis().FindBin(x)
#    if ibin < 1 or ibin > hist.GetNbinsX():
#        return None
#    return hist.GetBinContent(ibin)
#
#def _ratio_graph_from_hist_and_exp(hist, df_exp, min_my_y=1e-30):
#    x_vals = []
#    r_vals = []
#
#    for _, row in df_exp.iterrows():
#        x = float(row["W"])
#        y_exp = float(row["sigma_plot"])
#        y_me = _hist_value_at_x(hist, x)
#
#        if y_me is None:
#            continue
#        if y_me <= min_my_y:
#            continue
#
#        ratio = y_exp / y_me
#        x_vals.append(x)
#        r_vals.append(ratio)
#
#    gr = TGraph(
#        len(x_vals),
#        np.array(x_vals, dtype="float64"),
#        np.array(r_vals, dtype="float64")
#    )
#    gr.SetMarkerStyle(20)
#    gr.SetMarkerSize(1.0)
#    gr.SetMarkerColor(ROOT.kBlack)
#    gr.SetLineColor(ROOT.kBlack)
#    gr.SetLineWidth(2)
#    return gr, np.array(r_vals, dtype=float)
#
#def PlotRatio_ValeriiOverMine(hist_list, exp_dir,
#                              out_dir="ratio_valerii_over_mine_png",
#                              xmin=1.15, xmax=2.5,
#                              last_q2_xmax=2.25):
#    os.makedirs(out_dir, exist_ok=True)
#
#    exp_files = sorted(glob.glob(os.path.join(exp_dir, "*.dat")))
#    if not exp_files:
#        raise RuntimeError(f"No .dat files found in: {exp_dir}")
#
#    exp_map = {}
#    for f in exp_files:
#        q2 = _extract_q2_from_filename(f)
#        exp_map[q2] = f
#    exp_q2_vals = sorted(exp_map.keys())
#
#    # only Q2 bins 1..9
#    for iQ2 in range(1, 10):
#        h_in = hist_list[iQ2]
#        if not h_in:
#            print(f"Skip Q2 bin {iQ2}: histogram is missing")
#            continue
#
#        q2_center = _q2_center_from_iQ2(iQ2)
#        q2_match = min(exp_q2_vals, key=lambda q: abs(q - q2_center))
#        exp_file = exp_map[q2_match]
#
#        xmax_use = last_q2_xmax if iQ2 == 9 else xmax
#
#        df_exp = _load_exp_df(exp_file, xmin=xmin, xmax=xmax_use)
#        gr_ratio, ratios = _ratio_graph_from_hist_and_exp(h_in, df_exp)
#
#        if len(ratios) == 0:
#            print(f"Skip Q2 bin {iQ2}: no valid ratio points")
#            continue
#
#        ymin = min(0.9, 0.95 * np.min(ratios))
#        ymax = max(1.1, 1.05 * np.max(ratios))
#
#        c = TCanvas(f"c_ratio_Q2_{iQ2}", f"c_ratio_Q2_{iQ2}", 900, 700)
#        c.SetMargin(0.12, 0.05, 0.12, 0.08)
#        c.SetGrid()
#
#        frame = ROOT.TH1D(f"frame_ratio_Q2_{iQ2}", "", 100, xmin, xmax_use)
#        frame.SetMinimum(ymin)
#        frame.SetMaximum(ymax)
#        frame.SetTitle(GetTitle(iQ2))
#        frame.GetXaxis().SetTitle("W (GeV)")
#        frame.GetYaxis().SetTitle("Valerii / Mine")
#        frame.GetXaxis().CenterTitle()
#        frame.GetYaxis().CenterTitle()
#        frame.Draw()
#
#        # horizontal line at 1
#        line = ROOT.TLine(xmin, 1.0, xmax_use, 1.0)
#        line.SetLineStyle(2)
#        line.SetLineWidth(2)
#        line.Draw("SAME")
#
#        gr_ratio.Draw("P SAME")
#
#        leg = TLegend(0.14, 0.78, 0.44, 0.88)
#        leg.SetBorderSize(0)
#        leg.SetFillStyle(0)
#        leg.AddEntry(gr_ratio, "Valerii / Mine", "p")
#        leg.Draw()
#
#        out_name = os.path.join(out_dir, f"ratio_Q2bin_{iQ2}.png")
#        c.SaveAs(out_name)
#
#        print(f"Q2 bin {iQ2}: matched exp file Q2 = {q2_match:.4f}")
#        print(f"Saved: {out_name}")
#
## Choose which result you want to compare:
## xsec_integrated_1D_bayes_2
## xsec_integrated_LIN1D_bayes_2
#exp_data_dir = "../../HarryLeeDCC/paper_plots/exp_data"
#
#PlotRatio_ValeriiOverMine(
#    xsec_integrated_LIN1D_bayes_2,
#    exp_dir=exp_data_dir,
#    out_dir="NEW_attempt_ratio_valerii_over_mine_LIN1D_png",
#    xmin=1.15,
#    xmax=2.5,
#    last_q2_xmax=2.25
#)
#
#print("end cell 7")
#----------------------------------------End Cell 7 -------------------------------------------------------
