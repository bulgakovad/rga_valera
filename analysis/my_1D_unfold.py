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
# 1D unfolding RooUnfold by W in each fixed Q2 bin

def _must_get(root_file, key):
    obj = root_file.Get(key)
    if not obj:
        raise RuntimeError(f"Cannot find object: {key}")
    return obj

def _clone_detached(h, new_name):
    hc = h.Clone(new_name)
    hc.SetDirectory(0)
    return hc

# 1D unfolding response
def getResponse1D(fileSim, MSL, TSML, RSL, ORSL):
    nSec = 6
    nQ2bins = 11

    for sec in range(nSec):
        MSL.append([])
        TSML.append([])
        RSL.append([])
        ORSL.append([])

        for iQ2 in range(nQ2bins):
            fold = f'Deconv1D/S{sec}/Q2_{iQ2}/'

            h_meas  = _must_get(fileSim, fold + 'measured')
            h_truth = _must_get(fileSim, fold + 'truth')
            h_resp  = _must_get(fileSim, fold + 'response')

            MSL[-1].append(h_meas)
            TSML[-1].append(h_truth)
            RSL[-1].append(h_resp)
            ORSL[-1].append(ROOT.RooUnfoldResponse(MSL[-1][-1], TSML[-1][-1], RSL[-1][-1]))

# GetXSEC using 1D RooUnfold methods
def GetXSEC_roo1D(fileSim, fileData, outXSEC, addN, chargeFC_tmp, regParam, method='bayes'):
    # regParam is nIter for bayes or i_comp for SVD
    if (method != 'bayes' and method != 'bbb' and method != 'svd'):
        print('unsupported deconvolution method is used, please put bayes or bbb or svd')
        return

    # sectors
    nSec = 6
    # number of Q2 bins
    nQ2bins = 11

    # initialization of internal RooUnfold matrices and struct.
    tmp_measuredL = []
    tmp_truthL = []
    tmp_respMatL = []
    tmp_respObjL = []

    # read root file Sim
    getResponse1D(fileSim, tmp_measuredL, tmp_truthL, tmp_respMatL, tmp_respObjL)

    # measured from Data
    tmp_hMeas = []

    # read root file Data
    for sec in range(nSec):
        tmp_hMeas.append([])
        for iQ2 in range(nQ2bins):
            fold = f'Deconv1D/S{sec}/Q2_{iQ2}_my/'
            tmp_hMeas[-1].append(_must_get(fileData, fold + f'measured1D_DATA_S{sec}_Q2_{iQ2}'))

    # RooUnfoldObjs
    tmp_unfold = []
    tmp_hUnfold = []

    for sec in range(nSec):
        tmp_unfold.append([])
        tmp_hUnfold.append([])

        for iQ2 in range(nQ2bins):

            # method selection
            if (method == 'bayes'):
                unf = ROOT.RooUnfoldBayes(tmp_respObjL[sec][iQ2], tmp_hMeas[sec][iQ2], regParam)

            if (method == 'bbb'):
                unf = ROOT.RooUnfoldBinByBin(tmp_respObjL[sec][iQ2], tmp_hMeas[sec][iQ2])

            if (method == 'svd'):
                unf = ROOT.RooUnfoldSvd(tmp_respObjL[sec][iQ2], tmp_hMeas[sec][iQ2], regParam)

            #####################################
            # errors accounting
            unf.IncludeSystematics(0)
            unf.SetOverflow(0)
            # unf.SetVerbose(2)
            ######################################

            # Actual deconvolution
            h_unf = unf.Hunfold()
            h_unf = _clone_detached(h_unf, f'hUnfold1D_sec{sec}_Q2bin{iQ2}_{addN}')

            tmp_unfold[-1].append(unf)
            tmp_hUnfold[-1].append(h_unf)

    # unfolded 1D W-yields
    wYield1D = []

    for sec in range(nSec):
        wYield1D.append([])
        for iQ2 in range(nQ2bins):
            wYield1D[-1].append(_clone_detached(tmp_hUnfold[sec][iQ2], f'wYield1D_sec{sec}_Q2bin{iQ2}_{addN}'))

    # read file with all the corrections factors
    # corrDF = pd.read_csv('/home/valerii/Clas12/Inclusive/xsecs/rc_bcc_et_corr.dat') Not need?
    corrDF = pd.read_csv('rc_bcc_et_corr_targetFix.dat') 

    # W yield corrected for Lumin.
    yieldNorm = []

    for iQ2 in range(nQ2bins):

        yieldNorm.append([])
        outXSEC.append([])

        for sec in range(nSec):
            # Q2 convert binning to match binning in ana12.
            q = iQ2 + 5

            WcurrentBin = wYield1D[sec][iQ2]

            yieldNorm[-1].append(ROOT.TH1D('nrY' + str(iQ2) + ' ' + str(sec) + addN,
                                           'nrY' + str(iQ2) + ' ' + str(sec) + addN,
                                           WcurrentBin.GetNbinsX(),
                                           WcurrentBin.GetXaxis().GetXmin(),
                                           WcurrentBin.GetXaxis().GetXmax()))
            yieldNorm[-1][-1].SetDirectory(0)

            outXSEC[-1].append(ROOT.TH1D('xsec' + str(iQ2) + ' ' + str(sec) + addN,
                                         'xsec' + str(iQ2) + ' ' + str(sec) + addN,
                                         WcurrentBin.GetNbinsX(),
                                         WcurrentBin.GetXaxis().GetXmin(),
                                         WcurrentBin.GetXaxis().GetXmax()))
            outXSEC[-1][-1].SetDirectory(0)

            normalize(WcurrentBin, yieldNorm[-1][-1], q, chargeFC_tmp)

            ApplyCorr(yieldNorm[-1][-1], outXSEC[-1][-1], q, corrDF)    



    out_file = ROOT.TFile(addN + "_unfolded_output.root", "RECREATE")

    # Save unfolded 1D hists
    out_file.mkdir("Unfold1D")
    out_file.cd("Unfold1D")
    for sec in range(nSec):
        out_file.mkdir(f"Unfold1D/sec{sec}")
        out_file.cd(f"Unfold1D/sec{sec}")
        for iQ2 in range(nQ2bins):
            h1 = tmp_hUnfold[sec][iQ2]
            h1.SetName(f"hUnfold1D_sec{sec}_Q2bin{iQ2}")
            h1.Write()

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
    Q2bin_max = q2Min*math.exp((npad + 6)*deltaQ2)
    return  "             " + str(Q2bin_min)[0:4] + " < Q^{2} < " +  str(Q2bin_max)[0:4] + " GeV^{2}"
print("end cell 2")
#----------------------------------------End Cell 2 --------------------------------------------------------

#----------------------------------------Cell 3 -----------------------------------------------------------
epochs = [1,2,3]

xsec_allS_roo1D_bayes_1 = []
xsec_allS_roo1D_bayes_2 = []
xsec_allS_roo1D_bayes_3 = []

GetXSEC_roo1D(fSim, fData, xsec_allS_roo1D_bayes_2, 'bayes2_1D_2', fccConst, epochs[1], method = 'bayes')

xsec_integrated_1D_bayes_1 = []
xsec_integrated_1D_bayes_2 = []
xsec_integrated_1D_bayes_3 = []

IntOverSec(xsec_allS_roo1D_bayes_2, xsec_integrated_1D_bayes_2, 'xsec_integrated_1D_bayes_2')

"""
xsec_allS_roo1D_bayes_2_lastQ2 = []
GetXSEC_roo1D(fSim_lastQ2, fData_lastQ2, xsec_allS_roo1D_bayes_2_lastQ2, 'bayes2_1D_2_lastQ2', fccConst, epochs[1], method = 'bayes')
xsec_integrated_1D_bayes_2_lastQ2 = []
IntOverSec(xsec_allS_roo1D_bayes_2_lastQ2, xsec_integrated_1D_bayes_2_lastQ2, 'xsec_integrated_1D_bayes_2_lastQ2')
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


AddSysUNC(xsec_integrated_1D_bayes_2, statDF)
#AddSysUNC(xsec_integrated_1D_bayes_2_lastQ2, statDF)

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

GetSysHists(xsec_integrated_1D_bayes_2, sysUncHistL, sysDF, 'sys_Name')
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
from ROOT import TCanvas, TLegend, TGraph

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

    gr = TGraph(len(x), x, y)
    gr.SetMarkerStyle(20)              # red circles
    gr.SetMarkerSize(0.7)
    gr.SetMarkerColor(ROOT.kRed + 1)
    gr.SetLineColor(ROOT.kRed + 1)
    gr.SetLineWidth(2)
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
            y_vals.extend(df_exp["sigma"].to_numpy(dtype=float) * 1e-3)

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
        leg.AddEntry(gr_me, "My 1D unfolded", "p")
        leg.AddEntry(gr_exp, "Valerii published (2D unfolded)", "p")
        leg.Draw()

        out_name = os.path.join(out_dir, f"xsec_vs_W_Q2bin_{iQ2}.png")
        c.SaveAs(out_name)

        print(f"Q2 bin {iQ2}: theory Q2 center = {q2_center:.4f}, matched exp file Q2 = {q2_match:.4f}")
        print(f"Saved: {out_name}")

# Set this to the folder with Valera's .dat files
exp_data_dir = "../../HarryLeeDCC/paper_plots/exp_data"

PlotIntegratedXsecVsW_withExp(
    xsec_integrated_1D_bayes_2,
    exp_dir=exp_data_dir,
    out_dir='xsec_vs_W_overlay_png',
    xmin=1.1,
    xmax=2.6
)

print("end cell 6")
#----------------------------------------End Cell 6 -------------------------------------------------------


