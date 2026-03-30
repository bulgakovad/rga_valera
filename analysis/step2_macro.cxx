// source ~/.cshrc
// root
// gSystem->Load("libRooUnfold");
// .L step2_macro.cxx+
// fix the warnings!!!
// ana12_sys(1) data or ana12_sys(0) - simu

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <math.h>
#include "Riostream.h"
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TRandom3.h"
#include "TLorentzVector.h"
#include "TRandom.h"
#include "TMath.h"
#include <vector>
#include <set>
#include <TAxis.h>
#include <TLorentzRotation.h>
#include <vector>
#include <algorithm>
#include <functional>
#include <TChain.h>
#include <TCutG.h>
#include "TStyle.h"
#include "TROOT.h"
#include "TLatex.h"

#include "RooUnfoldResponse.h"
#include "RooUnfoldBayes.h"

using namespace std;

#include "Math/Vector3D.h"
#include "Math/Vector4D.h"

using namespace ROOT::Math;

//S0 : -19 38

struct binFold{
    double Wmin,Wmax, WbinSize,q2Min,deltaQ2;
};

int GetSectorByPhi(const double phi){
	if (phi > -20. && phi <= 40.) return 0;
	if (phi > 40. && phi <= 100.) return 1;
	if (phi > 100. && phi <= 160.) return 2;
	if (phi > -140. && phi <= -80.) return 4;
	if (phi > -80. && phi <= -20.) return 5;
	if (phi > 160. || phi <= -140) return 3;
	throw out_of_range("wrong phi");
	return -1;
}


void getBinFold2D(double wRec, double wGen, double q2Rec, double q2Gen, int& binWrec, int& binWgen, int& binQ2rec, int& binQ2gen, const binFold& binParam){
    
    int startQ2bin = 5;
    
    //if (wGen < binParam.Wmin || wGen > binParam.Wmax)
    //    throw out_of_range("Wgen is wrong folding 2D");
    
    binWgen =  int ((wGen - binParam.Wmin) / binParam.WbinSize);
    
    //if (wRec < binParam.Wmin || wRec > binParam.Wmax)
    //    throw out_of_range("wRec is wrong folding 2D");
    
    binWrec =  int ((wRec - binParam.Wmin) / binParam.WbinSize);
    
    binQ2gen = log(q2Gen/ binParam.q2Min)/ binParam.deltaQ2;
    binQ2gen -= startQ2bin;
    binQ2rec = log(q2Rec/ binParam.q2Min)/ binParam.deltaQ2;
    binQ2rec -= startQ2bin;
    
    return;
    
}

double kin_W(TLorentzVector ele, float Ebeam);
double kin_Q2(TLorentzVector ele, float Ebeam);
int PCALFidXY(float x, float y, int cutLevel);

float HTCCWeight(float x, float y, float HTCCEff[250][250],  int stepXHTCC, int stepYHTCC, int isData);
float newHTCCWeight(float x, float y, float HTCCEff[250][250],  int stepXHTCC, int stepYHTCC, int isData);
void readHTCCEff(float HTCCEff[250][250]);
bool SFTriangleCut(float ecinE, float pcalE, float cutParams[6][10][2], int sector, int pBin, float shift);
void readTriangleCut(float cutParams[6][10][2], int isData);
bool DCFidXY(float X, float Y, int region, int sector, int cutLevel);

const int SFCutSigma = 3;
const int fidCutLevel = 1;
const int qBinMax = 20;

//Valerii
bool SameBinW(float W1, float W2, size_t nWBins);
void fillVfromTree(float (&data)[100], const vector<float>* branch, const int icount);
void fillVfromTree(float (&data)[100], const vector<double>* branch, const int icount);

struct DCXY{
 double r1X,r1Y,r2X,r2Y,r3X,r3Y;
};

bool CutP(const float p, const int cutLevel);
bool CutPCALdepos(const float pcalE, const int cutLevel);
bool CutVz(const float _vz, const int cutLevel);
bool CutDCfid(const DCXY& dc, int sec, const int cutLevel);
bool SfCutValerii(const double sf,const double p,const int sec,const int cutLevel,const int isData);
bool SfCutValerii_Edepos(const double sf,const double Edep,const int sec,const int cutLevel,const int isData);

bool phiSpikeCut(const double phi, const double theta, const int sec, const int cutLevel){

	double shiftSys = 0;
	if (cutLevel == 0) shiftSys = 1.;
	if (cutLevel == 2) shiftSys = -1.;
    
    // Sec 1:
    if (sec == 0){
        if (theta > (13.5 + shiftSys) && theta < (27 - shiftSys))
            if (phi > (-5 + shiftSys) && phi < (0 - shiftSys))
                return false;
    }
    // Sec 4:
    if (sec == 3){
        if (theta > (13.5 + shiftSys) && theta < (29.5 - shiftSys))
            if (phi > (-6 + shiftSys)&& phi < (-1 - shiftSys))
                return false;
    }
    return true; 
    
}

bool thetaMinCut(const double theta, const int cutLevel){
    return theta > (4 + cutLevel);
}
bool thetaMaxCut(const double theta, const int cutLevel){
    return true;
    //return theta < (26 - cutLevel);
}

bool phiMaxCut(const double phi, const double cutLevel){
    return true;
    //return phi < (11. - cutLevel / 2);
}

bool phiMinCut(const double phi, const double cutLevel){
    return true;
    //return phi > (-8.5 + cutLevel / 2);
}


bool kinemCut(const double Q2, const double W){
    
    return true;
    
    //const double Q2_min = 2.3, Q2_max = 10.6, W_min = 0.7, W_max = 2.9; 
    //for last bin:
    //return W >= W_min && W <= W_max && Q2 >= Q2_min && Q2 <= Q2_max;
}

bool thetaVSphiDCcut(const double theta, const double phi, const int cutLevel){

    return true;
    
    // theta, phi, sec, cutLevel
    double shiftTheta = 0;
    double shiftPhi = 0;
    
    if (cutLevel == 0){
        shiftTheta = -0.5;
        shiftPhi = -1;
    }
    if (cutLevel == 2){
        shiftTheta = 0.5;
        shiftPhi = 1;
    }
    
    if (phi < -20 + shiftPhi || phi > 20 - shiftPhi)
        return false;
    
    //if (phi < -9 + shiftPhi) return false;
    //double p0 = 9.86752, p1 = -0.296311, p2 = 0.0285485;
    //double p0 = 10.519376239326663 , p1 = -0.3248872393298614 , p2 = 0.035445354266716925;
    //double p0 = 10.4278, p1 = -0.33, p2 = 0.03578;
    
    
    //9.105612343379079 ,  -0.06809857039000168 ,  0.03999793373541474
    double p0 = 9.1056, p1 = -0.06809857, p2 = 0.03999793;
    
    
    return theta > p0 + p1 *phi + p2*phi*phi + shiftTheta;
}

enum cutType{
    vzCut = 0, 
    dcCut = 1,
    sfCut = 2,
    momCut = 3,
    pcalCut = 4,
    triagCut = 5,
    badElemCut = 6,
    pcalDeposCut = 7,
    thetaMax = 8,
    thetaMin = 9,
    thetaVSphiDC = 10,
    phiMax = 11,
    phiMin = 12,
    phiSpike = 13,
    allCuts = 14
};

set<cutType> AddAlllCutsExcept(const cutType noCut){
    
    // List of all cuts:
    const vector<cutType> vCuts = {cutType::vzCut, cutType::dcCut, cutType::sfCut, 
                                   cutType::momCut, cutType::pcalCut, cutType::triagCut, 
                                   cutType::badElemCut, cutType::pcalDeposCut,cutType::thetaMax, 
                                   cutType::thetaMin, cutType::thetaVSphiDC, cutType::phiMax, cutType::phiMin, cutType::phiSpike};
    set<cutType> cutsToApply;
    for (const auto& cut : vCuts){
        if (cut != noCut){
            cutsToApply.insert(cut);
        }
    }
    return cutsToApply;
}
          
struct KinemE{
    double momentum;
    float vz;
    float calEnerg;
    float ecinE;
    float pcalE;
    int pBin;
    float pcalHx;
    float pcalHy;
    float theta;
    float phi;
    float phi_notCorrected;
    float ecinHx;
    float ecinHy;
    float ecoutHx;
    float ecoutHy;
    float pcalV;
    float pcalW;
    float pcalU;
};

void rotVect(TVector3& resVector3, const int sec){
	resVector3.RotateZ(-60*sec/57.2958);
	resVector3.RotateY(-25/57.2958);
}

double func(double x, double k, double b){
    return k * x + b;
}

struct line{
    double k;
    double b;
};

double isOutOfLines(double x, double y, line topLine, line botLine){
    return y > func(x, topLine.k, topLine.b) || y < func(x, botLine.k, botLine.b);
}

double isBetweenOfLines(double x, double y, line topLine, line botLine){
    return y < func(x, topLine.k, topLine.b) && y > func(x, botLine.k, botLine.b);
}

bool BadElementKnockOut(double Hx_pcal, double Hy_pcal, double Hx_ecin, double Hy_ecin, double Hx_ecout, double Hy_ecout, int sector, int cutLevel);

//bool ApplyAllCutExcept(TLorentzVector p4_electron, float (&triangleCutParams)[6][10][2], const TVector3& resVector3, const KinemE& kin, const DCXY& dc, const int sec, const cutType noCut,const int isData);

bool ApplyCuts(TLorentzVector p4_electron, float (&triangleCutParams)[6][10][2], const KinemE& kin, const TVector3& PCALvector3, const DCXY& dc, const int sec ,const int isData, const set<cutType>& cutsToApply);



void smear(TLorentzVector *V4, int q); 

TRandom3 *myMC;

void fastMC(void){
        myMC = new TRandom3(0);
}


vector<string> GetFilesPath(int typeData, int simQ2range){
    //simQ2range:
    //0 - fullQ2; 
    //1 - 4.08-10.6
    
    vector<string> result;
    if (typeData == 1){
        
        //Data:
        result.push_back("/lustre24/expphy/volatile/clas12/bulgakov/rga_valera/data/pass1_rga_inclusive_runs_valera.dat.root");
        result.push_back("/lustre24/expphy/volatile/clas12/bulgakov/rga_valera/data/pass1_rga_inclusive_runs_valera.dat_1.root");
        result.push_back("/lustre24/expphy/volatile/clas12/bulgakov/rga_valera/data/pass1_rga_inclusive_runs_valera.dat_2.root");
    
    }
    //Sim 1:
    if (typeData == 0 && simQ2range == 0){
    	result.push_back("/lustre24/expphy/cache/clas12/rg-a/production/montecarlo/inclusive_rates_calcom/SimIter3_1.datAna.root"); 
        result.push_back("/lustre24/expphy/cache/clas12/rg-a/production/montecarlo/inclusive_rates_calcom/SimIter3_2.datAna.root"); 
        result.push_back("/lustre24/expphy/cache/clas12/rg-a/production/montecarlo/inclusive_rates_calcom/SimIter3_3.datAna.root"); 
        result.push_back("/lustre24/expphy/cache/clas12/rg-a/production/montecarlo/inclusive_rates_calcom/SimIter3_4.datAna.root"); 
        result.push_back("/lustre24/expphy/cache/clas12/rg-a/production/montecarlo/inclusive_rates_calcom/SimIter3_5.datAna.root"); 
        result.push_back("/lustre24/expphy/cache/clas12/rg-a/production/montecarlo/inclusive_rates_calcom/SimIter3_6.datAna.root"); 
        result.push_back("/lustre24/expphy/cache/clas12/rg-a/production/montecarlo/inclusive_rates_calcom/SimIter3_7.datAna.root"); 
        result.push_back("/lustre24/expphy/cache/clas12/rg-a/production/montecarlo/inclusive_rates_calcom/SimIter3_8.datAna.root"); 
        result.push_back("/lustre24/expphy/cache/clas12/rg-a/production/montecarlo/inclusive_rates_calcom/SimIter3_9.datAna.root");
        result.push_back("/lustre24/expphy/cache/clas12/rg-a/production/montecarlo/inclusive_rates_calcom/SimIter3_10.datAna.root");
        result.push_back("/lustre24/expphy/cache/clas12/rg-a/production/montecarlo/inclusive_rates_calcom/SimIter3_11.datAna.root");
        result.push_back("/lustre24/expphy/cache/clas12/rg-a/production/montecarlo/inclusive_rates_calcom/SimIter3_12.datAna.root");
    
    	

    } 
    //sim 2, high Q2
    if (typeData == 0 && simQ2range == 1){
    }      
    return result;
}

void ana12_sys(int isData){
////////////// Simulation range:////////////////// 
    int simQ2range = 0;
    //other params:
    
    bool isSystematics = false;
    bool isUnfolding = true;
//////////////////////////////////////////////////////    
/////////// SMEARING: ////////////////////////////////  
//////////////////////////////////////////////////////  
    bool isFXsmearing = false;
    bool isRichardSmearing = false;
    bool isValeriiSmearing = true;
    
    //double smearingFactor = 1.7;
    
    if (isFXsmearing && isRichardSmearing) {
    	throw out_of_range("BOTH FX and RICHARD SMEARING ON, SELECT JUST ONE or NONE");
    	return;
    }

//////////////////////////////////////////////////////    
/////////// Mom Corr: ////////////////////////////////  
//////////////////////////////////////////////////////

    //Data:
    bool isMomCorrOld = false;
    bool isMomCorrMay2023 = true;
    //MC: 
    bool isMomCorrMC = false;
    
    if (isMomCorrOld && isMomCorrMay2023) {
    	throw out_of_range("BOTH OLD and NEW MC ON, SELECT JUST ONE or NONE");
    	return;
    }
    
//////////////////////////////////////////////////////    
/////////// END Mom Corr: ////////////////////////////////  
//////////////////////////////////////////////////////    
    bool isKinemCut = true;
    double maxWcut = 2.9;
    double minWcut = 0.8;
    
////////////// Q2 binning: ////////////////////////////////
	int nQ2 = 50;
	float q2Min = 1;
	float q2Max = 2500;
	float deltaQ2 = log(q2Max/q2Min)/nQ2;
//////////////////////////
    
    /////////////
    vector<string> file_list = GetFilesPath(isData, simQ2range);
    
////////////// Output Files names:////////////////
    
    //Place holder for output file:
	TFile *resultsF = new TFile("PH_0.root", "recreate");
    // Sim:
    // 4 interpolated weight
	if (isData == 0 && simQ2range == 0) 
		// 5 actually for 1.7
        //resultsF = new TFile("SIM_newSM_4f_f15_5.root", "recreate");
        //resultsF = new TFile("SIM_Apr4_iterPierre_N3_12f.root", "recreate");
        resultsF = new TFile("SIM_xsec_MCcorr_SMfix_Jul22.root", "recreate");
        
        //resultsF = new TFile("SIM_newSM_NOvzShift_1.root", "recreate");

	if (isData == 0 && simQ2range == 1) 
        resultsF = new TFile("noNeed.root", "recreate");

    //Data:
	if (isData == 1) 
        //resultsF = new TFile("Data_Apr4_iterPierre_N3_12f.root", "recreate");
        resultsF = new TFile("Data_xsec_MCcorr_SMfix_Jul22.root", "recreate");
    
    //Empty Target:
        //resultsF = new TFile("/w/hallb-scshelf2102/clas12/valerii/data/ana2/RGA/ET_Aug18_50MeV_shift_1_newFileList.root", "recreate");
    
    
////////////// W binning: ////////////////////////////////
    //BinSize = 50 MeV, 25 MeV Shift
    const float wBinSize = 0.05;
	const unsigned nWBins = 30;
	const float lowBorderW = 1.025;
	const float highBorderW = lowBorderW + wBinSize * nWBins;

    //original binning:
	//unsigned nWBins = 100;
	//float lowBorderW = 0;
	//float highBorderW = 5.0;

////////////// Kinematic params: ///////////////////////////
    
    float eBeam = 10.6041;
    
////////////// Output file structure: ///////////////////////////
    
	resultsF->mkdir("overview");

	resultsF->mkdir("HTCC");
	resultsF->mkdir("PCAL");
	resultsF->mkdir("Results");
	resultsF->mkdir("ResultsThetaPhi");
	resultsF->mkdir("EID");
	resultsF->mkdir("DCFid");
	resultsF->mkdir("Resolution");
	resultsF->mkdir("folding");
	resultsF->mkdir("Unfolding");

	resultsF->mkdir("UnfoldingS1");
	resultsF->mkdir("UnfoldingS2");
	resultsF->mkdir("UnfoldingS3");
	resultsF->mkdir("UnfoldingS4");
	resultsF->mkdir("UnfoldingS5");
	resultsF->mkdir("UnfoldingS6");
	
	const size_t nCutSys = 11; // the two last is for FTOF and HTCC 
	for (size_t iCutSys = 0; iCutSys < nCutSys; iCutSys ++){
		for (size_t iSec  = 0; iSec < 6; iSec ++){
			for (size_t iLevel  = 0; iLevel < 3; iLevel ++){
				string unfName = "UNF_CUT_" + to_string(iCutSys) + "_S_" + to_string(iSec)  + "_L_" + to_string(iLevel);
				resultsF->mkdir(unfName.c_str());
			}
		}
	}
	
	
	
// added by Valerii
	resultsF->mkdir("SectorDependences");
	resultsF->mkdir("Systematics");
    
    
/////////////////////////////////////////////////////////////////////////////////        
////////////// defining histograms: /////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////// 
    
/////////////////////////////////////////////////////////////////////////////////        
////////////// overview ///////////////////////////////////////////////////////// 
/////////////////////////////////////////////////////////////////////////////////        
	TH2F *wQ2Inclusive;
	TH2F *wQ2InclusiveGen;
	
	
	// Theta, P, W distributions 
	
	const size_t secN = 6;
    
    TH1F *thetaGen[secN];
    TH1F *momGen[secN];
    TH1F *thetaRec[secN];
    TH1F *momRec[secN];
    
    TH1F *thetaRec_noCut[secN];
    TH1F *momRec_noCut[secN];
    
    TH1F *Q2rec_noCut[secN];
    TH1F *Q2rec_Cut[secN];
    
    TH1F *phiRec_noCut[secN];
    TH1F *phiRec_Cut[secN];
    
    

    
    const double thMin =5., thMax = 36;
    const double pMin = 1.5, pMax = 11;
    
    const double q2MinNorm = 0.5, q2MaxNorm = 11;
    const double phiMinNorm = -80, phiMaxNorm = 80;
    
    const size_t nBinsNorm = 100;
    const size_t nBinsNormAlot = 200;
    
    
    
    
    for (size_t iSec = 0 ; iSec < secN; iSec++){
        string endName = "S" + to_string(iSec);
        thetaGen[iSec] = new TH1F(("thetaGen"+ endName).c_str(), ("thetaGen"+ endName).c_str(), nBinsNorm, thMin, thMax);
        thetaRec[iSec] = new TH1F(("thetaRec"+ endName).c_str(), ("thetaRec"+ endName).c_str(), nBinsNorm, thMin, thMax);
        thetaRec_noCut[iSec] = new TH1F(("thetaRec_noCut"+ endName).c_str(), ("thetaRec_noCut"+ endName).c_str(), 
                                        nBinsNorm, thMin, thMax);
        momGen[iSec] = new TH1F(("momGen"+ endName).c_str(), ("momGen"+ endName).c_str(), nBinsNorm, pMin, pMax);
        momRec[iSec] = new TH1F(("momRec"+ endName).c_str(), ("momRec"+ endName).c_str(), nBinsNorm, pMin, pMax);
        momRec_noCut[iSec] = new TH1F(("momRec_noCut"+ endName).c_str(), ("momRec_noCut"+ endName).c_str(), 
                                      nBinsNorm, pMin, pMax);
                                      
                                      
        Q2rec_noCut[iSec] = new TH1F(("Q2rec_noCut"+ endName).c_str(), ("Q2rec_noCut"+ endName).c_str(), nBinsNormAlot, q2MinNorm, q2MaxNorm);
        Q2rec_Cut[iSec] = new TH1F(("Q2rec_Cut"+ endName).c_str(), ("Q2rec_Cut"+ endName).c_str(), nBinsNormAlot, q2MinNorm, q2MaxNorm);
        
        phiRec_noCut[iSec] = new TH1F(("phiRec_noCut"+ endName).c_str(), ("phiRec_noCut"+ endName).c_str(), nBinsNormAlot, phiMinNorm, phiMaxNorm);
        phiRec_Cut[iSec] = new TH1F(("phiRec_Cut"+ endName).c_str(), ("phiRec_Cut"+ endName).c_str(), nBinsNormAlot, phiMinNorm, phiMaxNorm);
        
    }
    
   TH1D *wQ2_norm[6][qBinMax];
   TH1D *wQ2_norm_noCut[6][qBinMax];
   
   const double wMinNorm = 0.7, wMaxNorm = 3.3;
    
    for (int s = 0; s < 6; s++){
		for (int q = 0; q < qBinMax; q++){	
			string endName = "S" + to_string(s) + "Q" + to_string(q);
			wQ2_norm[s][q] = new TH1D(("wQ2_norm"+ endName).c_str(), ("wQ2_norm"+ endName).c_str(), nBinsNorm, wMinNorm, wMaxNorm);   
			wQ2_norm_noCut[s][q] = new TH1D(("wQ2_norm_noCut"+ endName).c_str(), ("wQ2_norm_noCut"+ endName).c_str(), nBinsNorm, wMinNorm, wMaxNorm);   
		}
	}
	
	
	TH2F *WvsQ2_noCut = new TH2F("WvsQ2_noCut", "WvsQ2_noCut", nBinsNormAlot, wMinNorm, wMaxNorm, nBinsNormAlot, q2MinNorm, q2MaxNorm);
	TH2F *WvsQ2_Cut = new TH2F("WvsQ2_Cut", "WvsQ2_Cut", nBinsNormAlot, wMinNorm, wMaxNorm, nBinsNormAlot, q2MinNorm, q2MaxNorm);
	///////////////////////
	
	
    
///// differen bin size effect and Bin purity: ///       
    vector<size_t> wNumBim = {size_t((highBorderW-lowBorderW)/0.03), nWBins, size_t((highBorderW-lowBorderW)/0.02), size_t((highBorderW-lowBorderW)/0.01)};
    const size_t wNumBinSize = wNumBim.size();
    TH1D *WQ2_WRecAndWGen[6][wNumBinSize][qBinMax];
    TH1D *wQ2BinSectorLogQ2StrictFidAllCutsVaryBinSize[6][wNumBinSize][qBinMax];
    TH1D *wGenVaryBinSize[6][wNumBinSize][qBinMax];
    
// for simulation comparasion plots:
	TH1D *CosthetaGen[6];
	TH1D *phiGen[6];
	TH1D *pGen[6];
    
    TH1F *simPecomp = new TH1F("simPecomp", "simPecomp", 200, 1.2,10.5);
    TH1F *simQ2comp = new TH1F("simQ2comp", "simQ2comp", 400, 0.5,10.5);
    TH1F *simWcomp = new TH1F("simWcomp", "simWcomp", 500, 0.8,4.1);
    TH1F *simPHIcomp = new TH1F("simPHIcomp", "simPHIcomp", 500, -200,200);
    TH1F *simTHETAcomp = new TH1F("simTHETAcomp", "simTHETAcomp", 500, 0,45);
    
    TH2F *simQ2Wcomp = new TH2F("simQ2Wcomp", "simQ2Wcomp", 400,0.8, 4.1, 800, 0.5, 11);
    TH2F *simPTHETHAcomp = new TH2F("simPTHETHAcomp", "simPTHETHAcomp", 200, 1.2, 10, 200, 0, 45);
    TH2F *simPTHETHAcompBefore = new TH2F("simPTHETHAcompBefore", "simPTHETHAcompBefore", 200, 1.2, 10, 200, 0, 45);
    TH2F *costhetaVSpGen = new TH2F("costhetaVSpGen", "costhetaVSpGen", 200, 1.2, 10, 100, 0.5, 1);
    
    TH2F *simPPHIcomp = new TH2F("simPPHIcomp", "simPPHIcomp", 200, 1.2, 10, 200, -200, 200);
    
// triag cut studies:
	TH2F *triagCut[6][12][3];
// sf cut studies:
	TH2F *sfCutProb[6][3];
    
// Empty Target Contribution runs studies:
    TH2F *SFempty[6][2];
    TH1F *VxEmpty[6][2];
    TH1F *VyEmpty[6][2];
    TH1F *VzEmpty[6][2];
    TH1F *NPEempty[6][2];
    TH1F *Wempty[6][2];
    
//FTOF:
	TH1F *ftofNoCuts = new TH1F("ftofNoCuts", "ftofNoCuts", 65, -0.5,64.5);
	TH1F *ftofAfterCuts = new TH1F("ftofAfterCuts", "ftofAfterCuts", 65, -0.5,64.5);
    
    
    
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////    Folding - Unfolding (It is fun :)     ///////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    const size_t foldBinNumber = 8;
    //const double wBinSize_fold = 0.002;// in GeV
    const double wBinSize_fold_1D = 0.01;// in GeV
    const double wBinSize_fold_2D = 0.05;// in GeV
    //const double wMin_fold = 1.025, wMax_fold = 2.725;
    const double wMin_fold = 0.875, wMax_fold = 2.875;
    
    //const size_t q2foldBinN = 10; // hardcoded  cause they are log, carefully change min max and mumber.
    const size_t q2foldBinN = 11; // hardcoded  cause they are log, carefully change min max and mumber.
    //const double q2min_Fold = 2.186724, q2max_Fold = 10.456396;
    const double q2min_Fold = 2.186724, q2max_Fold = 12.227588;//11
    //const double q2min_Fold = 1.869972, q2max_Fold = 12.227588;//12
    
    const size_t wBins_fold_1D = (wMax_fold - wMin_fold) / wBinSize_fold_1D;
    const size_t wBins_fold_2D = (wMax_fold - wMin_fold) / wBinSize_fold_2D;
    
    const size_t MatrixDimens2D = q2foldBinN * wBins_fold_2D;
    
    binFold binFoldingParams = {wMin_fold, wMax_fold, wBinSize_fold_2D, q2Min, deltaQ2};
    
    //log scale:
    const double q2BinSize_fold_1D = 1.;// in GeV
    const double q2Min_fold = 2., q2Max_fold = 12; //Q2 deconvolution
    const size_t q2Bins_fold_1D = 10;//(q2Max_fold - q2Min_fold) / q2BinSize_fold_1D;
    
    const size_t nQ2binsLog = 17;
    const double q2min_Fold_log = 0, q2max_Fold_log = nQ2binsLog * deltaQ2; 
    
    //////////////////////////////////// 1D case: ///////////////////////
	TH1D *wGEN_fold;
	TH1D *wREC_fold;// wil work ar a REC and a DATA depends on input file
	TH2D *responseMatrix;
    
	TH1D *xini_1D_gen[6][qBinMax];
	TH1D *bini_1D_rec[6][qBinMax];
	TH2D *Adet_1D_genVSrec[6][qBinMax];
    
	TH1D *xini_1D_genQ2[6][wBins_fold_2D];
	TH1D *bini_1D_recQ2[6][wBins_fold_2D];
	TH2D *Adet_1D_genVSrecQ2[6][wBins_fold_2D];
    
	TH1D *xini_1D_genQ2_orig[6][wBins_fold_2D];
	TH1D *bini_1D_recQ2_orig[6][wBins_fold_2D];
	TH2D *Adet_1D_genVSrecQ2_orig[6][wBins_fold_2D];
    
    wGEN_fold = new TH1D("wGEN_fold", "wGEN_fold", wBins_fold_1D, wMin_fold, wMax_fold);
    wREC_fold = new TH1D("wREC_fold", "wREC_fold", wBins_fold_1D, wMin_fold, wMax_fold);
    responseMatrix = new TH2D("responseMatrix", "responseMatrix", wBins_fold_1D, wMin_fold, wMax_fold, wBins_fold_1D, wMin_fold, wMax_fold);
    
	for (int s = 0; s < 6; s++){
        for (int q = 0; q < static_cast<int>(qBinMax); q++){
            string hist_name_end = "_S" + to_string(s) + "_Q2_" + to_string(q);
            xini_1D_gen[s][q] = new TH1D(("xini_1D_gen" + hist_name_end).c_str(),("xini_1D_gen" + hist_name_end).c_str(),
                                               wBins_fold_1D, wMin_fold, wMax_fold);
            bini_1D_rec[s][q] = new TH1D(("bini_1D_rec" + hist_name_end).c_str(),("bini_1D_rec" + hist_name_end).c_str(),
                                               wBins_fold_1D, wMin_fold, wMax_fold);
            
            Adet_1D_genVSrec[s][q] = new TH2D(("Adet_1D_genVSrec" + hist_name_end).c_str(), ("Adet_1D_genVSrec" + hist_name_end).c_str(), 
                                        wBins_fold_1D, wMin_fold, wMax_fold, wBins_fold_1D, wMin_fold, wMax_fold);
            
        }
        for (int iW = 0; iW < static_cast<int>(wBins_fold_2D); iW++){
            string hist_name_end = "_S" + to_string(s) + "_W_" + to_string(iW);
            xini_1D_genQ2[s][iW] = new TH1D(("xini_1D_genQ2" + hist_name_end).c_str(),("xini_1D_genQ2" + hist_name_end).c_str(),
                                               q2Bins_fold_1D, q2Min_fold, q2Max_fold);
            bini_1D_recQ2[s][iW] = new TH1D(("bini_1D_recQ2" + hist_name_end).c_str(),("bini_1D_recQ2" + hist_name_end).c_str(),
                                               q2Bins_fold_1D, q2Min_fold, q2Max_fold);
            
            Adet_1D_genVSrecQ2[s][iW] = new TH2D(("Adet_1D_genVSrecQ2" + hist_name_end).c_str(), ("Adet_1D_genVSrecQ2" + hist_name_end).c_str(), 
                                        q2Bins_fold_1D, q2Min_fold, q2Max_fold, q2Bins_fold_1D, q2Min_fold, q2Max_fold);
            
            
            
            xini_1D_genQ2_orig[s][iW] = new TH1D(("xini_1D_genQ2_orig" + hist_name_end).c_str(),("xini_1D_genQ2_orig" + hist_name_end).c_str(),
                                               nQ2binsLog , q2min_Fold_log, q2max_Fold_log);
            bini_1D_recQ2_orig[s][iW] = new TH1D(("bini_1D_recQ2_orig" + hist_name_end).c_str(),("bini_1D_recQ2_orig" + hist_name_end).c_str(),
                                               nQ2binsLog , q2min_Fold_log, q2max_Fold_log);
            
            Adet_1D_genVSrecQ2_orig[s][iW] = new TH2D(("Adet_1D_genVSrecQ2_orig" + hist_name_end).c_str(), ("Adet_1D_genVSrecQ2_orig" + hist_name_end).c_str(), 
                                        nQ2binsLog , q2min_Fold_log, q2max_Fold_log, nQ2binsLog , q2min_Fold_log, q2max_Fold_log);
            
        }
    }
    
    /////////////////////////////// 2D case: /////////////////////////////////////
    
    
    RooUnfoldResponse *hResponseRooSys[nCutSys][6][3];
    
    for (size_t iCutSys = 0; iCutSys < nCutSys; iCutSys ++){
    	for (size_t iSec = 0; iSec < 6 ;iSec ++){ 
    		for (size_t iLevel  = 0; iLevel < 3; iLevel ++){
				hResponseRooSys[iCutSys][iSec][iLevel] = new RooUnfoldResponse(MatrixDimens2D, -0.5, MatrixDimens2D - 0.5); 
			}
		}
	}
    
    
	TH2D *responseMatrix2Dcase;
	TH1D *wGEN_fold2Dcase;
	TH1D *wREC_fold2Dcase;// wil work ar a REC and a DATA depends on input file
    
	TH1D *xini_2D_gen[6];
	TH1D *bini_2D_rec[6];
	TH2D *Adet_2D_genVSrec[6];
	
	TH1D *hMeas[6];
	TH1D *hMeasSys[nCutSys][6][3];
	
    for (size_t iCutSys = 0; iCutSys < nCutSys; iCutSys ++){
    	for (size_t iSec = 0; iSec < 6 ;iSec ++){ 
    		for (size_t iLevel  = 0; iLevel < 3; iLevel ++){
				string hist_name_end =  "_C_" + to_string(iCutSys) + "_S_" + to_string(iSec)  + "_L_" + to_string(iLevel);
				hMeasSys[iCutSys][iSec][iLevel] = new TH1D(("hMeasSys" + hist_name_end).c_str(), ("hMeasSys" + hist_name_end).c_str(), MatrixDimens2D, -0.5, MatrixDimens2D - 0.5); 
			}
		}
	}
	
	
	TH1D *hMeasX[6];
	TH1D *hMeasY[6];
	
	RooUnfoldResponse *hResponseRoo[6];
	RooUnfoldResponse *hResponseRooWeightHalf[6];
	RooUnfoldResponse *hResponseRooNoRes[6];
	
	RooUnfoldResponse *hResponseRooOneAndHalfRes[6];
	RooUnfoldResponse *hResponseRooMisakWeight[6];
	
	for (size_t iSec = 0; iSec < 6 ;iSec ++){ 
		hResponseRoo[iSec] = new RooUnfoldResponse(MatrixDimens2D, -0.5, MatrixDimens2D - 0.5); 
		hResponseRooWeightHalf[iSec] = new RooUnfoldResponse(MatrixDimens2D, -0.5, MatrixDimens2D - 0.5); 
		hResponseRooNoRes[iSec] = new RooUnfoldResponse(MatrixDimens2D, -0.5, MatrixDimens2D - 0.5); 
		hResponseRooOneAndHalfRes[iSec] = new RooUnfoldResponse(MatrixDimens2D, -0.5, MatrixDimens2D - 0.5); 
		hResponseRooMisakWeight[iSec] = new RooUnfoldResponse(MatrixDimens2D, -0.5, MatrixDimens2D - 0.5); 
	}
	
	for (int s = 0; s < 6; s++){
        string hist_name_end = "_S" + to_string(s);
        bini_2D_rec[s] = new TH1D(("bini_2D_rec" + hist_name_end).c_str(), ("bini_2D_rec" + hist_name_end).c_str(), MatrixDimens2D, -0.5, MatrixDimens2D - 0.5);
        xini_2D_gen[s] = new TH1D(("xini_2D_gen" + hist_name_end).c_str(), ("xini_2D_gen" + hist_name_end).c_str(), MatrixDimens2D, -0.5, MatrixDimens2D - 0.5);
        Adet_2D_genVSrec[s] = new TH2D(("Adet_2D_genVSrec" + hist_name_end).c_str(), ("Adet_2D_genVSrec" + hist_name_end).c_str()
                                       , MatrixDimens2D, -0.5, MatrixDimens2D - 0.5, MatrixDimens2D, -0.5, MatrixDimens2D - 0.5);
                                       
                             
                                       
    }
    
    
    for (int s = 0; s < 6; s++){
    	string hist_name_end = "_S" + to_string(s);
        hMeas[s] = new TH1D(("measAll" + hist_name_end).c_str(), ("Test Measured" + hist_name_end).c_str(), MatrixDimens2D, -0.5, MatrixDimens2D - 0.5); 
        hMeasX[s] = new TH1D(("measX" + hist_name_end).c_str(), ("Test Measured" + hist_name_end).c_str(), MatrixDimens2D, -0.5, MatrixDimens2D - 0.5); 
        hMeasY[s] = new TH1D(("measY" + hist_name_end).c_str(), ("Test Measured" + hist_name_end).c_str(), MatrixDimens2D, -0.5, MatrixDimens2D - 0.5); 
    }
    // X and Y just bin Number no physics:
    responseMatrix2Dcase = new TH2D("responseMatrix2Dcase", "responseMatrix2Dcase", MatrixDimens2D, -0.5, MatrixDimens2D - 0.5, 
                                    MatrixDimens2D, -0.5, MatrixDimens2D - 0.5);
    wGEN_fold2Dcase = new TH1D("wGEN_fold2Dcase", "wGEN_fold2Dcase", MatrixDimens2D, -0.5, MatrixDimens2D - 0.5);
    wREC_fold2Dcase = new TH1D("wREC_fold2Dcase", "wREC_fold2Dcase", MatrixDimens2D, -0.5, MatrixDimens2D - 0.5);
    
    
    // Gen for test:
    TH1D *wGEN_testBefWeight[6];
    TH1D *wGEN_testAftWeight[6];
    TH1D *wGEN_testAftWeight_NoRes[6];
    
   TH1D * wGEN_testAftWeight_Misak[6];
    
    TH1D *wGEN_testBefWeightMissed[6];
    TH1D *wGEN_testAftWeightMissed[6];
    TH1D *wGEN_testAftWeight_NoResMissed[6];
    
    for (size_t iSec = 0; iSec < 6; iSec++){
    	const string endName = "_S" + to_string(iSec);
    	wGEN_testBefWeight[iSec] = new TH1D(("wGEN_testBefWeight" + endName).c_str(),("wGEN_testBefWeight" + endName).c_str(), 
    						10 * nWBins, lowBorderW, highBorderW); 
    	wGEN_testAftWeight[iSec] = new TH1D(("wGEN_testAftWeight" + endName).c_str(),("wGEN_testAftWeight" + endName).c_str(), 
    						10 * nWBins, lowBorderW, highBorderW); 
    						
    	wGEN_testAftWeight_NoRes[iSec] = new TH1D(("wGEN_testAftWeight_NoRes" + endName).c_str(),("wGEN_testAftWeight_NoRes" + endName).c_str(), 
    						10 * nWBins, lowBorderW, highBorderW);
    						
    	wGEN_testAftWeight_Misak[iSec] = new TH1D(("wGEN_testAftWeight_Misak" + endName).c_str(),("wGEN_testAftWeight_Misak" + endName).c_str(), 
    						10 * nWBins, lowBorderW, highBorderW); 					
    						
    	wGEN_testBefWeightMissed[iSec] = new TH1D(("wGEN_testBefWeightMissed" + endName).c_str(),("wGEN_testBefWeightMissed" + endName).c_str(), 
    						10 * nWBins, lowBorderW, highBorderW); 
    	wGEN_testAftWeightMissed[iSec] = new TH1D(("wGEN_testAftWeightMissed" + endName).c_str(),("wGEN_testAftWeightMissed" + endName).c_str(), 
    						10 * nWBins, lowBorderW, highBorderW); 
    						
    	wGEN_testAftWeight_NoResMissed[iSec] = new TH1D(("wGEN_testAftWeight_NoResMissed" + endName).c_str(),("wGEN_testAftWeight_NoResMissed" + endName).c_str(), 
    						10 * nWBins, lowBorderW, highBorderW); 
    }
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//2D and 1D using RooUnfold Deconvolutions:
	
	// measured D.
	TH2D *measured2D_secBySec[6];
	TH2D *measured2D_DATA[6];
	// Generated (not needed) only for binning scheme
	TH2D *truth2D_secBySec[6];
    //Valera asked to chck
    TH2D *for_Valera_Q2_VS_W[6];
    TH2D *for_Valera_Q2_VS_W_gen[6];
	
	//Actual 2D deconv. obj
	RooUnfoldResponse *rooResponse2D[6];
	
	//2D Q2 not equal binning:
	const Int_t nbins_q2Y_2D=11;
	const Double_t q2Ybins_2D[nbins_q2Y_2D+1] = {2.186724, 2.55713, 2.990278, 3.496797, 4.089114, 4.781762, 5.591738, 6.538914, 7.64653, 8.941764, 10.456396, 12.227588};
	
	for (size_t iSec = 0; iSec < 6 ;iSec ++){ 
	
		const string endName = "_S" + to_string(iSec);
		
		measured2D_DATA[iSec] = new TH2D(("measured2D_DATA" + endName).c_str(),("measured2D_DATA" + endName).c_str(), 
    						wBins_fold_2D, wMin_fold, wMax_fold, nbins_q2Y_2D, q2Ybins_2D); 
    						
		measured2D_secBySec[iSec] = new TH2D(("Sep_measured2D_my_secBySec" + endName).c_str(),("Sep_measured2D_my_secBySec" + endName).c_str(), 
    						wBins_fold_2D, wMin_fold, wMax_fold, nbins_q2Y_2D, q2Ybins_2D); 
    	truth2D_secBySec[iSec] = new TH2D(("Sep_truth2D_my_secBySec" + endName).c_str(),("Sep_truth2D_my_secBySec" + endName).c_str(), 
    						wBins_fold_2D, wMin_fold, wMax_fold, nbins_q2Y_2D, q2Ybins_2D); 
    						
            // НОВОЕ: что реально идёт в response-матрицу (Q2 vs W)
        for_Valera_Q2_VS_W[iSec] = new TH2D(
            ("for_Valera_Q2_VS_W" + endName).c_str(),
            ("for_Valera_Q2_VS_W" + endName).c_str(),
            wBins_fold_2D, wMin_fold, wMax_fold,        // X: W
            nbins_q2Y_2D, q2Ybins_2D                    // Y: Q2
        );

        for_Valera_Q2_VS_W_gen[iSec] = new TH2D(
            ("for_Valera_Q2_VS_W_gen_MISS" + endName).c_str(),
            ("for_Valera_Q2_VS_W_gen_MISS" + endName).c_str(),
            wBins_fold_2D, wMin_fold, wMax_fold,        // X: W_gen
            nbins_q2Y_2D, q2Ybins_2D                    // Y: Q2_gen
        );
        /////////////////////////////////////////////////////////////////////////////////////////////
        
    	rooResponse2D[iSec] = new RooUnfoldResponse(measured2D_secBySec[iSec], truth2D_secBySec[iSec], ("R2D" + endName).c_str()); 					
    	
	}
	
	//1D for SVD:
	TH1D *measured1D_secBySec[6][q2foldBinN];
	TH1D *measured1D_DATA[6][q2foldBinN];
    TH1D *truth1D_secBySec[6][q2foldBinN];
    TH1D *generated1D_DATA_sec[6][q2foldBinN];
    TH1D *generated1D_DATA[q2foldBinN];
    
    //Actual 1D deconv. obj
	RooUnfoldResponse *rooResponse1D[6][q2foldBinN];
	
	for (size_t iSec = 0; iSec < 6 ;iSec ++){ 
		for (size_t iQ2 = 0; iQ2 < q2foldBinN ;iQ2 ++){ 
	
			const string endName = "_S" + to_string(iSec) + "_Q2_" + to_string(iQ2);
			
			measured1D_secBySec[iSec][iQ2] = new TH1D(("Sep_measured1D_secBySec" + endName).c_str(),("Sep_measured1D_secBySec" + endName).c_str(), 
								wBins_fold_2D, wMin_fold, wMax_fold); 
								
			measured1D_DATA[iSec][iQ2] = new TH1D(("measured1D_DATA" + endName).c_str(),("measured1D_DATA" + endName).c_str(), 
								wBins_fold_2D, wMin_fold, wMax_fold); 
								
			generated1D_DATA_sec[iSec][iQ2] = new TH1D(("generated1D_DATA_sec" + endName).c_str(),("generated1D_DATA_sec" + endName).c_str(), 
								wBins_fold_2D, wMin_fold, wMax_fold); 
								
			if (iSec == 0){
			const string endName2 = "_Q2_" + to_string(iQ2);
				generated1D_DATA[iQ2] = new TH1D(("generated1D_DATA" + endName2).c_str(),("generated1D_DATA" + endName2).c_str(), 
								wBins_fold_2D, wMin_fold, wMax_fold); 
			}
								
			truth1D_secBySec[iSec][iQ2] = new TH1D(("Sep_truth1D_secBySec" + endName).c_str(),("Sep_truth1D_secBySec" + endName).c_str(), 
								wBins_fold_2D, wMin_fold, wMax_fold); 
								
			rooResponse1D[iSec][iQ2] = new RooUnfoldResponse(wBins_fold_2D, wMin_fold, wMax_fold, ("R1D" + endName).c_str());
		}
	}
	/// Folders creation:
	
	for (size_t iSec = 0; iSec < 6 ;iSec ++){ 
		const string unfName2D = "Deconv2D/S" + to_string(iSec);
		resultsF->mkdir(unfName2D.c_str());
		const string unfName2D_2 = "Deconv2D/S" + to_string(iSec) +"_my";
		resultsF->mkdir(unfName2D_2.c_str());
		for (size_t iQ2 = 0; iQ2 < q2foldBinN ;iQ2 ++){ 
			const string unfName1D = "Deconv1D/S" + to_string(iSec) + "/Q2_" + to_string(iQ2);
			resultsF->mkdir(unfName1D.c_str());
			const string unfName1D_2 = "Deconv1D/S" + to_string(iSec) + "/Q2_" + to_string(iQ2) +"_my";
			resultsF->mkdir(unfName1D_2.c_str());
		}
	}
    
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////          END Deconvolution           /////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    TH2F *sfCut_Sectors_bef[6];
	TH2F *sfCut_Sectors_aft[6];
    
    
    for (size_t iSec = 0 ; iSec < 6; iSec++){
        string hNameSF1 = "sfCut_Sectors_bef_S" + to_string(iSec);
        string hNameSF2 = "sfCut_Sectors_aft_S" + to_string(iSec);
        sfCut_Sectors_bef[iSec] = new TH2F(hNameSF1.c_str(),hNameSF1.c_str(), 80, 0., 5., 200, 0.1, 0.35);
        sfCut_Sectors_aft[iSec] = new TH2F(hNameSF2.c_str(),hNameSF2.c_str(), 80, 0., 5., 200, 0.1, 0.35);
     }
    
    
//detectors bad elements studies:
    
    const int nXbinDetec = 450, nYbinDetec = 450; // was 150
    
    // DC no fid cuts
	TH2F *xYDCR1detec = new TH2F("xYDCR1detec", "xYDCR1detec", nXbinDetec, -180, 180, nYbinDetec, -180, 180);
	TH2F *xYDCR2detec = new TH2F("xYDCR2detec", "xYDCR2detec", nXbinDetec, -250, 250, nYbinDetec, -350, 250);
	TH2F *xYDCR3detec = new TH2F("xYDCR3detec", "xYDCR3detec", nXbinDetec, -400, 400, nYbinDetec, -400, 400);
	
	TH2F *xYDpcalDetec = new TH2F("xYDpcalDetec", "xYDpcalDetec", nXbinDetec, -400, 400, nYbinDetec, -400, 400);
    
	TH2F *xYDhtccDetec = new TH2F("xYDhtccDetec", "xYDhtccDetec", nXbinDetec, -100, 100, nYbinDetec, -100, 100);
	TH2F *xYDftofDetec = new TH2F("xYDftofDetec", "xYDftofDetec", nXbinDetec, -400, 400, nYbinDetec, -400, 400);
	
	TH2F *xYDecinDetec = new TH2F("xYDecinDetec", "xYDecinDetec", nXbinDetec, -400, 400, nYbinDetec, -400, 400);
	TH2F *xYDecoutDetec = new TH2F("xYDecoutDetec", "xYDecoutDetec", nXbinDetec, -400, 400, nYbinDetec, -400, 400);
    
    // DC with fid cuts
	TH2F *xYDCR1detecKnock = new TH2F("xYDCR1detecKnock", "xYDCR1detecKnock", nXbinDetec, -180, 180, nYbinDetec, -180, 180);
	TH2F *xYDCR2detecKnock = new TH2F("xYDCR2detecKnock", "xYDCR2detecKnock", nXbinDetec, -250, 250, nYbinDetec, -250, 250);
	TH2F *xYDCR3detecKnock = new TH2F("xYDCR3detecKnock", "xYDCR3detecKnock", nXbinDetec, -400, 400, nYbinDetec, -400, 400);
	
	TH2F *xYDpcalDetecKnock = new TH2F("xYDpcalDetecKnock", "xYDpcalDetecKnock", nXbinDetec, -400, 400, nYbinDetec, -400, 400);
	
	TH2F *xYDhtccDetecKnock = new TH2F("xYDhtccDetecKnock", "xYDhtccDeteKnockc", nXbinDetec, -100, 100, nYbinDetec, -100, 100);
	TH2F *xYDftofDetecKnock = new TH2F("xYDftofDetecKnock", "xYDftofDetecKnock", nXbinDetec, -400, 400, nYbinDetec, -400, 400);
	
	TH2F *xYDecinDetecKnock = new TH2F("xYDecinDetecKnock", "xYDecinDetecKnock", nXbinDetec, -400, 400, nYbinDetec, -400, 400);
	TH2F *xYDecoutDetecKnock = new TH2F("xYDecoutDetecKnock", "xYDecoutDetecKnock", nXbinDetec, -400, 400, nYbinDetec, -400, 400);    
    
    
	TH2F *xYDCR1detecNo_BadEl_PhiSp = new TH2F("xYDCR1detecNo_BadEl_PhiSp", "xYDCR1detecNo_BadEl_PhiSp", nXbinDetec, -180, 180, nYbinDetec, -180, 180);
	TH2F *xYDCR2detecNo_BadEl_PhiSp = new TH2F("xYDCR2detecNo_BadEl_PhiSp", "xYDCR2detecNo_BadEl_PhiSp", nXbinDetec, -250, 250, nYbinDetec, -250, 250);
	TH2F *xYDCR3detecNo_BadEl_PhiSp = new TH2F("xYDCR3detecNo_BadEl_PhiSp", "xYDCR3detecNo_BadEl_PhiSp", nXbinDetec, -400, 400, nYbinDetec, -400, 400);
    
	TH2F *xYDhtccDetecNo_BadEl_PhiSp = new TH2F("xYDhtccDetecNo_BadEl_PhiSp", "xYDhtccDetecNo_BadEl_PhiSp", nXbinDetec, -100, 100, nYbinDetec, -100, 100);
	TH2F *xYDftofDetecNo_BadEl_PhiSp = new TH2F("xYDftofDetecNo_BadEl_PhiSp", "xYDftofDetecNo_BadEl_PhiSp", nXbinDetec, -400, 400, nYbinDetec, -400, 400);
	TH2F *xYDpcalDetecNo_BadEl_PhiSp = new TH2F("xYDpcalDetecNo_BadEl_PhiSp", "xYDpcalDetecNo_BadEl_PhiSp", nXbinDetec, -400, 400, nYbinDetec, -400, 400);
	TH2F *xYDecinDetecNo_BadEl_PhiSp = new TH2F("xYDecinDetecNo_BadEl_PhiSp", "xYDecinDetecNo_BadEl_PhiSp", nXbinDetec, -400, 400, nYbinDetec, -400, 400);
	TH2F *xYDecoutDetecNo_BadEl_PhiSp = new TH2F("xYDecoutDetecNo_BadEl_PhiSp", "xYDecoutDetecNo_BadEl_PhiSp", nXbinDetec, -400, 400, nYbinDetec, -400, 400);  
    
//DC phi distributions:
    
    const int nThetaBins =30;
    const int npBins = 3;
    TH1F * countDCthetaBF[6][nThetaBins][npBins];
    TH1F * countDCthetaAF[6][nThetaBins][npBins];    
    TH1F * countDCthetaBFnoP[6][nThetaBins];
    TH1F * countDCthetaAFnoP[6][nThetaBins];
    TH1F * countDCthetaBFstephan[6][nThetaBins];
    
///resolution:
    TH1D *wResQ2[qBinMax];
    TH1D *wResAllBins[qBinMax][nWBins];
    TH1D *q2ResAllBins[qBinMax][nWBins];  
    TH1D *q2Res;
    
/////////////////////////////////////////////////////////////////////////////////        
////////////// DCFid folder /////////////////////////////////////////////////////         
/////////////////////////////////////////////////////////////////////////////////   
    
// DC fid rotated:
// sec, layer,cut on/off
    TH2F *xyDC_noWeight[6][3][2];
    TH2F *xyDC_nChi2[6][3][2];
    TH2F *xyDC_ndf[6][3][2];
    
    TH2F *xyPCAL_noWeight[6][2];
    TH2F *xyECOUT_noWeight[6][2];

/////////////////////////////////////////////////////////////////////////////////        
////////////// results //////////////////////////////////////////////////////////    
/////////////////////////////////////////////////////////////////////////////////
    
//////////////  Nick's CS calculatuion, keep for self check
	TH1D *wQ2BinSectorLogQ2AllCutsSFSys[6][qBinMax][3];
	TH1D *wQ2BinSectorAccLogQ2AllCutsSFSys[6][qBinMax][3];
    
	TH1D *wQ2BinSectorLogQ2StrictFidAllCuts[6][qBinMax];
    
    
// acceptance
	TH1D *wQ2BinSectorGenLogQ2[6][qBinMax];
	TH1D *wQ2BinSectorAccLogQ2StrictFidAllCuts[6][qBinMax];
       
//systematics:    
    const int nCut = 14; // 0 - 12 cuts, 13 - all
	TH1D *wQ2BinSectorLogQ2AllCutsAllSys[nCut][6][qBinMax][3];
	TH1D *wQ2BinSectorAccLogQ2AllCutsAllSys[nCut][6][qBinMax][3];
/////////////////////////////////////////////////////////////////////////////////    
////////////// Creating histograms //////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

// w vs Q2 for sim and general
	float lowBorderQ2 = 0.0;
	float highBorderQ2 = 10.4;
	int nQ2Bins = 250;
    
	wQ2Inclusive = new TH2F("wQ2Inclusive", "wQ2Inclusive", 5*nWBins, lowBorderW, highBorderW, nQ2Bins, lowBorderQ2, highBorderQ2);
	wQ2InclusiveGen = new TH2F("wQ2InclusiveGen", "wQ2InclusiveGen", 5*nWBins, lowBorderW, highBorderW, 4 * nQ2Bins, lowBorderQ2, highBorderQ2);
    TH2F* pVSthetaMomCorr =  new TH2F("pVSthetaMomCorr", "pVSthetaMomCorr", 80, 2., 11., 80, 5, 35.);
    
	TH1D *wInclusiveBef[6];
	TH1D *wInclusiveAft[6];
	TH2D *wPInclusiveBef[6];
    
    
    for (size_t iSec = 0 ; iSec < 6; iSec++){
        string hName1 = "wInclusiveBefS" + to_string(iSec);
        string hName2 = "wInclusiveAftS" + to_string(iSec);
        wInclusiveBef[iSec] = new TH1D(hName1.c_str(),hName1.c_str(), 80, 0.8, 1.2);
        wInclusiveAft[iSec] = new TH1D(hName2.c_str(),hName2.c_str(), 80, 0.8, 1.2);
        
        string hName3 = "wPInclusiveBefS" + to_string(iSec);
        wPInclusiveBef[iSec] = new TH2D(hName3.c_str(),hName3.c_str(), 7, 4, 11, 30,  0.7, 1.4);
        
    }
// triag and SF cut studies:
    
    for (size_t iSec = 0 ; iSec < 6; iSec++){
        for (size_t iCut = 0; iCut < 3; iCut++){ 
            for  (size_t iPbin = 0; iPbin < 12; iPbin++){
                string hName = "triagS" + to_string(iSec) + "pB" + to_string(iPbin) + 'C' + to_string(iCut);
                triagCut[iSec][iPbin][iCut] = new TH2F(hName.c_str(),hName.c_str(), 200, 0, 0.35, 200, 0, 0.35);
            }
            string hNameSF = "sfS" + to_string(iSec) + 'C' + to_string(iCut);
            sfCutProb[iSec][iCut] = new TH2F(hNameSF.c_str(),hNameSF.c_str(), 100, 2., 10.4, 200, 0.1, 0.35);
        }
    }
    
    
// Empty Target Contribution runs studies:    
    for (size_t iSec = 0 ; iSec < 6; iSec++){
        for (size_t iCut = 0; iCut < 2; iCut++){
            string endName = "S" + to_string(iSec) + "C" + to_string(iCut);
            SFempty[iSec][iCut] = new TH2F(("SFempty" + endName).c_str(), ("SFempty" + endName).c_str(), 500, 0, 10, 250, 0, 0.5);
            Wempty[iSec][iCut] = new TH1F(("Wempty" + endName).c_str(), ("Wempty" + endName).c_str(),  nWBins, lowBorderW, highBorderW);     
            VxEmpty[iSec][iCut] = new TH1F(("VxEmpty" + endName).c_str(), ("VxEmpty" + endName).c_str(),  100, -15, 15);     
            VyEmpty[iSec][iCut] = new TH1F(("VyEmpty" + endName).c_str(), ("VyEmpty" + endName).c_str(),  100, -15, 15);     
            VzEmpty[iSec][iCut] = new TH1F(("VzEmpty" + endName).c_str(), ("VzEmpty" + endName).c_str(),  100, -15, 15);     
            NPEempty[iSec][iCut] = new TH1F(("NPEempty" + endName).c_str(), ("NPEempty" + endName).c_str(),  50, 0, 50);     
        }
    }
    
    

// DC + PCAL fid cuts:   
    int nBinsDCfid = 100;//was 1500
    int nBinsDCfidX = 120;
    int nBinsDCfidY = 120;
    double dcFidXmin = -180, dcFidXmax = 50;
    double dcFidYmin = -130, dcFidYmax = 130;
    
    TH2F* pcalFID_SF_V[2][6];
    TH2F* pcalFID_SF_W[2][6];
    TH2F* pcalFID_SF_U[2][6];
    
    for (size_t isCut = 0 ; isCut < 2; isCut++){
    	for (size_t iSec = 0 ; iSec < 6; iSec++){
    		string endName = "S" + to_string(iSec) + 'C' + to_string(isCut);
    	    pcalFID_SF_V[isCut][iSec] = new TH2F(("pcalFID_SF_V" + endName).c_str(), ("pcalFID_SF_V" + endName).c_str(), 
                                                              1000, 0, 500, 150, 0.1, 0.5);
            pcalFID_SF_W[isCut][iSec] = new TH2F(("pcalFID_SF_W" + endName).c_str(), ("pcalFID_SF_W" + endName).c_str(), 
                                                              1000, 0, 500, 150, 0.1, 0.5);
            pcalFID_SF_U[isCut][iSec] = new TH2F(("pcalFID_SF_U" + endName).c_str(), ("pcalFID_SF_U" + endName).c_str(), 
                                                              1000, 0, 500, 150, 0.1, 0.5);
    	}
    }
    
    for (size_t iSec = 0 ; iSec < 6; iSec++){
        for (size_t iLayer = 0; iLayer < 3; iLayer++){
            for (size_t isCut = 0; isCut < 2; isCut++){
                string endName = "S" + to_string(iSec) + 'R' + to_string(iLayer) + 'C' + to_string(isCut);

                xyDC_noWeight[iSec][iLayer][isCut] = new TH2F(("xyDC_noWeight" + endName).c_str(), ("xyDC_noWeight" + endName).c_str(), 
                                                              nBinsDCfidX, dcFidXmin, dcFidXmax, nBinsDCfidY, dcFidYmin, dcFidYmax);
                xyDC_nChi2[iSec][iLayer][isCut] = new TH2F(("xyDC_nChi2" + endName).c_str(), ("xyDC_nChi2" + endName).c_str(), 
                                                           nBinsDCfidX, dcFidXmin, dcFidXmax, nBinsDCfidY, dcFidYmin, dcFidYmax);
                xyDC_ndf[iSec][iLayer][isCut] = new TH2F(("xyDC_ndf" + endName).c_str(), ("xyDC_ndf" + endName).c_str(), 
                                                         nBinsDCfidX, dcFidXmin, dcFidXmax, nBinsDCfidY, dcFidYmin, dcFidYmax);
                if (iLayer == 0){
                    xyPCAL_noWeight[iSec][isCut] = new TH2F(("xyPCAL_noWeight" + endName).c_str(), ("xyPCAL_noWeight" + endName).c_str(), 
                                                            nBinsDCfidX, -270, 100, nBinsDCfidY, -150, 150);
                    xyECOUT_noWeight[iSec][isCut] = new TH2F(("xyECOUT_noWeight" + endName).c_str(), ("xyECOUT_noWeight" + endName).c_str(), 
                                                            nBinsDCfidX, -270, 100, nBinsDCfidY, -150, 150);                                        
                    
               }
            }
        }
    }

//DC phi distributions:
	for (int s = 0; s < 6; s++){
        for (int itheta = 0; itheta < nThetaBins; itheta++){
            string bfNameNoP = "countDCthetaBFnoPS" + to_string(s) + 'T' + to_string(itheta);
            string afNameNoP = "countDCthetaAFnoPS" + to_string(s) + 'T' + to_string(itheta);

            countDCthetaBFnoP[s][itheta] = new TH1F(bfNameNoP.c_str(), bfNameNoP.c_str(), 40, -40, 40);
            countDCthetaAFnoP[s][itheta] = new TH1F(afNameNoP.c_str(), afNameNoP.c_str(), 40, -40, 40);
            
            string bfNamestephan = "countDCthetaBFstefanS" + to_string(s) + 'T' + to_string(itheta);
            countDCthetaBFstephan[s][itheta] = new TH1F(bfNamestephan.c_str(), bfNamestephan.c_str(), 100, -50, 50);
            
            for (int ip = 0; ip < npBins; ip++){
                string bfName = "countDCthetaBFS" + to_string(s) + 'T' + to_string(itheta) + 'P' + to_string(ip);
                string afName = "countDCthetaAFS" + to_string(s) + 'T' + to_string(itheta) + 'P' + to_string(ip);

                countDCthetaBF[s][itheta][ip] = new TH1F(bfName.c_str(), bfName.c_str(), 100, -50, 50);
                countDCthetaAF[s][itheta][ip] = new TH1F(afName.c_str(), afName.c_str(), 100, -50, 50);
            }
        }
    }
//resolution in P 
    
	TH2F* resP_WQ2[qBinMax];
	TH2F* resW_WQ2[qBinMax];
	for (int q = 0; q < static_cast<int>(qBinMax); q++){
		string endName = "Q2_" + to_string(q);
		resP_WQ2[q] = new TH2F(("resP_WQ2" + endName).c_str(), ("resP_WQ2" + endName).c_str(), 34, 1.0, 2.7, 1000, -600., 600.);
		resW_WQ2[q] = new TH2F(("resW_WQ2" + endName).c_str(), ("resW_WQ2" + endName).c_str(), 34, 1.0, 2.7, 1000, -400., 400.);
	}
	
///resolution W and Q2:
    
    float dWabs = 0.3;
    float dQ2abs = 0.32;
    int nBinsDw = 80;
    
    for (size_t iq = 0; static_cast<int>(iq) < qBinMax; iq++){
        wResQ2[iq] = new TH1D(Form("wResQ2%d", static_cast<int>(iq)), Form("wResQ2%d", static_cast<int>(iq)), nBinsDw, 0 - dWabs, dWabs);
        for (size_t iW = 0; iW < nWBins; iW++){
            string endName = "Q2_" + to_string(iq) + "_W_" + to_string(iW); 
            string resQ2name = "q2ResAll" + endName;
            string resWname = "wAll" + endName;
            
            q2ResAllBins[iq][iW] = new TH1D(resQ2name.c_str(), resQ2name.c_str(), nBinsDw, 0 - dQ2abs, dQ2abs);
            wResAllBins[iq][iW] = new TH1D(resWname.c_str(), resWname.c_str(), nBinsDw, 0 - dWabs, dWabs);
        }
    }
    q2Res =  new TH1D("q2Res", "q2Res", nBinsDw, 0 - dQ2abs, dQ2abs); 
    
///bin purity
	for (int s = 0; s < 6; s++){
        for (int q = 0; q < static_cast<int>(qBinMax); q++){
            for (int t = 0; t < static_cast<int>(wNumBinSize); t++){
            
                WQ2_WRecAndWGen[s][t][q] = new TH1D(Form("WQ2_WRecAndWGenS%dQ%dB%d", s , q, t), Form("WQ2_WRecAndWGenS%dQ%d%d", s , q, t), wNumBim[t], lowBorderW, highBorderW);       
                wQ2BinSectorLogQ2StrictFidAllCutsVaryBinSize[s][t][q] = new TH1D(Form("wQ2BinSectorLogQ2StrictFidAllCutsVaryBinSizeS%dQ%dB%d", s , q, t),Form("wQ2BinSectorLogQ2StrictFidAllCutsVaryBinSizeS%dQ%dB%d", s , q, t), wNumBim[t], lowBorderW, highBorderW);  
                wGenVaryBinSize[s][t][q] = new TH1D(Form("wGenVaryBinSizeS%dQ%dB%d", s , q, t),Form("wGenVaryBinSizeS%dQ%dB%d", s , q, t), wNumBim[t], lowBorderW, highBorderW);
        
            }
        }
    }
    
//systematics:

    for (int iCut = 0; iCut < nCut; iCut++){	
        for (int iSec = 0; iSec < 6; iSec++){	
            for (int iQ2 = 0; iQ2 < qBinMax; iQ2++){
                for (int iCutLevel = 0; iCutLevel < 3; iCutLevel++){
                    string tmp_name = "wQ2BinSectorLogQ2AllCutsAllSysC" + to_string(iCut) + "S" + to_string(iSec) + "Q2" + to_string(iQ2) + "L" + to_string(iCutLevel);
                    string nameAcc = "wQ2BinSectorAccLogQ2AllCutsAllSysC" + to_string(iCut) + "S" + to_string(iSec) + "Q2" + to_string(iQ2) + "L" + to_string(iCutLevel);
                    wQ2BinSectorLogQ2AllCutsAllSys[iCut][iSec][iQ2][iCutLevel] = new TH1D(tmp_name.c_str(),tmp_name.c_str(), nWBins, lowBorderW, highBorderW);  
                    wQ2BinSectorAccLogQ2AllCutsAllSys[iCut][iSec][iQ2][iCutLevel] = new TH1D(nameAcc.c_str(),nameAcc.c_str(), nWBins, lowBorderW, highBorderW);  
                }
            }
        }
    }
///////////////////////////////////////////////////////////////////////////////////////////    
// Sector dependences:
//////////////////////////////////////////////////////////////////////////////////////////
    
    const size_t pNsd = 9, thNsd = 13, nCutCombination = 9; 
    const size_t wBins_ST = 8;
    const double minW_ST = 1.1, maxW_ST = 2.5;
    
    TH2F *phiVSw[nCutCombination][secN];
    TH2F *phiVSw_Q2bin[nCutCombination][secN][qBinMax];
    TH2F *phiVSq2_Wbin[nCutCombination][secN][wBins_ST];
    TH2F *phiVSq2[nCutCombination][secN];
    TH2F *wVSq2[nCutCombination][secN];
    
    TH2F *phiVSthFullRange[nCutCombination]; 
    TH2F *phiVSthFullRangeGen; 
    
    TH2F *thVSmom[nCutCombination][secN];
    TH2F *thVSphi[nCutCombination][secN];
    TH2F *phiVSmom[nCutCombination][secN];
    
    TH1F *thYield[nCutCombination][secN];
    TH1F *phiYield[nCutCombination][secN];
    TH1F *momYield[nCutCombination][secN];
    
    // theta distrib params:
    const int nThetaBinsSecStudies = 30, nThetaBinsSecStudies2D = 30;
    const double thetaMinSecStudies = 6., thetaMaxSecStudies = 36.;
    // phi distrib params:
    const int nPhiBinsSecStudies = 50, nPhiBinsSecStudies2D = 36;
    const double phiMinSecStudies = -30., phiMaxSecStudies = 42.;
    // mom distrib params:
    const int nMomBinsSecStudies = 120;
    const double momMinSecStudies = 2., momMaxSecStudies = 10.;
    // w and q2 params:
    const int nWBinsSecStudies = 60, nQ2BinsSecStudies2D = 70;
    const int nWBinsSecStudies_binsOVERq2 = 40, nQ2BinsSecStudies2D_binsOVERw = 40;
    const double wMinSecStudies = 1., wMaxSecStudies = 2.9;
    const double q2MinSecStudies = 2., q2MaxSecStudies = 10.4;
    
    
    //multipl. for 2D plots
    
    const float n2DtimesBinTheta = 4;
    const int n2DtimesBinPhi = 3;
    const int nPhiBinsSecStudies2D_noP = 60;
    const int nPhiBinsSecStudies2D_noP_wQ2bins = 25;
    
    TH1F *thYieldGen[secN];
    TH1F *phiYieldGen[secN];
    TH1F *momYieldGen[secN];
    
    TH2F *thVSphiGen[secN];
    TH2F *thVSmomGen[secN];
    TH2F *phiVSmomGen[secN];
    
    TH2F *phiVSwGen[secN];
    TH2F *phiVSw_Q2binGen[secN][qBinMax];
    TH2F *phiVSq2_WbinGen[secN][wBins_ST];
    TH2F *phiVSq2Gen[secN];
    TH2F *wVSq2Gen[secN];
    
    
    TH1F *phiYieldGenALL = new TH1F(("phiYieldGenALL"), ("phiYieldGenALL"), n2DtimesBinPhi, -200, 200);
    
    TH1F *wYieldPthet[nCutCombination][secN][pNsd][thNsd];
    TH1F *pYieldTh[nCutCombination][secN][thNsd];
///GEN:
    
    phiVSthFullRangeGen = new TH2F("phiVSthFullRangeGen", "phiVSthFullRangeGen", n2DtimesBinTheta*nThetaBinsSecStudies2D, thetaMinSecStudies, thetaMaxSecStudies, n2DtimesBinPhi * nPhiBinsSecStudies, phiMinSecStudies, phiMaxSecStudies);
    
    
    for(size_t iSec = 0 ; iSec < secN; iSec++){
        string endName = "S" + to_string(iSec);
        
        thYieldGen[iSec] = new TH1F(("thYieldGen"+ endName).c_str(), ("thYieldGen"+ endName).c_str(), nThetaBinsSecStudies, thetaMinSecStudies, thetaMaxSecStudies);
        phiYieldGen[iSec] = new TH1F(("phiYieldGen"+ endName).c_str(), ("phiYieldGen"+ endName).c_str(), nPhiBinsSecStudies, phiMinSecStudies, phiMaxSecStudies);
        momYieldGen[iSec] = new TH1F(("momYieldGen"+ endName).c_str(), ("momYieldGen"+ endName).c_str(), nMomBinsSecStudies, momMinSecStudies, momMaxSecStudies);
        
        thVSphiGen[iSec] = new TH2F(("thVSphiGen" + endName).c_str(), ("thVSphiGen" + endName).c_str(), nPhiBinsSecStudies2D_noP, phiMinSecStudies, phiMaxSecStudies, n2DtimesBinTheta*nThetaBinsSecStudies2D, thetaMinSecStudies, thetaMaxSecStudies);
        
        thVSmomGen[iSec] = new TH2F(("thVSmomGen" + endName).c_str(), ("thVSmomGen" + endName).c_str(),nMomBinsSecStudies, momMinSecStudies, momMaxSecStudies, n2DtimesBinTheta*nThetaBinsSecStudies2D, thetaMinSecStudies, thetaMaxSecStudies);
        
        phiVSmomGen[iSec] = new TH2F(("phiVSmomGen" + endName).c_str(), ("phiVSmomGen" + endName).c_str(), nMomBinsSecStudies, momMinSecStudies, momMaxSecStudies, nPhiBinsSecStudies2D_noP, phiMinSecStudies, phiMaxSecStudies);
        
        
        phiVSwGen[iSec] = new TH2F(("phiVSwGen" + endName).c_str(), ("phiVSwGen" + endName).c_str(),
                                                           nWBinsSecStudies, wMinSecStudies, wMaxSecStudies,
                                                           nPhiBinsSecStudies2D_noP, phiMinSecStudies, phiMaxSecStudies);
        phiVSq2Gen[iSec] = new TH2F(("phiVSq2Gen" + endName).c_str(), ("phiVSq2Gen" + endName).c_str(),
                                                           nQ2BinsSecStudies2D, q2MinSecStudies, q2MaxSecStudies,
                                                           nPhiBinsSecStudies2D_noP, phiMinSecStudies, phiMaxSecStudies);
        wVSq2Gen[iSec] = new TH2F(("wVSq2Gen" + endName).c_str(), ("wVSq2Gen" + endName).c_str(),
                                                           nWBinsSecStudies, wMinSecStudies, wMaxSecStudies,
                                                           nQ2BinsSecStudies2D, q2MinSecStudies, q2MaxSecStudies);
        for (int q = 0; q < static_cast<int>(qBinMax); q++){
            string endNameQ2 = "S" + to_string(iSec) + 'Q' + to_string(q);
            phiVSw_Q2binGen[iSec][q] = new TH2F(("phiVSw_Q2binGen" + endNameQ2).c_str(), ("phiVSw_Q2binGen" + endNameQ2).c_str(),
                                                       nWBinsSecStudies, wMinSecStudies, wMaxSecStudies,
                                                       nPhiBinsSecStudies2D_noP_wQ2bins, phiMinSecStudies, phiMaxSecStudies);
        }

        for (int wTMPbin = 0; wTMPbin < static_cast<int>(wBins_ST); wTMPbin++){
            string endNameW = "S" + to_string(iSec) + 'W' + to_string(wTMPbin);
            phiVSq2_WbinGen[iSec][wTMPbin] = new TH2F(("phiVSq2_WbinGen" + endNameW).c_str(), ("phiVSq2_WbinGen" + endNameW).c_str(),
                                                       nQ2BinsSecStudies2D, q2MinSecStudies, q2MaxSecStudies,
                                                       nPhiBinsSecStudies2D_noP_wQ2bins, phiMinSecStudies, phiMaxSecStudies);
        }
    }
    
    
//REC:    
    for (int iCut = 0; iCut < static_cast<int>(nCutCombination); iCut++){
        
        string endNameFull = 'C' + to_string(iCut);
        phiVSthFullRange[iCut]  = new TH2F(("phiVSthFullRange" + endNameFull).c_str(), ("phiVSthFullRange" + endNameFull).c_str(), n2DtimesBinTheta*nThetaBinsSecStudies2D, thetaMinSecStudies, thetaMaxSecStudies, n2DtimesBinPhi * nPhiBinsSecStudies, phiMinSecStudies, phiMaxSecStudies);
        
        
        for (size_t iSec = 0 ; iSec < secN; iSec++){
            string endName = 'C' + to_string(iCut) + "S" + to_string(iSec);
            
            
            phiVSw[iCut][iSec] = new TH2F(("phiVSw" + endName).c_str(), ("phiVSw" + endName).c_str(),
                                                           nWBinsSecStudies, wMinSecStudies, wMaxSecStudies,
                                                           nPhiBinsSecStudies2D_noP, phiMinSecStudies, phiMaxSecStudies);
            phiVSq2[iCut][iSec] = new TH2F(("phiVSwphiVSq2" + endName).c_str(), ("phiVSq2" + endName).c_str(),
                                                           nQ2BinsSecStudies2D, q2MinSecStudies, q2MaxSecStudies,
                                                           nPhiBinsSecStudies2D_noP, phiMinSecStudies, phiMaxSecStudies);
            wVSq2[iCut][iSec] = new TH2F(("wVSq2" + endName).c_str(), ("wVSq2" + endName).c_str(),
                                                           nWBinsSecStudies, wMinSecStudies, wMaxSecStudies,
                                                           nQ2BinsSecStudies2D, q2MinSecStudies, q2MaxSecStudies);
            
            for (int q = 0; q < static_cast<int>(qBinMax); q++){
                string endNameQ2 = 'C' + to_string(iCut) + "S" + to_string(iSec) + 'Q' + to_string(q);
                phiVSw_Q2bin[iCut][iSec][q] = new TH2F(("phiVSw_Q2bin" + endNameQ2).c_str(), ("phiVSw_Q2bin" + endNameQ2).c_str(),
                                                           nWBinsSecStudies, wMinSecStudies, wMaxSecStudies,
                                                           nPhiBinsSecStudies2D_noP_wQ2bins, phiMinSecStudies, phiMaxSecStudies);
            }
            
            for (int wTMPbin = 0; wTMPbin < static_cast<int>(wBins_ST); wTMPbin++){
                string endNameW = 'C' + to_string(iCut) + "S" + to_string(iSec) + 'W' + to_string(wTMPbin);
                phiVSq2_Wbin[iCut][iSec][wTMPbin] = new TH2F(("phiVSq2_Wbin" + endNameW).c_str(), ("phiVSq2_Wbin" + endNameW).c_str(),
                                                           nQ2BinsSecStudies2D, q2MinSecStudies, q2MaxSecStudies,
                                                           nPhiBinsSecStudies2D_noP_wQ2bins, phiMinSecStudies, phiMaxSecStudies);
            }
            
            
            thVSmom[iCut][iSec] = new TH2F(("thVSmom" + endName).c_str(), ("thVSmom" + endName).c_str(), nMomBinsSecStudies, momMinSecStudies, momMaxSecStudies, n2DtimesBinTheta*nThetaBinsSecStudies2D, thetaMinSecStudies, thetaMaxSecStudies);
            thVSphi[iCut][iSec] = new TH2F(("thVSphi" + endName).c_str(), ("thVSphi" + endName).c_str(), nPhiBinsSecStudies2D_noP, phiMinSecStudies, phiMaxSecStudies, n2DtimesBinTheta*nThetaBinsSecStudies2D, thetaMinSecStudies, thetaMaxSecStudies);
            
            phiVSmom[iCut][iSec] = new TH2F(("phiVSmom" + endName).c_str(), ("phiVSmom" + endName).c_str(), nMomBinsSecStudies, momMinSecStudies, momMaxSecStudies, nPhiBinsSecStudies2D_noP, phiMinSecStudies, phiMaxSecStudies);
            
            thYield[iCut][iSec] = new TH1F(("thYield" + endName).c_str(), ("thYield" + endName).c_str(), nThetaBinsSecStudies, thetaMinSecStudies, thetaMaxSecStudies);
            phiYield[iCut][iSec] = new TH1F(("phiYield" + endName).c_str(), ("phiYield" + endName).c_str(), nPhiBinsSecStudies, phiMinSecStudies, phiMaxSecStudies);
            momYield[iCut][iSec] = new TH1F(("momYield" + endName).c_str(), ("momYield" + endName).c_str(), nMomBinsSecStudies, momMinSecStudies, momMaxSecStudies);
            
            for (size_t iTh = 0; iTh < thNsd; iTh++){
                string endNameP = 'C' + to_string(iCut) + "S" + to_string(iSec) + "Th" + to_string(iTh);
                pYieldTh[iCut][iSec][iTh] = new TH1F(("pYieldTh" + endNameP).c_str(), ("pYieldTh" + endNameP).c_str(), 100, 2.5, 9);
            }
            for (size_t iP = 0; iP < pNsd; iP++){
                for (size_t iTh = 0; iTh < thNsd; iTh++){
                    string endNamePth = 'C' + to_string(iCut) + "S" + to_string(iSec) + 'P' + to_string(iP) + "Th" + to_string(iTh);
                    wYieldPthet[iCut][iSec][iP][iTh] = new TH1F(("wYieldPthet" + endNamePth).c_str(), ("wYieldPthet" + endNamePth).c_str(), 100, 1.15, 2.5);
                }
            }
        }
    }
// phi in theta and mom binning:
    const size_t iPbin_for_PhiYield = 12;
    const size_t iThetaBin_for_PhiYield = 20;
    const size_t iPhiBin_for_PhiYield = 30;
    
    TH1F *phiYield_P_Th[nCutCombination][secN][iPbin_for_PhiYield][iThetaBin_for_PhiYield];
    TH1F *phiYield_P_Th_Gen[secN][iPbin_for_PhiYield][iThetaBin_for_PhiYield];
    
    for (size_t iSec = 0 ; iSec < secN; iSec++){
        for (size_t iPbin = 0 ; iPbin < iPbin_for_PhiYield; iPbin++){
            for (size_t iThetaBin = 0 ; iThetaBin < iThetaBin_for_PhiYield; iThetaBin++){
                string endName = "S" + to_string(iSec) + 'P' + to_string(iPbin) + 'T' + to_string(iThetaBin);
                
                phiYield_P_Th_Gen[iSec][iPbin][iThetaBin] = new TH1F(("phiYield_P_Th_Gen" + endName).c_str(), ("phiYield_P_Th_Gen" + endName).c_str(), 
                                                   iPhiBin_for_PhiYield, phiMinSecStudies, phiMaxSecStudies);
                
                for (int iCut = 0; iCut < static_cast<int>(nCutCombination); iCut++){
                    string endNameC = 'C' + to_string(iCut) + "S" + to_string(iSec) + 'P' + to_string(iPbin) + 'T' + to_string(iThetaBin);
                    
                    phiYield_P_Th[iCut][iSec][iPbin][iThetaBin]  = new TH1F(("phiYield_P_Th" + endNameC).c_str(), 
                                                                            ("phiYield_P_Th" + endNameC).c_str(),
                                                                            iPhiBin_for_PhiYield, phiMinSecStudies, phiMaxSecStudies);                    
                }
            }
        }
    }
    //2d:
    TH2F *thVSphi_PbinGen[secN][iPbin_for_PhiYield];
    TH2F *thVSphi_Pbin[nCutCombination][secN][iPbin_for_PhiYield];
    for (size_t iSec = 0 ; iSec < secN; iSec++){
        for (size_t iPbin = 0 ; iPbin < iPbin_for_PhiYield; iPbin++){
            string endName = "S" + to_string(iSec) + 'P' + to_string(iPbin);
            thVSphi_PbinGen[iSec][iPbin] = new TH2F(("thVSphi_PbinGen" + endName).c_str(), ("thVSphi_PbinGen" + endName).c_str(), nPhiBinsSecStudies2D, phiMinSecStudies, phiMaxSecStudies, nThetaBinsSecStudies2D, thetaMinSecStudies, thetaMaxSecStudies);
                
            for (int iCut = 0; iCut < static_cast<int>(nCutCombination); iCut++){
                string endNameC = 'C' + to_string(iCut) + "S" + to_string(iSec) + 'P' + to_string(iPbin);
                thVSphi_Pbin[iCut][iSec][iPbin] = new TH2F(("thVSphi_Pbin" + endNameC).c_str(), ("thVSphi_Pbin" + endNameC).c_str(), nPhiBinsSecStudies2D, phiMinSecStudies, phiMaxSecStudies, nThetaBinsSecStudies2D, thetaMinSecStudies, thetaMaxSecStudies);
            }
        }
    }
// phi distrib in Q2 and W binning
    
    
// Nick's histograms:  
	for (int s = 0; s < 6; s++){
	  //overview
		CosthetaGen[s] = new TH1D(Form("CosthetaGenS%d", s + 1), Form("CosthetaGenS%d", s + 1), 100, 0.5, 1);
		pGen[s] = new TH1D(Form("pGenS%d", s + 1), Form("pGenS%d", s + 1), 200, 0, 10);
	  //inclusive

		for (int q = 0; q < qBinMax; q++){	
			wQ2BinSectorLogQ2StrictFidAllCuts[s][q] = new TH1D(Form("wQ2BinSectorLogQ2StrictFidAllCutsS%dQ%d", s , q), Form("wQ2BinSectorLogQ2StrictFidAllCutsS%dQ%d", s , q), nWBins, lowBorderW, highBorderW);            
			wQ2BinSectorAccLogQ2StrictFidAllCuts[s][q] = new TH1D(Form("wQ2BinSectorAccLogQ2StrictFidAllCutsS%dQ%d", s , q), Form("wQ2BinSectorAccLogQ2StrictFidAllCutsS%dQ%d", s , q), nWBins, lowBorderW, highBorderW);
            
			wQ2BinSectorGenLogQ2[s][q] = new TH1D(Form("wQ2BinSectorGenLogQ2S%dQ%d", s , q), Form("wQ2BinSectorGenLogQ2S%dQ%d", s , q), nWBins, lowBorderW, highBorderW);

// Legacy for self check:
			for (int k = 0; k < 3; k++){
				wQ2BinSectorLogQ2AllCutsSFSys[s][q][k] = new TH1D(Form("wQ2BinSectorLogQ2AllCutsSFSysS%dQ%dSys%d", s, q, k), Form("wQ2BinSectorLogQ2AllCutsSFSysS%dQ%dSys%d", s, q, k), nWBins, lowBorderW, highBorderW);
                
				wQ2BinSectorAccLogQ2AllCutsSFSys[s][q][k] = new TH1D(Form("wQ2BinSectorAccLogQ2AllCutsSFSysS%dQ%dSys%d", s, q, k), Form("wQ2BinSectorAccLogQ2AllCutsSFSysS%dQ%dSys%d", s, q, k), nWBins, lowBorderW, highBorderW);

			}
		}
	}
/////////////////////////////////////////////////////////////////////////////////    
////////////// Round 2 comments /////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
    TH1F *theta_ECAL[secN][2]; 
    TH1F *w_p_delta[secN][5][2]; 
    TH1F *w_q2_delta[secN][qBinMax][2]; 
    
    for (size_t iSec = 0 ; iSec < secN; iSec++){
		for (size_t iTh = 0 ; iTh < 2; iTh++){
			string endName = "S" + to_string(iSec) + "_Cuts_" + to_string(iTh);
			theta_ECAL[iSec][iTh] = new TH1F(("theta_ECAL" + endName).c_str(), ("theta_ECAL" + endName).c_str(), 
											200, 5., 37.);
			for (size_t iP = 0 ; iP < 5; iP++){
				string endName2 = "_P" + to_string(iP);
				w_p_delta[iSec][iP][iTh] = new TH1F(("w_p_delta" + endName + endName2).c_str(), ("w_p_delta" + endName + endName2).c_str(),40, 0.7, 1.4);
			}
			for (int q = 0; q < qBinMax; q++){	
				string endName2 = 'Q' + to_string(q);
				w_q2_delta[iSec][q][iTh] = new TH1F(("w_q2_delta" + endName + endName2).c_str(), ("w_q2_delta" + endName + endName2).c_str(),40, 0.7, 1.4);
			}							
		}
	}
 
    TH2F* WvsP_elastic[6];
    for (int iSec = 0; iSec < 6; iSec++){	
    	string endName = "S_" + to_string(iSec);
    	WvsP_elastic[iSec] = new TH2F(("WvsP_elastic_" + endName).c_str(), ("WvsP_elastic_" + endName).c_str(), 5, 2.5, 10, 45, 0.5, 1.4);
    }
    
    TH1F* WinQ2_zoom[6][qBinMax];
    
    for (int s = 0; s < 6; s++){
        for (int q = 0; q < static_cast<int>(qBinMax); q++){
            string hist_name_end = "_S" + to_string(s) + "_Q2_" + to_string(q);
            WinQ2_zoom[s][q] = new TH1F(("WinQ2_zoom" + hist_name_end).c_str(),("WinQ2_zoom" + hist_name_end).c_str(), 45, 0.5, 1.4);
        }
    }
   
/////////////////////////////////////////////////////////////////////////////////    
////////////// HTCC /////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

	const int nHTCCXBins = 250;
	const int nHTCCYBins = 250;
	int stepXHTCC = 1;
	int stepYHTCC = 1;

	TH1F *npe1D[nHTCCXBins][nHTCCYBins];

	for (int xC = 0; xC < nHTCCXBins; xC++){
		for (int yC = 0; yC < nHTCCYBins; yC++){
			npe1D[xC][yC] = new TH1F(Form("npe1DX%dY%d", xC, yC), Form("npe1DX%dY%d", xC, yC), 200, 0, 50);
		}
	}
	float HTCCEff[nHTCCXBins][nHTCCYBins];
	readHTCCEff(HTCCEff);
    
/////////////////////////////////////////////////////////////////////////////////    
////////////// triangle Cut Params //////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
    
	float triangleCutParams[6][10][2];
	readTriangleCut(triangleCutParams, isData);
/////////////////////////////////////////////////////////////////////////////////    
////////////// parameters ///////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
    
// Max event number:
	long long nEvents = 100000000000000;
// Counters:  
    long long count_GENERATED = 0;
    long long count_GENERATED2 = 0;
    long long count_RECONTRUCTED = 0;
    
	float toRD = 57.2958;
    

    
    const double pBinSizeTMP = (momMaxSecStudies - momMinSecStudies ) / iPbin_for_PhiYield;
    const double pThSizeTMP = (thetaMaxSecStudies - thetaMinSecStudies ) / iThetaBin_for_PhiYield;
    
/////////////////////////////////////////////////////////////////////////////////    
////////////// Mom Corr /////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////


// MC May 2023

auto dppC_may2023 = [&](float Px, float Py, float Pz, int sec, int ivec, int corON){
    // corON == 0 --> DOES NOT apply the momentum corrections (i.e., turns the corrections 'off')
    // corON == 1 --> Applies the momentum corrections for the experimental (real) data
    // corON == 2 --> Applies the momentum corrections for the Monte Carlo (simulated) data
    if(corON == 0){ // Momentum Corrections are OFF
        double dp = 0;
        return dp;
    }
    else{ // corON != 0 --> Applies the momentum corrections (i.e., turns the corrections 'on')
        // ivec = 0 --> Electron Corrections
        // ivec = 1 --> π+ Corrections
        // ivec = 2 --> π- Corrections
        // ivec = 3 --> Proton Corrections

        // Momentum Magnitude
        double pp = sqrt(Px*Px + Py*Py + Pz*Pz);

        // Initializing the correction factor
        double dp = 0;

        // Defining Phi Angle
        double Phi = (180/3.1415926)*atan2(Py, Px);

        // (Initial) Shift of the Phi Angle (done to realign sectors whose data is separated when plotted from ±180˚)
        if(((sec == 4 || sec == 3) && Phi < 0) || (sec > 4 && Phi < 90)){
            Phi += 360;
        }

        // Getting Local Phi Angle
        double PhiLocal = Phi - (sec - 1)*60;

        // Applying Shift Functions to Phi Angles (local shifted phi = phi)
        double phi = PhiLocal;

        // For Electron Shift
        if(ivec == 0){
            phi = PhiLocal - 30/pp;
        }

        if(corON == 2){ // Monte Carlo Simulated Corrections
            // Not Sector or Angle dependent (as of 3-21-2023)
            // Both particles were corrected at the same time using Extra_Name = "Multi_Dimension_Unfold_V1_"
            // Used ∆P = GEN - REC so the other particle does not affect how much the correction is needed
            if(ivec == 0){ // Electron Corrections
                // For MC REC (Unsmeared) ∆P(Electron) Vs Momentum Correction Equation:
                dp = (-8.2310e-04)*pp*pp + (9.0877e-03)*pp + (-1.5853e-02);
            }
            if(ivec == 1){ // Pi+ Pion Corrections
                // For MC REC (Unsmeared) ∆P(Pi+ Pion) Vs Momentum Correction Equation:
                dp = (-7.3067e-05)*pp*pp + (-8.1215e-06)*pp + (4.2144e-03);
            }
            return dp/pp;
        }
        else{
            //////////////////////////////////////////////////////////////////////////////////
            //==========//==========//     Electron Corrections     //==========//==========//
            //////////////////////////////////////////////////////////////////////////////////
            if(ivec == 0){
                if(sec == 1){
                    dp = ((-4.3303e-06)*phi*phi + (1.1006e-04)*phi + (-5.7235e-04))*pp*pp + ((3.2555e-05)*phi*phi + (-0.0014559)*phi + (0.0014878))*pp + ((-1.9577e-05)*phi*phi + (0.0017996)*phi + (0.025963));
                }
                if(sec == 2){
                    dp = ((-9.8045e-07)*phi*phi + (6.7395e-05)*phi + (-4.6757e-05))*pp*pp + ((-1.4958e-05)*phi*phi + (-0.0011191)*phi + (-0.0025143))*pp + ((1.2699e-04)*phi*phi + (0.0033121)*phi + (0.020819));
                }
                if(sec == 3){
                    dp = ((-5.9459e-07)*phi*phi + (-2.8289e-05)*phi + (-4.3541e-04))*pp*pp + ((-1.5025e-05)*phi*phi + (5.7730e-04)*phi + (-0.0077582))*pp + ((7.3348e-05)*phi*phi + (-0.001102)*phi + (0.057052));
                }
                if(sec == 4){
                    dp = ((-2.2714e-06)*phi*phi + (-3.0360e-05)*phi + (-8.9322e-04))*pp*pp + ((2.9737e-05)*phi*phi + (5.1142e-04)*phi + (0.0045641))*pp + ((-1.0582e-04)*phi*phi + (-5.6852e-04)*phi + (0.027506));
                }
                if(sec == 5){
                    dp = ((-1.1490e-06)*phi*phi + (-6.2147e-06)*phi + (-4.7235e-04))*pp*pp + ((3.7039e-06)*phi*phi + (-1.5943e-04)*phi + (-8.5238e-04))*pp + ((4.4069e-05)*phi*phi + (0.0014152)*phi + (0.031933));
                }
                if(sec == 6){
                    dp = ((1.1076e-06)*phi*phi + (4.0156e-05)*phi + (-1.6341e-04))*pp*pp + ((-2.8613e-05)*phi*phi + (-5.1861e-04)*phi + (-0.0056437))*pp + ((1.2419e-04)*phi*phi + (4.9084e-04)*phi + (0.049976));
                }
            }
            //////////////////////////////////////////////////////////////////////////////////
            //==========//==========//  Electron Corrections (End)  //==========//==========//
            //////////////////////////////////////////////////////////////////////////////////

            return dp/pp;
        }
    }
};

    
    /*
  double xx[] =
{
    0.0263375, 0.0158871, 0.0130852, -0.00366006, 0.00694866, 0.0197195,
    0.00767067, 0.00480921, -0.0175756, 0.0252757, 0.0156601, 0.00984872,
    0.00244435, 0.00681414, 0.0294068, 0.0059881, 0.00286992, 0.0179319,
    0.0171495, 0.00359637, -0.0046115, 0.00314739, 0.0136338, 0.0768753,
    0.00675454, -0.0118234, -0.0288654, 0.0189465, 0.0131816, 0.0262004,
    0.00375165, 0.00907457, 0.0486894, 0.00806305, 0.0006999, 0.00527513,
    0.0116485, 0.0105681, 0.0149848, 0.000318094, -0.00480124, 0.0395545,
    0.00824216, -0.00070659, -0.0057075, 0.0213057, 0.0112999, 0.0100216,
    0.000653685, 0.0093174, 0.0822385, 0.00808384, 0.000898799, -0.0172692,
};

double pars[6][3][3];
int ipar=0;

for(int isec=0;isec<6;isec++)
for(int ivec=0;ivec<3;ivec++)
{
    double dp1=xx[ipar++], dp5=xx[ipar++], dp9=xx[ipar++];

    pars[isec][ivec][0] = (dp1 - 2*dp5 + dp9)/32.;
    pars[isec][ivec][1] = (-7*dp1)/16. + (5*dp5)/8. - (3*dp9)/16.;
    pars[isec][ivec][2] = (45*dp1)/32. - (9*dp5)/16. + (5*dp9)/32.;
}

auto dpp = [&](float px, float py, float pz, int sec, int ivec)
{
    double pp = sqrt(px*px + py*py + pz*pz);

    double a=pars[sec-1][ivec][0],
           b=pars[sec-1][ivec][1],
           c=pars[sec-1][ivec][2];

    double dp = a*pp*pp + b*pp + c; //pol2 corr func
   
    //electron pol1 corr func for each sec and each phi bins
    if(ivec == 0)
    {
        if(sec == 1)
        {
            dp = 0.45*b*(pp-9)+0.1*c;
           
            //ep 3 phi bins
            //dp = -0.01*b*(pp-9)+1.35*c; //phi<-5
            //dp = 0.6*b*(pp-9)-0.3*c; //-5<phi<5
            //dp = 1.7*b*(pp-9)-1.5*c; //phi>5
        }
        if(sec == 2)
        {
            dp = -0.15*b*(pp-8.0)-0.3*c;

            //ep 3 phi bins
            //dp = -0.7*b*(pp-8.0)+0.4*c; //phi<-5
            //dp = -0.05*b*(pp-8.0)-0.4*c; //-5<phi<5
            //dp = 0.01*b*(pp-8.0)-1.5*c; //phi>5
        }
        if(sec == 3)
        {
            dp = 3.*b*(pp-5.4)-0.5*c;
         
            //ep 3 phi bins
            //dp = 0.04*b*(pp-5.4)-3.5*c; //phi<-5
            //dp = 0.06*b*(pp-5.4)-3.*c; //-5<phi<5
            //dp = 1.1*b*(pp-5.4)-0.7*c; //phi>5
        }
        if(sec == 4)
        {
            dp = 0.25*b*(pp-9.25)-0.3*c;
       
            //ep 3 phi bins
            //dp = 0.25*b*(pp-9.25)-0.7*c; //phi<-5
            //dp = 0.25*b*(pp-9.25)+0.05*c; //-5<phi<5
            //dp = 0.1*b*(pp-9.25)+1.1*c; //phi>5
        }
        if(sec == 5)
        {
            dp = 2.2*b*(pp-7.5)-0.5*c;
   
            //ep 3 phi bins
            //dp = 2.2*b*(pp-7.5)+0.5*c; //phi<-5
            //dp = 2.2*b*(pp-7.5)-0.1*c; //-5<phi<5
            //dp = 2.2*b*(pp-7.5)-0.6*c; //phi>5
        }
        if(sec == 6)
        {
            dp = 0.5*b*(pp-7)-0.6*c;
           
            //ep 3 phi bins
            //dp = 1.263*b*(pp-7)+0.5*c; //phi<-5
            //dp = 1.*b*(pp-7)-0.5*c; //-5<phi<5
            //dp = 0.5*b*(pp-7)-1.45*c; //phi>5
   
        }
    }
    return dp/pp;
};  
*/    
    
    
    


// MC from may 16 2022 with phi dependences based one one pion channel
 auto dppC = [&](float Px, float Py, float Pz, int sec, int ivec){

        // ivec = 0 --> Electron Corrections
        // ivec = 1 --> Pi+ Corrections (removed)
        // ivec = 2 --> Pi- Corrections (removed)
        // ivec = 3 --> Proton Corrections (removed)

        // Momentum Magnitude
        double pp = sqrt(Px*Px + Py*Py + Pz*Pz);

        // Initializing the correction factor
        double dp = 0;

        // Defining Phi Angle
        double Phi = (180/3.1415926)*atan2(Py, Px);

        // (Initial) Shift of the Phi Angle (done to realign sectors whose data is separated when plotted from ±180˚)
        if(((sec == 4 || sec == 3) && Phi < 0) || (sec > 4 && Phi < 90)){
            Phi += 360;
        }

        // Getting Local Phi Angle
        double PhiLocal = Phi - (sec - 1)*60;

        // Applying Shift Functions to Phi Angles (local shifted phi = phi)
        double phi = PhiLocal;

        // For Electron Shift
        if(ivec == 0){
            phi = PhiLocal - 30/pp;
        }

        // For Pi+ Pion/Proton Shift
        if(ivec == 1 || ivec == 3){
            phi = PhiLocal + (32/(pp-0.05));
        }

        // For Pi- Pion Shift
        if(ivec == 2){
            phi = PhiLocal - (32/(pp-0.05));
        }




        //==========//  PARTICLE = ELECTRON  //==========//
        
        
        if(ivec == 0){

            if(sec == 1){

                dp = ((1.57e-06)*phi*phi + (5.021e-05)*phi + (-1.74089e-03))*pp*pp + ((-2.192e-05)*phi*phi + (-1.12528e-03)*phi + (0.0146476))*pp + ((8.504e-05)*phi*phi + (2.08012e-03)*phi + (-0.0122501));

            }

            if(sec == 2){

                dp = ((-3.98e-06)*phi*phi + (1.66e-05)*phi + (-1.55918e-03))*pp*pp + ((2.136e-05)*phi*phi + (-5.7373e-04)*phi + (0.0143591))*pp + ((2.4e-06)*phi*phi + (1.6656e-03)*phi + (-0.0218711));

            }

            if(sec == 3){

                dp = ((5.57e-06)*phi*phi + (2.3e-07)*phi + (-2.26999e-03))*pp*pp + ((-7.761e-05)*phi*phi + (4.1437e-04)*phi + (0.0152985))*pp + ((2.2542e-04)*phi*phi + (-9.442e-04)*phi + (-0.0231432));

            }

            if(sec == 4){

                dp = ((3.48e-06)*phi*phi + (2.166e-05)*phi + (-2.29e-04))*pp*pp + ((-2.758e-05)*phi*phi + (7.226e-05)*phi + (-3.38e-03))*pp + ((3.166e-05)*phi*phi + (6.93e-05)*phi + (0.04767));

            }

            if(sec == 5){

                dp = ((1.19e-06)*phi*phi + (-2.286e-05)*phi + (-1.6332e-04))*pp*pp + ((-1.05e-06)*phi*phi + (7.04e-05)*phi + (-5.0754e-03))*pp + ((-7.22e-06)*phi*phi + (4.1748e-04)*phi + (0.04441));

            }

            if(sec == 6){

                dp = ((-5.97e-06)*phi*phi + (-3.689e-05)*phi + (5.782e-05))*pp*pp + ((6.573e-05)*phi*phi + (2.1376e-04)*phi + (-9.54576e-03))*pp + ((-1.7732e-04)*phi*phi + (-8.62e-04)*phi + (0.0618975));

            }

        }
        
        
        //==========//  PARTICLE = ELECTRON (END)  //==========//
        return dp/pp;
    };

/////////////////////////////////////////////////////////////////////////////////    
////////////// smearing /////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
    
    fastMC();
    
    
//===========================================================================//
//================//     (Modified) Smearing Function      //================//
//===========================================================================//
auto smear_func = [&](TLorentzVector V4, int ivec){
    // True generated values (i.e., values of the unsmeared TLorentzVector)
    double inM = V4.M();
    double smeared_P  = V4.P();
    double smeared_Th = V4.Theta();
    double smeared_Phi = V4.Phi();
    TLorentzVector V4_new(V4.X(), V4.Y(), V4.Z(), V4.E());
    // Calculate resolutions
    double smeared_ThD = TMath::RadToDeg()*smeared_Th;
    double momS1 = 0.0184291 - 0.0110083*smeared_ThD + 0.00227667*smeared_ThD*smeared_ThD - 0.000140152*smeared_ThD*smeared_ThD*smeared_ThD + (3.07424e-06)*smeared_ThD*smeared_ThD*smeared_ThD*smeared_ThD;
    double momS2 = 0.02*smeared_ThD;
    double momR  = 0.01 * TMath::Sqrt( TMath::Power(momS1*smeared_P,2) + TMath::Power(momS2,2));
    momR *= 2.0;
    if(ivec == 0){
        // From ∆P(Electron) Sigma Vs Momentum distributions:
        momR *= (2.0604e-02)*(V4.P())*(V4.P()) + (-1.1212e-01)*(V4.P()) + (7.1348e-01);
        momR *= (-1.1295e-02)*(V4.P())*(V4.P()) + (2.3416e-01)*(V4.P()) + (-1.0974e-01);
        // From ∆P(Electron) Sigma Vs Theta distributions:
        momR *= (-2.7657e-03)*(V4.Theta()*TMath::RadToDeg())*(V4.Theta()*TMath::RadToDeg()) +  (8.1714e-02)*(V4.Theta()*TMath::RadToDeg()) + (4.0196e-01);
        momR *= (-7.3974e-04)*(V4.Theta()*TMath::RadToDeg())*(V4.Theta()*TMath::RadToDeg()) +  (1.0908e-02)*(V4.Theta()*TMath::RadToDeg()) + (9.9876e-01);
    }
    if(ivec == 1){
        // From ∆P(Pi+ Pion) Sigma Vs Momentum distributions:
        momR *= (-5.2125e-02)*(V4.P())*(V4.P()) + (2.7110e-01)*(V4.P()) + (4.8534e-01);
        momR *= (-3.4607e-02)*(V4.P())*(V4.P()) + (1.5836e-01)*(V4.P()) + (6.8845e-01);
        // From ∆P(Pi+ Pion) Sigma Vs Theta distributions:
        momR *= (-5.7711e-04)*(V4.Theta()*TMath::RadToDeg())*(V4.Theta()*TMath::RadToDeg()) +  (2.4354e-02)*(V4.Theta()*TMath::RadToDeg()) + (6.6472e-01);
        momR *= (-1.3210e-03)*(V4.Theta()*TMath::RadToDeg())*(V4.Theta()*TMath::RadToDeg()) +  (3.5065e-02)*(V4.Theta()*TMath::RadToDeg()) + (8.3333e-01);
        // From ∆P(Pi+ Pion) Sigma Vs Theta distributions:
        momR *=  (6.9462e-03)*(V4.Theta()*TMath::RadToDeg())*(V4.Theta()*TMath::RadToDeg()) + (-4.8277e-01)*(V4.Theta()*TMath::RadToDeg()) + (9.8916e+00);
    }
    double theS1 = 0.004*smeared_ThD + 0.1;
    double theS2 = 0;
    double theR  = TMath::Sqrt(TMath::Power(theS1*TMath::Sqrt(smeared_P*smeared_P + 0.13957*0.13957)/(smeared_P*smeared_P),2) + TMath::Power(theS2,2) );
    theR *= 2.5;
    double phiS1 = 0.85 - 0.015*smeared_ThD;
    double phiS2 = 0.17 - 0.003*smeared_ThD;
    double phiR  = TMath::Sqrt(TMath::Power(phiS1*TMath::Sqrt(smeared_P*smeared_P + 0.13957*0.13957)/(smeared_P*smeared_P),2) + TMath::Power(phiS2,2) );
    phiR *= 3.5;
    // overwrite EB (i.e., applying the smear)
    smeared_Phi += TMath::DegToRad() * phiR * gRandom->Gaus(0,1);
    smeared_Th += TMath::DegToRad() * theR * gRandom->Gaus(0,1);
    smeared_P  += momR  * gRandom->Gaus(0,1) *  V4.P();
    V4_new.SetE( TMath::Sqrt( smeared_P*smeared_P + inM*inM )  );
    V4_new.SetRho( smeared_P );
    V4_new.SetTheta( smeared_Th );
    V4_new.SetPhi( smeared_Phi );
    return V4_new;
};
    
//New Smearing:
auto smear_oneFactor = [&](const TLorentzVector V4rec){
	// smear factor:
	
	double elTh_indeg = V4rec.Theta() * 57.2958;
	double sigmaRes = 0;
 	const double fp_4[5] = {     2.3762794461143009e+000,
    -2.8635604070815707e-001,
     3.8197546356590874e-002,
    -1.8427728047602948e-003,
     2.7446233737351622e-005};
     

 

	const double factor = fp_4[0] + fp_4[1]*elTh_indeg + fp_4[2]*elTh_indeg*elTh_indeg + fp_4[3]*elTh_indeg*elTh_indeg*elTh_indeg + fp_4[4]*elTh_indeg*elTh_indeg*elTh_indeg*elTh_indeg;;
	
	
	sigmaRes = 0.008007538378779988 -0.0010764451335380787* elTh_indeg + 6.943332514486485e-05* elTh_indeg * elTh_indeg  -1.1181109847470665e-06* elTh_indeg * elTh_indeg * elTh_indeg ;

    // True generated values (i.e., values of the unsmeared TLorentzVector)
    double inM = V4rec.M();
    double gausValue = myMC->Gaus(0,1);
    double smeared_P  = V4rec.P() +  V4rec.P() * sigmaRes * factor * gausValue;
    
//cout<<inM<<endl;
//    cout<<V4rec.P() * sigmaRes * factor * gausValue << " gaus:"<< gausValue <<endl;   
    double smeared_Th = V4rec.Theta();
    double smeared_Phi = V4rec.Phi();
    
    TLorentzVector V4_new(V4rec.X(), V4rec.Y(), V4rec.Z(), V4rec.E());
    
    V4_new.SetE( TMath::Sqrt( smeared_P*smeared_P + inM*inM )  );
    V4_new.SetRho( smeared_P );
    V4_new.SetTheta( smeared_Th );
    V4_new.SetPhi( smeared_Phi );
    
    return V4_new;
};
/*
auto smear_oneFactor = [&](const TLorentzVector V4rec, const double factor){
	// smear factor:

	
	double elTh_indeg = V4rec.Theta() * 57.2958;
	double sigmaRes = 0;

	if (elTh_indeg < 27){
		sigmaRes = 0.003871160894735281 -0.00019070051211886123* elTh_indeg  + 1.3021243288068022e-05* elTh_indeg* elTh_indeg;
		
	}else{
		sigmaRes =-0.0023611338970121534 + 0.0007469420215626479* elTh_indeg -1.322734449485198e-05* elTh_indeg* elTh_indeg;
	}		
    // True generated values (i.e., values of the unsmeared TLorentzVector)
    double inM = V4rec.M();
    double gausValue = myMC->Gaus(0,1);
    double smeared_P  = V4rec.P() +  V4rec.P() * sigmaRes * factor * gausValue;
    
//cout<<inM<<endl;
//    cout<<V4rec.P() * sigmaRes * factor * gausValue << " gaus:"<< gausValue <<endl;   
    double smeared_Th = V4rec.Theta();
    double smeared_Phi = V4rec.Phi();
    
    TLorentzVector V4_new(V4rec.X(), V4rec.Y(), V4rec.Z(), V4rec.E());
    
    V4_new.SetE( TMath::Sqrt( smeared_P*smeared_P + inM*inM )  );
    V4_new.SetRho( smeared_P );
    V4_new.SetTheta( smeared_Th );
    V4_new.SetPhi( smeared_Phi );
    
    return V4_new;
};
*/    
///////////////////////mom corr MC///////////////////////////////
auto dpp_Sim = [&](float mom){

	if (mom < 5.2 || mom > 10) return 0.;
	
	double dpp = 0.;
	if (mom < 7.75){
		dpp = 0.021292659130789704 - 0.0059428307653227735 * mom +  0.0004391295650694022 * mom * mom;
	}else{
		dpp = -0.015507777972874268 + 0.004419179617144075 * mom  -0.0002859299237885943 * mom * mom;
	}
	return dpp;
};
    
/////////////////////////////////////////////////////////////////////////////////    
////////////// TTree definition /////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
    
    size_t nfile = 1;
    for (string file : file_list){
        
    TFile *f = new TFile(file.c_str());
    TTree *anaTree=(TTree *) f->Get("out_tree");
        
	float part_px[100], part_py[100], part_pz[100], part_E[100];
	float part_vx[100], part_vy[100], part_vz[100];

	float part_px_gen[100], part_py_gen[100], part_pz_gen[100], part_E_gen[100];

	float ele_vx[100], ele_vy[100], ele_vz[100];
	float htcc_X[100], htcc_Y[100], htcc_NPE[100];
	float pcal_U[100], pcal_V[100], pcal_W[100];
	float pcal_LV[100], pcal_LW[100], pcal_LU[100];
	float pcal_X[100], pcal_Y[100], pcal_Z[100];
	float pcal_HX[100], pcal_HY[100], pcal_HZ[100];

	float ecin_U[100], ecin_V[100], ecin_W[100];
	float ecin_X[100], ecin_Y[100], ecin_Z[100];
	float ecin_HX[100], ecin_HY[100], ecin_HZ[100];

	float ecout_U[100], ecout_V[100], ecout_W[100];
	float ecout_X[100], ecout_Y[100], ecout_Z[100];
	float ecout_HX[100], ecout_HY[100], ecout_HZ[100];
        
	float ftof_HX[100],ftof_HY[100];

	float dc_XR1[100], dc_XR2[100], dc_XR3[100];
	float dc_YR1[100], dc_YR2[100], dc_YR3[100];
	float dc_ZR1[100], dc_ZR2[100], dc_ZR3[100];

	float pcal_E[100], ecin_E[100], ecout_E[100], ele_chi2[100], ele_ndf[100], ele_tr_chi2[100];
	int sectorE[100];
	int sector_Sci[100], layer_Sci[100], comp_Sci[100];

	vector<double>  *vpart_px      = 0;
	vector<double>  *vpart_py      = 0;
	vector<double>  *vpart_pz      = 0;
	
	vector<double>  *htccX      = 0;
	vector<double>  *htccY      = 0;
	vector<double>  *htccNPE      = 0;

	vector<double>  *pcalU      = 0;
	vector<double>  *pcalV      = 0;
	vector<double>  *pcalW      = 0;

	vector<double>  *pcalX      = 0;
	vector<double>  *pcalY      = 0;
	vector<double>  *pcalZ      = 0;
	vector<double>  *pcalHX      = 0;
	vector<double>  *pcalHY      = 0;
	vector<double>  *pcalHZ      = 0;


	vector<double>  *ecinU      = 0;
	vector<double>  *ecinV      = 0;
	vector<double>  *ecinW      = 0;
	vector<double>  *ecinX      = 0;
	vector<double>  *ecinY      = 0;
	vector<double>  *ecinZ      = 0;
	vector<double>  *ecinHX      = 0;
	vector<double>  *ecinHY      = 0;
	vector<double>  *ecinHZ      = 0;

	vector<double>  *ecoutU      = 0;
	vector<double>  *ecoutV      = 0;
	vector<double>  *ecoutW      = 0;
	vector<double>  *ecoutX      = 0;
	vector<double>  *ecoutY      = 0;
	vector<double>  *ecoutZ      = 0;
	vector<double>  *ecoutHX      = 0;
	vector<double>  *ecoutHY      = 0;
	vector<double>  *ecoutHZ      = 0;
        
	vector<double>  *ftofHX      = 0;
	vector<double>  *ftofHY      = 0;


	vector<double> *dcXR1       = 0;
	vector<double> *dcYR1       = 0;
	vector<double> *dcXR2       = 0;
	vector<double> *dcYR2       = 0;
	vector<double> *dcXR3       = 0;
	vector<double> *dcYR3       = 0;

	vector<double> *dcZR1       = 0;
	vector<double> *dcZR2       = 0;
	vector<double> *dcZR3       = 0;

	vector<double>  *gen_px      = 0;
	vector<double>  *gen_py      = 0;
	vector<double>  *gen_pz      = 0;

	vector<double>  *vpart_vx      = 0;
	vector<double>  *vpart_vy      = 0;
	vector<double>  *vpart_vz      = 0;
	vector<double>  *vpart_E       = 0;
	vector<int>    *v_sectorE     = 0;
	
	//Ftof:
	vector<int>    *sectorSci     = 0;
	vector<int>    *layerSci     = 0;
	vector<int>    *compSci     = 0;

	vector<float>  *pcalE         = 0;
	vector<float>  *ecinE         = 0;
	vector<float>  *ecoutE        = 0;
	vector<float>  *eleChi2      = 0;
	vector<float>  *eleNdf      = 0;
	vector<float>  *eletrchi2      = 0;
        
	TLorentzVector p4_ele[100], p4_ele_gen[100], p4_ele_final, p4_ele_corrected[100];

	anaTree->SetCacheSize(0);
        
	anaTree->SetBranchAddress("p4_ele_px", &vpart_px);
	anaTree->SetBranchAddress("p4_ele_py", &vpart_py);
	anaTree->SetBranchAddress("p4_ele_pz", &vpart_pz);
	anaTree->SetBranchAddress("p4_ele_E", &vpart_E);

	anaTree->SetBranchAddress("htccX", &htccX);
	anaTree->SetBranchAddress("htccY", &htccY);
	anaTree->SetBranchAddress("htccNPE", &htccNPE);

	anaTree->SetBranchAddress("pcalHX", &pcalHX);
	anaTree->SetBranchAddress("pcalHY", &pcalHY);
	anaTree->SetBranchAddress("pcalHZ", &pcalHZ);
	
	anaTree->SetBranchAddress("pcalLv", &pcalV);
	anaTree->SetBranchAddress("pcalLw", &pcalW);
	anaTree->SetBranchAddress("pcalLu", &pcalU);
        
    // added detectors:
        
	anaTree->SetBranchAddress("ecinHX", &ecinHX);
	anaTree->SetBranchAddress("ecinHY", &ecinHY);
	anaTree->SetBranchAddress("ecinHZ", &ecinHZ);

	anaTree->SetBranchAddress("ecoutHX", &ecoutHX);
	anaTree->SetBranchAddress("ecoutHY", &ecoutHY);
	anaTree->SetBranchAddress("ecoutHZ", &ecoutHZ);
        
	anaTree->SetBranchAddress("ftofHX", &ftofHX);
	anaTree->SetBranchAddress("ftofHY", &ftofHY);
    // end of added detectors        

	anaTree->SetBranchAddress("pcalE", &pcalE);
	anaTree->SetBranchAddress("ecinE", &ecinE);
	anaTree->SetBranchAddress("ecoutE", &ecoutE);
	anaTree->SetBranchAddress("ele_chi2", &eleChi2);
	anaTree->SetBranchAddress("ele_tr_chi2", &eletrchi2);
	anaTree->SetBranchAddress("ele_ndf", &eleNdf);
        
	anaTree->SetBranchAddress("p4_ele_vx", &vpart_vx);
	anaTree->SetBranchAddress("p4_ele_vy", &vpart_vy);
	anaTree->SetBranchAddress("p4_ele_vz", &vpart_vz);

	anaTree->SetBranchAddress("dcXR1", &dcXR1);
	anaTree->SetBranchAddress("dcYR1", &dcYR1);
	anaTree->SetBranchAddress("dcXR2", &dcXR2);
	anaTree->SetBranchAddress("dcYR2", &dcYR2);
	anaTree->SetBranchAddress("dcXR3", &dcXR3);
	anaTree->SetBranchAddress("dcYR3", &dcYR3);
        
    // MISSPRINT in ana12Short dcZR1 -> dcYZ1. I fixed it in the skimming for data, but Valerii sim file still has this misprint
     if (isData==1){
    anaTree->SetBranchAddress("dcZR1", &dcZR1);  // my fix.
     } else {
	anaTree->SetBranchAddress("dcYZ1", &dcZR1);
     }
	anaTree->SetBranchAddress("dcZR2", &dcZR2);
	anaTree->SetBranchAddress("dcZR3", &dcZR3);

	anaTree->SetBranchAddress("gen_px", &gen_px);
	anaTree->SetBranchAddress("gen_py", &gen_py);
	anaTree->SetBranchAddress("gen_pz", &gen_pz);

	anaTree->SetBranchAddress("sectorE", &v_sectorE);
	
	if (isData == 1){
		anaTree->SetBranchAddress("sectorSci", &sectorSci);
		anaTree->SetBranchAddress("layerSci", &layerSci);
		anaTree->SetBranchAddress("compSci", &compSci);
    }
    
    Long64_t n_entriesTree = anaTree->GetEntries();
    cout<<"ENTRIES:"<<n_entriesTree<<endl;
/////////////////////////////////////////////////////////////////////////////////    
////////////// Various cut combinations///////////////////
/////////////////////////////////////////////////////////////////////////////////
                            
    // one cut missing:
    const set<cutType> cuts_NoDCcut = AddAlllCutsExcept(cutType::dcCut);
    const set<cutType> cuts_NoBADelem = AddAlllCutsExcept(cutType::badElemCut);
    const set<cutType> cuts_NoThVSphi = AddAlllCutsExcept(cutType::thetaVSphiDC);
    const set<cutType> cuts_NoThMax = AddAlllCutsExcept(cutType::thetaMax);
    const set<cutType> cuts_NoPCALcut = AddAlllCutsExcept(cutType::pcalCut);
    const set<cutType> cuts_NoPhiMax = AddAlllCutsExcept(cutType::phiMax);
    const set<cutType> cuts_NoPhiMin = AddAlllCutsExcept(cutType::phiMin);
    const set<cutType> cuts_NoPhiSpike = AddAlllCutsExcept(cutType::phiSpike);
                              
                            
    //two cut missing
    // no theta max, no bad elem
    set<cutType> cuts_NoThetaMaxNoBadElem = AddAlllCutsExcept(cutType::allCuts);
    cuts_NoThetaMaxNoBadElem.erase(cutType::thetaMax);
    cuts_NoThetaMaxNoBadElem.erase(cutType::badElemCut);
        
                            
    // no theta max, phi max
    set<cutType> cuts_NoThetaMaxNoPhiMax = AddAlllCutsExcept(cutType::allCuts);
    cuts_NoThetaMaxNoPhiMax.erase(cutType::thetaMax);
    cuts_NoThetaMaxNoPhiMax.erase(cutType::phiMax);
                            
    // no theta max, phi max
    set<cutType> cuts_NoThetaMaxNoPhiMin = AddAlllCutsExcept(cutType::allCuts);
    cuts_NoThetaMaxNoPhiMin.erase(cutType::thetaMax);
    cuts_NoThetaMaxNoPhiMin.erase(cutType::phiMin);
                            
    // no theta max, thetaVSphi
    set<cutType> cuts_NoThetaMaxNoThetaVSphi= AddAlllCutsExcept(cutType::allCuts);
    cuts_NoThetaMaxNoThetaVSphi.erase(cutType::thetaMax);
    cuts_NoThetaMaxNoThetaVSphi.erase(cutType::thetaVSphiDC);
        
    // no phi spike, thetaVSphi
    set<cutType> cuts_NoSpikeNoThetaVSphi = AddAlllCutsExcept(cutType::allCuts);
    cuts_NoSpikeNoThetaVSphi.erase(cutType::phiSpike);
    cuts_NoSpikeNoThetaVSphi.erase(cutType::thetaVSphiDC);
                            
    // three cut missing:
    // no DC, PCAL and thetaVSphiDC cuts
    set<cutType> cuts_NoFidAtAll = AddAlllCutsExcept(cutType::allCuts);
    cuts_NoFidAtAll.erase(cutType::dcCut);
    cuts_NoFidAtAll.erase(cutType::pcalCut);
    cuts_NoFidAtAll.erase(cutType::thetaVSphiDC);
    cuts_NoFidAtAll.erase(cutType::phiSpike);
                            
    // no theta max, phi max, thetaVSphi
    set<cutType> cuts_NoThetaMaxNoThetaVSphiNoPhiMin= AddAlllCutsExcept(cutType::allCuts);
    cuts_NoThetaMaxNoThetaVSphiNoPhiMin.erase(cutType::thetaMax);
    cuts_NoThetaMaxNoThetaVSphiNoPhiMin.erase(cutType::phiMin);
    cuts_NoThetaMaxNoThetaVSphiNoPhiMin.erase(cutType::thetaVSphiDC);
                            
    // four cut missing:
    // no theta max, phi max, no phi min, thetaVSphi
    set<cutType> cuts_NoAngleCuts = AddAlllCutsExcept(cutType::allCuts);
    cuts_NoAngleCuts.erase(cutType::thetaMax);
    cuts_NoAngleCuts.erase(cutType::phiMax);
    cuts_NoAngleCuts.erase(cutType::phiMin);
    cuts_NoAngleCuts.erase(cutType::thetaVSphiDC);
        
    // five cut missing:
    // no theta max, phi max, no phi min, thetaVSphi
    set<cutType> cuts_NoAngleCutsNoBadElem = AddAlllCutsExcept(cutType::allCuts);
    cuts_NoAngleCutsNoBadElem.erase(cutType::thetaMax);
    cuts_NoAngleCutsNoBadElem.erase(cutType::phiMax);
    cuts_NoAngleCutsNoBadElem.erase(cutType::phiMin);
    cuts_NoAngleCutsNoBadElem.erase(cutType::thetaVSphiDC);
    cuts_NoAngleCutsNoBadElem.erase(cutType::badElemCut);
    cuts_NoAngleCutsNoBadElem.erase(cutType::phiSpike);
        
    // eight cut missing:
    // no theta max, phi max, no phi min, thetaVSphi
    set<cutType> cuts_NoAngleNoFid = AddAlllCutsExcept(cutType::allCuts);
    cuts_NoAngleNoFid.erase(cutType::thetaMax);
    cuts_NoAngleNoFid.erase(cutType::phiMax);
    cuts_NoAngleNoFid.erase(cutType::phiMin);
    cuts_NoAngleNoFid.erase(cutType::thetaVSphiDC);
    cuts_NoAngleNoFid.erase(cutType::badElemCut);
    cuts_NoAngleNoFid.erase(cutType::dcCut);
    cuts_NoAngleNoFid.erase(cutType::pcalCut);
    cuts_NoAngleNoFid.erase(cutType::phiSpike);
    
    set<cutType> cuts_NoAngleNoFidNoTrig = AddAlllCutsExcept(cutType::allCuts);
    cuts_NoAngleNoFidNoTrig.erase(cutType::thetaMax);
    cuts_NoAngleNoFidNoTrig.erase(cutType::phiMax);
    cuts_NoAngleNoFidNoTrig.erase(cutType::phiMin);
    cuts_NoAngleNoFidNoTrig.erase(cutType::thetaVSphiDC);
    cuts_NoAngleNoFidNoTrig.erase(cutType::badElemCut);
    cuts_NoAngleNoFidNoTrig.erase(cutType::dcCut);
    cuts_NoAngleNoFidNoTrig.erase(cutType::pcalCut);
    cuts_NoAngleNoFidNoTrig.erase(cutType::phiSpike);
    cuts_NoAngleNoFidNoTrig.erase(cutType::triagCut);
        
    // 
    set<cutType> cuts_NoCuts = AddAlllCutsExcept(cutType::allCuts);
    cuts_NoCuts.erase(cutType::thetaMax);
    cuts_NoCuts.erase(cutType::phiMax);
    cuts_NoCuts.erase(cutType::phiMin);
    cuts_NoCuts.erase(cutType::thetaVSphiDC);
    cuts_NoCuts.erase(cutType::badElemCut);
    cuts_NoCuts.erase(cutType::dcCut);
    cuts_NoCuts.erase(cutType::pcalCut);
    cuts_NoCuts.erase(cutType::phiSpike);
    cuts_NoCuts.erase(cutType::triagCut);
    cuts_NoCuts.erase(cutType::sfCut);
        
        
                            
    //all:
    const set<cutType> cuts_All = AddAlllCutsExcept(cutType::allCuts);    
        
                            
    //vzCut = 0
    const set<cutType> cuts_NoVz = AddAlllCutsExcept(cutType::vzCut);
    
    //dcCut = 1 
    //Defined before, const set<cutType> cuts_NoDCcut
                            
    //sfCut = 2
    const set<cutType> cuts_NoSFcut = AddAlllCutsExcept(cutType::sfCut);
                            
    //momCut = 3
    const set<cutType> cuts_NoMomCut = AddAlllCutsExcept(cutType::momCut);
    
    //pcalCut = 4
    //const set<cutType> cuts_NoPCALcut = AddAlllCutsExcept(cutType::pcalCut);
                            
    //triagCut = 5
    const set<cutType> cuts_NoTriag = AddAlllCutsExcept(cutType::triagCut);
                            
    //badElemCut = 6
    //Defined before const set<cutType> cuts_NoBADelem 
     
    //pcalDeposCut = 7
    const set<cutType> cuts_NoDeposE = AddAlllCutsExcept(cutType::pcalDeposCut);
                            
    //thetaMax = 8
    //Defined before, const set<cutType> cuts_NoThMax
    
    //thetaMin = 9
    const set<cutType> cuts_NoThetaMin = AddAlllCutsExcept(cutType::thetaMin);
                            
    //thetaVSphiDC = 10
    //Defined before, const set<cutType> cuts_NoT
        
        
    // No HTCC weight goes separately
    //const vector<pair<set<cutType>, int>> vCutComb = {{cuts_NoCuts,0}, {cuts_NoAngleNoFid, 1}, {cuts_NoPhiSpike,2},
    //                                                  {cuts_NoSpikeNoThetaVSphi,3}, {cuts_All, 4}};
                                                      
    const vector<pair<set<cutType>, int>> vCutComb = {{cuts_All, 0}};
        
        
    // 
        
/////////////////////////////////////////////////////////////////////////////////    
////////////// TTree reading ////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////        
	for(Long64_t k=0; k<n_entriesTree;k++){
		if (anaTree->GetEntry(k) == 0) continue;
        float WGen = -1, Q2Gen = -1, PGen_unf = -1;
        int Q2_bin_log = -1;
        
/////////////////////////////////////////////////////////////////////////////////    
////////////// Simulation GENERATED /////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////    
        
        bool isMissEvent = 1;
        bool isMissEvent_roo = 1;
        
        bool genEventExist = 0;
        int SectorGen = -1;
        
        int genBINmiss = 0;
        
		if (isData == 0){  
			int NPARTGen = gen_py->size();
            count_GENERATED++;
			for(Int_t i = 0; i < NPARTGen; i++){
                count_GENERATED2++;
				part_px_gen[i] = gen_px->at(i);
				part_py_gen[i] = gen_py->at(i);
				part_pz_gen[i] = gen_pz->at(i);
				p4_ele_gen[i].SetPxPyPzE(part_px_gen[i], part_py_gen[i], part_pz_gen[i], sqrt(part_px_gen[i]*part_px_gen[i] + part_py_gen[i]*part_py_gen[i] + part_pz_gen[i]*part_pz_gen[i] + 0.000511*0.000511));
				
				
				float localPhiGen = p4_ele_gen[i].Phi()*toRD;
                SectorGen = GetSectorByPhi(localPhiGen);
				
                
				WGen = kin_W(p4_ele_gen[i],  eBeam);
				Q2Gen = kin_Q2(p4_ele_gen[i],  eBeam);

                 // *** НОВЫЙ CUT на генерированные kinematics *** - аналогично тому, что на реконструированные и данные, от Валеры
                if (WGen < 0.9 || WGen > 2.8 || Q2Gen < 2.0 || Q2Gen > 11.0) continue;   // выкидываем этот сгенерированный e⁻
                // *********************************************

                genEventExist = 1;
                
				int q2NumberLogQ2 = log(Q2Gen/q2Min)/deltaQ2;
                Q2_bin_log = q2NumberLogQ2;
				//generated events
                if (Q2Gen < 11. && WGen > 0.8) wQ2BinSectorGenLogQ2[0][q2NumberLogQ2]->Fill(WGen, 1/6.0);
                
                
                if (q2NumberLogQ2 >= 5 && q2NumberLogQ2 < 16 &&  WGen >= binFoldingParams.Wmin  &&  WGen <= binFoldingParams.Wmax ){
                
                	
                    
                    if (q2NumberLogQ2 == foldBinNumber) wGEN_fold->Fill(WGen, 1/6.0);
                    
                    int wGenBinTMP = 0, q2GenBinTMP = 0, skipValue1 = 0, skipValue2 = 0;
                    double skipRECq2 = 3., skipRECw = 1.5;
                    
                    getBinFold2D(skipRECw, WGen, skipRECq2,  Q2Gen, skipValue1, wGenBinTMP, skipValue2, q2GenBinTMP, binFoldingParams);
                    
                    wGEN_fold2Dcase->Fill(wBins_fold_2D * q2GenBinTMP + wGenBinTMP, 1/6.0);
                    
                    genBINmiss = wBins_fold_2D * q2GenBinTMP + wGenBinTMP;
                    
                    for (size_t iSec = 0; iSec < 6; iSec++){
                        xini_1D_gen[iSec][q2NumberLogQ2]->Fill(WGen, 1/6.0);
                        xini_2D_gen[iSec]->Fill(wBins_fold_2D * q2GenBinTMP + wGenBinTMP, 1/6.0);

                        if (WGen > wMin_fold && WGen < wMax_fold){
                            int wGenBin = (int)((WGen - wMin_fold)/wBinSize_fold_2D);
                            xini_1D_genQ2[iSec][wGenBin]->Fill(Q2Gen, 1/6.0);
                            xini_1D_genQ2_orig[iSec][wGenBin]->Fill(log(Q2Gen/q2Min), 1/6.0);
                        }
                    }
                    
                }
                
                
                
                const double thetaGenRad = p4_ele_gen[i].Theta()*toRD;
                const double pGenRad = p4_ele_gen[i].P();
                PGen_unf  = p4_ele_gen[i].P();
				
                
                
				int phiBinRes = (localPhiGen + 180)/5;
				int thetaBinRes = thetaGenRad/3;
				
				if (Q2_bin_log >= 5 && Q2_bin_log < 16) {
				
					generated1D_DATA[Q2_bin_log - 5] -> Fill(WGen, 1./6.);
					generated1D_DATA_sec[SectorGen][Q2_bin_log - 5] -> Fill(WGen);
				}
                
				if (true){
                    pGen[0]->Fill(p4_ele_gen[i].P());
                    CosthetaGen[0]->Fill(cos(thetaGenRad / toRD));
                    simPTHETHAcompBefore->Fill(p4_ele_gen[i].P(),thetaGenRad);
                    costhetaVSpGen->Fill(p4_ele_gen[i].P(),cos( thetaGenRad/ toRD));
                    
                    //Gen Event:
                    for (int sectorElectronGen = 0; sectorElectronGen < 6 ; sectorElectronGen++){
                    
                    
                        thetaGen[sectorElectronGen] -> Fill(thetaGenRad);
                        momGen[sectorElectronGen] -> Fill(pGenRad); 
                    
                        
                        phiVSwGen[sectorElectronGen] -> Fill(WGen, localPhiGen);
                        phiVSw_Q2binGen[sectorElectronGen][Q2_bin_log] -> Fill(WGen, localPhiGen);
                        
                        if (WGen > minW_ST && WGen < maxW_ST){
                            float binSize_phiVSq2 = (maxW_ST - minW_ST) / wBins_ST;
                            float wBinBuff = (WGen - minW_ST) / binSize_phiVSq2;
                            int wBin_phiVSq2 = (int)wBinBuff;
                            
                            phiVSq2_WbinGen[sectorElectronGen][wBin_phiVSq2] -> Fill(Q2Gen, localPhiGen);
                        }
                        
                        phiVSq2Gen[sectorElectronGen] -> Fill(Q2Gen, localPhiGen);
                        wVSq2Gen[sectorElectronGen] -> Fill(WGen, Q2Gen);
                        
                        thYieldGen[sectorElectronGen] -> Fill(thetaGenRad,1/6.);
                        phiYieldGen[sectorElectronGen] -> Fill(localPhiGen);
                        momYieldGen[sectorElectronGen] -> Fill(pGenRad,1/6.); 
                        thVSphiGen[sectorElectronGen] -> Fill(localPhiGen, thetaGenRad);
                        
                        thVSmomGen[sectorElectronGen] -> Fill(pGenRad, thetaGenRad, 1/6.);
                        phiVSmomGen[sectorElectronGen] -> Fill(pGenRad, localPhiGen);
                        
                        
                        if ( pGenRad > momMinSecStudies && pGenRad < momMaxSecStudies &&  
                           thetaGenRad > thetaMinSecStudies && thetaGenRad < thetaMaxSecStudies){
                            
                            const size_t pBinTMpgen = (int) ((pGenRad - momMinSecStudies) / pBinSizeTMP);
                            const size_t thBinTMPgen = (int) ((thetaGenRad - thetaMinSecStudies) / pThSizeTMP);
                                
                                
                            phiYield_P_Th_Gen[sectorElectronGen][pBinTMpgen][thBinTMPgen] -> Fill(localPhiGen);
                            thVSphi_PbinGen[sectorElectronGen][pBinTMpgen] -> Fill(localPhiGen, thetaGenRad);
                        }
                        
                    }
                    phiVSthFullRangeGen -> Fill(thetaGenRad, localPhiGen);
                    phiYieldGenALL -> Fill(localPhiGen);


                    wQ2InclusiveGen->Fill(WGen, Q2Gen);
                    if (q2NumberLogQ2 > 0 && q2NumberLogQ2 < qBinMax && Q2Gen < 11. && WGen > 0.8){
                        for (auto& item :  wGenVaryBinSize[0]){
                            item[q2NumberLogQ2]->Fill(WGen);
                        }
                    }
                }
			}
		}
        
/////////////////////////////////////////////////////////////////////////////////    
////////////// Data and Reconstructed Filling vectors ///////////////////////////
/////////////////////////////////////////////////////////////////////////////////     

		 
        
		int NPart = vpart_px->size();
		if (NPart > 0){
			for(Int_t i = 0; i < 1; i++){
                
				part_px[i] = vpart_px->at(i);
				part_py[i] = vpart_py->at(i);
				part_pz[i] = vpart_pz->at(i);
				part_E[i] = vpart_E->at(i);
				sectorE[i] = v_sectorE->at(i);
				part_vz[i] = vpart_vz->at(i);
				part_vx[i] = vpart_vx->at(i);
				part_vy[i] = vpart_vy->at(i);
                if (isData == 1){
		            sector_Sci[i] = sectorSci->at(i);
		            layer_Sci[i] = layerSci->at(i);
		            comp_Sci[i] = compSci->at(i);
                }
                // not tested:
                fillVfromTree(ele_chi2, eleChi2, i);
                fillVfromTree(ele_tr_chi2, eletrchi2, i);
                
                fillVfromTree(pcal_E, pcalE, i);
                fillVfromTree(ecin_E, ecinE, i);
                fillVfromTree(ecout_E, ecoutE, i);
                
                fillVfromTree(htcc_X, htccX, i);
                fillVfromTree(htcc_Y, htccY, i);
                fillVfromTree(htcc_NPE, htccNPE, i);
                
                fillVfromTree(pcal_HX, pcalHX, i);
                fillVfromTree(pcal_HY, pcalHY, i);
                fillVfromTree(pcal_HZ, pcalHZ, i);
                
                fillVfromTree(pcal_LV, pcalV, i);
                fillVfromTree(pcal_LW, pcalW, i);
                fillVfromTree(pcal_LU, pcalU, i);
                
                fillVfromTree(ele_ndf, eleNdf, i);
                
                // verified
                fillVfromTree(ecin_HX, ecinHX, i);
                fillVfromTree(ecin_HY, ecinHY, i);
                fillVfromTree(ecin_HZ, ecinHZ, i);

                fillVfromTree(ecout_HX, ecoutHX, i);
                fillVfromTree(ecout_HY, ecoutHY, i);
                fillVfromTree(ecout_HZ, ecoutHZ, i);
                
                fillVfromTree(ftof_HX, ftofHX, i);
                fillVfromTree(ftof_HY, ftofHY, i);
                
                fillVfromTree(dc_ZR1, dcZR1, i);
                fillVfromTree(dc_ZR2, dcZR2, i);
                fillVfromTree(dc_ZR3, dcZR3, i);
                
                
                // not tested:
                fillVfromTree(dc_XR1, dcXR1, i);
                fillVfromTree(dc_YR1, dcYR1, i);
                
                fillVfromTree(dc_XR2, dcXR2, i);
                fillVfromTree(dc_YR2, dcYR2, i);

                fillVfromTree(dc_XR3, dcXR3, i);
                fillVfromTree(dc_YR3, dcYR3, i);
                
/////////////////////////////////////////////////////////////////////////////////    
////////////// Data and Reconstructed Analysis Preparation //////////////////////
/////////////////////////////////////////////////////////////////////////////////                 

				if (sectorE[i] - 1 > -1 &&  sectorE[i] - 1 < 7){
                    
					p4_ele[i].SetPxPyPzE(part_px[i], part_py[i], part_pz[i], part_E[i]); 
                    
                    //smearing, the second parameter is charge so for elec -1
                    if (isData == 0 && isFXsmearing) smear(&p4_ele[i], -1);  
                    if (isData == 0 && isRichardSmearing) {
                    
		                TLorentzVector ele_smeared  = smear_func(p4_ele[i],  0);
		                p4_ele[i] = ele_smeared;
                    }
                    
                   if (isData == 0 && isValeriiSmearing) {
                      //auto fe = dpp_Sim(p4_ele[i].P()) + 1;
                      //double newP = fe * p4_ele[i].P();
                      //TLorentzVector ele_smeared;
                      //ele_smeared.SetPxPyPzE(fe*p4_ele[i].Px(), fe*p4_ele[i].Py(), fe*p4_ele[i].Pz(), newP);
                      //p4_ele[i] = smear_oneFactor(ele_smeared);
                      
                      TLorentzVector ele_smeared  = smear_oneFactor(p4_ele[i]);
		      p4_ele[i] = ele_smeared;
                    }
                    
                    
                    

                    // converting 1:6 to 0:5
					int sectorElectron = sectorE[i] - 1;
					if (p4_ele[i].P() > 1.5 && p4_ele[i].P() < 11){
                        
                        //DoMomCorr
                        if (isData == 1 && isMomCorrMay2023){
                            
                            //may 2023
                            auto fe    = dppC_may2023(p4_ele[i].Px(),p4_ele[i].Py(),p4_ele[i].Pz(), sectorElectron+1,   0, 1) + 1;
                            
                            
                            double newP = fe * p4_ele[i].P();
                            p4_ele[i].SetPxPyPzE(fe*p4_ele[i].Px(), fe*p4_ele[i].Py(), fe*p4_ele[i].Pz(), newP);
                            
                        }
                        
                        if (isData == 1 && isMomCorrOld){
                            
                            // before may 2023
                            auto fe = dppC(p4_ele[i].Px(),p4_ele[i].Py(),p4_ele[i].Pz(), sectorElectron+1, 0) + 1;
                            
                            double newP = fe * p4_ele[i].P();
                            p4_ele[i].SetPxPyPzE(fe*p4_ele[i].Px(), fe*p4_ele[i].Py(), fe*p4_ele[i].Pz(), newP);
                            
                        }
                        
                        if (isData == 0 && isMomCorrMC){

				        	//may 2023
				            auto fe    = dppC_may2023(p4_ele[i].Px(),p4_ele[i].Py(),p4_ele[i].Pz(), sectorElectron+1,   0, 2) + 1;
				                    
				            double newP = fe * p4_ele[i].P();
				            p4_ele[i].SetPxPyPzE(fe*p4_ele[i].Px(), fe*p4_ele[i].Py(), fe*p4_ele[i].Pz(), newP);
				                    
				        }
                        
// Kinematic variables calculation:                       
						float W = kin_W(p4_ele[i],  eBeam);
						float Q2 = kin_Q2(p4_ele[i],  eBeam);
						int q2NumberLogQ2 = log(Q2/q2Min)/deltaQ2;
						int pBin = (int)(p4_ele[i].P()/1);
//Theta mesured/reconstructed:                       
						float thetaEMeasured = p4_ele[i].Theta()*toRD;
						float pEMeasured = p4_ele[i].P();
//Phi mesured/reconstructed:                       
						float localPhi = p4_ele[i].Phi()*toRD;
						
                        
//DC phi rotating
                        float phiDC = localPhi;
                        if (phiDC < 0 && sectorElectron > 0) 
                            phiDC = phiDC + 360. - sectorElectron * 60.;
                        else
                            phiDC=phiDC - sectorElectron * 60.;
                        float phiDCnotCorr = phiDC; 
                        phiDC = phiDC - 30./ pEMeasured;
                        
//phi in P and theta:
                        const size_t pBinTMP = (int) ((pEMeasured - momMinSecStudies) / pBinSizeTMP);
                        const size_t thBinTMP = (int) ((thetaEMeasured - thetaMinSecStudies) / pThSizeTMP);
                        
// Filling containers for cuts:
                        KinemE kin = {p4_ele[i].P(), part_vz[i], ecout_E[i] + ecin_E[i] + pcal_E[i], ecin_E[i], 
                                      pcal_E[i], pBin, pcal_HX[i], pcal_HY[i],thetaEMeasured, phiDC, phiDCnotCorr, ecin_HX[i],ecin_HY[i], ecout_HX[i], ecout_HY[i], pcal_LV[i], pcal_LW[i], pcal_LU[i]};
                        
// DC rotating                        
						TVector3 dcR1Vector3(dc_XR1[i], dc_YR1[i], dc_ZR1[i]);
						TVector3 dcR2Vector3(dc_XR2[i], dc_YR2[i], dc_ZR2[i]);
						TVector3 dcR3Vector3(dc_XR3[i], dc_YR3[i], dc_ZR3[i]);
                        
                        rotVect(dcR1Vector3, sectorElectron);
                        rotVect(dcR2Vector3, sectorElectron);
                        rotVect(dcR3Vector3, sectorElectron);
                        
                        DCXY dcParam = {dcR1Vector3.X(), dcR1Vector3.Y(),
                                       dcR2Vector3.X(), dcR2Vector3.Y(),
                                       dcR3Vector3.X(), dcR3Vector3.Y()
                                        };
// PCAL rotating
						TVector3 resVector3(pcal_HX[i], pcal_HY[i], pcal_HZ[i]);
						resVector3.RotateZ(-60*sectorElectron/57.2958);
						resVector3.RotateY(-25/57.2958);
						
						
						TVector3 ecoutResVector3(ecout_HX[i], ecout_HY[i], ecout_HZ[i]);
						ecoutResVector3.RotateZ(-60*sectorElectron/57.2958);
						ecoutResVector3.RotateY(-25/57.2958);
                        //const double htccWeightConstOLD =  HTCCWeight(htcc_X[i], htcc_Y[i], HTCCEff, stepXHTCC, stepYHTCC, isData);
                        
                        if  (isData && layer_Sci[i] == 2 && sector_Sci[i] == 6){
                        	ftofNoCuts->Fill(comp_Sci[i] );
                        }
                        
                        const double ftofCorr = (isData && sector_Sci[i] == 6 && (comp_Sci[i] > 32 && comp_Sci[i] < 44)) ? 1.013 : 1.;

                        if (W < 0.9 || W > 2.8 || Q2 < 2 || Q2 > 11) continue;  // Valera suggested this FIX  to make the last step work.
                        
                        const double htccWeightConst =  newHTCCWeight(htcc_X[i], htcc_Y[i], HTCCEff, stepXHTCC, stepYHTCC, isData) * ftofCorr;
                        const double htccOnlyConst =  newHTCCWeight(htcc_X[i], htcc_Y[i], HTCCEff, stepXHTCC, stepYHTCC, isData);
                        
/////////////////////////////////////////////////////////////////////////////////    
////////////// Data and Reconstructed Filling Histograms without cuts ///////////
/////////////////////////////////////////////////////////////////////////////////                        
                        
						wQ2Inclusive->Fill(W, Q2);
						
						if (W>0.7 && W < 1.5 && sectorElectron >=0 && sectorElectron <=6 && q2NumberLogQ2 >=4 && q2NumberLogQ2 <=18){
							WvsP_elastic[sectorElectron]->Fill(p4_ele[i].P(),W);
							WinQ2_zoom[sectorElectron][q2NumberLogQ2]->Fill(W);
							
							if (p4_ele[i].P() > 2.5 && p4_ele[i].P() < 10){
								int pBin_elas = (int)((p4_ele[i].P()-2.5)/1.5);
								w_p_delta[sectorElectron][pBin_elas][0]->Fill(W);
							}
							
							w_q2_delta[sectorElectron][q2NumberLogQ2][0]->Fill(W);
						}
						
						if ((ecout_E[i] + ecin_E[i] + pcal_E[i]) > 1.)
							theta_ECAL[sectorElectron][0]->Fill(thetaEMeasured);
						
                        
                        ///EmptyCheckNoCuts
                        SFempty[sectorElectron][0]->Fill(p4_ele[i].P(), (ecout_E[i] + ecin_E[i] + pcal_E[i])/p4_ele[i].P());
                        VxEmpty[sectorElectron][0]->Fill(part_vx[i]);
                        VyEmpty[sectorElectron][0]->Fill(part_vy[i]);
                        VzEmpty[sectorElectron][0]->Fill(part_vz[i]);
                        NPEempty[sectorElectron][0]->Fill(htcc_NPE[i]);
                        Wempty[sectorElectron][0]->Fill(W);
                        //////////////
                        
                        wInclusiveBef[sectorElectron]->Fill(W);
                        wPInclusiveBef[sectorElectron]->Fill(p4_ele[i].P(), W);
                        if (W > 0.8 && W < 1.05)
                            pVSthetaMomCorr->Fill(p4_ele[i].P(), thetaEMeasured);
                        
                      
                        wQ2_norm_noCut[sectorElectron][q2NumberLogQ2]->Fill(W);
                        thetaRec_noCut[sectorElectron] -> Fill(thetaEMeasured);
                        momRec_noCut[sectorElectron] -> Fill(kin.momentum);
                        
                        Q2rec_noCut[sectorElectron] -> Fill(Q2);
    					phiRec_noCut[sectorElectron] -> Fill(phiDC);
    
    					WvsQ2_noCut -> Fill(W,Q2, htccWeightConst);
    					
    					
                        
//check Q2 binning
						if (q2NumberLogQ2 < qBinMax && q2NumberLogQ2 > 0  && ((W >= minWcut && W <= maxWcut) || !isKinemCut) && Q2 >= 2.0){
                            
                            
						sfCut_Sectors_bef[sectorElectron] -> Fill((ecout_E[i] + ecin_E[i] + pcal_E[i]) , (ecout_E[i] + ecin_E[i] + pcal_E[i]) / p4_ele[i].P()); 

                            xYDhtccDetecKnock->Fill(htcc_X[i],htcc_Y[i]);
                            xYDftofDetecKnock->Fill(ftof_HX[i],ftof_HY[i]);

                            

                            
                            
                            int wBin_phiVSq2 = wBins_ST +1;
                            
                            if (W >= minW_ST && W <= maxW_ST){
                                float binSize_phiVSq2 = (maxW_ST - minW_ST) / wBins_ST;
                                float wBinBuff = (W - minW_ST) / binSize_phiVSq2;
                                wBin_phiVSq2 = (int)wBinBuff;
                            }

                            
                            
/////////////////////////////////////////////////////////////////////////////////    
////////////// Data and Reconstructed Filling Histograms all cuts but DC fid cut/
////////////// No W and Q2 cut///////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////     
                            
/////////////////////////////////////////////////////////////////////////////////    
////////////// Data and Reconstructed Filling Histograms ////////////////////////
////////////// ALL CUTS ON //////////////////////////////////////////////////////
////////////// W and Q2 cut OFF /////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////// 
         
 /////////////////////////////////////////////////////////////////////////////////    
////////////// Data and Reconstructed W and Q2 ON////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////                            
                            
                            if (kinemCut(Q2, W)){
                                
                                
                                size_t isPCALfidCut = 0;
                                xyDC_nChi2[sectorElectron][0][isPCALfidCut] -> Fill(dcR1Vector3.X(), dcR1Vector3.Y(),
                                                                                    abs(ele_tr_chi2[i]));
                                xyDC_nChi2[sectorElectron][1][isPCALfidCut] -> Fill(dcR2Vector3.X(), dcR2Vector3.Y(),
                                                                                    abs(ele_tr_chi2[i]));
                                xyDC_nChi2[sectorElectron][2][isPCALfidCut] -> Fill(dcR3Vector3.X(), dcR3Vector3.Y(),
                                                                                    abs(ele_tr_chi2[i]));
                                
                                xyDC_ndf[sectorElectron][0][isPCALfidCut] -> Fill(dcR1Vector3.X(), dcR1Vector3.Y(), abs(ele_ndf[i]));
                                xyDC_ndf[sectorElectron][1][isPCALfidCut] -> Fill(dcR2Vector3.X(), dcR2Vector3.Y(), abs(ele_ndf[i]));
                                xyDC_ndf[sectorElectron][2][isPCALfidCut] -> Fill(dcR3Vector3.X(), dcR3Vector3.Y(), abs(ele_ndf[i])); 
                                
                                
                                pcalFID_SF_V[0][sectorElectron]-> Fill(kin.pcalV,kin.calEnerg/kin.momentum);
                                pcalFID_SF_W[0][sectorElectron]-> Fill(kin.pcalW,kin.calEnerg/kin.momentum);
                                pcalFID_SF_U[0][sectorElectron]-> Fill(kin.pcalU,kin.calEnerg/kin.momentum);
                                
/////////////////////////////////////////////////////////////////////////////////    
////////////// Data and Reconstructed Filling Sector Depend Histograms //////////
////////////// W and Q2 cut ON///////////////////////////////////////////////////
////////////// 1st Filling of SD/////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
                                
                                
                                // check that we have enough histograms to store all cut combinations
                                if (nCutCombination < vCutComb.size() + 1) 
                                    throw out_of_range("Increase nCutCombination becuse it cannot store all vCutComb");
                                
                                // cut combination counter (for hist filling) 
                                for (const auto iCutComb : vCutComb){
                                    
                                    if (ApplyCuts(p4_ele[i], triangleCutParams, kin, resVector3, dcParam, sectorElectron, 
                                                  isData, iCutComb.first)){
                                        
                                    //if(ApplyAllCutExcept(p4_ele[i], triangleCutParams, 
                                    //                     resVector3, kin, dcParam, sectorElectron, ivCut, isData)){
                                        
                                        int cutN = iCutComb.second;
                                        thVSmom[cutN][sectorElectron] -> Fill(p4_ele[i].P(), thetaEMeasured, htccWeightConst);
                                        thVSphi[cutN][sectorElectron] -> Fill(phiDCnotCorr, thetaEMeasured, htccWeightConst);
                                        phiVSmom[cutN][sectorElectron] -> Fill(p4_ele[i].P(), phiDCnotCorr, htccWeightConst);
                                        
                                        thYield[cutN][sectorElectron] -> Fill(thetaEMeasured, htccWeightConst);
                                        phiYield[cutN][sectorElectron] -> Fill(phiDC, htccWeightConst);
                                        momYield[cutN][sectorElectron] -> Fill(p4_ele[i].P(), htccWeightConst);
                                        
                                        phiVSw[cutN][sectorElectron] -> Fill(W, phiDCnotCorr, htccWeightConst);
                                        phiVSw_Q2bin[cutN][sectorElectron][q2NumberLogQ2] -> Fill(W, phiDCnotCorr, htccWeightConst);
                                        if (wBin_phiVSq2 >= 0 && wBin_phiVSq2 < static_cast<int>(wBins_ST)){
                                            phiVSq2_Wbin[cutN][sectorElectron][wBin_phiVSq2] -> Fill(Q2, phiDCnotCorr, htccWeightConst);
                                        }
                                        
                                        
                                        
                                        phiVSq2[cutN][sectorElectron] -> Fill(Q2, phiDCnotCorr, htccWeightConst);
                                        wVSq2[cutN][sectorElectron] -> Fill(W, Q2, htccWeightConst);
                                        
                                        if (thetaEMeasured > thetaMinSecStudies && thetaEMeasured < thetaMaxSecStudies
                                           && pEMeasured > momMinSecStudies && pEMeasured < momMaxSecStudies
                                           ){
                                            phiYield_P_Th[cutN][sectorElectron][pBinTMP][thBinTMP] -> Fill(phiDC, htccWeightConst);
                                            thVSphi_Pbin[cutN][sectorElectron][pBinTMP] -> Fill(phiDC, thetaEMeasured, htccWeightConst);
                                        }
                                        
                                        phiVSthFullRange[cutN] -> Fill(thetaEMeasured, phiDC, htccWeightConst);


/////////////////////////////////////////////////////////////////////////////////    
////////////// Data and Reconstructed Filling Histograms ////////////////////////
////////////// ALL CUTS /////////////////////////////////////////////////////////
////////////// W and Q2 cut On///////////////////////////////////////////////////
////////////// if 'if' is true then we use all CUt combination///////////////////
/////////////////////////////////////////////////////////////////////////////////  
                                        if (iCutComb.first == cuts_All){
                                        
                                        	if ((ecout_E[i] + ecin_E[i] + pcal_E[i]) > 1.)
												theta_ECAL[sectorElectron][1]->Fill(thetaEMeasured);
											
											if (W>0.7 && W < 1.5 && sectorElectron >=0 && sectorElectron <=6 && q2NumberLogQ2 >=4 && q2NumberLogQ2 <=18){
												if (p4_ele[i].P() > 2.5 && p4_ele[i].P() < 10){
													int pBin_elas = (int)((p4_ele[i].P()-2.5)/1.5);
													w_p_delta[sectorElectron][pBin_elas][1]->Fill(W);
												}
												w_q2_delta[sectorElectron][q2NumberLogQ2][1]->Fill(W);
											}
                                        
                                        	sfCut_Sectors_aft[sectorElectron] -> Fill((ecout_E[i] + ecin_E[i] + pcal_E[i]) , (ecout_E[i] + ecin_E[i] + pcal_E[i]) / p4_ele[i].P()); 
                                        
                                            if (W > wMin_fold && W < wMax_fold && Q2 > q2min_Fold && Q2 < q2max_Fold){
                                                thetaRec[sectorElectron] -> Fill(thetaEMeasured, htccWeightConst);
                                            	momRec[sectorElectron] -> Fill(kin.momentum, htccWeightConst);
		                                        Q2rec_Cut[sectorElectron] -> Fill(Q2, htccWeightConst);
												phiRec_Cut[sectorElectron] -> Fill(phiDC, htccWeightConst);
                                            }
                                            wQ2_norm[sectorElectron][q2NumberLogQ2]->Fill(W, htccWeightConst);

    
                                            wInclusiveAft[sectorElectron]->Fill(W);
                                            WvsQ2_Cut -> Fill(W,Q2, htccWeightConst);
                                            
                                            if (isData == 0) {
                                            	resP_WQ2[q2NumberLogQ2]->Fill(PGen_unf - p4_ele[i].P(),W);
                                            	resW_WQ2[q2NumberLogQ2]->Fill(WGen - W,W);
                                            }
                                            
                                            
                                            if (isUnfolding == false) continue;
                          

                                            // FID plots with PCAL and Theta vs Phi plot cut
                                            isPCALfidCut = 1;
                                            xyDC_noWeight[sectorElectron][0][isPCALfidCut] -> Fill(dcR1Vector3.X(), dcR1Vector3.Y());
                                            xyDC_noWeight[sectorElectron][1][isPCALfidCut] -> Fill(dcR2Vector3.X(), dcR2Vector3.Y());
                                            xyDC_noWeight[sectorElectron][2][isPCALfidCut] -> Fill(dcR3Vector3.X(), dcR3Vector3.Y());

                                            xyDC_nChi2[sectorElectron][0][isPCALfidCut] -> Fill(dcR1Vector3.X(), 
                                                                                                dcR1Vector3.Y(), abs(ele_tr_chi2[i]));
                                            xyDC_nChi2[sectorElectron][1][isPCALfidCut] -> Fill(dcR2Vector3.X(),
                                                                                                dcR2Vector3.Y(), abs(ele_tr_chi2[i]));
                                            xyDC_nChi2[sectorElectron][2][isPCALfidCut] -> Fill(dcR3Vector3.X(),
                                                                                                dcR3Vector3.Y(), abs(ele_tr_chi2[i]));

                                            xyDC_ndf[sectorElectron][0][isPCALfidCut] -> Fill(dcR1Vector3.X(), dcR1Vector3.Y(), 
                                                                                              abs(ele_ndf[i]));
                                            xyDC_ndf[sectorElectron][1][isPCALfidCut] -> Fill(dcR2Vector3.X(), dcR2Vector3.Y(), 
                                                                                              abs(ele_ndf[i]));
                                            xyDC_ndf[sectorElectron][2][isPCALfidCut] -> Fill(dcR3Vector3.X(), dcR3Vector3.Y(), 
                                                                                              abs(ele_ndf[i]));   

                                            xyPCAL_noWeight[sectorElectron][isPCALfidCut] -> Fill(resVector3.X(), resVector3.Y());
                                            xyECOUT_noWeight[sectorElectron][isPCALfidCut] -> Fill(ecoutResVector3.X(), ecoutResVector3.Y());
                                            
                                            
                                            pcalFID_SF_V[1][sectorElectron]-> Fill(kin.pcalV,kin.calEnerg/kin.momentum);
                                			pcalFID_SF_W[1][sectorElectron]-> Fill(kin.pcalW,kin.calEnerg/kin.momentum);
                                			pcalFID_SF_U[1][sectorElectron]-> Fill(kin.pcalU,kin.calEnerg/kin.momentum);

/////////////////////////////////////////////////////////////////////////////////    
////////////// ALL CUTS + KINEM /////////////////////////////////////////////////    
/////////////////////////////////////////////////////////////////////////////////    


										    if  (isData && layer_Sci[i] == 2 && sector_Sci[i] == 6){
										    	ftofAfterCuts->Fill(comp_Sci[i], htccWeightConst);
										    }


                                            count_RECONTRUCTED++;

                                            if (thetaEMeasured > 9 && kin.momentum < 11){
                                                  countDCthetaAF[sectorElectron][(int)thetaEMeasured - 9][((int)kin.momentum - 2)/3] -> Fill(phiDC, htccWeightConst);
                                                  countDCthetaAFnoP[sectorElectron][(int)thetaEMeasured - 9] -> Fill(phiDC, htccWeightConst);
                                            }

                                            //Resolution:
                                            wResQ2[q2NumberLogQ2]->Fill(WGen - W, htccWeightConst);
                                            q2Res->Fill(Q2Gen - Q2, htccWeightConst);


                                            int wBinRes = int((W - lowBorderW) / wBinSize);

                                            if (wBinRes > 0 && wBinRes < (int)nWBins && WGen > 0.9 && Q2Gen > 2.){
                                                q2ResAllBins[q2NumberLogQ2][wBinRes] -> Fill(Q2Gen - Q2, htccWeightConst);
                                                wResAllBins[q2NumberLogQ2][wBinRes] -> Fill(WGen - W, htccWeightConst);
                                            }
                                            //Bad Detectors:

                                            xYDCR1detec->Fill(dc_XR1[i],dc_YR1[i], htccWeightConst);
                                            xYDCR2detec->Fill(dc_XR2[i],dc_YR2[i], htccWeightConst);
                                            xYDCR3detec->Fill(dc_XR3[i],dc_YR3[i], htccWeightConst);
                                            
                                            xYDhtccDetec->Fill(htcc_X[i],htcc_Y[i], htccWeightConst);
                                            xYDftofDetec->Fill(ftof_HX[i],ftof_HY[i], htccWeightConst);
                                            
                                            xYDpcalDetec->Fill(pcal_HX[i],pcal_HY[i], htccWeightConst);
                                            xYDecinDetec->Fill(ecin_HX[i],ecin_HY[i], htccWeightConst);
                                            xYDecoutDetec->Fill(ecout_HX[i],ecout_HY[i], htccWeightConst);  

                                            wQ2BinSectorLogQ2StrictFidAllCuts[sectorElectron][q2NumberLogQ2]->Fill(W, htccWeightConst);
                                            
                                            
                                            //Unfolding filling:




                                          
                                          	if (q2NumberLogQ2 >= 5 && q2NumberLogQ2 < 16){
                                                
                                                //measured2D_secBySec[sectorElectron]->Fill(W, Q2, htccWeightConst);
                                                measured2D_DATA[sectorElectron]->Fill(W, Q2, htccWeightConst);
                                                measured1D_DATA[sectorElectron][q2NumberLogQ2 - 5]->Fill(W, htccWeightConst);
												measured1D_secBySec[sectorElectron][q2NumberLogQ2 - 5]->Fill(W, htccWeightConst);
												
												
												if (isData == 0){
													isMissEvent_roo = 0;
                                                	//Actual 2D deconv. obj
													rooResponse2D[sectorElectron]->Fill(W, Q2, WGen, Q2Gen, htccWeightConst);

                                                    ////////////////////////////////////////////////////НОВОЕ: Valera asked to check this
                                                    // Moved to Miss->
                                                    //////////////////////////////////////////////////

													//if (Q2_bin_log == q2NumberLogQ2 &&  WGen >= binFoldingParams.Wmin  &&  WGen <= binFoldingParams.Wmax ) 
													if (Q2_bin_log == q2NumberLogQ2) 
														rooResponse1D[sectorElectron][q2NumberLogQ2 - 5]->Fill(W, WGen, htccWeightConst);
													else
														rooResponse1D[sectorElectron][q2NumberLogQ2 - 5]->Fake(W);
                                                }
                                            }
                                            
                                            
                                            if (q2NumberLogQ2 >= 5 && q2NumberLogQ2 < 16 && W >= binFoldingParams.Wmin  &&  W < binFoldingParams.Wmax){
                                                   //&& (isData==1 || (Q2_bin_log >= 5 && Q2_bin_log < 15 &&  WGen >= binFoldingParams.Wmin  &&  WGen <= binFoldingParams.Wmax ) ){
                                                    
                                                    

                                                

                                                /// END RooUnfold Build-in functions
                                                    
                    
                                                int wGenBinTMP = 0, q2GenBinTMP = 0, wRecBinP = 0, q2RecBinTMP = 0;
                    
                                                getBinFold2D(W, WGen, Q2,  Q2Gen, wRecBinP, wGenBinTMP, 
                                                             q2RecBinTMP, q2GenBinTMP, binFoldingParams);
                                                
                                                
                                                float folding2D_yBin = wBins_fold_2D * q2RecBinTMP + wRecBinP;
                                                float folding2D_xBin = wBins_fold_2D * q2GenBinTMP + wGenBinTMP;
                                                //cout<< folding2D_yBin << ' ';
                                                
                                                bini_1D_rec[sectorElectron][q2NumberLogQ2]->Fill(W, htccWeightConst);
                                                bini_2D_rec[sectorElectron]->Fill(wBins_fold_2D * q2RecBinTMP + wRecBinP, htccWeightConst);
                                                //cout<< folding2D_yBin << ' ';
                                                if (W > wMin_fold && W < wMax_fold){
                                                    int wBinUnfold = (int)((W - wMin_fold)/wBinSize_fold_2D);
                                                    bini_1D_recQ2[sectorElectron][wBinUnfold]->Fill(Q2, htccWeightConst);
                                                    bini_1D_recQ2_orig[sectorElectron][wBinUnfold]->Fill(log(Q2/q2Min), htccWeightConst);
                                                }

                    
                                                //it is not W anymore
                                                wREC_fold2Dcase->Fill(folding2D_yBin, htccWeightConst);
                                                //cout<< folding2D_yBin << ' ';
                                                
                                            	hMeas[sectorElectron]->Fill(folding2D_yBin, htccWeightConst);
                                            	
                                            	hMeasX[sectorElectron]->Fill(folding2D_xBin, htccWeightConst);
                                            	hMeasY[sectorElectron]->Fill(folding2D_yBin, htccWeightConst);
                                            		
                                            	//cout<<'\n';
                    
                                                if (isData == 0){
                                                
                                                	isMissEvent = 0;
                                                    responseMatrix2Dcase->Fill(folding2D_xBin, folding2D_yBin, htccWeightConst);
                                                    Adet_2D_genVSrec[sectorElectron]->Fill(folding2D_xBin, folding2D_yBin, htccWeightConst);
                                                    Adet_1D_genVSrec[sectorElectron][q2NumberLogQ2]->Fill(W, WGen, htccWeightConst);
                                                    
                                                                                                
				                                    
				                       
				                                    wGEN_testBefWeight[sectorElectron]->Fill(WGen, htccWeightConst);
				                                    
				                                    if (folding2D_xBin < MatrixDimens2D && folding2D_xBin > -0.1)
				                                    	hResponseRoo[sectorElectron]->Fill(folding2D_yBin, folding2D_xBin, htccWeightConst);
				                                    else
				                                    	hResponseRoo[sectorElectron]->Fake(folding2D_yBin, htccWeightConst);
				                                    
                                                    if (W > wMin_fold && W < wMax_fold){
                                                        int wBinUnfold = (int)((W - wMin_fold)/wBinSize_fold_2D);
                                                        Adet_1D_genVSrecQ2[sectorElectron][wBinUnfold]->Fill(Q2, Q2Gen, htccWeightConst);
                                                        Adet_1D_genVSrecQ2_orig[sectorElectron][wBinUnfold]->Fill(log(Q2/q2Min), log(Q2Gen/q2Min), htccWeightConst);
                                                    }
                                                
                                                
                                                if (q2NumberLogQ2 == foldBinNumber) {
                                                    wREC_fold->Fill(W, htccWeightConst);
                                                    if (isData == 0){
                                                        responseMatrix->Fill(WGen, W, htccWeightConst);
                                                    }
                                                }
                                                
                                                }
                                            }//unfolding fitting
                                            
                                        }
/////////////////////////////////////////////////////////////////////////////////    
///////////////////////////////// No triag and Fid //////////////
/////////////////////////////////////////////////////////////////////////////////
										if (iCutComb.first == cuts_NoAngleNoFidNoTrig){
											triagCut[sectorElectron][pBin][0]->Fill(ecin_E[i]/p4_ele[i].P(), pcal_E[i]/p4_ele[i].P());
										}
/////////////////////////////////////////////////////////////////////////////////    
///////////////////////////////// cuts_NoAngleNoFid  //////////////
/////////////////////////////////////////////////////////////////////////////////             
                                        if (iCutComb.first == cuts_NoAngleNoFid){
                                        	xYDpcalDetecKnock->Fill(pcal_HX[i],pcal_HY[i], htccWeightConst);
                                        	xYDCR1detecKnock->Fill(dc_XR1[i],dc_YR1[i], htccWeightConst);
										    xYDCR2detecKnock->Fill(dc_XR2[i],dc_YR2[i], htccWeightConst);
										    xYDCR3detecKnock->Fill(dc_XR3[i],dc_YR3[i], htccWeightConst);
										    
										    xYDecinDetecKnock->Fill(ecin_HX[i],ecin_HY[i], htccWeightConst);
                            				xYDecoutDetecKnock->Fill(ecout_HX[i],ecout_HY[i], htccWeightConst); 
										    
										    
										    size_t isPCALfidCut = 0;
						                    xyDC_noWeight[sectorElectron][0][isPCALfidCut] -> Fill(dcR1Vector3.X(), dcR1Vector3.Y(), htccWeightConst);
						                    xyDC_noWeight[sectorElectron][1][isPCALfidCut] -> Fill(dcR2Vector3.X(), dcR2Vector3.Y(), htccWeightConst);
						                    xyDC_noWeight[sectorElectron][2][isPCALfidCut] -> Fill(dcR3Vector3.X(), dcR3Vector3.Y(), htccWeightConst);
						                    xyPCAL_noWeight[sectorElectron][isPCALfidCut] -> Fill(resVector3.X(), resVector3.Y(), htccWeightConst);
						                    xyECOUT_noWeight[sectorElectron][isPCALfidCut] -> Fill(ecoutResVector3.X(), ecoutResVector3.Y(), htccWeightConst);
						                    
						                    
                                        }
                                        
/////////////////////////////////////////////////////////////////////////////////    
///////////////////////////////// all Cuts no HTCC ////////////////////////////////////////////////
/////////////////////////////////////////////////////////////    
                                        
/////////////////////////////////////////////////////////////////////////////////    
////////////// END ALL CUTS + KINEM /////////////////////////////////////////////    
/////////////////////////////////////////////////////////////////////////////////                                 

                                        
                                    }// loop over cut combination
                                }
                            } // kinem
                            
/////////////////////////////////////////////////////////////////////////////////    
////////////// Systematics //////////////////////////////////////////////////////    
////////////// Prepare all cut combinations /////////////////////////////////////   
/////////////////////////////////////////////////////////////////////////////////

                            
/////////////////////////////////////////////////////////////////////////////////    
////////////// Systematics Calculation //////////////////////////////////////////
////////////// PCAL /////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////  
                        if (isSystematics){
                        
                        	int wGenBinTMP = 0, q2GenBinTMP = 0, wRecBinP = 0, q2RecBinTMP = 0;
                        	float folding2D_yBin = 0, folding2D_xBin = 0;
                        	if (q2NumberLogQ2 >= 5 && q2NumberLogQ2 < 16 &&
                               	(isData==1 || (Q2_bin_log >= 5 && Q2_bin_log < 16 && 
                               	WGen >= binFoldingParams.Wmin  &&  WGen <= binFoldingParams.Wmax ) ) && 
                                W >= binFoldingParams.Wmin  &&  W < binFoldingParams.Wmax
                               ){
                    
                                
                    			getBinFold2D(W, WGen, Q2,  Q2Gen, wRecBinP, wGenBinTMP, q2RecBinTMP, q2GenBinTMP, binFoldingParams);
                                                
                                folding2D_yBin = wBins_fold_2D * q2RecBinTMP + wRecBinP;
                                folding2D_xBin = wBins_fold_2D * q2GenBinTMP + wGenBinTMP;
                        		}
                        		else continue;
                        
                            
                            if (ApplyCuts(p4_ele[i], triangleCutParams, kin, resVector3, dcParam, sectorElectron, isData, cuts_NoPCALcut)){
                            
                            //if(ApplyAllCutExcept(p4_ele[i], triangleCutParams, resVector3, kin, dcParam, sectorElectron, cutType::pcalCut, isData)){
                                for (int k = 0; k < 3; k++){
									if (PCALFidXY(resVector3.X(), resVector3.Y(), k) == 1){
                                        //pcal fiducial
                                        hResponseRooSys[0][sectorElectron][k]->Fill(folding2D_yBin, folding2D_xBin, htccWeightConst);
                                        hMeasSys[0][sectorElectron][k]->Fill(folding2D_yBin, htccWeightConst);
                                    }
                                }
                            }
                            
/////////////////////////////////////////////////////////////////////////////////    
////////////// Systematics Calculation //////////////////////////////////////////
////////////// SF /////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////// 
                            
                            if (ApplyCuts(p4_ele[i], triangleCutParams, kin, resVector3, dcParam, sectorElectron, isData, cuts_NoSFcut)){
                            //if(ApplyAllCutExcept(p4_ele[i], triangleCutParams, resVector3, kin, dcParam, sectorElectron, cutType::sfCut, isData)){
       sfCutProb[sectorElectron][0]->Fill(p4_ele[i].P(), (ecout_E[i] + ecin_E[i] + pcal_E[i]) / p4_ele[i].P());                         
                                for (int k = 0; k < 3; k++){
                                    //samgling fraction systematics
									if (SfCutValerii_Edepos( (ecout_E[i] + ecin_E[i] + pcal_E[i]) /  p4_ele[i].P(),  ecout_E[i] + ecin_E[i] + pcal_E[i], sectorElectron, k, isData)){
										hResponseRooSys[1][sectorElectron][k]->Fill(folding2D_yBin, folding2D_xBin, htccWeightConst);
										hMeasSys[1][sectorElectron][k]->Fill(folding2D_yBin, htccWeightConst);
									}
                                }
                            }
                            
/////////////////////////////////////////////////////////////////////////////////    
////////////// Systematics Calculation //////////////////////////////////////////
////////////// TRIAG /////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
                            
                            if (ApplyCuts(p4_ele[i], triangleCutParams, kin, resVector3, dcParam, sectorElectron, isData, cuts_NoTriag)){
                            //if(ApplyAllCutExcept(p4_ele[i], triangleCutParams, resVector3, kin, dcParam, sectorElectron, cutType::triagCut, isData)){
                                triagCut[sectorElectron][pBin][0]->Fill(ecin_E[i]/p4_ele[i].P(), pcal_E[i]/p4_ele[i].P());
								if (SFTriangleCut(ecin_E[i]/p4_ele[i].P(), pcal_E[i]/p4_ele[i].P(), triangleCutParams, sectorElectron, pBin, -0.005)) {
										hResponseRooSys[2][sectorElectron][0]->Fill(folding2D_yBin, folding2D_xBin, htccWeightConst);
										hMeasSys[2][sectorElectron][0]->Fill(folding2D_yBin, htccWeightConst);
                                    
								}
                                if (SFTriangleCut(ecin_E[i]/p4_ele[i].P(), pcal_E[i]/p4_ele[i].P(), triangleCutParams, sectorElectron, pBin, 0)) {
										hResponseRooSys[2][sectorElectron][1]->Fill(folding2D_yBin, folding2D_xBin, htccWeightConst);
										hMeasSys[2][sectorElectron][1]->Fill(folding2D_yBin, htccWeightConst);
                                    
								}
								if (SFTriangleCut(ecin_E[i]/p4_ele[i].P(), pcal_E[i]/p4_ele[i].P(), triangleCutParams, sectorElectron, pBin, 0.005)) {
										hResponseRooSys[2][sectorElectron][2]->Fill(folding2D_yBin, folding2D_xBin, htccWeightConst);
										hMeasSys[2][sectorElectron][2]->Fill(folding2D_yBin, htccWeightConst);
								}
                            }
                            
/////////////////////////////////////////////////////////////////////////////////    
////////////// Systematics Calculation //////////////////////////////////////////
////////////// VZ /////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
                            
                            if (ApplyCuts(p4_ele[i], triangleCutParams, kin, resVector3, dcParam, sectorElectron, isData, cuts_NoVz)){
                            //if(ApplyAllCutExcept(p4_ele[i], triangleCutParams, resVector3, kin, dcParam, sectorElectron, cutType::vzCut, isData)){
                                for (int k = 0; k < 3; k++){
                                    if (CutVz(kin.vz,k)){
										hResponseRooSys[3][sectorElectron][k]->Fill(folding2D_yBin, folding2D_xBin, htccWeightConst);
										hMeasSys[3][sectorElectron][k]->Fill(folding2D_yBin, htccWeightConst);
										}
                                }
                            }
                            
/////////////////////////////////////////////////////////////////////////////////    
////////////// Systematics Calculation //////////////////////////////////////////
////////////// Momentum /////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/*
    vzCut = 0,  // +
    dcCut = 1, // +
    sfCut = 2, //+
    momCut = 3, // no need
    pcalCut = 4, //+
    triagCut = 5, //+
    badElemCut = 6, //+
    pcalDeposCut = 7,//+
    thetaMax = 8,// no need
    thetaMin = 9,// no need
    thetaVSphiDC = 10,// no need
    phiMax = 11,// no need
    phiMin = 12,// no need
    phiSpike = 13,
    allCuts = 14   //+                
    */           
/////////////////////////////////////////////////////////////////////////////////    
////////////// Systematics Calculation //////////////////////////////////////////
////////////// DC fid ///////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
                            
                            if (ApplyCuts(p4_ele[i], triangleCutParams, kin, resVector3, dcParam, sectorElectron, isData, cuts_NoDCcut)){
                            //if(ApplyAllCutExcept(p4_ele[i], triangleCutParams, resVector3, kin, dcParam, sectorElectron, cutType::dcCut, isData)){
                                for (int k = 0; k < 3; k++){
                                    if (CutDCfid(dcParam,sectorElectron,k)){
										hResponseRooSys[4][sectorElectron][k]->Fill(folding2D_yBin, folding2D_xBin, htccWeightConst);
										hMeasSys[4][sectorElectron][k]->Fill(folding2D_yBin, htccWeightConst);
										}
                                }
                            }
                            
/////////////////////////////////////////////////////////////////////////////////    
////////////// Systematics Calculation //////////////////////////////////////////
////////////// PCAL DEPOSITED ENERGY ////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
                            
                            if (ApplyCuts(p4_ele[i], triangleCutParams, kin, resVector3, dcParam, sectorElectron, isData, cuts_NoDeposE)){
                            //if(ApplyAllCutExcept(p4_ele[i], triangleCutParams, resVector3, kin, dcParam, sectorElectron, cutType::pcalDeposCut, isData)){
                                for (int k = 0; k < 3; k++){
                                    if (CutPCALdepos(kin.pcalE,k)){
										hResponseRooSys[5][sectorElectron][k]->Fill(folding2D_yBin, folding2D_xBin, htccWeightConst);
										hMeasSys[5][sectorElectron][k]->Fill(folding2D_yBin, htccWeightConst);
										}
                                }
                            }
                            
/////////////////////////////////////////////////////////////////////////////////    
////////////// Systematics Calculation //////////////////////////////////////////
////////////// BAD STRIPS PCAL ////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
                            
                            if (ApplyCuts(p4_ele[i], triangleCutParams, kin, resVector3, dcParam, sectorElectron, isData, cuts_NoBADelem)){
                            //if(ApplyAllCutExcept(p4_ele[i], triangleCutParams, resVector3, kin, dcParam, sectorElectron, cutType::badElemCut, isData)){
                                for (int k = 0; k < 3; k++){
                                    if (BadElementKnockOut(pcal_HX[i],pcal_HY[i], ecin_HX[i],ecin_HY[i], ecout_HX[i], ecout_HY[i], sectorElectron + 1, k)){
										hResponseRooSys[6][sectorElectron][k]->Fill(folding2D_yBin, folding2D_xBin, htccWeightConst);
										hMeasSys[6][sectorElectron][k]->Fill(folding2D_yBin, htccWeightConst);
									}
                                }
                            }                            
                            
/////////////////////////////////////////////////////////////////////////////////    
////////////// Systematics Calculation //////////////////////////////////////////
////////////// Theta Min ////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////////    
////////////// Systematics Calculation //////////////////////////////////////////
////////////// Theta Max ////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

                            
/////////////////////////////////////////////////////////////////////////////////    
////////////// Systematics Calculation //////////////////////////////////////////
////////////// Theta vs Phi DC fid //////////////////////////////////////////// 
/////////////////////////////////////////////////////////////////////////////////
                      /*      
                            if (ApplyCuts(p4_ele[i], triangleCutParams, kin, resVector3, dcParam, sectorElectron, isData, cuts_NoThVSphi)){
                                for (int k = 0; k < 3; k++){
                                    if (thetaVSphiDCcut(kin.theta, kin.phi, k))
                                        wQ2BinSectorLogQ2AllCutsAllSys[static_cast<int>(cutType::phiMax)][sectorElectron][q2NumberLogQ2][k]->Fill(W, htccWeightConst);
                                }
                            }
                        */    
/////////////////////////////////////////////////////////////////////////////////    
////////////// Systematics Calculation //////////////////////////////////////////
////////////// Phi max DC cut //////////////////////////////////////////// 
/////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////////    
////////////// Systematics Calculation //////////////////////////////////////////
////////////// Phi spike //////////////////////////////////////////// 
/////////////////////////////////////////////////////////////////////////////////
                            if (ApplyCuts(p4_ele[i], triangleCutParams, kin, resVector3, dcParam, sectorElectron, isData, cuts_NoPhiSpike)){
                            	for (int k = 0; k < 3; k++){
		                        	if (phiSpikeCut(kin.phi_notCorrected, kin.theta, sectorElectron, k)){
										hResponseRooSys[7][sectorElectron][k]->Fill(folding2D_yBin, folding2D_xBin, htccWeightConst);
										hMeasSys[7][sectorElectron][k]->Fill(folding2D_yBin, htccWeightConst);
									}
								}

                            } 
/////////////////////////////////////////////////////////////////////////////////    
////////////// ALL Cuts CS for check //////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
                            
                            if (ApplyCuts(p4_ele[i], triangleCutParams, kin, resVector3, dcParam, sectorElectron, isData, cuts_All)){
                                for (int k = 0; k < 3; k++){
										hResponseRooSys[8][sectorElectron][k]->Fill(folding2D_yBin, folding2D_xBin, htccOnlyConst);
										hMeasSys[8][sectorElectron][k]->Fill(folding2D_yBin, htccOnlyConst);
										
										hResponseRooSys[9][sectorElectron][k]->Fill(folding2D_yBin, folding2D_xBin);
										hMeasSys[9][sectorElectron][k]->Fill(folding2D_yBin);
										
										hResponseRooSys[10][sectorElectron][k]->Fill(folding2D_yBin, folding2D_xBin, htccWeightConst);
										hMeasSys[10][sectorElectron][k]->Fill(folding2D_yBin, htccWeightConst);
										
										
										
                                }
                            }
                        }// end if sys
                            
/////////////////////////////////////////////////////////////////////////////////    
////////////// Systematics Calculation //////////////////////////////////////////
////////////// END ///////////////// ////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
						//Q
                        }
                        //p
					}
                    //sector
				}
                //particle
			}
            //Npart
		}
		if (k%3000000 == 0) {
            std::cout.precision(2);
            cout << "File:" << nfile << "/" << file_list.size() << " EventsFile:" <<  k/1000000 << "/" << n_entriesTree/1000000 << "M";
            if(!isData) cout << ", total GEN:" << count_GENERATED2/1000000 <<"M, REC:" << count_RECONTRUCTED/1000000 << "M, REC/GEN:" << (float)count_RECONTRUCTED/count_GENERATED2;
            cout <<", Est. Done:" << 100*((float)k/(n_entriesTree * file_list.size()) +((float)nfile-1) / file_list.size()) << "%" << endl;
        }
		if(k == nEvents) break;
		//if(k > 10000000) break;
		
		
		
		if (genEventExist){ 
		
			if (isMissEvent){
				hResponseRoo[SectorGen]->Miss(genBINmiss);
				
				for (size_t iCutSys = 0; iCutSys < nCutSys; iCutSys ++){
					for (size_t iLevel  = 0; iLevel < 3; iLevel ++){
						hResponseRooSys[iCutSys][SectorGen][iLevel]->Miss(genBINmiss);
					}
				}
			}
			
			if (isMissEvent_roo){
				rooResponse2D[SectorGen]->Miss(WGen, Q2Gen);
                for_Valera_Q2_VS_W_gen[SectorGen]->Fill(WGen, Q2Gen);
                //for_Valera_Q2_VS_W[sectorElectron]->Fill(W, Q2);

				if (Q2_bin_log>=5 && Q2_bin_log<16) rooResponse1D[SectorGen][Q2_bin_log - 5]->Miss(WGen);
			}
		 }
		
	}
        f->Close();
        nfile++;
       // if (nfile==2) break;
    }
    cout<<" GENERATED: "<< count_GENERATED <<endl;
    cout<<" GENERATED2: "<< count_GENERATED2 <<endl;
    cout<<" RECONSTRUCTED: "<< count_RECONTRUCTED <<endl;
    
	cout << "writing"<<endl;
	float g, r;
	

	
	
	resultsF->cd("overview");
	
	
	for (size_t iSec = 0 ; iSec < secN; iSec++){
		for (size_t iTh = 0 ; iTh < 2; iTh++){
			theta_ECAL[iSec][iTh]->Write();
			for (size_t iP = 0 ; iP < 5; iP++){
				w_p_delta[iSec][iP][iTh]->Write();
			}
			for (int q = 0; q < qBinMax; q++){	
				w_q2_delta[iSec][q][iTh]->Write();
			}							
		}
	}
	
	///EmptyRuns:
    for (size_t iSec = 0; iSec < 6; iSec++){
    
    	WvsP_elastic[iSec]-> Write();
    
        for (size_t iCut = 0; iCut < 3; iCut++){
            for (size_t iPbin = 0; iPbin <12 ; iPbin++){
                    triagCut[iSec][iPbin][iCut]->Write();
                }
                sfCutProb[iSec][iCut]->Write();
        }
        
        for (size_t iCut = 0; iCut < 2; iCut++){
            SFempty[iSec][iCut]->Write();
            VxEmpty[iSec][iCut]->Write();
            VyEmpty[iSec][iCut]->Write();
            VzEmpty[iSec][iCut]->Write();;
            NPEempty[iSec][iCut]->Write();
            Wempty[iSec][iCut]->Write();
        }
    }
	
	//FTOF:
	ftofNoCuts->Write();
	ftofAfterCuts->Write();
	
	
	WvsQ2_noCut -> Write();
	WvsQ2_Cut -> Write();
	
	xYDCR1detecKnock->Write();
    xYDCR2detecKnock->Write();
    xYDCR3detecKnock->Write();
    xYDCR1detec->Write();
    xYDCR2detec->Write();
    xYDCR3detec->Write();
    xYDpcalDetec->Write();
    xYDpcalDetecKnock->Write();
    xYDecinDetecKnock->Write();
    xYDecoutDetecKnock->Write();
    xYDecinDetec->Write();
    xYDecoutDetec->Write();
    
    
    for (size_t iSec = 0 ; iSec < 6; iSec++){
        for (size_t iLayer = 0; iLayer < 3; iLayer++){
            for (size_t isCut = 0; isCut < 2; isCut++){
                
                xyDC_noWeight[iSec][iLayer][isCut] -> Write();
                if (iLayer == 0)
                    xyPCAL_noWeight[iSec][isCut] -> Write();
                    xyECOUT_noWeight[iSec][isCut] -> Write();
            }
        }
    }
	
	for (size_t iSec = 0; iSec < 6; iSec++){
        thetaRec[iSec]-> Write();
        momRec[iSec]-> Write();

        thetaRec_noCut[iSec]-> Write();
        momRec_noCut[iSec]-> Write();
        
        Q2rec_noCut[iSec]-> Write();
    	Q2rec_Cut[iSec]-> Write();
    
    	phiRec_noCut[iSec]-> Write();
    	phiRec_Cut[iSec]-> Write();
    
        
        for (int iQ2 = 0; iQ2 < qBinMax; iQ2++){
        	WinQ2_zoom[iSec][iQ2]-> Write();
		    wQ2_norm[iSec][iQ2]-> Write();
	   		wQ2_norm_noCut[iSec][iQ2]-> Write();
        }
    }
    
       for (size_t iSec = 0; iSec < 6; iSec++){
        for (size_t iCut = 0; iCut < 3; iCut++){
            for (size_t iPbin = 0; iPbin <12 ; iPbin++){
                    triagCut[iSec][iPbin][iCut]->Write();
                }
                sfCutProb[iSec][iCut]->Write();
        }
    }
    
    
/////////////////////////
	if (isData == 0){
        for (int t = 0; t < 6; t++){
            thetaGen[t]-> Write();
            momGen[t]-> Write();
        }
    }
    
    
	
	
	/*
	
	wQ2Inclusive->Write();
	wQ2InclusiveGen->Write();
    //resolution:
    for (int iplot = 0; iplot < qBinMax; iplot++){
        wResQ2[iplot]->Write();
        for (int iW = 0; iW < (int)nWBins; iW++){
            q2ResAllBins[iplot][iW] -> Write();
            wResAllBins[iplot][iW]  -> Write();
        }
    }
    q2Res->Write();
    
    
   
    
    //Detectors:
    

    xYDhtccDetec->Write();
    xYDftofDetec->Write();
    
    xYDecinDetec->Write();
    xYDecoutDetec->Write();
    

    
    
    xYDhtccDetecKnock->Write();
    xYDftofDetecKnock->Write();
    
    xYDecinDetecKnock->Write();
    xYDecoutDetecKnock->Write();
    
    xYDCR1detecNo_BadEl_PhiSp->Write();
    xYDCR2detecNo_BadEl_PhiSp->Write();
    xYDCR3detecNo_BadEl_PhiSp->Write();
    xYDhtccDetecNo_BadEl_PhiSp->Write();
    xYDftofDetecNo_BadEl_PhiSp->Write();
    xYDpcalDetecNo_BadEl_PhiSp->Write();
    xYDecinDetecNo_BadEl_PhiSp->Write();
    xYDecoutDetecNo_BadEl_PhiSp->Write();    
    
    
    
    //DC fid:
    pVSthetaMomCorr->Write();
    for (size_t iSec = 0; iSec < 6; iSec++){
        wInclusiveBef[iSec]->Write();
        wInclusiveAft[iSec]->Write();
        wPInclusiveBef[iSec]->Write();
        
        for (size_t iTheta = 0; iTheta < nThetaBins; iTheta++){
            countDCthetaBFnoP[iSec][iTheta]->Write();
            countDCthetaAFnoP[iSec][iTheta]->Write();
            countDCthetaBFstephan[iSec][iTheta]->Write();
            for (int ip = 0; ip < npBins; ip++){
                countDCthetaAF[iSec][iTheta][ip] ->Write();
                countDCthetaBF[iSec][iTheta][ip] ->Write();
            }
        }
    }
        

    //////
    
	cout << "writing 2"<<endl; 
    for (int t = 0; t < 6; t++){
        for(size_t n_Wdiv = 0; n_Wdiv < wNumBinSize; n_Wdiv++){
            for (int q = 0; q < qBinMax; q++){
                WQ2_WRecAndWGen[t][n_Wdiv][q]->Write();
                wQ2BinSectorLogQ2StrictFidAllCutsVaryBinSize[t][n_Wdiv][q]->Write();
                wGenVaryBinSize[t][n_Wdiv][q]->Write();
            }
        }
    }
	cout << "writing 3"<<endl;
    ///////////////////for sim comp//////////
    simPecomp->Write();
    simQ2comp->Write();
    simWcomp->Write();
    simPHIcomp->Write();
    simTHETAcomp->Write();
    
    simQ2Wcomp->Write();
    simPTHETHAcomp->Write();
    simPTHETHAcompBefore->Write();
    costhetaVSpGen->Write();
    simPPHIcomp->Write();
    for (int isec = 0; isec < 6; isec++){
        CosthetaGen[isec]->Write();
        
    }
    
/////////////////////////


	if (isData == 0){
        for (int iCut = 0; iCut < nCut; iCut++){	
            for (int iSec = 0; iSec < 6; iSec++){	
                for (int iQ2 = 0; iQ2 < qBinMax; iQ2++){
                    for (int iCutLevel = 0; iCutLevel < 3; iCutLevel++){
                        for (int ibin = 0; ibin < wQ2BinSectorLogQ2AllCutsAllSys[0][0][0][0]->GetNbinsX(); ibin++){
                            g = wQ2BinSectorGenLogQ2[0][iQ2]->GetBinContent(ibin + 1, 1);
                            r = wQ2BinSectorLogQ2AllCutsAllSys[iCut][iSec][iQ2][iCutLevel]->GetBinContent(ibin + 1);
                            if (g > 0){
                                wQ2BinSectorAccLogQ2AllCutsAllSys[iCut][iSec][iQ2][iCutLevel]->SetBinContent(ibin + 1,r/g);
                                wQ2BinSectorAccLogQ2AllCutsAllSys[iCut][iSec][iQ2][iCutLevel]->SetBinError(ibin + 1,1/g*sqrt(r*r/g + r));
                            }
                        }
                    }
                }
            }
        } 
        
        
		for (int t = 0; t < 6; t++){
			for (int q = 0; q < qBinMax; q++){
				for (int d = 0; d < 250; d++){
					g = wQ2BinSectorGenLogQ2[0][q]->GetBinContent(d + 1, 1);
					r = wQ2BinSectorLogQ2StrictFidAllCuts[t][q]->GetBinContent(d + 1);
					if (g > 0){
						wQ2BinSectorAccLogQ2StrictFidAllCuts[t][q]->SetBinContent(d + 1, r/g);
						wQ2BinSectorAccLogQ2StrictFidAllCuts[t][q]->SetBinError(d + 1,1/g*sqrt(r*r/g + r));
					}
				}
			}
		}        
	cout << "writing 4"<<endl;

		for (int k = 0; k < 3; k++){
			for (int t = 0; t < 6; t++){
				for (int q = 0; q < qBinMax; q++){
					for (int d = 0; d < 250; d++){
						g = wQ2BinSectorGenLogQ2[0][q]->GetBinContent(d + 1, 1);
						r = wQ2BinSectorLogQ2AllCutsSFSys[t][q][k]->GetBinContent(d + 1);
						if (g > 0 ){
							wQ2BinSectorAccLogQ2AllCutsSFSys[t][q][k]->SetBinContent(d + 1, r/g);
							wQ2BinSectorAccLogQ2AllCutsSFSys[t][q][k]->SetBinError(d + 1,1/g*sqrt(r*r/g + r));
						}
					}
				}
			}
		}
    }
    
    
	cout << "writing 5"<<endl;

	for (int t = 0; t < 6; t++){
		pGen[t]->Write();

	}

	resultsF->cd("HTCC");

	for (int xC = 0; xC < nHTCCXBins; xC++){
		for (int yC = 0; yC < nHTCCXBins; yC++){
			npe1D[xC][yC]->Write();
		}
	}
	cout << "writing 6"<<endl;

	resultsF->cd("ResultsThetaPhi");
    
    phiVSthFullRangeGen -> Write();
    
    for (int t = 0; t < 6; t++){
        
        thYieldGen[t] -> Write(); 
        phiYieldGen[t]-> Write();
        momYieldGen[t] -> Write();
        thVSphiGen[t] -> Write();
        thVSmomGen[t] -> Write();
        phiVSmomGen[t] -> Write();
        
        phiVSwGen[t] -> Write();
        phiVSq2Gen[t] -> Write();
        wVSq2Gen[t] -> Write();
        
        for (int iQ2 = 0; iQ2 < qBinMax; iQ2++)
            phiVSw_Q2binGen[t][iQ2] -> Write();
        
        for (int iW = 0; iW < static_cast<int>(wBins_ST); iW++)
            phiVSq2_WbinGen[t][iW] -> Write();
        
    }  
    phiYieldGenALL -> Write();
    
    for (auto item : phiVSthFullRange){
        item -> Write();
    }
    
	for (int s = 0; s < 6; s++){
        for (int iCut = 0; iCut < static_cast<int>(nCutCombination); iCut++){
            thVSmom[iCut][s] -> Write(); 
            thVSphi[iCut][s] -> Write(); 
            phiVSmom[iCut][s] -> Write(); 
            
            phiVSw[iCut][s] -> Write(); 
            phiVSq2[iCut][s] -> Write(); 
            wVSq2[iCut][s] -> Write(); 
            
            for (int iQ2 = 0; iQ2 < qBinMax; iQ2++){
                phiVSw_Q2bin[iCut][s][iQ2] -> Write(); 
            }
            
            for (int iW = 0; iW < static_cast<int>(wBins_ST); iW++)
                phiVSq2_Wbin[iCut][s][iW] -> Write();
            
            
            thYield[iCut][s] -> Write(); 
            phiYield[iCut][s] -> Write(); 
            momYield[iCut][s] -> Write(); 

            for (size_t iTh = 0; iTh < thNsd; iTh++)
                pYieldTh[iCut][s][iTh] -> Write(); 

            for (size_t iP = 0; iP < pNsd; iP++){
                for (size_t iTh = 0; iTh < thNsd; iTh++){
                    wYieldPthet[iCut][s][iP][iTh] -> Write(); 
                }
            }
        }
    }
        
*/

	resultsF->cd("Results");
    
        for (int iCut = 0; iCut < nCut; iCut++){	
            for (int iSec = 0; iSec < 6; iSec++){	
                for (int iQ2 = 0; iQ2 < qBinMax; iQ2++){
                    for (int iCutLevel = 0; iCutLevel < 3; iCutLevel++){
    
                        wQ2BinSectorLogQ2AllCutsAllSys[iCut][iSec][iQ2][iCutLevel]->Write();           
                        wQ2BinSectorAccLogQ2AllCutsAllSys[iCut][iSec][iQ2][iCutLevel]->Write();
        }}}}

	for (int t = 0; t < 6; t++){
		for (int q = 0; q < qBinMax; q++){
			wQ2BinSectorGenLogQ2[t][q]->Write();
			wQ2BinSectorLogQ2StrictFidAllCuts[t][q]->Write();
			wQ2BinSectorAccLogQ2StrictFidAllCuts[t][q]->Write();
            
			for (int k = 0; k < 3; k++){	
				wQ2BinSectorLogQ2AllCutsSFSys[t][q][k]->Write();
				wQ2BinSectorAccLogQ2AllCutsSFSys[t][q][k]->Write();


			}
		}
	}
	
	
	resultsF->cd("Resolution");
	for (int iQ2 = 0; iQ2 < qBinMax; iQ2++){
		resP_WQ2[iQ2]->Write();
		resW_WQ2[iQ2]->Write();
	}




    /*
	cout << "writing 7"<<endl;

	//resultsF->cd("EID");
	*/
	resultsF->cd("DCFid");
    
    for (size_t iSec = 0 ; iSec < 6; iSec++){
    	sfCut_Sectors_bef[iSec]->Write();
    	sfCut_Sectors_aft[iSec]->Write();
        for (size_t iLayer = 0; iLayer < 3; iLayer++){
            for (size_t isCut = 0; isCut < 2; isCut++){
            
            	
                
                xyDC_noWeight[iSec][iLayer][isCut] -> Write();
                xyDC_nChi2[iSec][iLayer][isCut] -> Write();
                xyDC_ndf[iSec][iLayer][isCut] -> Write();
                if (iLayer == 0)
                    xyPCAL_noWeight[iSec][isCut] -> Write();
                    pcalFID_SF_V[isCut][iSec] -> Write();
                    pcalFID_SF_W[isCut][iSec] -> Write();
                    pcalFID_SF_U[isCut][iSec] -> Write();
            }
        }
    }
    
	cout << "writing 9"<<endl;
    resultsF->cd("SectorDependences");
    for (size_t iSec = 0 ; iSec < secN; iSec++){
        for (size_t iPbin = 0 ; iPbin < iPbin_for_PhiYield; iPbin++){
            thVSphi_PbinGen[iSec][iPbin] -> Write();
            
            for (int iCut = 0; iCut < static_cast<int>(nCutCombination); iCut++)
                thVSphi_Pbin[iCut][iSec][iPbin] -> Write();
            
            for (size_t iThetaBin = 0 ; iThetaBin < iThetaBin_for_PhiYield; iThetaBin++){
                phiYield_P_Th_Gen[iSec][iPbin][iThetaBin] -> Write();
                for (int iCut = 0; iCut < static_cast<int>(nCutCombination); iCut++){
                    phiYield_P_Th[iCut][iSec][iPbin][iThetaBin] -> Write();
                }
            }
        }
    }
    
    
    resultsF->cd("folding");
    
    wGEN_fold-> Write();
    wREC_fold-> Write();
    responseMatrix-> Write();
    
    
    responseMatrix2Dcase-> Write();
    wGEN_fold2Dcase-> Write();
    wREC_fold2Dcase-> Write();
    
    for (size_t iSec = 0 ; iSec < secN; iSec++){
        xini_2D_gen[iSec]-> Write();
        bini_2D_rec[iSec]-> Write();
        Adet_2D_genVSrec[iSec]-> Write();
        for (int q = 0; q < qBinMax; q++){
            xini_1D_gen[iSec][q]-> Write();
            bini_1D_rec[iSec][q]-> Write();
            Adet_1D_genVSrec[iSec][q]-> Write();
        }
        for (int iW = 0; iW < static_cast<int>(wBins_fold_2D); iW++){
            xini_1D_genQ2[iSec][iW]-> Write();
            bini_1D_recQ2[iSec][iW]-> Write();
            Adet_1D_genVSrecQ2[iSec][iW]-> Write();   
            
            xini_1D_genQ2_orig[iSec][iW]-> Write();
            bini_1D_recQ2_orig[iSec][iW]-> Write();
            Adet_1D_genVSrecQ2_orig[iSec][iW]-> Write();   
            
        }
        

    }
	///
	
	resultsF->cd("Unfolding");
		
		for (size_t iSec = 0; iSec < 6; iSec++){
    		hMeas[iSec]-> Write(); 
			hMeasX[iSec]-> Write(); 
			hMeasY[iSec]-> Write();
		}
		
		for (size_t iCutSys = 0; iCutSys < nCutSys; iCutSys ++){
			for (size_t iSec = 0; iSec < 6 ;iSec ++){ 
				for (size_t iLevel  = 0; iLevel < 3; iLevel ++){
					hMeasSys[iCutSys][iSec][iLevel]->Write();
				}
			}
		}
	
	for (size_t iSec  = 0; iSec < 6; iSec ++){
			const string unfName2D = "Deconv2D/S" + to_string(iSec) +"_my";
			resultsF->cd(unfName2D.c_str());
			
			measured2D_secBySec[iSec]->Write();
			measured2D_DATA[iSec]->Write();

                // НОВОЕ: гисты для Валеры
            for_Valera_Q2_VS_W[iSec]->Write();
            for_Valera_Q2_VS_W_gen[iSec]->Write();
            ////////////////////////////
                                                
			for (size_t iQ2 = 0; iQ2 < q2foldBinN ;iQ2 ++){ 
				const string unfName1D = "Deconv1D/S" + to_string(iSec) + "/Q2_" + to_string(iQ2) +"_my";
				resultsF->cd(unfName1D.c_str());
				
				measured1D_secBySec[iSec][iQ2]->Write();
				measured1D_DATA[iSec][iQ2]->Write();
			}
	}
	

	if (isData == 0){
		resultsF->cd("Unfolding");
		
		for (size_t iSec = 0; iSec < 6; iSec++){
			wGEN_testBefWeight[iSec]->Write();
			wGEN_testAftWeight[iSec]->Write();
			wGEN_testBefWeightMissed[iSec]->Write();
			wGEN_testAftWeightMissed[iSec]->Write();
  			wGEN_testAftWeight_NoRes[iSec]->Write();
  			wGEN_testAftWeight_Misak[iSec]->Write();
    		wGEN_testAftWeight_NoResMissed[iSec]->Write();
    		
		}
		
		
		for (size_t iSec  = 0; iSec < 6; iSec ++){
			const string unfName2D = "Deconv2D/S" + to_string(iSec);
			resultsF->cd(unfName2D.c_str());
			
			rooResponse2D[iSec]->Hmeasured()->Write();
			rooResponse2D[iSec]->Htruth()->Write(); 
			rooResponse2D[iSec]->Hresponse()->Write();
			rooResponse2D[iSec]->Hfakes()->Write();
			
			
			for (size_t iQ2 = 0; iQ2 < q2foldBinN ;iQ2 ++){ 
				const string unfName1D = "Deconv1D/S" + to_string(iSec) + "/Q2_" + to_string(iQ2);
				resultsF->cd(unfName1D.c_str());
				
				rooResponse1D[iSec][iQ2]->Hmeasured()->Write();
				rooResponse1D[iSec][iQ2]->Htruth()->Write(); 
				rooResponse1D[iSec][iQ2]->Hresponse()->Write();
				rooResponse1D[iSec][iQ2]->Hfakes()->Write();
				
				generated1D_DATA_sec[iSec][iQ2]->Write();
				if (iSec == 0) generated1D_DATA[iQ2]->Write();
			}
		}
		
		
	for (size_t iCutSys = 0; iCutSys < nCutSys; iCutSys ++){
		for (size_t iSec  = 0; iSec < 6; iSec ++){
			for (size_t iLevel  = 0; iLevel < 3; iLevel ++){
				string unfName = "UNF_CUT_" + to_string(iCutSys) + "_S_" + to_string(iSec)  + "_L_" + to_string(iLevel);
				resultsF->cd(unfName.c_str());
				
				hResponseRooSys[iCutSys][iSec][iLevel]->Hmeasured()->Write();
				hResponseRooSys[iCutSys][iSec][iLevel]->Htruth()->Write(); 
				hResponseRooSys[iCutSys][iSec][iLevel]->Hresponse()->Write();
				hResponseRooSys[iCutSys][iSec][iLevel]->Hfakes()->Write();
			}
		}
	}
		
		
		
		resultsF->cd("UnfoldingS1");
		
		hResponseRoo[0]->Hmeasured()->Write();
		hResponseRoo[0]->Htruth()->Write(); 
		hResponseRoo[0]->Hresponse()->Write();
		hResponseRoo[0]->Hfakes()->Write();
		
		resultsF->cd("UnfoldingS2");
		
		hResponseRoo[1]->Hmeasured()->Write();
		hResponseRoo[1]->Htruth()->Write(); 
		hResponseRoo[1]->Hresponse()->Write();
		hResponseRoo[1]->Hfakes()->Write();
		
		resultsF->cd("UnfoldingS3");
		
		hResponseRoo[2]->Hmeasured()->Write();
		hResponseRoo[2]->Htruth()->Write(); 
		hResponseRoo[2]->Hresponse()->Write();
		hResponseRoo[2]->Hfakes()->Write();
		
		resultsF->cd("UnfoldingS4");
		
		hResponseRoo[3]->Hmeasured()->Write();
		hResponseRoo[3]->Htruth()->Write(); 
		hResponseRoo[3]->Hresponse()->Write();
		hResponseRoo[3]->Hfakes()->Write();
		
		resultsF->cd("UnfoldingS5");
		
		hResponseRoo[4]->Hmeasured()->Write();
		hResponseRoo[4]->Htruth()->Write(); 
		hResponseRoo[4]->Hresponse()->Write();
		hResponseRoo[4]->Hfakes()->Write();
		
		resultsF->cd("UnfoldingS6");
		
		hResponseRoo[5]->Hmeasured()->Write();
		hResponseRoo[5]->Htruth()->Write(); 
		hResponseRoo[5]->Hresponse()->Write();
		hResponseRoo[5]->Hfakes()->Write();
	
	
	}
    cout << "writing 10"<<endl;

	resultsF->Close();
	cout << "writing 11"<<endl;
}



double kin_W(TLorentzVector ele, float Ebeam){
	TLorentzVector beam(0,0,Ebeam,Ebeam);
	TLorentzVector target(0,0,0,0.93827);
	TLorentzVector fGamma = beam - ele;
	TLorentzVector fCM = fGamma + target;
	return fCM.M();
}

double kin_Q2(TLorentzVector ele, float Ebeam){
	TLorentzVector beam(0,0,Ebeam,Ebeam);
	TLorentzVector target(0,0,0,0.93827);
	TLorentzVector fGamma = beam - ele;
	return -fGamma.M2();
}

bool DCFidXY(float X, float Y, int region, int sector, int cutLevel){
    
    
    float cutLimit = 0;
    float angle = 0.495;
    
    if (region == 1){
        cutLimit = 72;
        angle = 0.50;
       }
    if (region == 2){
        cutLimit = 114;
        angle = 0.505;
       }
    if (region == 3){
        cutLimit = 180;
	}
	
	cutLimit -= 2*(cutLevel - 1) * region * 0.6;
	
    
 /*   
    float cutLimit = 0;
    float angle = 0.55;
    
    if (region == 1){
        cutLimit = 72;
        angle = 0.6;
    }
    if (region == 2){
        cutLimit = 114;
        angle = 0.575;
    }
    if (region == 3) cutLimit = 180;
     
    cutLimit -= (cutLevel - 0.5) * region * 0.6;
    */
    return (Y >= -angle*(X + cutLimit) && Y <= angle*(X + cutLimit));

		}
		
		int PCALFid_VW(float V, float W, float U, int cutLevel){
		
			if (cutLevel < 0 || cutLevel > 2) throw out_of_range("wrong CutLevel pcal FID, VW");
		
			float cutValues[3] = {17,19,22};
			float cutValuesU[3] = {390,395,400};
			
			return ((V > cutValues[cutLevel]) && (W > cutValues[cutLevel]) && (U < cutValuesU[cutLevel]));
		
		}

		int PCALFidXY(float x, float y, int cutLevel){
			int fidCut = 0;
			float cutLimit = 0;
//			REgular cuts picture
//			if (cutLevel == 0) cutLimit = 255;
//			if (cutLevel == 1) cutLimit = 250;
//			if (cutLevel == 2) cutLimit = 245;

//Nick cut:
//			if (cutLevel == 0) cutLimit = 255;
//			if (cutLevel == 1) cutLimit = 252.5;
//			if (cutLevel == 2) cutLimit = 250;

//tighter cut:
			if (cutLevel == 0) cutLimit = 255;
			if (cutLevel == 1) cutLimit = 252.5;
			if (cutLevel == 2) cutLimit = 250;            
            
            float shift = 2.5;

			if (y > -0.5*(x + cutLimit - shift) && y < +0.5*(x + cutLimit - shift)) fidCut = 1;
			return fidCut;
		}



		float HTCCWeight(float x, float y, float HTCCEff[250][250], int stepXHTCC, int stepYHTCC, int isData){
			int htccXBin = (x + 125)/stepXHTCC;
			int htccYBin = (y + 125)/stepYHTCC;
			float eff = HTCCEff[htccXBin][htccYBin];
			if (isData == 1)
				if (eff < 0.95) eff = 0;
			if (isData == 0){
				if (eff < 0.95) eff = 0;
				else eff = 1;
			}
			if (eff > 0) eff = 1/eff;
			return eff;
		}

		float newHTCCWeight(float x, float y, float HTCCEff[250][250], int stepXHTCC, int stepYHTCC, int isData){
            
			int htccXBin = (x + 125)/stepXHTCC;
			int htccYBin = (y + 125)/stepYHTCC;
			float eff = HTCCEff[htccXBin][htccYBin];
			if (isData == 1)
				if (eff < 0.70) eff = 0;
			if (isData == 0){
				if (eff < 0.70) eff = 0;
				else eff = 1;
			}
			if (eff > 0) eff = 1/eff;
			return eff;
            
		}



		void readHTCCEff(float HTCCEff[250][250]){
			ifstream htccEffFile("../params/HTCCEfficiencyData.dat");
			for (int x = 0; x < 250; x++){
				for (int y = 0; y < 250; y++){
					htccEffFile >> HTCCEff[x][y];
				}
			}
			htccEffFile.close();
		}
		void readTriangleCut(float cutParams[6][10][2], int isData){
			ifstream fData("../params/dataTriangleCut.dat");
			ifstream fSim("../params/simTriangleCut.dat");

			
			for (int s = 0; s < 6; s++){
				for (int t = 0; t < 10; t++){
					for (int l = 0; l < 2; l++){
						if (isData == 1) fData >> cutParams[s][t][l];
						if (isData == 0) fSim >> cutParams[s][t][l];
						
					}
				}
			}
		}

		bool SFTriangleCut(float ecinE, float pcalE, float cutParams[6][10][2], int sector, int pBin, float shift){
		
			if (pBin < 5) return true;
            //for system:
            //if (abs(shift)>0.0001) return true;
                
			float zero = cutParams[sector][pBin][0];
			float slope = cutParams[sector][pBin][1];
			if (zero + slope*ecinE + shift< pcalE) return true;
			else return false;
		}

bool SameBinW(float W1, float W2, size_t nWBins){
	float lowBorderW = 0;
	float highBorderW = 5.0;
    size_t binW1 = nWBins * W1/(highBorderW - lowBorderW);
    size_t binW2 = nWBins * W2/(highBorderW - lowBorderW);
    return binW1 == binW2;
}

void fillVfromTree(float (&data)[100], const vector<double>* branch, const int icount){
    if (static_cast<int>(branch->size()) > 0)
        data[icount] = branch->at(icount);
    else
        data[icount] = 0;
}
void fillVfromTree(float (&data)[100], const vector<float>* branch, const int icount){
    if (static_cast<int>(branch->size()) > 0)
        data[icount] = branch->at(icount);
    else
        data[icount] = 0;
}

bool CutP(const float p, const int cutLevel){
    float pCut[3] = {1.5,2.,2.5};
    return p > pCut[cutLevel];
}


bool CutPCALdepos(const float pcalE, const int cutLevel){
    float pcalEcut[3] = {0.06, 0.07, 0.08}; 
    return pcalE > pcalEcut[cutLevel];
}


bool CutVz(const float _vz, const int cutLevel){
    float minVz[3] = {-8.5, -8, -7.5};
    float maxVz[3] = {2.5,2,1.5};

    return (_vz > minVz[cutLevel] && _vz < maxVz[cutLevel]);
}

bool CutDCfid(const DCXY& dc, int sec, const int cutLevel){
    return (
         DCFidXY(dc.r1X, dc.r1Y, 1, sec + 1, cutLevel) && 
         DCFidXY(dc.r2X, dc.r2Y, 2, sec + 1, cutLevel) &&
         DCFidXY(dc.r3X, dc.r3Y, 3, sec + 1, cutLevel)
         );
}

bool BadElementKnockOut(double Hx_pcal, double Hy_pcal, double Hx_ecin, double Hy_ecin, double Hx_ecout, double Hy_ecout, int sector, int cutLevel){
    
    double widthChange = 0.25;
    if (cutLevel == 0)  widthChange = -0.25;
    if (cutLevel == 2)  widthChange = 0.75;
    
    
    //MANUAL LINES NOT FROM CCDB, I LOOKED DATA/SIM Plots and found abnormal behaviuor
    
    float pcalManualParams[6][2] = {{107.2766, -10602.9779}, //S2
                                    {98.9667, -10262.0167}, //S2
                                    
        
                                    {98.0644, 5825.4023}, //S5
                                    {99.9337, 5098.3456}, //S5
                                    
                                     //S6
                                    {0.4547, -275.9317},
                                    {0.4547, -285.9317}}; //S6
    
    if (sector == 6){
        double k = tan(-30.6*3.1415/180);
        double b = -185;
        // wire swap remove
        if ((isOutOfLines(Hx_pcal, Hy_pcal, {k, b + widthChange} , {k, b - widthChange - 2}) &&
               isOutOfLines(Hx_pcal, Hy_pcal, {k, b + widthChange - 8.3} , {k, b - widthChange - 8.3 - 2.2})
               ) == false ) return false;
        
        
        float k1 = pcalManualParams[4][0];
        float b1 = pcalManualParams[4][1];
        float k2 = pcalManualParams[5][0];
        float b2 = pcalManualParams[5][1];
        
        // manual line
        if (isBetweenOfLines(Hx_pcal, Hy_pcal,  {k1, b1 + widthChange}, {k2, b2 - widthChange})) return false;
        
        return true;
    }
    
    
    //PCAL
    //Sec 2, Lay=1, Comp = 24
    //Sec 2, Lay=2,Comp = 43
    //Sec 3, Lay=2,Comp = 6
    //Sec 3, Lay=3, Comp = 41
    //Sec 4, Lay=1, Comp = 60
    //Sec 4, Lay=2, Comp = 20

    float linesParams[8][2] = {{-0.5774, 156.9137}, // REMOVED MAY 2023 Sec 2, Lay=1, Comp = 24
                                {-0.5878, 152.4101}, // REMOVED MAY 2023 Sec 2,
                                {0.5897, 120.7937}, // KEEP MAY 2023 
                                {0.5913, 114.3872}, // KEEP MAY 2023
                                
                                //{-82.3095, -24692.8538}, // -302.38 - -313.71
                                //{114.2353, 36036.7541}, //
                                {-0.5934, 129.9189},// REMOVED MAY 2023 Sec 3, Lay=3, Comp = 41
                                {-0.5899, 125.2137},// REMOVED MAY 2023
                                
                                //{1484.1452, 455714.6153}, // -306.95 - 315.36
                                //{-2214.6082, -698249.0095}, // 
                                {0.7971, 288.6657},// REMOVED MAY 2023 Sec 4, Lay=2, Comp = 20
                                {0.7848, 277.9871}};// REMOVED MAY 2023 Sec 4, Lay=2, Comp = 20
    

    
    //ECIN
    // Sec 1, Lay=5,Comp = 8
    // Sec 1, Lay=6,Comp = 19
    // Sec 1, Lay=6,Comp = 21
    // Sec 4, Lay=5,Comp = 13
    // Sec 4, Lay=5,Comp = 15
    // Sec 5, Lay=5,Comp = 16
    
    /*
    float ecinParams[12][2] = {{-0.5782, 345.7499},
                                {-0.6039, 337.3586},
                                {0.5679, -212.2414},
                                {0.5997, -231.8112},
                                {0.573, -193.1586},
                                {0.5765, -205.2132},
                                {-0.5599, -273.4762},
                                {-0.5593, -283.4213},
                                {-0.5667, -254.0252},
                                {-0.5804, -269.582},
                                {0.5901, -246.3225},
                                {0.6137, -254.2039}};
    */
    // ECOUT
    // Sec 1, Lay=7,Comp = 32
    // Sec 1, Lay=9,Comp = 19
    // Sec 1, Lay=9,Comp = 21
    // Sec 1, Lay=9,Comp = 24
    // Sec 5, Lay=7,Comp = 18
    
    float ecoutParams[8][2] = {//{351.542, -129188.6636}, // replace by 357.14 - 367.99
                                //{-1656.4878, 591767.5446}, // reverse
                                {0.5765, -222.5612},
                                {0.5618, -233.3443},
                                {0.572, -199.4151},
                                {0.5688, -210.0406},
                                {0.5684, -165.1428},
                                {0.5669, -175.8632},
                                {-0.5841, -252.1105}, // KEEP
                                {-0.5775, -263.2072}}; // KEEP
                                
    
    // SECTOR 1
    if (sector == 1){
        
        double k = tan(29.5*3.1415/180);
        double b = -92;
        if ( (isOutOfLines(Hx_pcal, Hy_pcal, {k, b + widthChange} , {k, b - widthChange - 2.4}) &&
                isOutOfLines(Hx_pcal, Hy_pcal, {k, b + widthChange - 9.1} , {k, b - widthChange - 9.1 - 2.4}) &&
                isOutOfLines(Hx_pcal, Hy_pcal, {k, b + widthChange - 127} , {k, b - widthChange - 127 - 2.4}) &&
                isOutOfLines(Hx_pcal, Hy_pcal, {k, b + widthChange - 127 - 8} , {k, b - widthChange -127 - 8 - 2.4}) 
               ) == false) return false;
        
        
        // ALL ECIN REMOVED ON MAY 2023, I do not see those structures in data
        // ALL ECOUT REMOVED FROM SEC 1 ON MAY 2023, I do not see those structures in data
        /*
        float k1 = ecinParams[0][0];
        float b1 = ecinParams[0][1];
        
        float k2 = ecinParams[1][0];
        float b2 = ecinParams[1][1];
        
        if (isBetweenOfLines(Hx_ecin, Hy_ecin,  {k1, b1 + widthChange}, {k2, b2 - widthChange})) return false;
        
        k1 = ecinParams[2][0];
        b1 = ecinParams[2][1];
        k2 = ecinParams[3][0];
        b2 = ecinParams[3][1];
        
        if (isBetweenOfLines(Hx_ecin, Hy_ecin,  {k1, b1 + widthChange}, {k2, b2 - widthChange})) return false;
        
        k1 = ecinParams[4][0];
        b1 = ecinParams[4][1];
        k2 = ecinParams[5][0];
        b2 = ecinParams[5][1];
        
        if (isBetweenOfLines(Hx_ecin, Hy_ecin,  {k1, b1 + widthChange}, {k2, b2 - widthChange})) return false;
        
        
        if (Hx_ecout > 357.14 && Hx_ecout < 367.99) return false;
        
        k1 = ecoutParams[0][0];
        b1 = ecoutParams[0][1];
        k2 = ecoutParams[1][0];
        b2 = ecoutParams[1][1];

        if (isBetweenOfLines(Hx_ecout, Hy_ecout,  {k1, b1 + widthChange}, {k2, b2 - widthChange})) return false;
        
        k1 = ecoutParams[2][0];
        b1 = ecoutParams[2][1];
        k2 = ecoutParams[3][0];
        b2 = ecoutParams[3][1];
        
        if (isBetweenOfLines(Hx_ecout, Hy_ecout,  {k1, b1 + widthChange}, {k2, b2 - widthChange})) return false;
        
        k1 = ecoutParams[4][0];
        b1 = ecoutParams[4][1];
        k2 = ecoutParams[5][0];
        b2 = ecoutParams[5][1];
        
        if (isBetweenOfLines(Hx_ecout, Hy_ecout,  {k1, b1 + widthChange}, {k2, b2 - widthChange})) return false;
        */
        
        return true;
        
    }

    
    
    // SECTOR 2
    if (sector == 2){
    
    /*
    //REMOVED MAY 2023
        float k1 = linesParams[0][0];
        float b1 = linesParams[0][1];
        
        float k2 = linesParams[1][0];
        float b2 = linesParams[1][1];
        
        if (isBetweenOfLines(Hx_pcal, Hy_pcal,  {k1, b1 + widthChange}, {k2, b2 - widthChange})) return false;
      */  
        float k1 = linesParams[2][0];
        float b1 = linesParams[2][1];
        float k2 = linesParams[3][0];
        float b2 = linesParams[3][1];
        
        if (isBetweenOfLines(Hx_pcal, Hy_pcal,  {k1, b1 + widthChange}, {k2, b2 - widthChange})) return false;
        
        k1 = pcalManualParams[0][0];
        b1 = pcalManualParams[0][1];
        k2 = pcalManualParams[1][0];
        b2 = pcalManualParams[1][1];
        
        if (isBetweenOfLines(Hx_pcal, Hy_pcal,  {k1, b1 + widthChange}, {k2, b2 - widthChange})) return false;
        
        return true;
        
    }
    
    // SECTOR 3
    if (sector == 3){
    
    	if (Hx_pcal > -313.71 && Hx_pcal <  -302.38) return false;
        
        /*
        // was removed before may 2023
        float k1 = linesParams[4][0];
        float b1 = linesParams[4][1];
        
        float k2 = linesParams[5][0];
        float b2 = linesParams[5][1];
        
        if (isBetweenOfLines(Hx_pcal, Hy_pcal,  {k1, b1 + widthChange}, {k2, b2 - widthChange})) return false;
        
        
        //REMOVED MAY 2023
        float k1 = linesParams[4][0];
        float b1 = linesParams[4][1];
        float k2 = linesParams[5][0];
        float b2 = linesParams[5][1];
        
        if (isBetweenOfLines(Hx_pcal, Hy_pcal,  {k1, b1 + widthChange}, {k2, b2 - widthChange})) return false;
        */
        return true;
        
    }
    
    // SECTOR 4
    if (sector == 4){
    
        double k = tan(-29.6*3.1415/180);
        double b = -232.8;
        
        if (Hx_pcal > -127.5 && Hx_pcal < -122.5) return false;
        return (isOutOfLines(Hx_pcal, Hy_pcal, {k, b + widthChange} , {k, b - widthChange - 3.5}));
    
    
		//REMOVED MAY 2023, did not work amyway, and not needed
        //if (Hx_pcal > -315.36 && Hx_pcal <  -306.95) return false;
        // KEEP MAY 2023
        //if (Hx_pcal > -127.50 && Hx_pcal < -122.5) return false;
        
        
        /*
        //REMOVED MAY 2023,
        float k1 = linesParams[6][0];
        float b1 = linesParams[6][1];
        float k2 = linesParams[7][0];
        float b2 = linesParams[7][1];
        
        
        if (isBetweenOfLines(Hx_pcal, Hy_pcal,  {k1, b1 + widthChange}, {k2, b2 - widthChange})) return false;
        
        //ALL ECIN REMOVED ON MAY 2023
        k1 = ecinParams[6][0];
        b1 = ecinParams[6][1];
        k2 = ecinParams[7][0];
        b2 = ecinParams[7][1];
        
        if (isBetweenOfLines(Hx_ecin, Hy_ecin,  {k1, b1 + widthChange}, {k2, b2 - widthChange})) return false;
        
        k1 = ecinParams[8][0];
        b1 = ecinParams[8][1];
        k2 = ecinParams[9][0];
        b2 = ecinParams[9][1];
        
        if (isBetweenOfLines(Hx_ecin, Hy_ecin,  {k1, b1 + widthChange}, {k2, b2 - widthChange})) return false;
        */
        
        return true;
        
    }
    
    // SECTOR 5
    if (sector == 5){
    	//ALL ECIN REMOVED ON MAY 2023
    	/*
        float k1 = ecinParams[10][0];
        float b1 = ecinParams[10][1];
        
        float k2 = ecinParams[11][0];
        float b2 = ecinParams[11][1];
        
        if (isBetweenOfLines(Hx_ecin, Hy_ecin,  {k1, b1 + widthChange}, {k2, b2 - widthChange})) return false;
        
         */
        
        float k1 = ecoutParams[6][0];
        float b1 = ecoutParams[6][1];
        float k2 = ecoutParams[7][0];
        float b2 = ecoutParams[7][1];
        
        //ecout July 2023
        if (isBetweenOfLines(Hx_ecout, Hy_ecout,  {k1, b1 + widthChange}, {k2, b2 - widthChange})) return false;
       
        k1 = pcalManualParams[2][0];
        b1 = pcalManualParams[2][1];
        k2 = pcalManualParams[3][0];
        b2 = pcalManualParams[3][1];
        
        if (isBetweenOfLines(Hx_pcal, Hy_pcal,  {k1, b1 + widthChange}, {k2, b2 - widthChange})) return false;
        
        return true;
        
    }
    
    return true;
}
        
        
bool ApplyCuts(TLorentzVector p4_electron, float (&triangleCutParams)[6][10][2], const KinemE& kin, 
               const TVector3& PCALvector3, const DCXY& dc, const int sec, const int isData, const set<cutType>& cutsToApply){
    
    int cutLevel = 1;
    
    return 
        (
        (cutsToApply.count(cutType::vzCut) == 0 || CutVz(kin.vz, cutLevel)) &&
        (cutsToApply.count(cutType::dcCut) == 0 || CutDCfid(dc,sec,cutLevel)) &&
        //(cutsToApply.count(cutType::sfCut) == 0 || SfCutValerii(kin.calEnerg / p4_electron.P(), p4_electron.P(), sec, cutLevel, isData)) &&
        (cutsToApply.count(cutType::sfCut) == 0 || SfCutValerii_Edepos(kin.calEnerg / p4_electron.P(), kin.calEnerg, sec, cutLevel, isData)) &&
        
        (cutsToApply.count(cutType::momCut) == 0 || CutP(kin.momentum, cutLevel)) &&
        //(cutsToApply.count(cutType::pcalCut) == 0 || PCALFidXY(PCALvector3.X(), PCALvector3.Y(), fidCutLevel) == 1) &&
        (cutsToApply.count(cutType::pcalCut) == 0 || PCALFid_VW(kin.pcalV, kin.pcalW, kin.pcalU, cutLevel)) &&
        
        (cutsToApply.count(cutType::triagCut) == 0 || SFTriangleCut(kin.ecinE/kin.momentum, kin.pcalE/kin.momentum, triangleCutParams, sec, kin.pBin, 0)) &&
        (cutsToApply.count(cutType::badElemCut) == 0 || BadElementKnockOut(kin.pcalHx, kin.pcalHy, kin.ecinHx, kin.ecinHy, kin.ecoutHx, kin.ecoutHy, sec + 1, cutLevel)) &&
        (cutsToApply.count(cutType::pcalDeposCut) == 0 || CutPCALdepos(kin.pcalE, cutLevel)) &&
        (cutsToApply.count(cutType::thetaMin) == 0 || thetaMinCut(kin.theta, cutLevel)) &&
        (cutsToApply.count(cutType::thetaMax) == 0 || thetaMaxCut(kin.theta, cutLevel)) &&
        (cutsToApply.count(cutType::thetaVSphiDC) == 0 || thetaVSphiDCcut(kin.theta, kin.phi, cutLevel)) &&
        (cutsToApply.count(cutType::phiMax) == 0 || phiMaxCut(kin.phi, cutLevel)) &&
        (cutsToApply.count(cutType::phiMin) == 0 || phiMinCut(kin.phi, cutLevel)) &&
        (cutsToApply.count(cutType::phiSpike) == 0 || phiSpikeCut(kin.phi_notCorrected, kin.theta, sec, cutLevel))
        
        
        );
}
        

void smear(TLorentzVector *V4, int q){
        //smearing degree
        float coeff = 1;
        ///
    
        double inM = V4->M();
	//true generated values
        double sP  = V4->P();
        double sTh = V4->Theta();
        double sPh = V4->Phi();

	//calculate resolutions
        double sThD = TMath::RadToDeg()*sTh;
        double momS1 = 0.0184291 -0.0110083*sThD + 0.00227667*sThD*sThD -0.000140152*sThD*sThD*sThD + 3.07424e-06*sThD*sThD*sThD*sThD;
        double momS2 = 0.02*sThD;
        double momR  = 0.01 * TMath::Sqrt( TMath::Power(momS1*sP,2) + TMath::Power(momS2,2) );
        momR *= (2.0 * coeff);

        double theS1 = 0.004*sThD + 0.1;
        double theS2 = 0;
        double theR  = TMath::Sqrt(TMath::Power(theS1*TMath::Sqrt(sP*sP+0.13957*0.13957)/(sP*sP),2) + TMath::Power(theS2,2) );
        theR *= (2.5 * coeff);

        double phiS1 = 0.85-0.015*sThD;
        double phiS2 = 0.17-0.003*sThD;
        double phiR  = TMath::Sqrt(TMath::Power(phiS1*TMath::Sqrt(sP*sP+0.13957*0.13957)/(sP*sP),2) + TMath::Power(phiS2,2) );
        phiR *= (3.5 * coeff);

	//overwrite EB
        sPh += TMath::DegToRad() * phiR * myMC->Gaus();
        sTh += TMath::DegToRad() * theR * myMC->Gaus();
        sP  += momR  * myMC->Gaus() *  V4->P() ;
	// EB_rec_mom = GEN_mom + resolution_momentum x gaussian x GEN_mom
	// EB_rec_ang = GEN_ang + resolution_angle x gaussian 

        V4->SetE( TMath::Sqrt( sP*sP + inM*inM )  );
        V4->SetRho( sP );
        V4->SetTheta( sTh );
        V4->SetPhi( sPh );
}

bool SfCutValerii_Edepos(const double sf,const double Edep,const int sec,const int cutLevel,const int isData){
    double SigmaRange = 3.5;
    double sigmaRangeChange = 0.5;
    
    if (cutLevel == 0) SigmaRange += sigmaRangeChange;
        
    if (cutLevel == 2) SigmaRange -= sigmaRangeChange;    
    
    double meanAll[3][7] = {{0.28617, 0.27974, 0.27476, 0.27264, 0.27074, 0.27575, 0.29045},
    						{-0.04047, -0.03786, -0.03409, -0.03304, -0.03211, -0.03394, -0.03974},
    						{-0.00296, -0.00118, -0.00141, -0.00067, 0.00051, -0.00143, -0.00286}};
    						
    double sigmaAll[3][7] = {{0.0172, 0.01878, 0.017, 0.01566, 0.01586, 0.01684, 0.01485},
    						 {-0.00122, -0.00296, -0.00224, 0.00028, -0.00135, -0.00185, -0.00053},
    						 {-0.0012, -0.00135, -0.00129, -0.00131, -0.00092, -0.00116, -0.00137 }};
    
    
    if (isData){

        
        double mean = meanAll[0][sec] + meanAll[1][sec] / Edep + meanAll[2][sec]/ ( Edep * Edep);
        double sigma = sigmaAll[0][sec] + sigmaAll[1][sec] / Edep + sigmaAll[2][sec] / ( Edep * Edep);
        
        double lowSFcut = mean - sigma * SigmaRange;
        
        return sf > lowSFcut;
    }
    else{

        double mean = meanAll[0][6]  + meanAll[1][6]  / Edep + meanAll[2][6] / ( Edep * Edep);
        double sigma = sigmaAll[0][6] + sigmaAll[1][6]  / Edep + sigmaAll[2][6]  / ( Edep * Edep);
        
        double lowSFcut = mean - sigma * SigmaRange;
        
        return sf > lowSFcut;
    }
    return 0;
}


bool SfCutValerii(const double sf,const double p,const int sec,const int cutLevel,const int isData){
    double SigmaRange = 3.5;
    double sigmaRangeChange = 0.5;
    
    if (cutLevel == 0) SigmaRange += sigmaRangeChange;
        
    if (cutLevel == 2) SigmaRange -= sigmaRangeChange;    
    
    if (isData){
        double meanData[6][3] = {{0.2585689824186437 ,  -0.019148913279231475 ,  -2.445738318524526e-05},
                                 {0.2585194202814996 ,  -0.019183919510484277 ,  -2.0767047278959854e-05},
                                 {0.25927775279658183 ,  -0.022744847949129454 ,  -3.0330696118736726e-05},
                                 {0.25853685002945714 ,  -0.018813335671120535 ,  -2.3681492564347605e-05},
                                 {0.25901666933641815 ,  -0.02129640869658156 ,  -2.8957129490435763e-05},
                                 {0.2584975622013138 ,  -0.018814414894893453 ,  -2.4677033027390794e-05}};
        
        double mean = meanData[sec][0] + meanData[sec][1] / p + meanData[sec][2] * p * p;
        
        double sigmaData[6][3] = {{0.0066405375788197996 ,  0.023908728433871818 ,  -6.138313458231744e-06},
                                  {0.006580254602726699 ,  0.024062242116243006 ,  -6.151014797249325e-06},
                                  {0.006542354326952354 ,  0.024212245324333637 ,  -5.9866080894261774e-06},
                                  {0.00655152117631963 ,  0.02414747030066422 ,  -6.048535470838711e-06},
                                  {0.0066180312283823095 ,  0.02411984106712406 ,  -6.951079063770554e-06},
                                  {0.006817841778346441 ,  0.023245092249594285 ,  -7.528693712532906e-06}};
        
        
        double sigma = sigmaData[sec][0] + sigmaData[sec][1] / p + sigmaData[sec][2] * p * p;
        
        double lowSFcut = mean - sigma * SigmaRange;
        
        return sf > lowSFcut;
    }
    else{
        double meanSim[3] = {0.26610015434872497 ,  -0.039307652943123106 ,  -0.0002455737726300279};
        double sigmaSim[3] = {0.013087032567661805 ,  0.007742099148408121 ,  -4.589415842038303e-05};
        
        double mean = meanSim[0] + meanSim[1] / p + meanSim[2] * p * p;
        double sigma = sigmaSim[0] + sigmaSim[1] / p + sigmaSim[2] * p * p;
        
        double lowSFcut = mean - sigma * SigmaRange;
        
        return sf > lowSFcut;
    }
    return 0;
}
