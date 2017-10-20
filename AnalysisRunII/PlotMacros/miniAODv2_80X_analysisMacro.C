/* to do -- add in signal curves with multiplyier option -- formatting needed */
/* to do -- add in option to not draw data on signal region plots */


///////////////////////////////////////////////////////////
// root -l
// .L miniAODv2_80X_analysisMacro.C+
// miniAODv2_80X_analysisMacro()
////////////////////////////////////////////////////////////

#include <iostream>
#include <map>
#include "TChain.h"
#include "THStack.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TVector2.h"
#include "TMath.h"
#include <fstream>

#include "TLine.h"
#include "TH1F.h"
#include "TCut.h"
#include <iomanip>
#include "TFile.h"
#include <stdio.h>
#include "TApplication.h"

#include "tdrstyle.C"
#include "loadTChains.C"

///////////////////////////////////////////////////////////
// global vars. & function decs. 
///////////////////////////////////////////////////////////

/* channel bools */
bool do_mTau = 0;
bool do_eTau = 0;
bool do_eMu = 0;
bool do_TauTau = 0;
bool do_All = 0;

bool printSyst = 0;

/* interactive choices*/

int choice = -1;
int sigChoice = -1;
int nonMVAChoice = -1;
bool doControlTMVA;
bool doSyst = 0;
int choice2 = -1;
int choice3 = -1;
int choice4 = -1;
int allChanBothDR = -1;

std::string metCutVal = "";
std::string mtCutVal = "";
bool cRegion = 0;
bool signalRegionMETonly = 0;
bool evalBaselineSys = 0;
bool evalEnergyShiftSys = 0;
int drCut = 0;

// BKG error init

double vvtError = 0.;
double vvjError = 0.;
double wError = 0.;
double zttError = 0.;
double zlError = 0.;
double zjError = 0.;
double zvvError = 0.;
double ewkError = 0.;
double tttError = 0.;
double ttjError = 0.;
double smhError = 0.;
double qcdError = 0.;

/* DATA maximum factor on hists y-axis */

float DATA_MAX = 1.7;
float DATA_MIN = 1e-1;
float LOG_MAX = 1000.;

/* set log y */
bool setLog  = 0;

/* plot signals */

bool shouldIplotSignals = 1;
bool splitPlotSMH = 0;

/* for safety blind met (in data) above this value, and m_vis between these */
double met_blind = 200;
std::string control_met_cutoff = "105.";

/* plotting region key requirement -- usually a cut on the met!
   for control region data usually choose met < 100
   if no cut is applied, make sure to set signalRegionMETonly to 1
   (no other dists. will be shown in that case)*/

std::map <int,std::string> drCutMap;
std::map <int,std::string> drCutNameMap;
std::map <std::string,std::string> met_options;
std::map <std::string,std::string> mt_options;
std::map <std::string,TCut> jetCutMap;
std::map <std::string,TCut> weightMap;
std::map <std::string,TCut> wNorm_cut_options;
std::map <std::string,TCut> TMVACutMap;
std::map <std::string,TCut> wControl_cut_options;

double lumi_sf = 1.0;

/* data marker style */

int dataMarker = 8;

/* some hist binnings */

float mtBinning[3] = {20.,0.,250.}; /* for mt in all regions */
float mt_totBinning[3] = {60.,0.,1500.};
float metBinning[3] = {50.,0.,1000.};  /* for met in all regions */
float m_visBinning[3] = {25.,0.,125.}; /* for m_vis in all regions */
float m_visConBinning[3] = {26.,70.,200.};
float mBinning[3] = {30,0.,3.};
float ptBinning[3] = {20.,0.,200.}; /* for pT in all regions */
float nbtagBinning[3] = {8.,-0.5,7.5}; /* for nbtag in all regions */
float tmvaBinning[3] = {200,0.0,1.0};
float phiBinning[3] = {15,-3.2,3.2};
float drBinning[3] = {25,0.,2.5};
float etaBinning[3] = {15,-3.,3.};
float pt_ttBinning[3] = {20,0.,200.};
float mt_totWconBinning[3] = {6,0.,300.};
float lptBinning[3] = {35.,0.,3500.};
float m_minBinning[3] = {20,0.,1000.};
float p_chiBinning[3] = {30,-150.,150.};
float mt_totTestBinning[3] = {2000,0.,2000.};

/* some hist x-axis labels */
std::map<std::string,std::string > x_axisLabels;

/* stack plot colors */

std::map<std::string,int> colors;

/* channel label */

std::string chan_label;

/* canvas dim */

int canDim[2] = {900,750};

/* lumi and sqrt(s) */

std::string lumiAndRtS = "[35.870/fb @ 13 TeV]";

/* sample weights (sample, weight) */

std::map<std::string,std::string> weights;

/* QCD OS to SS scale factor -- to be measured */

double qcdOStoSS = 1.00;


//txt file with sys
std::string saveNameValues = "YieldsAug16.output.txt";
ofstream out_data(saveNameValues.c_str());

/* string to hold main title for hists */

std::string global_title;

/* string for canvas save prefix */

std::string saveName;

/* given a TCut obtain a [data-(non-W)]/W scale factor */
/* args are cut, canvas name */

double shapeScaleVarBin(TCut, TCut, TChain *, std::string, int, float[]);
double shapeScale(TCut, TCut, TChain *, std::string, std::string, float[3]);
double wjetsNorm(TCut, TCut, TCut, TCut, TCut, TCut, bool, std::string, int, float[], std::string, std::string, std::string, bool);

/* given a TCut and a w mc scale factor fill a TH1F with a qcd shape under the TCut */
/* args are cut, w sf, qcd hist to fill, string expression to fill (pt_1 m_vis etc.), binning array,
canvas name string */

void fillQCD_Shape(TCut, TCut, double, TH1F *, std::string, float[3], std::string, std::string, std::string, bool, bool);
void fillQCD_ShapeVarBin(TCut, TCut, double, TCut, bool, TH1F *, std::string, int, float[], std::string, std::string, std::string, bool, bool);

/* given a TCut and a w mc scale factor and a qcd shape draw string under float[3] binning */
/* args are cut, w sf, qcd hist to draw, string expression to fill (pt_1 m_vis etc.), 
binning array, canvas name */

void drawSignalRegion(TCut, TCut, TCut, TCut, double, double, std::string, float[3], std::string, std::string, std::string, bool, bool);

void drawSignalRegionVarBin(TCut, TCut, TCut, TCut, double, TCut, TCut, bool, double, std::string, int, float[], std::string, std::string, std::string, bool, bool);

void countTotalsVarBin(TCut, TCut, TCut, TCut, double, double, std::string, int, float[], std::string, std::string, std::string);
/* determine BKG fraction distribution */

void findBkgFractions(TCut, double, std::string, float[3]);

/* find systematic for baseline variations */

void evalSysBaseline(std::vector<std::vector<double>>, TCut, TCut, TCut, TCut, TCut, TCut, TCut, TCut, std::string, float[3], std::string);
//for  channels with negligable W contribution
void evalSysBaselineNoW(std::vector<std::vector<double>>, TCut, TCut, TCut, TCut, std::string, float[3], std::string);

/* run channels */

void doMuTau();
void doEleTau();
void doEleMu();
void doTauTau();

/* draw a standard TCanvas and TLatex for this analysis */
/* as well as a ratio plot (canvas must have 2 dirs) */
void setBkgErr(std::string);
void drawTcanvasAndLegendAndRatio(TCanvas *, TH1F*, TH1F*, bool); /* data, sumBkg ,  bool drawSignals?*/

///////////////////////////////////////////////////////////
// main function
///////////////////////////////////////////////////////////

int miniAODv2_80X_analysisMacro()
{

	while(choice <1 || choice > 5)
	{
		std::cout<<" Please enter a channel choice \n";
		std::cout<<" 			enter 1 for muTau \n";
		std::cout<<" 			enter 2 for eTau \n";
		std::cout<<" 			enter 3 for eMu \n";
		std::cout<<" 			enter 4 for TauTau \n";
        std::cout<<" 			enter 5 for all channels \n";
	    cin>>choice;
	}

	/* string for canvas saving */
	saveName = "";

	if(choice==1) {do_mTau = 1; saveName += "MuTau_"; }
	if(choice==2) {do_eTau = 1; saveName += "EleTau_"; }
	if(choice==3) {do_eMu = 1; saveName += "EleMu_"; }
	if(choice==4) {do_TauTau = 1; saveName += "TauTau_"; }
    if(choice==5) {do_All = 1; saveName += "AllChanRun_"; }

	while(choice2!=1 && choice2!=2 && choice2!=3 && choice2!=4)
	{
		std::cout<<" Choose a plot set : \n";
		std::cout<<" 			enter 1 signal region dist. (above a cut) with blind tail \n";
		std::cout<<" 			enter 2 for many basic kinematic plots in control region \n";
		cin>>choice2;
	}
    
	if(choice2==1)
	{
        std::cout<< "Choose a parameter: \n";
        std::cout<<" 			enter 1 for MET dis. with blinded tail \n";
		std::cout<<" 			enter 2 for TMVA variable blinded above 0.5 (for EleTau and MuTau ONLY, EleMu and TauTau will remain in MET)\n";
        std::cout<<" 			enter 3 for Mt dis. with blinded tail\n";
        
        cin>>sigChoice;
        
        if (sigChoice==1)
        {
            metBinning[0] = 50;
            metBinning[1] = 0;
            metBinning[2] = 1000;
            std::cout<<" enter a met cut value: ";
            cin >> metCutVal;
            std::cout<<" enter a Mt cut value: ";
            cin >> mtCutVal;
            saveName += (("signalRegion_met_over_" + metCutVal + "_").c_str());
        }
        else if (sigChoice==2 )
        {
            
            std::cout<<" enter a met cut value: ";
            cin >> metCutVal;
            std::cout<<" enter a Mt_vis cut value: ";
            cin >> mtCutVal;
            if (choice==5 || choice==3 )
            {
            std::cout<<" enter a non-MVA variable: 1 for MET, 2 for Mt";
            cin>>nonMVAChoice;
            saveName += (("signalRegion_METChannelsOver" + metCutVal + "_TMVA").c_str());
            }
        }
        else if (sigChoice==3)
        {
            std::cout<<" enter a met cut value: ";
            cin >> metCutVal;
            std::cout<<" enter a max m_vis cut value: ";
            cin >> mtCutVal;
            saveName += (("signalRegion_Mt_over_" + mtCutVal + "_").c_str());
        }
		std::cout<<"    enter 1 to include systematics shapes \n";
        std::cout<<"    enter 0 to exclude systematics shapes \n";
        cin>>doSyst;
	}
	if(choice2==2)
	{
        cRegion = 1;
        std::cout<<" enter a Mt cut value: ";
        cin >> mtCutVal;
        if (choice==5||choice==1||choice==2||choice==4)
        {
            std::cout<< "Enter 1 for TMVA control plots: \n";
            std::cout<<" 			enter 1 for TMVA \n";
            std::cout<<" 			enter 0 for many kinematic plots \n";
            cin>>doControlTMVA;
        }
		metBinning[0] = 10;
		metBinning[1] = 100;
		metBinning[2] = 200;
        mt_totBinning[0] = 10;
		mt_totBinning[1] = 0;
		mt_totBinning[2] = 200;
		saveName += ("controlRegion_");
	}
	while(choice4!=1 && choice4!=2 && choice4!=3 && choice!=5)
	{
		std::cout<<" enter DR range : \n";
		std::cout<<"	   enter 1 for low  [ dr(tau,tau)<1 ] \n";
		std::cout<<"	   enter 2 for medium  [ dr(tau,tau)>=1 && dr(tau,tau) < 2.0] \n";
        std::cout<<"	   enter 3 for all regions above minimum\n";
		cin>>choice4;

		if(choice4==1)
        {
         drCut = 1;
        }
		if(choice4==2)
        {
         drCut = 2;
        }
        if(choice4==3)
        {
         drCut = 3;
        }
	}
	while(choice3!=1 && choice3!=2)
	{
		std::cout<<" Choose a plot style : \n";
		std::cout<<" 			enter 1 logarithmic \n";
		std::cout<<" 			enter 2 normal \n";
		cin>>choice3;
	}
	if(choice3==1) { setLog = 1; saveName += "log_"; }
	else { setLog = 0; saveName += "norm_";}

	if(choice2==1) {DATA_MIN = 1e-2; }

	/////////////////

	/* some formatting */
	setTDRStyle();
	/* int colors */
	colors["ZTT"] = 796;
    colors["ZL"] = 851;
    colors["ZJ"] = 851;
    colors["TT"] = 592;
	colors["TTT"] = 592;
    colors["TTJ"] = 592;
    colors["VV"] = 618;
	colors["W"] = 634;
	colors["QCD"] = 606;
    colors["EWK"] = 626;
    colors["ZVV"] = 623;
    colors["GluGluHTauTau"] = 40;
    colors["VBFHTauTau"] = 48;
	colors["ZHTauTau"] = 50;
	colors["MONO"] = 4;
    colors["SMH"] = 857;

	/* x - axis labels */

    x_axisLabels["mt_tot"] = "M_{T} Total [GeV]";
	x_axisLabels["m_vis"] = "Visible di-#tau Mass [GeV]";
	x_axisLabels["nbtag"] = "Number of b-tags";
	x_axisLabels["mvamet"] = "Missing Transverse Energy [GeV]";
    x_axisLabels["met"] = "PF Missing Transverse Energy [GeV]";
    x_axisLabels["LPT"] = "Likelihood Transverse Momentum of Pair Mother [GeV]";
    x_axisLabels["DeltaR_leg1_leg2"] = "DeltaR between legs";
    x_axisLabels["DeltaPhi_leg1_leg2"] = "DeltaPhi between legs";
    x_axisLabels["DeltaPhi_PFMET_Higgs"] = "DeltaPhi between PFMET and pair";
    x_axisLabels["DeltaPhi_MVAMET_Higgs"] = "DeltaPhi between MVAMET and pair";
    x_axisLabels["njetspt20"] = "Number of Jets (pt > 20 GeV)";
    
    x_axisLabels["mvaVar_et_MZP600A0400"] = "mvaVar_et_MZP600A0400";
    x_axisLabels["mvaVar_et_MZP800A0400"] = "mvaVar_et_MZP800A0400";
    x_axisLabels["mvaVar_et_MZP1000A0400"] = "mvaVar_et_MZP1000A0400";
    x_axisLabels["mvaVar_et_MZP1200A0400"] = "mvaVar_et_MZP1200A0400";
    x_axisLabels["mvaVar_mt_MZP600A0400"] = "mvaVar_mt_MZP600A0400";
    x_axisLabels["mvaVar_mt_MZP800A0400"] = "mvaVar_mt_MZP800A0400";
    x_axisLabels["mvaVar_mt_MZP1000A0400"] = "mvaVar_mt_MZP1000A0400";
    x_axisLabels["mvaVar_mt_MZP1200A0400"] = "mvaVar_mt_MZP1200A0400";
    x_axisLabels["mvaVar_tt_MZP600A0400"] = "mvaVar_tt_MZP600A0400";
    x_axisLabels["mvaVar_tt_MZP800A0400"] = "mvaVar_tt_MZP800A0400";
    x_axisLabels["mvaVar_tt_MZP1000A0400"] = "mvaVar_tt_MZP1000A0400";
    x_axisLabels["mvaVar_tt_MZP1200A0400"] = "mvaVar_tt_MZP1200A0400";

	if(do_eTau)
	{
		x_axisLabels["mt_1"] = "Transverse Mass (e,MET) [GeV]";
        x_axisLabels["mt_2"] = "Transverse Mass (#tau_{h},MET) [GeV]";
		x_axisLabels["pt_1"] = "Transverse Momentum of Electron [GeV]";
		x_axisLabels["pt_2"] = "Transverse Momentum of #tau_{h} [GeV]";
        x_axisLabels["M_min"] = "Function M_min [GeV]";
        x_axisLabels["P_chi"] = "Function P_chi [GeV]";
        x_axisLabels["m_1"] = "electron mass Mass [GeV]";
        x_axisLabels["m_2"] = "#tau_{h} Mass [GeV]";
        x_axisLabels["phi_1"] = "electron Phi";
        x_axisLabels["phi_2"] = "#tau_{h} Phi";
        x_axisLabels["eta_1"] = "electron eta";
        x_axisLabels["eta_2"] = "#tau_{h} eta";
        
	}

	if(do_eMu) 
	{
		x_axisLabels["mt_1"] = "Transverse Mass (e,MET) [GeV]";
		x_axisLabels["pt_1"] = "Transverse Momentum of Electron [GeV]";
		x_axisLabels["pt_2"] = "Transverse Momentum of Muon [GeV]";
	}

	if(do_mTau) 
	{
		x_axisLabels["mt_1"] = "Transverse Mass (#mu,MET) [GeV]";
        x_axisLabels["mt_2"] = "Transverse Mass (#tau_{h},MET) [GeV]";
		x_axisLabels["pt_1"] = "Transverse Momentum of Muon [GeV]";
		x_axisLabels["pt_2"] = "Transverse Momentum of #tau_{h} [GeV]";
        x_axisLabels["M_min"] = "Function M_min [GeV]";
        x_axisLabels["P_chi"] = "Function P_chi [GeV]";
        x_axisLabels["m_1"] = "muon Mass [GeV]";
        x_axisLabels["m_2"] = "#tau_{h} Mass [GeV]";
        x_axisLabels["phi_1"] = "muon Phi";
        x_axisLabels["phi_2"] = "#tau_{h} Phi";
        x_axisLabels["eta_1"] = "muon eta";
        x_axisLabels["eta_2"] = "#tau_{h} eta";
	}

	if(do_TauTau) 
	{
		x_axisLabels["mt_1"] = "Transverse Mass (lead-#tau_{h},MET) [GeV]";
        x_axisLabels["mt_2"] = "Transverse Mass (leg 2 #tau_{h},MET) [GeV]";
		x_axisLabels["pt_1"] = "Transverse Momentum of lead-#tau_{h} [GeV]";
		x_axisLabels["pt_2"] = "Transverse Momentum of 2nd-#tau_{h} [GeV]";
        x_axisLabels["m_1"] = "leg 1 #tau_{h} Mass [GeV]";
        x_axisLabels["m_2"] = "leg 2 #tau_{h} Mass [GeV]";
        x_axisLabels["phi_1"] = "leg 1 #tau_{h} Phi";
        x_axisLabels["phi_2"] = "leg 2 #tau_{h} Phi";
        x_axisLabels["eta_1"] = "leg 1 #tau_{h} eta";
        x_axisLabels["eta_2"] = "leg 2 #tau_{h} eta";
	}
    
    if(do_All)
	{
		x_axisLabels["mt_1"] = "Transverse Mass [GeV]";
		x_axisLabels["pt_1"] = "Transverse Momentum of Leg 1 [GeV]";
		x_axisLabels["pt_2"] = "Transverse Momentum of Leg 2 [GeV]";
	}

    ///CutMaps///
    drCutMap[1] = "DeltaR_leg1_leg2 < 1.";
    drCutMap[2] = "DeltaR_leg1_leg2 >= 1. && DeltaR_leg1_leg2 < 2.0";
    drCutMap[3] = "DeltaR_leg1_leg2 > 0.3 && DeltaR_leg1_leg2 < 2.0";
    drCutNameMap[1] = "low";
    drCutNameMap[2] = "high";
    drCutNameMap[3] = "both";
/*
    TMVACutMap["ET_600"] = TCut("mvaVar_et_MZP600A0400 > -0.5 && met > 125.0 && m_vis < 125");
    TMVACutMap["ET_800"] =  TCut("mvaVar_et_MZP800A0400 > -0.5 && met > 125.0 && m_vis < 125");
    TMVACutMap["ET_1000"] = TCut("mvaVar_et_MZP1000A0400 > -0.5 && met > 125.0 && m_vis < 125");
    TMVACutMap["ET_1200"] = TCut("mvaVar_et_MZP1200A0400 > -0.5 && met > 125.0 && m_vis < 125");
    
    TMVACutMap["MT_600"] = TCut("mvaVar_mt_MZP600A0400 > -0.5 && met > 125.0 && m_vis < 125");
    TMVACutMap["MT_800"] = TCut("mvaVar_mt_MZP800A0400 > -0.5 && met > 125.0 && m_vis < 125");
    TMVACutMap["MT_1000"] = TCut("mvaVar_mt_MZP1000A0400 > -0.5 && met > 125.0 && m_vis < 125");
    TMVACutMap["MT_1200"] = TCut("mvaVar_mt_MZP1200A0400 > -0.5 && met > 125.0 && m_vis < 125");
    
    TMVACutMap["TT_600"] = TCut("mvaVar_tt_MZP600A0400 > -0.5 && met > 125.0 && m_vis < 125");
    TMVACutMap["TT_800"] = TCut("mvaVar_tt_MZP800A0400 > -0.5 && met > 125.0 && m_vis < 125");
    TMVACutMap["TT_1000"] = TCut("mvaVar_tt_MZP1000A0400 > -0.5 && met > 125.0 && m_vis < 125");
    TMVACutMap["TT_1200"] = TCut("mvaVar_tt_MZP1200A0400 > -0.5 && met > 125.0 && m_vis < 125");
*/

    TMVACutMap["ET_600"] = TCut("mvaVar_et_MZP600A0400 > -0.5");
    TMVACutMap["ET_800"] =  TCut("mvaVar_et_MZP800A0400 > -0.5");
    TMVACutMap["ET_1000"] = TCut("mvaVar_et_MZP1000A0400 > -0.5");
    TMVACutMap["ET_1200"] = TCut("mvaVar_et_MZP1200A0400 > -0.5");
    
    TMVACutMap["MT_600"] = TCut("mvaVar_mt_MZP600A0400 > -0.5");
    TMVACutMap["MT_800"] = TCut("mvaVar_mt_MZP800A0400 > -0.5");
    TMVACutMap["MT_1000"] = TCut("mvaVar_mt_MZP1000A0400 > -0.5");
    TMVACutMap["MT_1200"] = TCut("mvaVar_mt_MZP1200A0400 > -0.5");
    
    TMVACutMap["TT_600"] = TCut("mvaVar_tt_MZP600A0400 > -0.5");
    TMVACutMap["TT_800"] = TCut("mvaVar_tt_MZP800A0400 > -0.5");
    TMVACutMap["TT_1000"] = TCut("mvaVar_tt_MZP1000A0400 > -0.5");
    TMVACutMap["TT_1200"] = TCut("mvaVar_tt_MZP1200A0400 > -0.5");

    met_options["NOM"] = "met > " + metCutVal;
    met_options["MET_UES_UP"] = "pfmet_type1_UnclusteredEnUp_Pt > " + metCutVal;
    met_options["MET_UES_DOWN"] = "pfmet_type1_UnclusteredEnDown_Pt > " + metCutVal;
    met_options["JECUP"]= "pfmet_type1_JetEnUp_Pt > " + metCutVal;
    met_options["JECDOWN"]= "pfmet_type1_JetEnDown_Pt > " + metCutVal;
    
    met_options["TES_UP"]= "pfmet_type1_TESUp_Pt > " + metCutVal;
    met_options["TES_DOWN"]= "pfmet_type1_TESDown_Pt > " + metCutVal;
    
    met_options["TES_DM0UP"]= "pfmet_type1_dm0TESUp_Pt > " + metCutVal;
    met_options["TES_DM0DOWN"]= "pfmet_type1_dm0TESDown_Pt > " + metCutVal;
    met_options["TES_DM1UP"]= "pfmet_type1_dm1TESUp_Pt > " + metCutVal;
    met_options["TES_DM1DOWN"]= "pfmet_type1_dm1TESDown_Pt > " + metCutVal;
    met_options["TES_DM10UP"]= "pfmet_type1_dm10TESUp_Pt > " + metCutVal;
    met_options["TES_DM10DOWN"]= "pfmet_type1_dm10TESDown_Pt > " + metCutVal;

    mt_options["NOM"] = "m_vis < " + mtCutVal;
    mt_options["MET_UES_UP"] = "m_vis < " + mtCutVal;
    mt_options["MET_UES_DOWN"] = "m_vis < " + mtCutVal;
    
    mt_options["TES_UP"] = "m_vis_TESUp < " + mtCutVal;
    mt_options["TES_DOWN"] = "m_vis_TESDown < " + mtCutVal;
    
    mt_options["TES_DM0UP"] = "m_vis_dm0TESUp < " + mtCutVal;
    mt_options["TES_DM0DOWN"] = "m_vis_dm0TESDown < " + mtCutVal;
    mt_options["TES_DM1UP"] = "m_vis_dm1TESUp < " + mtCutVal;
    mt_options["TES_DM1DOWN"] = "m_vis_dm1TESDown < " + mtCutVal;
    mt_options["TES_DM10UP"] = "m_vis_dm10TESUp < " + mtCutVal;
    mt_options["TES_DM10DOWN"] = "m_vis_dm10TESDown < " + mtCutVal;
    
    weightMap["ZUP"] = TCut("ZReWeight_WeightUp");
    weightMap["ZDOWN"] = TCut("ZReWeight_WeightDown");
    weightMap["KUP"] = TCut("KReWeight_WeightUp");
    weightMap["KDOWN"] = TCut("KReWeight_WeightDown");
    weightMap["JTFUP"] = TCut("JTF_WeightUp");
    weightMap["JTFDOWN"] = TCut("JTF_WeightDown");
    weightMap["HPTUP"] = TCut("highPtTauEff_WeightUp");
    weightMap["HPTDOWN"] = TCut("highPtTauEff_WeightDown");
    weightMap["TTPTUP"] = TCut("weight_ttPtUp");
    weightMap["TTPTDOWN"] = TCut("weight_ttPtDown");
    
    weightMap["ZZUP"] = TCut("ZZReWeight_WeightUp");
    weightMap["ZZDOWN"] = TCut("ZZReWeight_WeightDown");
    
    weightMap["TOPPTUP"] = TCut("TopQuarkPtWeight_Weight");

    weightMap["WWUP"] = TCut("WWReWeight_WeightUp");
    weightMap["WWDOWN"] = TCut("(1/WWReWeight_WeightUp)");
    
    jetCutMap["NOM"] = TCut("(nbtag==0 && bcsv_1_LooseWp < .8484 && bcsv_2_LooseWp < .8484)");
    jetCutMap["BSFJECUP"] = TCut("(nbtag_JECshiftedUp==0 && bcsv_1_LooseWp_JECshiftedUp < .8484 && bcsv_2_LooseWp_JECshiftedUp < .8484)");
    jetCutMap["BSFJECDOWN"] = TCut("(nbtag_JECshiftedDown==0 && bcsv_1_LooseWp_JECshiftedDown < .8484 && bcsv_2_LooseWp_JECshiftedDown < .8484)");
    
    //If using SF for bTag
    //jetCutMap["NOM"] = TCut("BtagEventSFproduct_or_DataTag_Central");
    //jetCutMap["BSFJECUP"] = TCut("BtagEventSFproduct_or_DataTag_Central_JECshiftedUp");
    //jetCutMap["BSFJECDOWN"] = TCut("BtagEventSFproduct_or_DataTag_Central_JECshiftedUp");
    
    jetCutMap["BSFUP"] = TCut("BtagEventSFproduct_or_DataTag_Up");
    jetCutMap["BSFDOWN"] = TCut("BtagEventSFproduct_or_DataTag_Down");

    jetCutMap["BCON"] = TCut("nbtag==1");

    wControl_cut_options["NOM"] = TCut("mt_1 > 80");
    wControl_cut_options["MET_UES_UP"] = TCut("pfmet_type1_UnclusteredEnUp_MT1 > 80");
    wControl_cut_options["MET_UES_DOWN"] = TCut("pfmet_type1_UnclusteredEnDown_MT1 > 80");
    
	/* max */
	if(setLog == 1) DATA_MAX = LOG_MAX;

	/* mc sample */
    weights["DY"] = " (35.870 * final_weight) ";
    weights["ZTT"] = " (35.870 * final_weight) ";
    weights["ZL"] = " (35.870 * final_weight) ";
    weights["ZJ"] = " (35.870 * final_weight) ";
    weights["EWK"] = " (35.870 * final_weight) ";
	weights["TTT"] = " (35.870 * final_weight) ";
    weights["TTJ"] = " (35.870 * final_weight) ";
    weights["VV"] = " (35.870 * final_weight) ";
	weights["W"] = " (35.870 * final_weight) ";
	weights["QCD"] = " (35.870 * final_weight) ";
    weights["ZVV"] = " (35.870 * final_weight) ";
    weights["GluGluHTauTau"] = " (35.870 * final_weight) ";
    weights["VBFHTauTau"] = " (35.870 * final_weight) ";
	weights["ZHTauTau"] = " (35.870 * final_weight) ";
	weights["MONO"] = " (35.870 * 0.06272 * final_weight) ";

	std::cout<<" NOTE ************** mono-Higgs normalized to "<<weights["MONO"]<<"\n";
	std::cout<<" event yields shown with extra (lumi) factor of "<<lumi_sf<<"\n";

	/* do channels */

	//std::cout<<saveName<<"\n";

	if(do_mTau) doMuTau();
	if(do_eTau) doEleTau();
	if(do_eMu) doEleMu();
	if(do_TauTau) doTauTau();
    if(do_All)
    {
        std::cout << "Enter 1 to separate DR Regions" << std::endl;
        std::cout << "Enter 2 to combine" << std::endl;
        cin>>allChanBothDR;
        if (allChanBothDR==1)
        {
            drCut=1;
            std::cout << "Low Region" << std::endl;
            doMuTau();
            doEleTau();
            //doEleMu();
            doTauTau();
            
            drCut=2;
            std::cout << "High Region" << std::endl;
            doMuTau();
            doEleTau();
            //doEleMu();
            doTauTau();
        }
        else
        {
            drCut=3;
            std::cout << "Combined DR Region" << std::endl;
            doMuTau();
            doEleTau();
            //doEleMu();
            doTauTau();
        }
    }
	
	std::cout<<" press any key to exit root\n";
	std::string k_;
	cin>>k_;
	gApplication->Terminate();
    
	return 0;
}
////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////
// function imps. 
////////////////////////////////////////////////////////////

void doTauTau()
{
    setBkgErr("tt");

    //float mt_totVarBinning[] = {0,20,40,60,90,120,150,210,260,325,400,2000};
    float mt_totVarBinning[] = {0,20,40,60,90,120,150,180,210,235,260,285,325,400,2000};
    int binnum_mt_tot = sizeof(mt_totVarBinning)/sizeof(*mt_totVarBinning) - 1;

	/* set channel label */
    
    /*
    TCut SScut_TauTau_base_RelTauIso("isOsPair==0 && pairGoodForTrigger==1 && extramuon_veto==0  && extraelec_veto==0  && dilepton_veto==0 && againstElectronVLooseMVA6_1 > 0.5 && againstMuonLoose3_1 > 0.5 && againstElectronVLooseMVA6_2 > 0.5 && againstMuonLoose3_2 > 0.5 && decayModeFinding_1 > 0.5 && decayModeFinding_2 > 0.5 && byLooseIsolationMVArun2v1DBdR03oldDMwLT_1 > 0.5 && byLooseIsolationMVArun2v1DBdR03oldDMwLT_2 > 0.5");
    
    TCut SScut_TauTau_base("isOsPair==0 && pairGoodForTrigger==1 && extramuon_veto==0  && extraelec_veto==0  && dilepton_veto==0 && againstElectronVLooseMVA6_1 > 0.5 && againstMuonLoose3_1 > 0.5 && againstElectronVLooseMVA6_2 > 0.5 && againstMuonLoose3_2 > 0.5 && decayModeFinding_1 > 0.5 && decayModeFinding_2 > 0.5 && byTightIsolationMVArun2v1DBdR03oldDMwLT_1 > 0.5 && byTightIsolationMVArun2v1DBdR03oldDMwLT_2 > 0.5");
    
    TCut signalCut_base("isOsPair==1 && pairGoodForTrigger==1 && extramuon_veto==0  && extraelec_veto==0  && dilepton_veto==0 && againstElectronVLooseMVA6_1 > 0.5 && againstMuonLoose3_1 > 0.5 && againstElectronVLooseMVA6_2 > 0.5 && againstMuonLoose3_2 > 0.5 && decayModeFinding_1 > 0.5 && decayModeFinding_2 > 0.5 && byTightIsolationMVArun2v1DBdR03oldDMwLT_1 > 0.5 && byTightIsolationMVArun2v1DBdR03oldDMwLT_2 > 0.5");
    
    TCut rel_signalCut_base("isOsPair==1 && pairGoodForTrigger==1 && extramuon_veto==0  && extraelec_veto==0  && dilepton_veto==0 && againstElectronVLooseMVA6_1 > 0.5 && againstMuonLoose3_1 > 0.5 && againstElectronVLooseMVA6_2 > 0.5 && againstMuonLoose3_2 > 0.5 && decayModeFinding_1 > 0.5 && decayModeFinding_2 > 0.5 && byLooseIsolationMVArun2v1DBdR03oldDMwLT_1 > 0.5 && byLooseIsolationMVArun2v1DBdR03oldDMwLT_2 > 0.5");
    
    TCut w_Control_base("isOsPair==1 && pairGoodForTrigger==1 && extramuon_veto==0  && extraelec_veto==0  && dilepton_veto==0 && againstElectronVLooseMVA6_1 > 0.5 && againstMuonLoose3_1 > 0.5 && againstElectronVLooseMVA6_2 > 0.5 && againstMuonLoose3_2 > 0.5 && decayModeFinding_1 > 0.5 && decayModeFinding_2 > 0.5 && byLooseIsolationMVArun2v1DBdR03oldDMwLT_1 > 0.5 && byLooseIsolationMVArun2v1DBdR03oldDMwLT_2 > 0.5 && byTightIsolationMVArun2v1DBdR03oldDMwLT_1 < 0.5 && byTightIsolationMVArun2v1DBdR03oldDMwLT_2 < 0.5");
    
    TCut wSS_Control_base("isOsPair==0 && pairGoodForTrigger==1 && extramuon_veto==0  && extraelec_veto==0  && dilepton_veto==0 && againstElectronVLooseMVA6_1 > 0.5 && againstMuonLoose3_1 > 0.5 && againstElectronVLooseMVA6_2 > 0.5 && againstMuonLoose3_2 > 0.5 && decayModeFinding_1 > 0.5 && decayModeFinding_2 > 0.5 && byLooseIsolationMVArun2v1DBdR03oldDMwLT_1 > 0.5 && byLooseIsolationMVArun2v1DBdR03oldDMwLT_2 > 0.5 && byTightIsolationMVArun2v1DBdR03oldDMwLT_1 < 0.5 && byTightIsolationMVArun2v1DBdR03oldDMwLT_2 < 0.5");
    */
    
    TCut SScut_TauTau_base_RelTauIso("isOsPair==0 && pairGoodForTrigger==1 && extramuon_veto==0  && extraelec_veto==0  && dilepton_veto==0 && againstElectronVLooseMVA6_1 > 0.5 && againstMuonLoose3_1 > 0.5 && againstElectronVLooseMVA6_2 > 0.5 && againstMuonLoose3_2 > 0.5 && decayModeFinding_1 > 0.5 && decayModeFinding_2 > 0.5 && byLooseIsolationMVArun2v1DBdR03oldDMwLT_1 > 0.5 && byLooseIsolationMVArun2v1DBdR03oldDMwLT_2 > 0.5");
    
    TCut SScut_TauTau_base("isOsPair==0 && pairGoodForTrigger==1 && extramuon_veto==0  && extraelec_veto==0  && dilepton_veto==0 && againstElectronVLooseMVA6_1 > 0.5 && againstMuonLoose3_1 > 0.5 && againstElectronVLooseMVA6_2 > 0.5 && againstMuonLoose3_2 > 0.5 && decayModeFinding_1 > 0.5 && decayModeFinding_2 > 0.5 && byLooseIsolationMVArun2v1DBdR03oldDMwLT_1 > 0.5 && byLooseIsolationMVArun2v1DBdR03oldDMwLT_2 > 0.5");
    
    TCut signalCut_base("isOsPair==1 && pairGoodForTrigger==1 && extramuon_veto==0  && extraelec_veto==0  && dilepton_veto==0 && againstElectronVLooseMVA6_1 > 0.5 && againstMuonLoose3_1 > 0.5 && againstElectronVLooseMVA6_2 > 0.5 && againstMuonLoose3_2 > 0.5 && decayModeFinding_1 > 0.5 && decayModeFinding_2 > 0.5 && byLooseIsolationMVArun2v1DBdR03oldDMwLT_1 > 0.5 && byLooseIsolationMVArun2v1DBdR03oldDMwLT_2 > 0.5");
    
    TCut rel_signalCut_base("isOsPair==1 && pairGoodForTrigger==1 && extramuon_veto==0  && extraelec_veto==0  && dilepton_veto==0 && againstElectronVLooseMVA6_1 > 0.5 && againstMuonLoose3_1 > 0.5 && againstElectronVLooseMVA6_2 > 0.5 && againstMuonLoose3_2 > 0.5 && decayModeFinding_1 > 0.5 && decayModeFinding_2 > 0.5 && byLooseIsolationMVArun2v1DBdR03oldDMwLT_1 > 0.5 && byLooseIsolationMVArun2v1DBdR03oldDMwLT_2 > 0.5");
    
    TCut w_Control_base("isOsPair==1 && pairGoodForTrigger==1 && extramuon_veto==0  && extraelec_veto==0  && dilepton_veto==0 && againstElectronVLooseMVA6_1 > 0.5 && againstMuonLoose3_1 > 0.5 && againstElectronVLooseMVA6_2 > 0.5 && againstMuonLoose3_2 > 0.5 && decayModeFinding_1 > 0.5 && decayModeFinding_2 > 0.5 && byLooseIsolationMVArun2v1DBdR03oldDMwLT_1 < 0.5 && byLooseIsolationMVArun2v1DBdR03oldDMwLT_2 < 0.5");
    
    TCut wSS_Control_base("isOsPair==0 && pairGoodForTrigger==1 && extramuon_veto==0  && extraelec_veto==0  && dilepton_veto==0 && againstElectronVLooseMVA6_1 > 0.5 && againstMuonLoose3_1 > 0.5 && againstElectronVLooseMVA6_2 > 0.5 && againstMuonLoose3_2 > 0.5 && decayModeFinding_1 > 0.5 && decayModeFinding_2 > 0.5 && byLooseIsolationMVArun2v1DBdR03oldDMwLT_1 < 0.5 && byLooseIsolationMVArun2v1DBdR03oldDMwLT_2 < 0.5");
    
    std::string metCutTmp;
    if (cRegion) {metCutTmp = "met > " + control_met_cutoff;}
    else {metCutTmp = met_options["NOM"];}
    
    TCut cut_options_nom(("pt_1 > 55. && pt_2 > 40. && " + drCutMap[drCut] + " && " + metCutTmp + " && " + mt_options["NOM"] + " && pt_tt > 65. && mt_tot <2000.").c_str());
    TCut cut_weights_nom = jetCutMap["NOM"];

    TCut mt_totBlind("mt_tot < 200");
    
    //Systematic Cuts
    TCut cut_options_UESUp(("pt_1 > 55. && pt_2 > 40. && " + drCutMap[drCut] + " && " + met_options["MET_UES_UP"] + " && " + mt_options["NOM"] + " && pt_tt > 65. && mt_tot <2000.").c_str());
    TCut cut_options_UESDown(("pt_1 > 55. && pt_2 > 40. && " + drCutMap[drCut] + " && " + met_options["MET_UES_DOWN"] + " && " + mt_options["NOM"] + " && pt_tt > 65. && mt_tot <2000.").c_str());
    
    TCut cut_options_TESUp(("pt_1_TESUp > 55. && pt_2_TESUp > 40. && " + drCutMap[drCut] + " && " + met_options["TES_UP"] + " && " + mt_options["TES_UP"] + " && pt_tt_TESUp > 65. && mt_tot <2000.").c_str());
    TCut cut_options_TESDown(("pt_1_TESDown > 55. && pt_2_TESDown > 40. && " + drCutMap[drCut] + " && " + met_options["TES_DOWN"] + " && " + mt_options["TES_DOWN"] + " && pt_tt_TESDown > 65. && mt_tot <2000.").c_str());
    
    TCut cut_options_dm0TESUp(("pt_1_dm0TESUp > 55. && pt_2_dm0TESUp > 40. && " + drCutMap[drCut] + " && " + met_options["TES_DM0UP"] + " && " + mt_options["TES_DM0UP"] + " && pt_tt_dm0TESUp > 65. && mt_tot <2000.").c_str());
    TCut cut_options_dm0TESDown(("pt_1_dm0TESDown > 55. && pt_2_dm0TESDown > 40. && " + drCutMap[drCut] + " && " + met_options["TES_DM0DOWN"] + " && " + mt_options["TES_DM0DOWN"] + " && pt_tt_dm0TESDown > 65. && mt_tot <2000.").c_str());
    
    TCut cut_options_dm1TESUp(("pt_1_dm1TESUp > 55. && pt_2_dm1TESUp > 40. && " + drCutMap[drCut] + " && " + met_options["TES_DM1UP"] + " && " + mt_options["TES_DM1UP"] + " && pt_tt_dm1TESUp > 65. && mt_tot <2000.").c_str());
    TCut cut_options_dm1TESDown(("pt_1_dm1TESDown > 55. && pt_2_dm1TESDown > 40. && " + drCutMap[drCut] + " && " + met_options["TES_DM1DOWN"] + " && " + mt_options["TES_DM1DOWN"] + " && pt_tt_dm1TESDown > 65. && mt_tot <2000.").c_str());
    
    TCut cut_options_dm10TESUp(("pt_1_dm10TESUp > 55. && pt_2_dm10TESUp > 40. && " + drCutMap[drCut] + " && " + met_options["TES_DM10UP"] + " && " + mt_options["TES_DM10UP"] + " && pt_tt_dm10TESUp > 65. && mt_tot <2000.").c_str());
    TCut cut_options_dm10TESDown(("pt_1_dm10TESDown > 55. && pt_2_dm10TESDown > 40. && " + drCutMap[drCut] + " && " + met_options["TES_DM10DOWN"] + " && " + mt_options["TES_DM10DOWN"] + " && pt_tt_dm10TESDown > 65. && mt_tot <2000.").c_str());
    
    TCut cut_options_TauFakeUp(("pt_1 > 55. && pt_2 > 40. && " + drCutMap[drCut] + " && " + met_options["NOM"] + " && " + mt_options["NOM"] + " && pt_tt > 65. && mt_tot <2000.").c_str());
    TCut cut_options_TauFakeDown(("pt_1 > 55. && pt_2 > 40. && " + drCutMap[drCut] + " && " + met_options["NOM"] + " && " + mt_options["NOM"] + " && pt_tt > 65. && mt_tot <2000.").c_str());
    
    TCut cut_options_JECUp(("pt_1 > 55. && pt_2 > 40. && " + drCutMap[drCut] + " && " + met_options["JECUP"] + " && " + mt_options["NOM"] + " && pt_tt > 65. && mt_tot <2000.").c_str());
    TCut cut_options_JECDown(("pt_1 > 55. && pt_2 > 40. && " + drCutMap[drCut] + " && " + met_options["JECDOWN"] + " && " + mt_options["NOM"] + " && pt_tt > 65. && mt_tot <2000.").c_str());
    
    TCut cut_options_HPTUp(("pt_1 > 55. && pt_2 > 40. && " + drCutMap[drCut] + " && " + met_options["NOM"] + " && " + mt_options["NOM"] + " && pt_tt > 65. && mt_tot <2000.").c_str());
    TCut cut_options_HPTDown(("pt_1 > 55. && pt_2 > 40. && " + drCutMap[drCut] + " && " + met_options["NOM"] + " && " + mt_options["NOM"] + " && pt_tt > 65. && mt_tot <2000.").c_str());
    
    TCut cut_options_bTagUp(("pt_1 > 55. && pt_2 > 40. && " + drCutMap[drCut] + " && " + met_options["NOM"] + " && " + mt_options["NOM"] + " && pt_tt > 65. && mt_tot <2000.").c_str());
    TCut cut_options_bTagDown(("pt_1 > 55. && pt_2 > 40. && " + drCutMap[drCut] + " && " + met_options["NOM"] + " && " + mt_options["NOM"] + " && pt_tt > 65. && mt_tot <2000.").c_str());
    
    //Systematic cuts for WCon

    TCut cut_options_wCon(("pt_1 > 55. && pt_2 > 40. && DeltaR_leg1_leg2 > 0.3 && met > 105. && mt_1 > 80. && pt_tt > 65. && mt_tot <2000."));
    
    TCut cut_options_wCon_UESUp(("pt_1 > 55. && pt_2 > 40. && DeltaR_leg1_leg2 > 0.3 && pfmet_type1_UnclusteredEnUp_Pt > 105. && pfmet_type1_UnclusteredEnUp_MT1 > 80. && pt_tt > 65. && mt_tot <2000."));
    TCut cut_options_wCon_UESDown(("pt_1 > 55. && pt_2 > 40. && DeltaR_leg1_leg2 > 0.3 && pfmet_type1_UnclusteredEnDown_Pt > 105. && pfmet_type1_UnclusteredEnDown_MT1 > 80. && pt_tt > 65. && mt_tot <2000."));
    
    TCut cut_options_wCon_dm0TESUp(("pt_1_dm0TESUp > 55. && pt_2_dm0TESUp > 40. && DeltaR_leg1_leg2 > 0.3 && pfmet_type1_dm0TESUp_Pt > 105. && mt_1_dm0TESUp > 80. && pt_tt_dm0TESUp > 65. && mt_tot <2000."));
    TCut cut_options_wCon_dm0TESDown(("pt_1_dm0TESDown > 55. && pt_2_dm0TESDown > 40. && DeltaR_leg1_leg2 > 0.3 && pfmet_type1_dm0TESDown_Pt > 105. && mt_1_dm0TESDown > 80. && pt_tt_dm0TESDown > 65. && mt_tot <2000."));
    
    TCut cut_options_wCon_dm1TESUp(("pt_1_dm1TESUp > 55. && pt_2_dm1TESUp > 40. && DeltaR_leg1_leg2 > 0.3 && pfmet_type1_dm1TESUp_Pt > 105. && mt_1_dm1TESUp > 80. && pt_tt_dm1TESUp > 65. && mt_tot <2000."));
    TCut cut_options_wCon_dm1TESDown(("pt_1_dm1TESDown > 55. && pt_2_dm1TESDown > 40. && DeltaR_leg1_leg2 > 0.3 && pfmet_type1_dm1TESDown_Pt > 105. && mt_1_dm1TESDown > 80. && pt_tt_dm1TESDown > 65. && mt_tot <2000."));
    
    TCut cut_options_wCon_dm10TESUp(("pt_1_dm10TESUp > 55. && pt_2_dm10TESUp > 40. && DeltaR_leg1_leg2 > 0.3 && pfmet_type1_dm10TESUp_Pt > 105. && mt_1_dm10TESUp > 80. && pt_tt_dm10TESUp > 65. && mt_tot <2000."));
    TCut cut_options_wCon_dm10TESDown(("pt_1_dm10TESDown > 55. && pt_2_dm10TESDown > 40. && DeltaR_leg1_leg2 > 0.3 && pfmet_type1_dm10TESDown_Pt > 105. && mt_1_dm10TESDown > 80. && pt_tt_dm10TESDown > 65. && mt_tot <2000."));
    
    TCut cut_options_wCon_JECUp(("pt_1 > 55. && pt_2 > 40. && DeltaR_leg1_leg2 > 0.3 && pfmet_type1_JetEnUp_Pt > 105. && pfmt_1_JEnUp > 80. && pt_tt > 65. && mt_tot <2000."));
    TCut cut_options_wCon_JECDown(("pt_1 > 55. && pt_2 > 40. && DeltaR_leg1_leg2 > 0.3 && pfmet_type1_JetEnDown_Pt > 105. && pfmt_1_JEnDown > 80. && pt_tt > 65. && mt_tot <2000."));
    
    //Weights
    TCut cut_weights_UESUp = jetCutMap["NOM"];
    TCut cut_weights_UESDown = jetCutMap["NOM"];
    
    TCut cut_weights_ZPtUp = jetCutMap["NOM"]; cut_weights_ZPtUp *= weightMap["ZUP"];
    TCut cut_weights_ZPtDown = jetCutMap["NOM"]; cut_weights_ZPtDown *= weightMap["ZDOWN"];
    
    TCut cut_weights_KPtUp = jetCutMap["NOM"]; cut_weights_KPtUp *= weightMap["KUP"];
    TCut cut_weights_KPtDown = jetCutMap["NOM"]; cut_weights_KPtDown *= weightMap["KDOWN"];
    
    TCut cut_weights_TauFakeUp = jetCutMap["NOM"]; cut_weights_TauFakeUp *= weightMap["JTFUP"];
    TCut cut_weights_TauFakeDown = jetCutMap["NOM"]; cut_weights_TauFakeDown *= weightMap["JTFDOWN"];
    
    TCut cut_weights_ttPtUp = jetCutMap["NOM"]; cut_weights_ttPtUp *= weightMap["TTPTUP"];
    TCut cut_weights_ttPtDown = jetCutMap["NOM"]; cut_weights_ttPtDown *= weightMap["TTPTDOWN"];
    
    TCut cut_weights_JECUp = jetCutMap["BSFJECUP"];
    TCut cut_weights_JECDown = jetCutMap["BSFJECDOWN"];
    
    TCut cut_weights_HPTUp = jetCutMap["NOM"]; cut_weights_HPTUp *= weightMap["HPTUP"];
    TCut cut_weights_HPTDown = jetCutMap["NOM"]; cut_weights_HPTDown *= weightMap["HPTDOWN"];
    
    TCut cut_weights_bTagUp = jetCutMap["BSFUP"];
    TCut cut_weights_bTagDown = jetCutMap["BSFDOWN"];
    
    TCut cut_weights_TopPtUp = jetCutMap["NOM"]; cut_weights_TopPtUp *= weightMap["TOPPTUP"];
    TCut cut_weights_TopPtDown = jetCutMap["NOM"];
    
    TCut cut_weights_WWUp = jetCutMap["NOM"]; cut_weights_WWUp *= weightMap["WWUP"];
    TCut cut_weights_WWDown = jetCutMap["NOM"]; cut_weights_WWDown *= weightMap["WWDOWN"];
    
    TCut cut_weights_ZZUp = jetCutMap["NOM"]; cut_weights_ZZUp *= weightMap["ZZUP"];
    TCut cut_weights_ZZDown = jetCutMap["NOM"]; cut_weights_ZZDown *= weightMap["ZZDown"];
    
	chan_label = "#tau_{h} + #tau_{h}";
    
    /* setup files for TauTau channel */
    setup_files_TauTau();
	
    double wSF_SS_TauTau = 1.0; /* temp */

    //double wSF_OS_TauTau = wjetsNorm((w_Control_base + cut_options_wCon)*cut_weights_nom, (w_Control_base + cut_options_wCon)*cut_weights_nom, (wSS_Control_base + cut_options_wCon)*cut_weights_nom, (SScut_TauTau_base_RelTauIso + cut_options_wCon)*cut_weights_nom, "TauTauWjetsOS", "tt");
    
    double wSF_OS_TauTau = 1.0;
    
    if (sigChoice==1)
    {
        
        global_title = "QCD (Same Sign) Estimate Region";

        drawSignalRegion((signalCut_base + cut_options_nom)*cut_weights_nom, (rel_signalCut_base + cut_options_nom)*cut_weights_nom, (SScut_TauTau_base + cut_options_nom)*cut_weights_nom, (SScut_TauTau_base_RelTauIso + cut_options_nom)*cut_weights_nom, wSF_SS_TauTau, wSF_OS_TauTau, "met", metBinning, "sig region met", "tt", "", 1, 1);
        
        reset_files();

    }
    else if (sigChoice==2)
    {
        
        global_title = "QCD (Same Sign) Estimate Region";
        
        drawSignalRegion((signalCut_base + TMVACutMap["TT_600"])*cut_weights_nom, (signalCut_base + TMVACutMap["TT_600"])*cut_weights_nom, (SScut_TauTau_base + TMVACutMap["TT_600"])*cut_weights_nom, (SScut_TauTau_base_RelTauIso + TMVACutMap["TT_600"])*cut_weights_nom, wSF_SS_TauTau, 1.0 , "mvaVar_tt_MZP600A0400", tmvaBinning, "sig region Function TMVA variable", "tt", "600", 1, 1);
        drawSignalRegion((signalCut_base + TMVACutMap["TT_800"])*cut_weights_nom, (signalCut_base + TMVACutMap["TT_800"])*cut_weights_nom, (SScut_TauTau_base + TMVACutMap["TT_800"])*cut_weights_nom, (SScut_TauTau_base_RelTauIso + TMVACutMap["TT_800"])*cut_weights_nom, wSF_SS_TauTau, 1.0 , "mvaVar_tt_MZP800A0400", tmvaBinning, "sig region Function TMVA variable", "tt", "800", 1, 1);
        drawSignalRegion((signalCut_base + TMVACutMap["TT_1000"])*cut_weights_nom, (signalCut_base + TMVACutMap["TT_1000"])*cut_weights_nom, (SScut_TauTau_base + TMVACutMap["TT_1000"])*cut_weights_nom, (SScut_TauTau_base_RelTauIso + TMVACutMap["TT_1000"])*cut_weights_nom, wSF_SS_TauTau, 1.0 , "mvaVar_tt_MZP1000A0400", tmvaBinning, "sig region Function TMVA variable", "tt", "1000", 1, 1);
        drawSignalRegion((signalCut_base + TMVACutMap["TT_1200"])*cut_weights_nom, (signalCut_base + TMVACutMap["TT_1200"])*cut_weights_nom, (SScut_TauTau_base + TMVACutMap["TT_1200"])*cut_weights_nom, (SScut_TauTau_base_RelTauIso + TMVACutMap["TT_1200"])*cut_weights_nom, wSF_SS_TauTau, 1.0 , "mvaVar_tt_MZP1200A0400", tmvaBinning, "sig region Function TMVA variable", "tt", "1200", 1, 1);
        reset_files();
        
        if (doSyst)
        {
            std::cout << "syst for tt started" << std::endl;
            global_title = "QCD (Same Sign) Estimate Region";
        }
        reset_files();
        
    }
    else if (sigChoice==3)
    {
        global_title = "QCD (Same Sign) Estimate Region";
        
        drawSignalRegionVarBin((signalCut_base + cut_options_nom)*cut_weights_nom, (rel_signalCut_base + cut_options_nom)*cut_weights_nom, (SScut_TauTau_base + cut_options_nom)*cut_weights_nom, (SScut_TauTau_base_RelTauIso + cut_options_nom)*cut_weights_nom, wSF_SS_TauTau, (signalCut_base + cut_options_nom)*cut_weights_nom, (SScut_TauTau_base + cut_options_nom)*cut_weights_nom, 0, wSF_OS_TauTau, "mt_tot", binnum_mt_tot, mt_totVarBinning, "sig region Mt_total", "tt", "", 1, 1);
        if (printSyst) countTotalsVarBin((signalCut_base + cut_options_nom)*cut_weights_nom, (rel_signalCut_base + cut_options_nom)*cut_weights_nom, (SScut_TauTau_base + cut_options_nom)*cut_weights_nom, (SScut_TauTau_base_RelTauIso + cut_options_nom)*cut_weights_nom, wSF_SS_TauTau, wSF_OS_TauTau, "mt_tot", binnum_mt_tot, mt_totVarBinning, "control region Mt_total", "tt", "NOM");
        
        if (doSyst)
        {
            std::cout << "syst for tt started" << std::endl;
            //syst
            
            std::cout << "QCD Shapes for syst filled" << std::endl;
            
            drawSignalRegionVarBin((signalCut_base + cut_options_UESUp)*cut_weights_UESUp, (rel_signalCut_base + cut_options_UESUp)*cut_weights_UESUp, (SScut_TauTau_base + cut_options_UESUp)*cut_weights_UESUp, (SScut_TauTau_base_RelTauIso + cut_options_UESUp)*cut_weights_UESUp, wSF_SS_TauTau, (signalCut_base + cut_options_nom)*cut_weights_nom, (SScut_TauTau_base + cut_options_nom)*cut_weights_nom, 0, wSF_OS_TauTau, "mt_tot_UESUp", binnum_mt_tot, mt_totVarBinning,  "", "tt", "_CMS_scale_m_13TeVUp", 0, 1);
            drawSignalRegionVarBin((signalCut_base + cut_options_UESDown)*cut_weights_UESDown, (rel_signalCut_base + cut_options_UESDown)*cut_weights_UESDown, (SScut_TauTau_base + cut_options_UESDown)*cut_weights_UESDown, (SScut_TauTau_base_RelTauIso + cut_options_UESDown)*cut_weights_UESDown, wSF_SS_TauTau, (signalCut_base + cut_options_nom)*cut_weights_nom, (SScut_TauTau_base + cut_options_nom)*cut_weights_nom, 0, wSF_OS_TauTau, "mt_tot_UESDown", binnum_mt_tot, mt_totVarBinning,  "", "tt", "_CMS_scale_m_13TeVDown", 0, 1);
            drawSignalRegionVarBin((signalCut_base + cut_options_nom)*cut_weights_ZPtUp, (rel_signalCut_base + cut_options_nom)*cut_weights_ZPtUp, (SScut_TauTau_base + cut_options_nom)*cut_weights_ZPtUp, (SScut_TauTau_base_RelTauIso + cut_options_nom)*cut_weights_ZPtUp, wSF_SS_TauTau, (signalCut_base + cut_options_nom)*cut_weights_nom, (SScut_TauTau_base + cut_options_nom)*cut_weights_nom, 0, wSF_OS_TauTau, "mt_tot", binnum_mt_tot, mt_totVarBinning,  "", "tt", "_CMS_xtt_dyShape_13TeVUp", 0, 1);
            drawSignalRegionVarBin((signalCut_base + cut_options_nom)*cut_weights_ZPtDown, (rel_signalCut_base + cut_options_nom)*cut_weights_ZPtDown, (SScut_TauTau_base + cut_options_nom)*cut_weights_ZPtDown, (SScut_TauTau_base_RelTauIso + cut_options_nom)*cut_weights_ZPtDown, wSF_SS_TauTau, (signalCut_base + cut_options_nom)*cut_weights_nom, (SScut_TauTau_base + cut_options_nom)*cut_weights_nom, 0, wSF_OS_TauTau, "mt_tot", binnum_mt_tot, mt_totVarBinning,  "", "tt", "_CMS_xtt_dyShape_13TeVDown", 0, 1);
            drawSignalRegionVarBin((signalCut_base + cut_options_nom)*cut_weights_KPtUp, (rel_signalCut_base + cut_options_nom)*cut_weights_KPtUp, (SScut_TauTau_base + cut_options_nom)*cut_weights_KPtUp, (SScut_TauTau_base_RelTauIso + cut_options_nom)*cut_weights_KPtUp, wSF_SS_TauTau, (signalCut_base + cut_options_nom)*cut_weights_nom, (SScut_TauTau_base + cut_options_nom)*cut_weights_nom, 0, wSF_OS_TauTau, "mt_tot", binnum_mt_tot, mt_totVarBinning,  "", "tt", "_CMS_xtt_wShape_13TeVUp", 0, 1);
            drawSignalRegionVarBin((signalCut_base + cut_options_nom)*cut_weights_KPtDown, (rel_signalCut_base + cut_options_nom)*cut_weights_KPtDown, (SScut_TauTau_base + cut_options_nom)*cut_weights_KPtDown, (SScut_TauTau_base_RelTauIso + cut_options_nom)*cut_weights_KPtDown, wSF_SS_TauTau, (signalCut_base + cut_options_nom)*cut_weights_nom, (SScut_TauTau_base + cut_options_nom)*cut_weights_nom, 0, wSF_OS_TauTau, "mt_tot", binnum_mt_tot, mt_totVarBinning,  "", "tt", "_CMS_xtt_wShape_13TeVDown", 0, 1);
            drawSignalRegionVarBin((signalCut_base + cut_options_TauFakeUp)*cut_weights_TauFakeUp, (rel_signalCut_base + cut_options_TauFakeUp)*cut_weights_TauFakeUp, (SScut_TauTau_base + cut_options_TauFakeUp)*cut_weights_TauFakeUp, (SScut_TauTau_base_RelTauIso + cut_options_TauFakeUp)*cut_weights_TauFakeUp, wSF_SS_TauTau, (signalCut_base + cut_options_nom)*cut_weights_nom, (SScut_TauTau_base + cut_options_nom)*cut_weights_nom, 0, wSF_OS_TauTau, "mt_tot", binnum_mt_tot, mt_totVarBinning,  "", "tt", "_CMS_xtt_jetToTauFake_13TeVUp", 0, 1);
            drawSignalRegionVarBin((signalCut_base + cut_options_TauFakeDown)*cut_weights_TauFakeDown, (rel_signalCut_base + cut_options_TauFakeDown)*cut_weights_TauFakeDown, (SScut_TauTau_base + cut_options_TauFakeDown)*cut_weights_TauFakeDown, (SScut_TauTau_base_RelTauIso + cut_options_TauFakeDown)*cut_weights_TauFakeDown, wSF_SS_TauTau, (signalCut_base + cut_options_nom)*cut_weights_nom, (SScut_TauTau_base + cut_options_nom)*cut_weights_nom, 0, wSF_OS_TauTau, "mt_tot", binnum_mt_tot, mt_totVarBinning,  "", "tt", "_CMS_xtt_jetToTauFake_13TeVDown", 0, 1);
            drawSignalRegionVarBin((signalCut_base + cut_options_bTagUp)*cut_weights_bTagUp, (rel_signalCut_base + cut_options_bTagUp)*cut_weights_bTagUp, (SScut_TauTau_base + cut_options_bTagUp)*cut_weights_bTagUp, (SScut_TauTau_base_RelTauIso + cut_options_bTagUp)*cut_weights_bTagUp, wSF_SS_TauTau, (signalCut_base + cut_options_nom)*cut_weights_nom, (SScut_TauTau_base + cut_options_nom)*cut_weights_nom, 0, wSF_OS_TauTau, "mt_tot", binnum_mt_tot, mt_totVarBinning,  "", "tt", "_CMS_xtt_bTag_13TeVUp", 0, 1);
            drawSignalRegionVarBin((signalCut_base + cut_options_bTagDown)*cut_weights_bTagDown, (rel_signalCut_base + cut_options_bTagDown)*cut_weights_bTagDown, (SScut_TauTau_base + cut_options_bTagDown)*cut_weights_bTagDown, (SScut_TauTau_base_RelTauIso + cut_options_bTagDown)*cut_weights_bTagDown, wSF_SS_TauTau, (signalCut_base + cut_options_nom)*cut_weights_nom, (SScut_TauTau_base + cut_options_nom)*cut_weights_nom, 0, wSF_OS_TauTau, "mt_tot", binnum_mt_tot, mt_totVarBinning,  "", "tt", "_CMS_xtt_bTag_13TeVDown", 0, 1);
            drawSignalRegionVarBin((signalCut_base + cut_options_JECUp)*cut_weights_JECUp, (rel_signalCut_base + cut_options_JECUp)*cut_weights_JECUp, (SScut_TauTau_base + cut_options_JECUp)*cut_weights_JECUp, (SScut_TauTau_base_RelTauIso + cut_options_JECUp)*cut_weights_JECUp, wSF_SS_TauTau, (signalCut_base + cut_options_nom)*cut_weights_nom, (SScut_TauTau_base + cut_options_nom)*cut_weights_nom, 0, wSF_OS_TauTau, "mt_tot_JEnUp", binnum_mt_tot, mt_totVarBinning,  "", "tt", "_CMS_scale_j_13TeVUp", 0, 1);
            drawSignalRegionVarBin((signalCut_base + cut_options_JECDown)*cut_weights_JECDown, (rel_signalCut_base + cut_options_JECDown)*cut_weights_JECDown, (SScut_TauTau_base + cut_options_JECDown)*cut_weights_JECDown, (SScut_TauTau_base_RelTauIso + cut_options_JECDown)*cut_weights_JECDown, wSF_SS_TauTau, (signalCut_base + cut_options_nom)*cut_weights_nom, (SScut_TauTau_base + cut_options_nom)*cut_weights_nom, 0, wSF_OS_TauTau, "mt_tot_JEnDown", binnum_mt_tot, mt_totVarBinning,  "", "tt","_CMS_scale_j_13TeVDown", 0, 1);
            drawSignalRegionVarBin((signalCut_base + cut_options_HPTUp)*cut_weights_HPTUp, (rel_signalCut_base + cut_options_HPTUp)*cut_weights_HPTUp, (SScut_TauTau_base + cut_options_HPTUp)*cut_weights_HPTUp, (SScut_TauTau_base_RelTauIso + cut_options_HPTUp)*cut_weights_HPTUp, wSF_SS_TauTau, (signalCut_base + cut_options_nom)*cut_weights_nom, (SScut_TauTau_base + cut_options_nom)*cut_weights_nom, 0, wSF_OS_TauTau, "mt_tot", binnum_mt_tot, mt_totVarBinning,  "", "tt", "_CMS_xtt_highTauEffi_13TeVUp", 0, 1);
            drawSignalRegionVarBin((signalCut_base + cut_options_HPTDown)*cut_weights_HPTDown, (rel_signalCut_base + cut_options_HPTDown)*cut_weights_HPTDown, (SScut_TauTau_base + cut_options_HPTDown)*cut_weights_HPTDown, (SScut_TauTau_base_RelTauIso + cut_options_HPTDown)*cut_weights_HPTDown, wSF_SS_TauTau, (signalCut_base + cut_options_nom)*cut_weights_nom, (SScut_TauTau_base + cut_options_nom)*cut_weights_nom, 0, wSF_OS_TauTau, "mt_tot", binnum_mt_tot, mt_totVarBinning,  "", "tt", "_CMS_xtt_highTauEffi_13TeVDown", 0, 1);
            drawSignalRegionVarBin((signalCut_base + cut_options_nom)*cut_weights_ttPtUp, (rel_signalCut_base + cut_options_nom)*cut_weights_ttPtUp, (SScut_TauTau_base + cut_options_nom)*cut_weights_ttPtUp, (SScut_TauTau_base_RelTauIso + cut_options_nom)*cut_weights_ttPtUp, wSF_SS_TauTau, (signalCut_base + cut_options_nom)*cut_weights_nom, (SScut_TauTau_base + cut_options_nom)*cut_weights_nom, 0, wSF_OS_TauTau, "mt_tot", binnum_mt_tot, mt_totVarBinning,  "", "tt", "_CMS_xtt_tt_trigger_13TeVUp", 0, 1);
            drawSignalRegionVarBin((signalCut_base + cut_options_nom)*cut_weights_ttPtDown, (rel_signalCut_base + cut_options_HPTDown)*cut_weights_ttPtDown, (SScut_TauTau_base + cut_options_nom)*cut_weights_ttPtDown, (SScut_TauTau_base_RelTauIso + cut_options_nom)*cut_weights_ttPtDown, wSF_SS_TauTau, (signalCut_base + cut_options_nom)*cut_weights_nom, (SScut_TauTau_base + cut_options_nom)*cut_weights_nom, 0, wSF_OS_TauTau, "mt_tot", binnum_mt_tot, mt_totVarBinning,  "", "tt", "_CMS_xtt_tt_trigger_13TeVDown", 0, 1);
            drawSignalRegionVarBin((signalCut_base + cut_options_TESUp)*cut_weights_nom, (rel_signalCut_base + cut_options_TESUp)*cut_weights_nom, (SScut_TauTau_base + cut_options_TESUp)*cut_weights_nom, (SScut_TauTau_base_RelTauIso + cut_options_TESUp)*cut_weights_nom, wSF_SS_TauTau, (signalCut_base + cut_options_nom)*cut_weights_nom, (SScut_TauTau_base + cut_options_nom)*cut_weights_nom, 0, wSF_OS_TauTau, "mt_tot_TESUp", binnum_mt_tot, mt_totVarBinning,  "", "tt", "_CMS_scale_t_13TeVUp", 0, 1);
            drawSignalRegionVarBin((signalCut_base + cut_options_TESDown)*cut_weights_nom, (rel_signalCut_base + cut_options_TESUp)*cut_weights_nom, (SScut_TauTau_base + cut_options_TESUp)*cut_weights_nom, (SScut_TauTau_base_RelTauIso + cut_options_TESUp)*cut_weights_nom, wSF_SS_TauTau, (signalCut_base + cut_options_nom)*cut_weights_nom, (SScut_TauTau_base + cut_options_nom)*cut_weights_nom, 0, wSF_OS_TauTau, "mt_tot_TESDown", binnum_mt_tot, mt_totVarBinning,  "", "tt","_CMS_scale_t_13TeVDown", 0, 1);
            drawSignalRegionVarBin((signalCut_base + cut_options_dm0TESUp)*cut_weights_nom, (rel_signalCut_base + cut_options_dm0TESUp)*cut_weights_nom, (SScut_TauTau_base + cut_options_dm0TESUp)*cut_weights_nom, (SScut_TauTau_base_RelTauIso + cut_options_dm0TESUp)*cut_weights_nom, wSF_SS_TauTau, (signalCut_base + cut_options_nom)*cut_weights_nom, (SScut_TauTau_base + cut_options_nom)*cut_weights_nom, 0, wSF_OS_TauTau, "mt_tot_dm0TESUp", binnum_mt_tot, mt_totVarBinning,  "", "tt", "_CMS_scale_t_dm0_13TeVUp", 0, 1);
            drawSignalRegionVarBin((signalCut_base + cut_options_dm0TESDown)*cut_weights_nom, (rel_signalCut_base + cut_options_dm0TESUp)*cut_weights_nom, (SScut_TauTau_base + cut_options_dm0TESUp)*cut_weights_nom, (SScut_TauTau_base_RelTauIso + cut_options_dm0TESUp)*cut_weights_nom, wSF_SS_TauTau, (signalCut_base + cut_options_nom)*cut_weights_nom, (SScut_TauTau_base + cut_options_nom)*cut_weights_nom, 0, wSF_OS_TauTau, "mt_tot_dm0TESDown", binnum_mt_tot, mt_totVarBinning,  "", "tt","_CMS_scale_t_dm0_13TeVDown", 0, 1);
            drawSignalRegionVarBin((signalCut_base + cut_options_dm1TESUp)*cut_weights_nom, (rel_signalCut_base + cut_options_dm1TESUp)*cut_weights_nom, (SScut_TauTau_base + cut_options_dm1TESUp)*cut_weights_nom, (SScut_TauTau_base_RelTauIso + cut_options_dm1TESUp)*cut_weights_nom, wSF_SS_TauTau, (signalCut_base + cut_options_nom)*cut_weights_nom, (SScut_TauTau_base + cut_options_nom)*cut_weights_nom, 0, wSF_OS_TauTau, "mt_tot_dm1TESUp", binnum_mt_tot, mt_totVarBinning,  "", "tt", "_CMS_scale_t_dm1_13TeVUp", 0, 1);
            drawSignalRegionVarBin((signalCut_base + cut_options_dm1TESDown)*cut_weights_nom, (rel_signalCut_base + cut_options_dm1TESUp)*cut_weights_nom, (SScut_TauTau_base + cut_options_dm1TESUp)*cut_weights_nom, (SScut_TauTau_base_RelTauIso + cut_options_dm1TESUp)*cut_weights_nom, wSF_SS_TauTau, (signalCut_base + cut_options_nom)*cut_weights_nom, (SScut_TauTau_base + cut_options_nom)*cut_weights_nom, 0, wSF_OS_TauTau, "mt_tot_dm1TESDown", binnum_mt_tot, mt_totVarBinning,  "", "tt","_CMS_scale_t_dm1_13TeVDown", 0, 1);
            drawSignalRegionVarBin((signalCut_base + cut_options_dm10TESUp)*cut_weights_nom, (rel_signalCut_base + cut_options_dm10TESUp)*cut_weights_nom, (SScut_TauTau_base + cut_options_dm10TESUp)*cut_weights_nom, (SScut_TauTau_base_RelTauIso + cut_options_dm10TESUp)*cut_weights_nom, wSF_SS_TauTau, (signalCut_base + cut_options_nom)*cut_weights_nom, (SScut_TauTau_base + cut_options_nom)*cut_weights_nom, 0, wSF_OS_TauTau, "mt_tot_dm10TESUp", binnum_mt_tot, mt_totVarBinning,  "", "tt", "_CMS_scale_t_dm10_13TeVUp", 0, 1);
            drawSignalRegionVarBin((signalCut_base + cut_options_dm10TESDown)*cut_weights_nom, (rel_signalCut_base + cut_options_dm10TESUp)*cut_weights_nom, (SScut_TauTau_base + cut_options_dm10TESUp)*cut_weights_nom, (SScut_TauTau_base_RelTauIso + cut_options_dm10TESUp)*cut_weights_nom, wSF_SS_TauTau, (signalCut_base + cut_options_nom)*cut_weights_nom, (SScut_TauTau_base + cut_options_nom)*cut_weights_nom, 0, wSF_OS_TauTau, "mt_tot_dm10TESDown", binnum_mt_tot, mt_totVarBinning,  "", "tt","_CMS_scale_t_dm10_13TeVDown", 0, 1);
            drawSignalRegionVarBin((signalCut_base + cut_options_nom)*cut_weights_TopPtUp, (rel_signalCut_base + cut_options_nom)*cut_weights_TopPtUp, (SScut_TauTau_base + cut_options_nom)*cut_weights_TopPtUp, (SScut_TauTau_base_RelTauIso + cut_options_nom)*cut_weights_TopPtUp, wSF_SS_TauTau, (signalCut_base + cut_options_nom)*cut_weights_nom, (SScut_TauTau_base + cut_options_nom)*cut_weights_nom, 0, wSF_OS_TauTau, "mt_tot", binnum_mt_tot, mt_totVarBinning,  "", "tt", "_CMS_xtt_ttbarShape_13TeVUp", 0, 1);
            drawSignalRegionVarBin((signalCut_base + cut_options_nom)*cut_weights_TopPtDown, (rel_signalCut_base + cut_options_nom)*cut_weights_TopPtDown, (SScut_TauTau_base + cut_options_nom)*cut_weights_nom, (SScut_TauTau_base_RelTauIso + cut_options_nom)*cut_weights_nom, wSF_SS_TauTau,(signalCut_base + cut_options_nom)*cut_weights_nom, (SScut_TauTau_base + cut_options_nom)*cut_weights_nom, 0, wSF_OS_TauTau, "mt_tot", binnum_mt_tot, mt_totVarBinning,  "", "tt","_CMS_xtt_ttbarShape_13TeVDown", 0, 1);
            drawSignalRegionVarBin((signalCut_base + cut_options_nom)*cut_weights_WWUp, (rel_signalCut_base + cut_options_nom)*cut_weights_WWUp, (SScut_TauTau_base + cut_options_nom)*cut_weights_WWUp, (SScut_TauTau_base_RelTauIso + cut_options_nom)*cut_weights_WWUp, wSF_SS_TauTau, (signalCut_base + cut_options_nom)*cut_weights_WWUp, (rel_signalCut_base + cut_options_nom)*cut_weights_WWUp, 0, wSF_OS_TauTau, "mt_tot", binnum_mt_tot, mt_totVarBinning,  "", "tt", "_CMS_xtt_WWNLOewk_13TeVUp", 0, 1);
            drawSignalRegionVarBin((signalCut_base + cut_options_nom)*cut_weights_WWDown, (rel_signalCut_base + cut_options_nom)*cut_weights_WWDown, (SScut_TauTau_base + cut_options_nom)*cut_weights_WWDown, (SScut_TauTau_base_RelTauIso + cut_options_nom)*cut_weights_WWDown, wSF_SS_TauTau, (signalCut_base + cut_options_nom)*cut_weights_nom, (SScut_TauTau_base + cut_options_nom)*cut_weights_nom, 0, wSF_OS_TauTau, "mt_tot", binnum_mt_tot, mt_totVarBinning,  "", "tt","_CMS_xtt_WWNLOewk_13TeVDown", 0, 1);
            drawSignalRegionVarBin((signalCut_base + cut_options_nom)*cut_weights_ZZUp, (rel_signalCut_base + cut_options_nom)*cut_weights_ZZUp, (SScut_TauTau_base + cut_options_nom)*cut_weights_ZZUp, (SScut_TauTau_base_RelTauIso + cut_options_nom)*cut_weights_ZZUp, wSF_SS_TauTau, (signalCut_base + cut_options_nom)*cut_weights_nom, (SScut_TauTau_base + cut_options_nom)*cut_weights_nom, 0, wSF_OS_TauTau, "mt_tot", binnum_mt_tot, mt_totVarBinning,  "", "tt", "_CMS_xtt_ZZNLOewk_13TeVUp", 0, 1);
            drawSignalRegionVarBin((signalCut_base + cut_options_nom)*cut_weights_ZZDown, (rel_signalCut_base + cut_options_nom)*cut_weights_ZZDown, (SScut_TauTau_base + cut_options_nom)*cut_weights_ZZDown, (SScut_TauTau_base_RelTauIso + cut_options_nom)*cut_weights_ZZDown, wSF_SS_TauTau, (signalCut_base + cut_options_nom)*cut_weights_nom, (SScut_TauTau_base + cut_options_nom)*cut_weights_nom, 0, wSF_OS_TauTau, "mt_tot", binnum_mt_tot, mt_totVarBinning,  "", "tt","_CMS_xtt_ZZNLOewk_13TeVDown", 0, 1);
            
            if(printSyst)
            {
            
                countTotalsVarBin((signalCut_base + cut_options_UESUp)*cut_weights_UESUp, (rel_signalCut_base + cut_options_UESUp)*cut_weights_UESUp, (SScut_TauTau_base + cut_options_UESUp)*cut_weights_UESUp, (SScut_TauTau_base_RelTauIso + cut_options_UESUp)*cut_weights_UESUp, wSF_SS_TauTau, wSF_OS_TauTau, "mt_tot_UESUp", binnum_mt_tot, mt_totVarBinning, "Mt_total", "tt", "_CMS_scale_m_13TeVUp");
                countTotalsVarBin((signalCut_base + cut_options_UESDown)*cut_weights_UESDown, (rel_signalCut_base + cut_options_UESDown)*cut_weights_UESDown, (SScut_TauTau_base + cut_options_UESDown)*cut_weights_UESDown, (SScut_TauTau_base_RelTauIso + cut_options_UESDown)*cut_weights_UESDown, wSF_SS_TauTau, wSF_OS_TauTau, "mt_tot_UESDown", binnum_mt_tot, mt_totVarBinning, "Mt_total", "tt", "_CMS_scale_m_13TeVDown");
                countTotalsVarBin((signalCut_base + cut_options_nom)*cut_weights_ZPtUp, (rel_signalCut_base + cut_options_nom)*cut_weights_ZPtUp, (SScut_TauTau_base + cut_options_nom)*cut_weights_ZPtUp, (SScut_TauTau_base_RelTauIso + cut_options_nom)*cut_weights_ZPtUp, wSF_SS_TauTau, wSF_OS_TauTau, "mt_tot", binnum_mt_tot, mt_totVarBinning, "Mt_total", "tt", "_CMS_xtt_dyShape_13TeVUp");
                countTotalsVarBin((signalCut_base + cut_options_nom)*cut_weights_ZPtDown, (rel_signalCut_base + cut_options_nom)*cut_weights_ZPtDown, (SScut_TauTau_base + cut_options_nom)*cut_weights_ZPtDown, (SScut_TauTau_base_RelTauIso + cut_options_nom)*cut_weights_ZPtDown, wSF_SS_TauTau, wSF_OS_TauTau, "mt_tot", binnum_mt_tot, mt_totVarBinning, "Mt_total", "tt", "_CMS_xtt_dyShape_13TeVDown");
                countTotalsVarBin((signalCut_base + cut_options_nom)*cut_weights_KPtUp, (rel_signalCut_base + cut_options_nom)*cut_weights_KPtUp, (SScut_TauTau_base + cut_options_nom)*cut_weights_KPtUp, (SScut_TauTau_base_RelTauIso + cut_options_nom)*cut_weights_KPtUp, wSF_SS_TauTau, wSF_OS_TauTau, "mt_tot", binnum_mt_tot, mt_totVarBinning, "Mt_total", "tt", "_CMS_xtt_wShape_13TeVUp");
                countTotalsVarBin((signalCut_base + cut_options_nom)*cut_weights_KPtDown, (rel_signalCut_base + cut_options_nom)*cut_weights_KPtDown, (SScut_TauTau_base + cut_options_nom)*cut_weights_KPtDown, (SScut_TauTau_base_RelTauIso + cut_options_nom)*cut_weights_KPtDown, wSF_SS_TauTau, wSF_OS_TauTau, "mt_tot", binnum_mt_tot, mt_totVarBinning, "Mt_total", "tt", "_CMS_xtt_wShape_13TeVDown");
                countTotalsVarBin((signalCut_base + cut_options_TauFakeUp)*cut_weights_TauFakeUp, (rel_signalCut_base + cut_options_TauFakeUp)*cut_weights_TauFakeUp, (SScut_TauTau_base + cut_options_TauFakeUp)*cut_weights_TauFakeUp, (SScut_TauTau_base_RelTauIso + cut_options_TauFakeUp)*cut_weights_TauFakeUp, wSF_SS_TauTau, wSF_OS_TauTau, "mt_tot", binnum_mt_tot, mt_totVarBinning, "Mt_total", "tt", "_CMS_xtt_jetToTauFake_13TeVUp");
                countTotalsVarBin((signalCut_base + cut_options_TauFakeDown)*cut_weights_TauFakeDown, (rel_signalCut_base + cut_options_TauFakeDown)*cut_weights_TauFakeDown, (SScut_TauTau_base + cut_options_TauFakeDown)*cut_weights_TauFakeDown, (SScut_TauTau_base_RelTauIso + cut_options_TauFakeDown)*cut_weights_TauFakeDown, wSF_SS_TauTau, wSF_OS_TauTau, "mt_tot", binnum_mt_tot, mt_totVarBinning, "Mt_total", "tt", "_CMS_xtt_jetToTauFake_13TeVDown");
                countTotalsVarBin((signalCut_base + cut_options_bTagUp)*cut_weights_bTagUp, (rel_signalCut_base + cut_options_bTagUp)*cut_weights_bTagUp, (SScut_TauTau_base + cut_options_bTagUp)*cut_weights_bTagUp, (SScut_TauTau_base_RelTauIso + cut_options_bTagUp)*cut_weights_bTagUp, wSF_SS_TauTau, wSF_OS_TauTau, "mt_tot", binnum_mt_tot, mt_totVarBinning, "Mt_total", "tt", "_CMS_xtt_bTag_13TeVUp");
                countTotalsVarBin((signalCut_base + cut_options_bTagDown)*cut_weights_bTagDown, (rel_signalCut_base + cut_options_bTagDown)*cut_weights_bTagDown, (SScut_TauTau_base + cut_options_bTagDown)*cut_weights_bTagDown, (SScut_TauTau_base_RelTauIso + cut_options_bTagDown)*cut_weights_bTagDown, wSF_SS_TauTau, wSF_OS_TauTau, "mt_tot", binnum_mt_tot, mt_totVarBinning, "Mt_total", "tt", "_CMS_xtt_bTag_13TeVDown");
                countTotalsVarBin((signalCut_base + cut_options_JECUp)*cut_weights_JECUp, (rel_signalCut_base + cut_options_JECUp)*cut_weights_JECUp, (SScut_TauTau_base + cut_options_JECUp)*cut_weights_JECUp, (SScut_TauTau_base_RelTauIso + cut_options_JECUp)*cut_weights_JECUp, wSF_SS_TauTau, wSF_OS_TauTau, "mt_tot_JEnUp", binnum_mt_tot, mt_totVarBinning, "Mt_total", "tt", "_CMS_scale_j_13TeVUp");
                countTotalsVarBin((signalCut_base + cut_options_JECDown)*cut_weights_JECDown, (rel_signalCut_base + cut_options_JECDown)*cut_weights_JECDown, (SScut_TauTau_base + cut_options_JECDown)*cut_weights_JECDown, (SScut_TauTau_base_RelTauIso + cut_options_JECDown)*cut_weights_JECDown, wSF_SS_TauTau, wSF_OS_TauTau, "mt_tot_JEnDown", binnum_mt_tot, mt_totVarBinning, "Mt_total", "tt", "_CMS_scale_j_13TeVDown");
                countTotalsVarBin((signalCut_base + cut_options_HPTUp)*cut_weights_HPTUp, (rel_signalCut_base + cut_options_HPTUp)*cut_weights_HPTUp, (SScut_TauTau_base + cut_options_HPTUp)*cut_weights_HPTUp, (SScut_TauTau_base_RelTauIso + cut_options_HPTUp)*cut_weights_HPTUp, wSF_SS_TauTau, wSF_OS_TauTau, "mt_tot", binnum_mt_tot, mt_totVarBinning, "Mt_total", "tt", "_CMS_xtt_highTauEffi_13TeVUp");
                countTotalsVarBin((signalCut_base + cut_options_HPTDown)*cut_weights_HPTDown, (rel_signalCut_base + cut_options_HPTDown)*cut_weights_HPTDown, (SScut_TauTau_base + cut_options_HPTDown)*cut_weights_HPTDown, (SScut_TauTau_base_RelTauIso + cut_options_HPTDown)*cut_weights_HPTDown, wSF_SS_TauTau, wSF_OS_TauTau, "mt_tot", binnum_mt_tot, mt_totVarBinning, "Mt_total", "tt", "_CMS_xtt_highTauEffi_13TeVDown");
                countTotalsVarBin((signalCut_base + cut_options_nom)*cut_weights_ttPtUp, (rel_signalCut_base + cut_options_nom)*cut_weights_ttPtUp, (SScut_TauTau_base + cut_options_nom)*cut_weights_ttPtUp, (SScut_TauTau_base_RelTauIso + cut_options_nom)*cut_weights_ttPtUp, wSF_SS_TauTau, wSF_OS_TauTau, "mt_tot", binnum_mt_tot, mt_totVarBinning, "Mt_total", "tt", "_CMS_xtt_tt_trigger_13TeVUp");
                countTotalsVarBin((signalCut_base + cut_options_HPTDown)*cut_weights_ttPtDown, (rel_signalCut_base + cut_options_HPTDown)*cut_weights_ttPtDown, (SScut_TauTau_base + cut_options_nom)*cut_weights_ttPtDown, (SScut_TauTau_base_RelTauIso + cut_options_nom)*cut_weights_ttPtDown, wSF_SS_TauTau, wSF_OS_TauTau, "mt_tot", binnum_mt_tot, mt_totVarBinning, "Mt_total", "tt", "_CMS_xtt_tt_trigger_13TeVDown");
                countTotalsVarBin((signalCut_base + cut_options_TESUp)*cut_weights_nom, (rel_signalCut_base + cut_options_TESUp)*cut_weights_nom, (SScut_TauTau_base + cut_options_TESUp)*cut_weights_nom, (SScut_TauTau_base_RelTauIso + cut_options_TESUp)*cut_weights_nom, wSF_SS_TauTau, wSF_OS_TauTau, "mt_tot", binnum_mt_tot, mt_totVarBinning, "Mt_total", "tt", "_CMS_scale_t_13TeVUp");
                countTotalsVarBin((signalCut_base + cut_options_TESDown)*cut_weights_nom, (rel_signalCut_base + cut_options_TESDown)*cut_weights_nom, (SScut_TauTau_base + cut_options_TESDown)*cut_weights_nom, (SScut_TauTau_base_RelTauIso + cut_options_TESDown)*cut_weights_nom, wSF_SS_TauTau, wSF_OS_TauTau, "mt_tot", binnum_mt_tot, mt_totVarBinning, "Mt_total", "tt", "_CMS_scale_t_13TeVDown");
                countTotalsVarBin((signalCut_base + cut_options_nom)*cut_weights_TopPtUp, (rel_signalCut_base + cut_options_nom)*cut_weights_TopPtUp, (SScut_TauTau_base + cut_options_nom)*cut_weights_TopPtUp, (SScut_TauTau_base_RelTauIso + cut_options_nom)*cut_weights_TopPtUp, wSF_SS_TauTau, wSF_OS_TauTau, "mt_tot", binnum_mt_tot, mt_totVarBinning, "Mt_total", "tt", "_CMS_xtt_ttbarShape_13TeVUp");
                countTotalsVarBin((signalCut_base + cut_options_nom)*cut_weights_TopPtDown, (rel_signalCut_base + cut_options_nom)*cut_weights_TopPtDown, (SScut_TauTau_base + cut_options_nom)*cut_weights_TopPtDown, (SScut_TauTau_base_RelTauIso + cut_options_nom)*cut_weights_TopPtDown, wSF_SS_TauTau, wSF_OS_TauTau, "mt_tot", binnum_mt_tot, mt_totVarBinning, "Mt_total", "tt", "_CMS_xtt_ttbarShape_13TeVDown");
                countTotalsVarBin((signalCut_base + cut_options_nom)*cut_weights_WWUp, (rel_signalCut_base + cut_options_nom)*cut_weights_WWUp, (SScut_TauTau_base + cut_options_nom)*cut_weights_WWUp, (SScut_TauTau_base_RelTauIso + cut_options_nom)*cut_weights_WWUp, wSF_SS_TauTau, wSF_OS_TauTau, "mt_tot", binnum_mt_tot, mt_totVarBinning, "Mt_total", "tt", "_CMS_xtt_WWNLOewk_13TeVUp");
                countTotalsVarBin((signalCut_base + cut_options_nom)*cut_weights_WWDown, (rel_signalCut_base + cut_options_nom)*cut_weights_WWDown, (SScut_TauTau_base + cut_options_nom)*cut_weights_WWDown, (SScut_TauTau_base_RelTauIso + cut_options_nom)*cut_weights_WWDown, wSF_SS_TauTau, wSF_OS_TauTau, "mt_tot", binnum_mt_tot, mt_totVarBinning, "Mt_total", "tt", "_CMS_xtt_WWNLOewk_13TeVDown");
                countTotalsVarBin((signalCut_base + cut_options_nom)*cut_weights_ZZUp, (rel_signalCut_base + cut_options_nom)*cut_weights_ZZUp, (SScut_TauTau_base + cut_options_nom)*cut_weights_ZZUp, (SScut_TauTau_base_RelTauIso + cut_options_nom)*cut_weights_ZZUp, wSF_SS_TauTau, wSF_OS_TauTau, "mt_tot", binnum_mt_tot, mt_totVarBinning, "Mt_total", "tt", "_CMS_xtt_ZZNLOewk_13TeVUp");
                countTotalsVarBin((signalCut_base + cut_options_nom)*cut_weights_ZZDown, (rel_signalCut_base + cut_options_nom)*cut_weights_ZZDown, (SScut_TauTau_base + cut_options_nom)*cut_weights_ZZDown, (SScut_TauTau_base_RelTauIso + cut_options_nom)*cut_weights_ZZDown, wSF_SS_TauTau, wSF_OS_TauTau, "mt_tot", binnum_mt_tot, mt_totVarBinning, "Mt_total", "tt", "_CMS_xtt_ZZNLOewk_13TeVDown");
            }
        }
        reset_files();
        
    }
    else if (choice2==2)
    {
        
        if(doControlTMVA)
        {
                drawSignalRegion((signalCut_base + TMVACutMap["TT_600"])*cut_weights_nom, (signalCut_base + TMVACutMap["TT_600"])*cut_weights_nom, (SScut_TauTau_base + TMVACutMap["TT_600"])*cut_weights_nom, (SScut_TauTau_base_RelTauIso + TMVACutMap["TT_600"])*cut_weights_nom, wSF_SS_TauTau, 1.0 , "mvaVar_tt_MZP600A0400", tmvaBinning, "sig region Function TMVA variable", "tt", "", 1, 0);
                drawSignalRegion((signalCut_base + TMVACutMap["TT_800"])*cut_weights_nom, (signalCut_base + TMVACutMap["TT_800"])*cut_weights_nom, (SScut_TauTau_base + TMVACutMap["TT_800"])*cut_weights_nom, (SScut_TauTau_base_RelTauIso + TMVACutMap["TT_800"])*cut_weights_nom, wSF_SS_TauTau, 1.0 , "mvaVar_tt_MZP800A0400", tmvaBinning, "sig region Function TMVA variable", "tt", "", 1, 0);
                drawSignalRegion((signalCut_base + TMVACutMap["TT_1000"])*cut_weights_nom, (signalCut_base + TMVACutMap["TT_1000"])*cut_weights_nom, (SScut_TauTau_base + TMVACutMap["TT_1000"])*cut_weights_nom, (SScut_TauTau_base_RelTauIso + TMVACutMap["TT_1000"])*cut_weights_nom, wSF_SS_TauTau, 1.0 , "mvaVar_tt_MZP1000A0400", tmvaBinning, "sig region Function TMVA variable", "tt", "", 1, 0);
                drawSignalRegion((signalCut_base + TMVACutMap["TT_1200"])*cut_weights_nom, (signalCut_base + TMVACutMap["TT_1200"])*cut_weights_nom, (SScut_TauTau_base + TMVACutMap["TT_1200"])*cut_weights_nom, (SScut_TauTau_base_RelTauIso + TMVACutMap["TT_1200"])*cut_weights_nom, wSF_SS_TauTau, 1.0 , "mvaVar_tt_MZP1200A0400", tmvaBinning, "sig region Function TMVA variable", "tt", "", 1, 0);
        }
        else
        {

            countTotalsVarBin((signalCut_base + cut_options_nom)*cut_weights_nom, (rel_signalCut_base + cut_options_nom)*cut_weights_nom, (SScut_TauTau_base + cut_options_nom)*cut_weights_nom, (SScut_TauTau_base_RelTauIso + cut_options_nom)*cut_weights_nom, wSF_SS_TauTau, wSF_OS_TauTau, "mt_tot", binnum_mt_tot, mt_totVarBinning, "control region Mt_total", "tt", "");
        
            drawSignalRegion((signalCut_base + cut_options_nom)*cut_weights_nom, (rel_signalCut_base + cut_options_nom)*cut_weights_nom, (SScut_TauTau_base + cut_options_nom)*cut_weights_nom, (SScut_TauTau_base_RelTauIso + cut_options_nom)*cut_weights_nom, wSF_SS_TauTau, wSF_OS_TauTau, "mt_1", mtBinning, " (met<100) mt", "tt", "", 0, 1);
            drawSignalRegion((signalCut_base + cut_options_nom)*cut_weights_nom, (rel_signalCut_base + cut_options_nom)*cut_weights_nom, (SScut_TauTau_base + cut_options_nom)*cut_weights_nom, (SScut_TauTau_base_RelTauIso + cut_options_nom)*cut_weights_nom, wSF_SS_TauTau, wSF_OS_TauTau, "m_vis", m_visBinning, "Control mvis", "tt", "", 0, 1);
            drawSignalRegion((signalCut_base + cut_options_nom)*cut_weights_nom, (rel_signalCut_base + cut_options_nom)*cut_weights_nom, (SScut_TauTau_base + cut_options_nom)*cut_weights_nom, (SScut_TauTau_base_RelTauIso + cut_options_nom)*cut_weights_nom, wSF_SS_TauTau, wSF_OS_TauTau, "nbtag", nbtagBinning, "Control nbtag", "tt", "", 0, 1);
            drawSignalRegion((signalCut_base + cut_options_nom)*cut_weights_nom, (rel_signalCut_base + cut_options_nom)*cut_weights_nom, (SScut_TauTau_base + cut_options_nom)*cut_weights_nom, (SScut_TauTau_base_RelTauIso + cut_options_nom)*cut_weights_nom, wSF_SS_TauTau, wSF_OS_TauTau, "met", metBinning, "Control met", "tt", "", 0, 1);
            drawSignalRegion((signalCut_base + cut_options_nom)*cut_weights_nom, (rel_signalCut_base + cut_options_nom)*cut_weights_nom, (SScut_TauTau_base + cut_options_nom)*cut_weights_nom, (SScut_TauTau_base_RelTauIso + cut_options_nom)*cut_weights_nom, wSF_SS_TauTau, wSF_OS_TauTau, "pt_1", ptBinning, "Control pt_1", "tt", "", 0, 1);
            drawSignalRegion((signalCut_base + cut_options_nom)*cut_weights_nom, (rel_signalCut_base + cut_options_nom)*cut_weights_nom, (SScut_TauTau_base + cut_options_nom)*cut_weights_nom, (SScut_TauTau_base_RelTauIso + cut_options_nom)*cut_weights_nom, wSF_SS_TauTau, wSF_OS_TauTau, "pt_2", ptBinning, "Control pt_2", "tt", "", 0, 1);
            drawSignalRegion((signalCut_base + cut_options_nom)*cut_weights_nom, (rel_signalCut_base + cut_options_nom)*cut_weights_nom, (SScut_TauTau_base + cut_options_nom)*cut_weights_nom, (SScut_TauTau_base_RelTauIso + cut_options_nom)*cut_weights_nom, wSF_SS_TauTau, wSF_OS_TauTau, "pt_tt", ptBinning, "Control pt_tt", "tt", "", 0, 1);
            
            drawSignalRegion((signalCut_base + cut_options_nom)*cut_weights_nom, (rel_signalCut_base + cut_options_nom)*cut_weights_nom, (SScut_TauTau_base + cut_options_nom)*cut_weights_nom, (SScut_TauTau_base_RelTauIso + cut_options_nom)*cut_weights_nom, wSF_SS_TauTau, wSF_SS_TauTau, "DeltaR_leg1_leg2", drBinning, "Control DeltaR_leg1_leg2", "tt", "", 0, 1);
            drawSignalRegion((signalCut_base + cut_options_nom)*cut_weights_nom, (rel_signalCut_base + cut_options_nom)*cut_weights_nom, (SScut_TauTau_base + cut_options_nom)*cut_weights_nom, (SScut_TauTau_base_RelTauIso + cut_options_nom)*cut_weights_nom, wSF_SS_TauTau, wSF_OS_TauTau, "m_1", mBinning, "Control m_1", "tt", "", 0, 1);
            drawSignalRegion((signalCut_base + cut_options_nom)*cut_weights_nom, (rel_signalCut_base + cut_options_nom)*cut_weights_nom, (SScut_TauTau_base + cut_options_nom)*cut_weights_nom, (SScut_TauTau_base_RelTauIso + cut_options_nom)*cut_weights_nom, wSF_SS_TauTau, wSF_OS_TauTau, "m_2", mBinning, "Control m_2", "tt", "", 0, 1);
            drawSignalRegion((signalCut_base + cut_options_nom)*cut_weights_nom, (rel_signalCut_base + cut_options_nom)*cut_weights_nom, (SScut_TauTau_base + cut_options_nom)*cut_weights_nom, (SScut_TauTau_base_RelTauIso + cut_options_nom)*cut_weights_nom, wSF_SS_TauTau, wSF_OS_TauTau,  "mt_2", mtBinning, "Control mt_2", "tt", "", 0, 1);
            drawSignalRegion((signalCut_base + cut_options_nom)*cut_weights_nom, (rel_signalCut_base + cut_options_nom)*cut_weights_nom, (SScut_TauTau_base + cut_options_nom)*cut_weights_nom, (SScut_TauTau_base_RelTauIso + cut_options_nom)*cut_weights_nom, wSF_SS_TauTau, wSF_OS_TauTau, "phi_1", phiBinning, "Control phi_1", "tt", "", 0, 1);
            drawSignalRegion((signalCut_base + cut_options_nom)*cut_weights_nom, (rel_signalCut_base + cut_options_nom)*cut_weights_nom, (SScut_TauTau_base + cut_options_nom)*cut_weights_nom, (SScut_TauTau_base_RelTauIso + cut_options_nom)*cut_weights_nom, wSF_SS_TauTau, wSF_OS_TauTau, "phi_2", phiBinning, "Control phi_2", "tt", "", 0, 1);
            drawSignalRegion((signalCut_base + cut_options_nom)*cut_weights_nom, (rel_signalCut_base + cut_options_nom)*cut_weights_nom, (SScut_TauTau_base + cut_options_nom)*cut_weights_nom, (SScut_TauTau_base_RelTauIso + cut_options_nom)*cut_weights_nom, wSF_SS_TauTau, wSF_OS_TauTau, "eta_1", etaBinning, "Control eta_1", "tt", "", 0, 1);
            drawSignalRegion((signalCut_base + cut_options_nom)*cut_weights_nom, (rel_signalCut_base + cut_options_nom)*cut_weights_nom, (SScut_TauTau_base + cut_options_nom)*cut_weights_nom, (SScut_TauTau_base_RelTauIso + cut_options_nom)*cut_weights_nom, wSF_SS_TauTau, wSF_OS_TauTau, "eta_2", etaBinning, "Control eta_2", "tt", "", 0, 1);
            drawSignalRegion((signalCut_base + cut_options_nom)*cut_weights_nom, (rel_signalCut_base + cut_options_nom)*cut_weights_nom, (SScut_TauTau_base + cut_options_nom)*cut_weights_nom, (SScut_TauTau_base_RelTauIso + cut_options_nom)*cut_weights_nom, wSF_SS_TauTau, wSF_OS_TauTau, "njetspt20", nbtagBinning, "Control njetspt20", "tt", "", 0, 1);
        }
        reset_files();
    }
}


void doEleMu()
{
    // 15% (using EleTau as placeholder)
    float mt_totVarBinning[] = {0, 4, 16, 23, 30, 35, 39, 44, 48, 52, 55, 58, 61, 65, 69, 73, 77, 81, 86, 92, 98, 105, 114, 122, 133, 143, 155, 168, 177, 185, 196, 212, 229, 246, 1200};
    int binnum_mt_tot = sizeof(mt_totVarBinning)/sizeof(*mt_totVarBinning) - 1;

	/* set channel label */

    TCut SScut_eleMu_base_RelTauIso("pairGoodForTrigger==1 && isOsPair==0 && extramuon_veto==0  && extraelec_veto==0  && dilepton_veto==0  && iso_1 < 0.10 && iso_2 < 0.15");
    TCut SScut_eleMu_base("pairGoodForTrigger==1 && isOsPair==0 && extramuon_veto==0  && extraelec_veto==0  && dilepton_veto==0  && iso_1 < 0.10 && iso_2 < 0.15");
    TCut signalCut_base("pairGoodForTrigger==1 && isOsPair==1 && extramuon_veto==0  && extraelec_veto==0  && dilepton_veto==0  && iso_1 < 0.10 && iso_2 < 0.15");
    TCut rel_signalCut_base("pairGoodForTrigger==1 && isOsPair==1 && extramuon_veto==0  && extraelec_veto==0  && dilepton_veto==0  && iso_1 < 0.10 && iso_2 < 0.15");
    TCut w_Control_base("pairGoodForTrigger==1 && isOsPair==1 && extramuon_veto==0  && extraelec_veto==0  && dilepton_veto==0  && iso_1 < 0.10 && iso_2 < 0.15");

    std::string metCutTmp;
    if (cRegion) {metCutTmp = "met > " + control_met_cutoff;}
    else {metCutTmp = met_options["NOM"];}
    TCut cut_options_nom((drCutMap[drCut] + " && " + metCutTmp + " && " + mt_options["NOM"] + " && pt_tt > 65").c_str());
    TCut cut_weights_nom = jetCutMap["NOM"];
    
    TCut cut_options_bCon((metCutTmp + " && pt_tt > 65 && m_vis > 70.").c_str());
    TCut cut_weights_bCon = jetCutMap["BCON"];
    
    TCut mt_totBlind("mt_tot < 200");
    
    //Systematic Cuts
    TCut cut_options_UESUp((drCutMap[drCut] + " && " + met_options["MET_UES_UP"] + " && " + mt_options["NOM"] + " && pt_tt > 65").c_str());
    TCut cut_options_UESDown((drCutMap[drCut] + " && " + met_options["MET_UES_DOWN"] + " && " + mt_options["NOM"] + " && pt_tt > 65").c_str());
    
    TCut cut_options_TauFakeUp((drCutMap[drCut] + " && " + met_options["NOM"] + " && " + mt_options["NOM"] + " && pt_tt > 65").c_str());
    TCut cut_options_TauFakeDown((drCutMap[drCut] + " && " + met_options["NOM"] + " && " + mt_options["NOM"] + " && pt_tt > 65").c_str());
    
    TCut cut_options_JECUp((drCutMap[drCut] + " && " + met_options["JECUP"] + " && " + mt_options["NOM"] + " && pt_tt > 65").c_str());
    TCut cut_options_JECDown((drCutMap[drCut] + " && " + met_options["JECDOWN"] + " && " + mt_options["NOM"] + " && pt_tt > 65").c_str());
    
    TCut cut_options_HPTUp((drCutMap[drCut] + " && " + met_options["NOM"] + " && " + mt_options["NOM"] + " && pt_tt > 65").c_str());
    TCut cut_options_HPTDown((drCutMap[drCut] + " && " + met_options["NOM"] + " && " + mt_options["NOM"] + " && pt_tt > 65").c_str());
    
    TCut cut_options_bTagUp((drCutMap[drCut] + " && " + met_options["NOM"] + " && " + mt_options["NOM"] + " && pt_tt > 65").c_str());
    TCut cut_options_bTagDown((drCutMap[drCut] + " && " + met_options["NOM"] + " && " + mt_options["NOM"] + " && pt_tt > 65").c_str());
    
    //Weights
    TCut cut_weights_UESUp = jetCutMap["NOM"];
    TCut cut_weights_UESDown = jetCutMap["NOM"];
    
    TCut cut_weights_ZPtUp = jetCutMap["NOM"]; cut_weights_ZPtUp *= weightMap["ZUP"];
    TCut cut_weights_ZPtDown = jetCutMap["NOM"]; cut_weights_ZPtDown *= weightMap["ZDOWN"];
    
    TCut cut_weights_KPtUp = jetCutMap["NOM"]; cut_weights_KPtUp *= weightMap["KUP"];
    TCut cut_weights_KPtDown = jetCutMap["NOM"]; cut_weights_KPtDown *= weightMap["KDOWN"];
    
    TCut cut_weights_TauFakeUp = jetCutMap["NOM"]; cut_weights_TauFakeUp *= weightMap["JTFUP"];
    TCut cut_weights_TauFakeDown = jetCutMap["NOM"]; cut_weights_TauFakeDown *= weightMap["JTFDOWN"];
    
    TCut cut_weights_JECUp = jetCutMap["BSFJECUP"];
    TCut cut_weights_JECDown = jetCutMap["BSFJECDOWN"];
    
    TCut cut_weights_HPTUp = weightMap["HPTUP"];
    TCut cut_weights_HPTDown = weightMap["HPTDOWN"];
    
    TCut cut_weights_bTagUp = jetCutMap["BSFUP"];
    TCut cut_weights_bTagDown = jetCutMap["BSFDOWN"];

	chan_label = "e + #mu";

    
    /* need to make this function more general for any dist. */
    
    double wSF_SS_eleMu = 1.0; /* in eMu case W bkg is negligible */
    double wSF_OS_eleMu = 1.0; /* in eMu case W bkg is negligible */

    if (choice2==1)
    {
        if (sigChoice==1 || nonMVAChoice==1)
        {
            setup_files_eleMu();
            
            global_title = "QCD (Same Sign) Estimate Region";
            drawSignalRegion((signalCut_base + cut_options_nom)*cut_weights_nom, (rel_signalCut_base + cut_options_nom)*cut_weights_nom, (SScut_eleMu_base + cut_options_nom)*cut_weights_nom, (SScut_eleMu_base_RelTauIso + cut_options_nom)*cut_weights_nom, wSF_SS_eleMu, wSF_OS_eleMu, "met", metBinning, "met", "em", "", 1, 1);
            reset_files();
            
        }
        else if (sigChoice==3 || nonMVAChoice==2)
        {
            setup_files_eleMu();
            global_title = "QCD (Same Sign) Estimate Region";
            drawSignalRegion((signalCut_base + cut_options_nom)*cut_weights_nom, (rel_signalCut_base + cut_options_nom)*cut_weights_nom, (SScut_eleMu_base + cut_options_nom)*cut_weights_nom, (SScut_eleMu_base_RelTauIso + cut_options_nom)*cut_weights_nom, wSF_SS_eleMu, wSF_OS_eleMu, "mt_tot", mt_totVarBinning, "sig region Mt_total", "em", "", 1, 1);
            
            if (doSyst)
            {
                std::cout << "syst for em started" << std::endl;
                //syst
                
            }
            reset_files();

        }
    }
    else if (choice2==2)
    {
        setup_files_eleMu();
        global_title = "QCD (Same Sign) Estimate Region";

        /* now draw the signal region in mt */

        global_title = "MET > " + control_met_cutoff + " GeV";

        drawSignalRegion((signalCut_base + cut_options_nom)*cut_weights_nom, (rel_signalCut_base + cut_options_nom)*cut_weights_nom, (SScut_eleMu_base + cut_options_nom)*cut_weights_nom, (SScut_eleMu_base_RelTauIso + cut_options_nom)*cut_weights_nom, wSF_SS_eleMu, wSF_OS_eleMu, "mt_1", mtBinning, "Control mt", "em", "", 1, 0);
        drawSignalRegion((signalCut_base + cut_options_nom)*cut_weights_nom, (rel_signalCut_base + cut_options_nom)*cut_weights_nom, (SScut_eleMu_base + cut_options_nom)*cut_weights_nom, (SScut_eleMu_base_RelTauIso + cut_options_nom)*cut_weights_nom, wSF_SS_eleMu, wSF_OS_eleMu, "m_vis", m_visBinning, "Control mvis", "em", "", 1, 0);
        drawSignalRegion((signalCut_base + cut_options_nom)*cut_weights_nom, (rel_signalCut_base + cut_options_nom)*cut_weights_nom, (SScut_eleMu_base + cut_options_nom)*cut_weights_nom, (SScut_eleMu_base_RelTauIso + cut_options_nom)*cut_weights_nom, wSF_SS_eleMu, wSF_OS_eleMu, "nbtag", nbtagBinning, "Control nbtag", "em", "", 1, 0);
        drawSignalRegion((signalCut_base + cut_options_nom)*cut_weights_nom, (rel_signalCut_base + cut_options_nom)*cut_weights_nom, (SScut_eleMu_base + cut_options_nom)*cut_weights_nom, (SScut_eleMu_base_RelTauIso + cut_options_nom)*cut_weights_nom, wSF_SS_eleMu, wSF_OS_eleMu, "met", metBinning, "Control met", "em", "", 1, 0);
        drawSignalRegion((signalCut_base + cut_options_nom)*cut_weights_nom, (rel_signalCut_base + cut_options_nom)*cut_weights_nom, (SScut_eleMu_base + cut_options_nom)*cut_weights_nom, (SScut_eleMu_base_RelTauIso + cut_options_nom)*cut_weights_nom, wSF_SS_eleMu, wSF_OS_eleMu, "pt_1", ptBinning, "Control pt_1", "em", "", 1, 0);
        drawSignalRegion((signalCut_base + cut_options_nom)*cut_weights_nom, (rel_signalCut_base + cut_options_nom)*cut_weights_nom, (SScut_eleMu_base + cut_options_nom)*cut_weights_nom, (SScut_eleMu_base_RelTauIso + cut_options_nom)*cut_weights_nom, wSF_SS_eleMu, wSF_OS_eleMu, "pt_2", ptBinning, "Control pt_2", "em", "", 1, 0);
        drawSignalRegion((signalCut_base + cut_options_nom)*cut_weights_nom, (rel_signalCut_base + cut_options_nom)*cut_weights_nom, (SScut_eleMu_base + cut_options_nom)*cut_weights_nom, (SScut_eleMu_base_RelTauIso + cut_options_nom)*cut_weights_nom, wSF_SS_eleMu, wSF_OS_eleMu, "LPT", lptBinning, "Control LPT", "em", "", 1, 0);
        drawSignalRegion((signalCut_base + cut_options_nom)*cut_weights_nom, (rel_signalCut_base + cut_options_nom)*cut_weights_nom, (SScut_eleMu_base + cut_options_nom)*cut_weights_nom, (SScut_eleMu_base_RelTauIso + cut_options_nom)*cut_weights_nom, wSF_SS_eleMu, wSF_SS_eleMu, "DeltaPhi_leg1_leg2", phiBinning, "Control Function DeltaPhi_leg1_leg2", "em", "", 1, 0);
        drawSignalRegion((signalCut_base + cut_options_nom)*cut_weights_nom, (rel_signalCut_base + cut_options_nom)*cut_weights_nom, (SScut_eleMu_base + cut_options_nom)*cut_weights_nom, (SScut_eleMu_base_RelTauIso + cut_options_nom)*cut_weights_nom, wSF_SS_eleMu, wSF_SS_eleMu, "DeltaPhi_MVAMET_Higgs", phiBinning, "Control Function DeltaPhi_MVAMET_Higgs", "em", "", 1, 0);
        drawSignalRegion((signalCut_base + cut_options_nom)*cut_weights_nom, (rel_signalCut_base + cut_options_nom)*cut_weights_nom, (SScut_eleMu_base + cut_options_nom)*cut_weights_nom, (SScut_eleMu_base_RelTauIso + cut_options_nom)*cut_weights_nom, wSF_SS_eleMu, wSF_SS_eleMu, "DeltaPhi_PFMET_Higgs", phiBinning, "Control Function DeltaPhi_PFMET_Higgs", "em", "", 1, 0);

        reset_files();
    }
    

}

void doEleTau()
{
    setBkgErr("et");

    float mt_totVarBinning[] = {0,20,40,60,90,120,150,180,210,235,260,285,325,400,2000};
    int binnum_mt_tot = sizeof(mt_totVarBinning)/sizeof(*mt_totVarBinning) - 1;
    
    /*
    TCut SScut_eleTau_base_RelTauIso("isOsPair==0 && pairGoodForTrigger==1 && extramuon_veto==0  && extraelec_veto==0  && dilepton_veto==0 && iso_1 < 0.1 && againstElectronTightMVA6_2 > 0.5 && againstMuonLoose3_2 > 0.5 && decayModeFinding_2 > 0.5 && byLooseIsolationMVArun2v1DBoldDMwLT_2 > 0.5");
    
    TCut rel_signalCut_base("isOsPair==1 && pairGoodForTrigger==1 && extramuon_veto==0  && extraelec_veto==0  && dilepton_veto==0 && iso_1 < 0.1 && againstElectronTightMVA6_2 > 0.5 && againstMuonLoose3_2 > 0.5 && decayModeFinding_2 > 0.5 && byLooseIsolationMVArun2v1DBoldDMwLT_2 > 0.5");
    */
    
    TCut SScut_eleTau_base_RelTauIso("isOsPair==0 && pairGoodForTrigger==1 && extramuon_veto==0  && extraelec_veto==0  && dilepton_veto==0 && iso_1 < 0.3 && againstElectronTightMVA6_2 > 0.5 && againstMuonLoose3_2 > 0.5 && decayModeFinding_2 > 0.5");
    
    TCut rel_signalCut_base("isOsPair==1 && pairGoodForTrigger==1 && extramuon_veto==0  && extraelec_veto==0  && dilepton_veto==0 && iso_1 < 0.3 && againstElectronTightMVA6_2 > 0.5 && againstMuonLoose3_2 > 0.5 && decayModeFinding_2 > 0.5");
    
    TCut wSS_Control_base_RelTauIso("isOsPair==0 && pairGoodForTrigger==1 && extramuon_veto==0  && extraelec_veto==0  && dilepton_veto==0 && iso_1 < 0.3 && againstElectronTightMVA6_2 > 0.5 && againstMuonLoose3_2 > 0.5 && byTightIsolationMVArun2v1DBoldDMwLT_2 < 0.5 && decayModeFinding_2 > 0.5");
    
    TCut SScut_eleTau_base("isOsPair==0 && pairGoodForTrigger==1 && extramuon_veto==0  && extraelec_veto==0  && dilepton_veto==0 && iso_1 < 0.1 && againstElectronTightMVA6_2 > 0.5 && againstMuonLoose3_2 > 0.5 && byTightIsolationMVArun2v1DBoldDMwLT_2 > 0.5 && decayModeFinding_2 > 0.5");
    
    TCut signalCut_base("isOsPair==1 && pairGoodForTrigger==1 && extramuon_veto==0  && extraelec_veto==0  && dilepton_veto==0 && iso_1 < 0.1 && againstElectronTightMVA6_2 > 0.5 && againstMuonLoose3_2 > 0.5 && byTightIsolationMVArun2v1DBoldDMwLT_2 > 0.5 && decayModeFinding_2 > 0.5");
    
    TCut w_Control_base("isOsPair==1 && pairGoodForTrigger==1 && extramuon_veto==0  && extraelec_veto==0  && dilepton_veto==0 && iso_1 < 0.1 && againstElectronTightMVA6_2 > 0.5 && againstMuonLoose3_2 > 0.5 && byTightIsolationMVArun2v1DBoldDMwLT_2 < 0.5 && decayModeFinding_2 > 0.5");
    
    TCut wSS_Control_base("isOsPair==0 && pairGoodForTrigger==1 && extramuon_veto==0  && extraelec_veto==0  && dilepton_veto==0 && iso_1 < 0.1 && againstElectronTightMVA6_2 > 0.5 && againstMuonLoose3_2 > 0.5 && byTightIsolationMVArun2v1DBoldDMwLT_2 < 0.5 && decayModeFinding_2 > 0.5");
    
    std::string metCutTmp;
    if (cRegion)
    {
        if(doControlTMVA){metCutTmp = "met > 0.";}
        else{metCutTmp = "met > " + control_met_cutoff;}
    }
    else
    {
        metCutTmp = met_options["NOM"];
    }
    
    TCut cut_options_nom(("pt_1 > 26. && pt_2 > 20. && " + drCutMap[drCut] + " && " + metCutTmp + " && " + mt_options["NOM"] + " && pt_tt > 65. && mt_tot <2000.").c_str());
    TCut cut_weights_nom = jetCutMap["NOM"];
    
    TCut mt_totBlind("mt_tot < 200");
    
    //Systematic Cuts
    TCut cut_options_UESUp(("pt_1 > 26. && pt_2 > 20. && " + drCutMap[drCut] + " && " + met_options["MET_UES_UP"] + " && " + mt_options["NOM"] + " && pt_tt > 65. && mt_tot <2000.").c_str());
    TCut cut_options_UESDown(("pt_1 > 26. && pt_2 > 20. && " + drCutMap[drCut] + " && " + met_options["MET_UES_DOWN"] + " && " + mt_options["NOM"] + " && pt_tt > 65. && mt_tot <2000.").c_str());
    
    TCut cut_options_TESUp(("pt_1_TESUp > 26. && pt_2_TESUp > 20. && " + drCutMap[drCut] + " && " + met_options["TES_UP"] + " && " + mt_options["TES_UP"] + " && pt_tt_TESUp > 65. && mt_tot <2000.").c_str());
    TCut cut_options_TESDown(("pt_1_TESDown > 26. && pt_2_TESDown > 20. && " + drCutMap[drCut] + " && " + met_options["TES_DOWN"] + " && " + mt_options["TES_DOWN"] + " && pt_tt_TESDown > 65. && mt_tot <2000.").c_str());
    
    TCut cut_options_dm0TESUp(("pt_1_dm0TESUp > 26. && pt_2_dm0TESUp > 20. && " + drCutMap[drCut] + " && " + met_options["TES_DM0UP"] + " && " + mt_options["TES_DM0UP"] + " && pt_tt_dm0TESUp > 65. && mt_tot <2000.").c_str());
    TCut cut_options_dm0TESDown(("pt_1_dm0TESDown > 26. && pt_2_dm0TESDown > 20. && " + drCutMap[drCut] + " && " + met_options["TES_DM0DOWN"] + " && " + mt_options["TES_DM0DOWN"] + " && pt_tt_dm0TESDown > 65. && mt_tot <2000.").c_str());
    
    TCut cut_options_dm1TESUp(("pt_1_dm1TESUp > 26. && pt_2_dm1TESUp > 20. && " + drCutMap[drCut] + " && " + met_options["TES_DM1UP"] + " && " + mt_options["TES_DM1UP"] + " && pt_tt_dm1TESUp > 65. && mt_tot <2000.").c_str());
    TCut cut_options_dm1TESDown(("pt_1_dm1TESDown > 26. && pt_2_dm1TESDown > 20. && " + drCutMap[drCut] + " && " + met_options["TES_DM1DOWN"] + " && " + mt_options["TES_DM1DOWN"] + " && pt_tt_dm1TESDown > 65. && mt_tot <2000.").c_str());
    
    TCut cut_options_dm10TESUp(("pt_1_dm10TESUp > 26. && pt_2_dm10TESUp > 20. && " + drCutMap[drCut] + " && " + met_options["TES_DM10UP"] + " && " + mt_options["TES_DM10UP"] + " && pt_tt_dm10TESUp > 65. && mt_tot <2000.").c_str());
    TCut cut_options_dm10TESDown(("pt_1_dm10TESDown > 26. && pt_2_dm10TESDown > 20. && " + drCutMap[drCut] + " && " + met_options["TES_DM10DOWN"] + " && " + mt_options["TES_DM10DOWN"] + " && pt_tt_dm10TESDown > 65. && mt_tot <2000.").c_str());
    
    TCut cut_options_TauFakeUp(("pt_1 > 26. && pt_2 > 20. && " + drCutMap[drCut] + " && " + met_options["NOM"] + " && " + mt_options["NOM"] + " && pt_tt > 65. && mt_tot <2000.").c_str());
    TCut cut_options_TauFakeDown(("pt_1 > 26. && pt_2 > 20. && " + drCutMap[drCut] + " && " + met_options["NOM"] + " && " + mt_options["NOM"] + " && pt_tt > 65. && mt_tot <2000.").c_str());
    
    TCut cut_options_JECUp(("pt_1 > 26. && pt_2 > 20. && " + drCutMap[drCut] + " && " + met_options["JECUP"] + " && " + mt_options["NOM"] + " && pt_tt > 65. && mt_tot <2000.").c_str());
    TCut cut_options_JECDown(("pt_1 > 26. && pt_2 > 20. && " + drCutMap[drCut] + " && " + met_options["JECDOWN"] + " && " + mt_options["NOM"] + " && pt_tt > 65. && mt_tot <2000.").c_str());
    
    TCut cut_options_HPTUp(("pt_1 > 26. && pt_2 > 20. && " + drCutMap[drCut] + " && " + met_options["NOM"] + " && " + mt_options["NOM"] + " && pt_tt > 65. && mt_tot <2000.").c_str());
    TCut cut_options_HPTDown(("pt_1 > 26. && pt_2 > 20. && " + drCutMap[drCut] + " && " + met_options["NOM"] + " && " + mt_options["NOM"] + " && pt_tt > 65. && mt_tot <2000.").c_str());
    
    TCut cut_options_bTagUp(("pt_1 > 26. && pt_2 > 20. && " + drCutMap[drCut] + " && " + met_options["NOM"] + " && " + mt_options["NOM"] + " && pt_tt > 65. && mt_tot <2000.").c_str());
    TCut cut_options_bTagDown(("pt_1 > 26. && pt_2 > 20. && " + drCutMap[drCut] + " && " + met_options["NOM"] + " && " + mt_options["NOM"] + " && pt_tt > 65. && mt_tot <2000.").c_str());
    
    //Systematic cuts for WCon

    TCut cut_options_wCon(("pt_1 > 26. && pt_2 > 20. && DeltaR_leg1_leg2 > 0.3 && met > 105. && mt_1 > 80. && pt_tt > 65. && mt_tot <2000."));
    
    TCut cut_options_wCon_UESUp(("pt_1 > 26. && pt_2 > 20. && DeltaR_leg1_leg2 > 0.3 && pfmet_type1_UnclusteredEnUp_Pt > 105. && pfmet_type1_UnclusteredEnUp_MT1 > 80. && pt_tt > 65. && mt_tot <2000."));
    TCut cut_options_wCon_UESDown(("pt_1 > 26. && pt_2 > 20. && DeltaR_leg1_leg2 > 0.3 && pfmet_type1_UnclusteredEnDown_Pt > 105. && pfmet_type1_UnclusteredEnDown_MT1 > 80. && pt_tt > 65. && mt_tot <2000."));
    
    TCut cut_options_wCon_TESUp(("pt_1_TESUp > 26. && pt_2_TESUp > 20. && DeltaR_leg1_leg2 > 0.3 && pfmet_type1_TESUp_Pt > 105. && mt_1_TESUp > 80. && pt_tt_TESUp > 65. && mt_tot <2000."));
    TCut cut_options_wCon_TESDown(("pt_1_TESDown > 26. && pt_2_TESDown > 20. && DeltaR_leg1_leg2 > 0.3 && pfmet_type1_TESDown_Pt > 105. && mt_1_TESDown > 80. && pt_tt_TESDown > 65. && mt_tot <2000."));
    
    TCut cut_options_wCon_dm0TESUp(("pt_1_dm0TESUp > 26. && pt_2_dm0TESUp > 20. && DeltaR_leg1_leg2 > 0.3 && pfmet_type1_dm0TESUp_Pt > 105. && mt_1_dm0TESUp > 80. && pt_tt_dm0TESUp > 65. && mt_tot <2000."));
    TCut cut_options_wCon_dm0TESDown(("pt_1_dm0TESDown > 26. && pt_2_dm0TESDown > 20. && DeltaR_leg1_leg2 > 0.3 && pfmet_type1_dm0TESDown_Pt > 105. && mt_1_dm0TESDown > 80. && pt_tt_dm0TESDown > 65. && mt_tot <2000."));
    
    TCut cut_options_wCon_dm1TESUp(("pt_1_dm1TESUp > 26. && pt_2_dm1TESUp > 20. && DeltaR_leg1_leg2 > 0.3 && pfmet_type1_dm1TESUp_Pt > 105. && mt_1_dm1TESUp > 80. && pt_tt_dm1TESUp > 65. && mt_tot <2000."));
    TCut cut_options_wCon_dm1TESDown(("pt_1_dm1TESDown > 26. && pt_2_dm1TESDown > 20. && DeltaR_leg1_leg2 > 0.3 && pfmet_type1_dm1TESDown_Pt > 105. && mt_1_dm1TESDown > 80. && pt_tt_dm1TESDown > 65. && mt_tot <2000."));
    
    TCut cut_options_wCon_dm10TESUp(("pt_1_dm10TESUp > 26. && pt_2_dm10TESUp > 20. && DeltaR_leg1_leg2 > 0.3 && pfmet_type1_dm10TESUp_Pt > 105. && mt_1_dm10TESUp > 80. && pt_tt_dm10TESUp > 65. && mt_tot <2000."));
    TCut cut_options_wCon_dm10TESDown(("pt_1_dm10TESDown > 26. && pt_2_dm10TESDown > 20. && DeltaR_leg1_leg2 > 0.3 && pfmet_type1_dm10TESDown_Pt > 105. && mt_1_dm10TESDown > 80. && pt_tt_dm10TESDown > 65. && mt_tot <2000."));
    
    TCut cut_options_wCon_JECUp(("pt_1 > 26. && pt_2 > 20. && DeltaR_leg1_leg2 > 0.3 && pfmet_type1_JetEnUp_Pt > 105. && pfmt_1_JEnUp > 80. && pt_tt > 65. && mt_tot <2000."));
    TCut cut_options_wCon_JECDown(("pt_1 > 26. && pt_2 > 20. && DeltaR_leg1_leg2 > 0.3 && pfmet_type1_JetEnDown_Pt > 105. && pfmt_1_JEnDown > 80. && pt_tt > 65. && mt_tot <2000."));
    
    
    //Weights
    TCut cut_weights_UESUp = jetCutMap["NOM"];
    TCut cut_weights_UESDown = jetCutMap["NOM"];
    
    TCut cut_weights_ZPtUp = jetCutMap["NOM"]; cut_weights_ZPtUp *= weightMap["ZUP"];
    TCut cut_weights_ZPtDown = jetCutMap["NOM"]; cut_weights_ZPtDown *= weightMap["ZDOWN"];
    
    TCut cut_weights_KPtUp = jetCutMap["NOM"]; cut_weights_KPtUp *= weightMap["KUP"];
    TCut cut_weights_KPtDown = jetCutMap["NOM"]; cut_weights_KPtDown *= weightMap["KDOWN"];
    
    TCut cut_weights_TauFakeUp = jetCutMap["NOM"]; cut_weights_TauFakeUp *= weightMap["JTFUP"];
    TCut cut_weights_TauFakeDown = jetCutMap["NOM"]; cut_weights_TauFakeDown *= weightMap["JTFDOWN"];
    
    TCut cut_weights_JECUp = jetCutMap["BSFJECUP"];
    TCut cut_weights_JECDown = jetCutMap["BSFJECDOWN"];
    
    TCut cut_weights_HPTUp = jetCutMap["NOM"]; cut_weights_HPTUp *= weightMap["HPTUP"];
    TCut cut_weights_HPTDown = jetCutMap["NOM"]; cut_weights_HPTDown *= weightMap["HPTDOWN"];
    
    TCut cut_weights_bTagUp = jetCutMap["BSFUP"];
    TCut cut_weights_bTagDown = jetCutMap["BSFDOWN"];
    
    TCut cut_weights_TopPtUp = jetCutMap["NOM"]; cut_weights_TopPtUp *= weightMap["TOPPTUP"];
    TCut cut_weights_TopPtDown = jetCutMap["NOM"];
    
    TCut cut_weights_WWUp = jetCutMap["NOM"]; cut_weights_WWUp *= weightMap["WWUP"];
    TCut cut_weights_WWDown = jetCutMap["NOM"]; cut_weights_WWDown *= weightMap["WWDOWN"];
    
    TCut cut_weights_ZZUp = jetCutMap["NOM"]; cut_weights_ZZUp *= weightMap["ZZUP"];
    TCut cut_weights_ZZDown = jetCutMap["NOM"]; cut_weights_ZZDown *= weightMap["ZZDown"];
    
	chan_label = "e + #tau_{h}";
    
    /* setup files for eleTau channel */
    setup_files_eleTau();
    
    double wSF_SS_eleTau = 1.0; /* temp */

    double wSF_OS_eleTau = wjetsNorm((w_Control_base + cut_options_wCon)*cut_weights_nom, (w_Control_base + cut_options_wCon)*cut_weights_nom, (wSS_Control_base + cut_options_wCon)*cut_weights_nom, (wSS_Control_base_RelTauIso + cut_options_wCon)*cut_weights_nom, (w_Control_base + cut_options_wCon)*cut_weights_nom, (wSS_Control_base + cut_options_wCon)*cut_weights_nom, 0, "mt_tot", binnum_mt_tot, mt_totVarBinning , "eleTauWjetsOS", "et", "", 1);
    
    double wSF_OS_eleTau_UESUp = 1.0;
    double wSF_OS_eleTau_UESDown = 1.0;
    double wSF_OS_eleTau_ZPtUp = 1.0;
    double wSF_OS_eleTau_ZPtDown = 1.0;
    double wSF_OS_eleTau_KPtUp = 1.0;
    double wSF_OS_eleTau_KPtDown = 1.0;
    double wSF_OS_eleTau_TauFakeUp = 1.0;
    double wSF_OS_eleTau_TauFakeDown = 1.0;
    double wSF_OS_eleTau_ZLUp = 1.0;
    double wSF_OS_eleTau_ZLDown = 1.0;
    double wSF_OS_eleTau_bTagUp = 1.0;
    double wSF_OS_eleTau_bTagDown = 1.0;
    double wSF_OS_eleTau_JECUp = 1.0;
    double wSF_OS_eleTau_JECDown = 1.0;
    double wSF_OS_eleTau_HPTUp = 1.0;
    double wSF_OS_eleTau_HPTDown = 1.0;
    double wSF_OS_eleTau_TESUp = 1.0;
    double wSF_OS_eleTau_TESDown = 1.0;
    double wSF_OS_eleTau_dm0TESUp = 1.0;
    double wSF_OS_eleTau_dm0TESDown = 1.0;
    double wSF_OS_eleTau_dm1TESUp = 1.0;
    double wSF_OS_eleTau_dm1TESDown = 1.0;
    double wSF_OS_eleTau_dm10TESUp = 1.0;
    double wSF_OS_eleTau_dm10TESDown = 1.0;
    double wSF_OS_eleTau_TopPtUp = 1.0;
    double wSF_OS_eleTau_TopPtDown = 1.0;
    double wSF_OS_eleTau_ZZnloUp = 1.0;
    double wSF_OS_eleTau_ZZnloDown = 1.0;
    double wSF_OS_eleTau_WWnloUp = 1.0;
    double wSF_OS_eleTau_WWnloDown = 1.0;
    
    std::cout << " WJets SF EleTau: " << wSF_OS_eleTau<< std::endl;
    
    if (doSyst)
    {

        wSF_OS_eleTau_UESUp = wjetsNorm((w_Control_base + cut_options_wCon_UESUp)*cut_weights_nom, (w_Control_base + cut_options_wCon_UESUp)*cut_weights_nom, (wSS_Control_base + cut_options_wCon_UESUp)*cut_weights_nom, (wSS_Control_base_RelTauIso + cut_options_wCon_UESUp)*cut_weights_nom, (w_Control_base + cut_options_wCon)*cut_weights_nom, (wSS_Control_base + cut_options_wCon)*cut_weights_nom, 0, "mt_tot_UESUp", binnum_mt_tot, mt_totVarBinning , "eleTauWjetsOS_UESUp", "et", "_CMS_scale_m_13TeVUp", 1);
        wSF_OS_eleTau_UESDown = wjetsNorm((w_Control_base + cut_options_wCon_UESDown)*cut_weights_nom, (w_Control_base + cut_options_wCon_UESDown)*cut_weights_nom, (wSS_Control_base + cut_options_wCon_UESDown)*cut_weights_nom, (wSS_Control_base_RelTauIso + cut_options_wCon_UESDown)*cut_weights_nom, (w_Control_base + cut_options_wCon)*cut_weights_nom, (wSS_Control_base + cut_options_wCon)*cut_weights_nom, 0, "mt_tot_UESDown", binnum_mt_tot, mt_totVarBinning , "eleTauWjetsOS_UESDown", "et", "_CMS_scale_m_13TeVDown", 1);
        
        wSF_OS_eleTau_ZPtUp = wjetsNorm((w_Control_base + cut_options_wCon)*cut_weights_ZPtUp, (w_Control_base + cut_options_wCon)*cut_weights_ZPtUp, (wSS_Control_base + cut_options_wCon)*cut_weights_ZPtUp, (wSS_Control_base_RelTauIso + cut_options_wCon)*cut_weights_ZPtUp, (w_Control_base + cut_options_wCon)*cut_weights_nom, (wSS_Control_base + cut_options_wCon)*cut_weights_nom, 0, "mt_tot", binnum_mt_tot, mt_totVarBinning , "eleTauWjetsOS", "et", "_CMS_xtt_dyShape_13TeVUp", 1);
        wSF_OS_eleTau_ZPtDown = wjetsNorm((w_Control_base + cut_options_wCon)*cut_weights_ZPtDown, (w_Control_base + cut_options_wCon)*cut_weights_ZPtDown, (wSS_Control_base + cut_options_wCon)*cut_weights_ZPtDown, (wSS_Control_base_RelTauIso + cut_options_wCon)*cut_weights_ZPtDown, (w_Control_base + cut_options_wCon)*cut_weights_nom, (wSS_Control_base + cut_options_wCon)*cut_weights_nom, 0, "mt_tot", binnum_mt_tot, mt_totVarBinning , "eleTauWjetsOS", "et", "_CMS_xtt_dyShape_13TeVDown", 1);
        
        wSF_OS_eleTau_KPtUp = wjetsNorm((w_Control_base + cut_options_wCon)*cut_weights_KPtUp, (w_Control_base + cut_options_wCon)*cut_weights_KPtUp, (wSS_Control_base + cut_options_wCon)*cut_weights_ZPtUp, (wSS_Control_base_RelTauIso + cut_options_wCon)*cut_weights_ZPtUp, (w_Control_base + cut_options_wCon)*cut_weights_nom, (wSS_Control_base + cut_options_wCon)*cut_weights_nom, 0, "mt_tot", binnum_mt_tot, mt_totVarBinning , "eleTauWjetsOS", "et", "_CMS_xtt_wShape_13TeVUp", 1);
        wSF_OS_eleTau_KPtDown = wjetsNorm((w_Control_base + cut_options_wCon)*cut_weights_KPtDown, (w_Control_base + cut_options_wCon)*cut_weights_KPtDown, (wSS_Control_base + cut_options_wCon)*cut_weights_KPtDown, (wSS_Control_base_RelTauIso + cut_options_wCon)*cut_weights_KPtDown, (w_Control_base + cut_options_wCon)*cut_weights_nom, (wSS_Control_base + cut_options_wCon)*cut_weights_nom, 0, "mt_tot", binnum_mt_tot, mt_totVarBinning , "eleTauWjetsOS", "et", "_CMS_xtt_wShape_13TeVDown", 1);
        
        wSF_OS_eleTau_TauFakeUp = wjetsNorm((w_Control_base + cut_options_wCon)*cut_weights_TauFakeUp, (w_Control_base + cut_options_wCon)*cut_weights_TauFakeUp, (wSS_Control_base + cut_options_wCon)*cut_weights_ZPtUp, (wSS_Control_base_RelTauIso + cut_options_wCon)*cut_weights_ZPtUp, (w_Control_base + cut_options_wCon)*cut_weights_nom, (wSS_Control_base_RelTauIso + cut_options_wCon)*cut_weights_nom, 1, "mt_tot", binnum_mt_tot, mt_totVarBinning , "eleTauWjetsOS", "et", "_CMS_xtt_jetToTauFake_13TeVUp", 1);
        wSF_OS_eleTau_TauFakeDown = wjetsNorm((w_Control_base + cut_options_wCon)*cut_weights_TauFakeDown, (w_Control_base + cut_options_wCon)*cut_weights_TauFakeDown, (wSS_Control_base + cut_options_wCon)*cut_weights_TauFakeDown, (wSS_Control_base_RelTauIso + cut_options_wCon)*cut_weights_TauFakeDown, (w_Control_base + cut_options_wCon)*cut_weights_nom, (wSS_Control_base_RelTauIso + cut_options_wCon)*cut_weights_nom, 1, "mt_tot", binnum_mt_tot, mt_totVarBinning , "eleTauWjetsOS", "et", "_CMS_xtt_jetToTauFake_13TeVDown", 1);
        
        wSF_OS_eleTau_ZLUp = wjetsNorm((w_Control_base + cut_options_wCon)*cut_weights_nom, (w_Control_base + cut_options_wCon)*cut_weights_nom, (wSS_Control_base + cut_options_wCon)*cut_weights_nom, (wSS_Control_base_RelTauIso + cut_options_wCon)*cut_weights_nom, (w_Control_base + cut_options_wCon)*cut_weights_nom, (wSS_Control_base + cut_options_wCon)*cut_weights_nom, 0, "mt_tot*0.98", binnum_mt_tot, mt_totVarBinning , "eleTauWjetsOS", "et", "_CMS_xtt_ZLScale_etau_13TeVUp", 1);
        wSF_OS_eleTau_ZLDown = wjetsNorm((w_Control_base + cut_options_wCon)*cut_weights_nom, (w_Control_base + cut_options_wCon)*cut_weights_nom, (wSS_Control_base + cut_options_wCon)*cut_weights_nom, (wSS_Control_base_RelTauIso + cut_options_wCon)*cut_weights_nom, (w_Control_base + cut_options_wCon)*cut_weights_nom, (wSS_Control_base + cut_options_wCon)*cut_weights_nom, 0, "mt_tot*1.02", binnum_mt_tot, mt_totVarBinning , "eleTauWjetsOS", "et", "_CMS_xtt_ZLScale_etau_13TeVDown", 1);
        
        wSF_OS_eleTau_bTagUp = wjetsNorm((w_Control_base + cut_options_wCon)*cut_weights_bTagUp, (w_Control_base + cut_options_wCon)*cut_weights_bTagUp, (wSS_Control_base + cut_options_wCon)*cut_weights_bTagUp, (wSS_Control_base_RelTauIso + cut_options_wCon)*cut_weights_bTagUp, (w_Control_base + cut_options_wCon)*cut_weights_nom, (wSS_Control_base + cut_options_wCon)*cut_weights_nom, 0, "mt_tot", binnum_mt_tot, mt_totVarBinning , "eleTauWjetsOS", "et", "_CMS_xtt_bTag_13TeVUp", 1);
        wSF_OS_eleTau_bTagDown = wjetsNorm((w_Control_base + cut_options_wCon)*cut_weights_bTagDown, (w_Control_base + cut_options_wCon)*cut_weights_bTagDown, (wSS_Control_base + cut_options_wCon)*cut_weights_bTagDown, (wSS_Control_base_RelTauIso + cut_options_wCon)*cut_weights_bTagDown, (w_Control_base + cut_options_wCon)*cut_weights_nom, (wSS_Control_base + cut_options_wCon)*cut_weights_nom, 0, "mt_tot", binnum_mt_tot, mt_totVarBinning , "eleTauWjetsOS", "et", "_CMS_xtt_bTag_13TeVDown", 1);
        
        wSF_OS_eleTau_JECUp = wjetsNorm((w_Control_base + cut_options_wCon_JECUp)*cut_weights_JECUp, (w_Control_base + cut_options_wCon_JECUp)*cut_weights_JECUp, (wSS_Control_base + cut_options_wCon_JECUp)*cut_weights_JECUp, (wSS_Control_base_RelTauIso + cut_options_wCon_JECUp)*cut_weights_JECUp, (w_Control_base + cut_options_wCon)*cut_weights_nom, (wSS_Control_base + cut_options_wCon)*cut_weights_nom, 0, "mt_tot_JEnUp", binnum_mt_tot, mt_totVarBinning , "eleTauWjetsOS", "et", "_CMS_scale_j_13TeVUp", 1);
        wSF_OS_eleTau_JECDown = wjetsNorm((w_Control_base + cut_options_wCon_JECDown)*cut_weights_JECDown, (w_Control_base + cut_options_wCon_JECDown)*cut_weights_JECDown, (wSS_Control_base + cut_options_wCon_JECDown)*cut_weights_JECDown, (wSS_Control_base_RelTauIso + cut_options_wCon_JECDown)*cut_weights_JECDown, (w_Control_base + cut_options_wCon)*cut_weights_nom, (wSS_Control_base + cut_options_wCon)*cut_weights_nom, 0, "mt_tot_JEnDown", binnum_mt_tot, mt_totVarBinning , "eleTauWjetsOS", "et", "_CMS_scale_j_13TeVDown", 1);
        
        wSF_OS_eleTau_HPTUp = wjetsNorm((w_Control_base + cut_options_wCon)*cut_weights_HPTUp, (w_Control_base + cut_options_wCon)*cut_weights_HPTUp, (wSS_Control_base + cut_options_wCon)*cut_weights_HPTUp, (wSS_Control_base_RelTauIso + cut_options_wCon)*cut_weights_HPTUp, (w_Control_base + cut_options_wCon)*cut_weights_nom, (wSS_Control_base + cut_options_wCon)*cut_weights_nom, 0, "mt_tot", binnum_mt_tot, mt_totVarBinning , "eleTauWjetsOS", "et", "_CMS_xtt_highTauEffi_13TeVUp", 1);
        wSF_OS_eleTau_HPTDown = wjetsNorm((w_Control_base + cut_options_wCon)*cut_weights_HPTDown, (w_Control_base + cut_options_wCon)*cut_weights_HPTDown, (wSS_Control_base + cut_options_wCon)*cut_weights_HPTDown, (wSS_Control_base_RelTauIso + cut_options_wCon)*cut_weights_HPTDown, (w_Control_base + cut_options_wCon)*cut_weights_nom, (wSS_Control_base + cut_options_wCon)*cut_weights_nom, 0, "mt_tot", binnum_mt_tot, mt_totVarBinning , "eleTauWjetsOS", "et", "_CMS_xtt_highTauEffi_13TeVDown", 1);
        
        wSF_OS_eleTau_TESUp = wjetsNorm((w_Control_base + cut_options_wCon_TESUp)*cut_weights_nom, (w_Control_base + cut_options_wCon_TESUp)*cut_weights_nom, (wSS_Control_base + cut_options_wCon_TESUp)*cut_weights_nom, (wSS_Control_base_RelTauIso + cut_options_wCon_TESUp)*cut_weights_nom, (w_Control_base + cut_options_wCon)*cut_weights_nom, (wSS_Control_base + cut_options_wCon)*cut_weights_nom, 0, "mt_tot_TESUp", binnum_mt_tot, mt_totVarBinning , "eleTauWjetsOS", "et", "_CMS_scale_t_13TeVUp", 1);
        wSF_OS_eleTau_TESDown = wjetsNorm((w_Control_base + cut_options_wCon_TESDown)*cut_weights_nom, (w_Control_base + cut_options_wCon_TESDown)*cut_weights_nom, (wSS_Control_base + cut_options_wCon_TESDown)*cut_weights_nom, (wSS_Control_base_RelTauIso + cut_options_wCon_TESDown)*cut_weights_nom, (w_Control_base + cut_options_wCon)*cut_weights_nom, (wSS_Control_base + cut_options_wCon)*cut_weights_nom, 0, "mt_tot_TESDown", binnum_mt_tot, mt_totVarBinning , "eleTauWjetsOS", "et", "_CMS_scale_t_13TeVDown", 1);
        
        wSF_OS_eleTau_dm0TESUp = wjetsNorm((w_Control_base + cut_options_wCon_dm0TESUp)*cut_weights_nom, (w_Control_base + cut_options_wCon_dm0TESUp)*cut_weights_nom, (wSS_Control_base + cut_options_wCon_dm0TESUp)*cut_weights_nom, (wSS_Control_base_RelTauIso + cut_options_wCon_dm0TESUp)*cut_weights_nom, (w_Control_base + cut_options_wCon)*cut_weights_nom, (wSS_Control_base + cut_options_wCon)*cut_weights_nom, 0, "mt_tot_dm0TESUp", binnum_mt_tot, mt_totVarBinning , "eleTauWjetsOS", "et", "_CMS_scale_t_dm0_13TeVUp", 1);
        wSF_OS_eleTau_dm0TESDown = wjetsNorm((w_Control_base + cut_options_wCon_dm0TESDown)*cut_weights_nom, (w_Control_base + cut_options_wCon_dm0TESDown)*cut_weights_nom, (wSS_Control_base + cut_options_wCon_dm0TESDown)*cut_weights_nom, (wSS_Control_base_RelTauIso + cut_options_wCon_dm0TESDown)*cut_weights_nom, (w_Control_base + cut_options_wCon)*cut_weights_nom, (wSS_Control_base + cut_options_wCon)*cut_weights_nom, 0, "mt_tot_dm0TESDown", binnum_mt_tot, mt_totVarBinning , "eleTauWjetsOS", "et", "_CMS_scale_t_dm0_13TeVDown", 1);
        
        wSF_OS_eleTau_dm1TESUp = wjetsNorm((w_Control_base + cut_options_wCon_dm1TESUp)*cut_weights_nom, (w_Control_base + cut_options_wCon_dm1TESUp)*cut_weights_nom, (wSS_Control_base + cut_options_wCon_dm1TESUp)*cut_weights_nom, (wSS_Control_base_RelTauIso + cut_options_wCon_dm1TESUp)*cut_weights_nom, (w_Control_base + cut_options_wCon)*cut_weights_nom, (wSS_Control_base + cut_options_wCon)*cut_weights_nom, 0, "mt_tot_dm1TESUp", binnum_mt_tot, mt_totVarBinning , "eleTauWjetsOS", "et", "_CMS_scale_t_dm1_13TeVUp", 1);
        wSF_OS_eleTau_dm1TESDown = wjetsNorm((w_Control_base + cut_options_wCon_dm1TESDown)*cut_weights_nom, (w_Control_base + cut_options_wCon_dm1TESDown)*cut_weights_nom, (wSS_Control_base + cut_options_wCon_dm1TESDown)*cut_weights_nom, (wSS_Control_base_RelTauIso + cut_options_wCon_dm1TESDown)*cut_weights_nom, (w_Control_base + cut_options_wCon)*cut_weights_nom, (wSS_Control_base + cut_options_wCon)*cut_weights_nom, 0, "mt_tot_dm1TESDown", binnum_mt_tot, mt_totVarBinning , "eleTauWjetsOS", "et", "_CMS_scale_t_dm1_13TeVDown", 1);
        
        wSF_OS_eleTau_dm10TESUp = wjetsNorm((w_Control_base + cut_options_wCon_dm10TESUp)*cut_weights_nom, (w_Control_base + cut_options_wCon_dm10TESUp)*cut_weights_nom, (wSS_Control_base + cut_options_wCon_dm10TESUp)*cut_weights_nom, (wSS_Control_base_RelTauIso + cut_options_wCon_dm10TESUp)*cut_weights_nom, (w_Control_base + cut_options_wCon)*cut_weights_nom, (wSS_Control_base + cut_options_wCon)*cut_weights_nom, 0, "mt_tot_dm10TESUp", binnum_mt_tot, mt_totVarBinning , "eleTauWjetsOS", "et", "_CMS_scale_t_dm10_13TeVUp", 1);
        wSF_OS_eleTau_dm10TESDown = wjetsNorm((w_Control_base + cut_options_wCon_dm10TESDown)*cut_weights_nom, (w_Control_base + cut_options_wCon_dm10TESDown)*cut_weights_nom, (wSS_Control_base + cut_options_wCon_dm10TESDown)*cut_weights_nom, (wSS_Control_base_RelTauIso + cut_options_wCon_dm10TESDown)*cut_weights_nom, (w_Control_base + cut_options_wCon)*cut_weights_nom, (wSS_Control_base + cut_options_wCon)*cut_weights_nom, 0, "mt_tot_dm10TESDown", binnum_mt_tot, mt_totVarBinning , "eleTauWjetsOS", "et", "_CMS_scale_t_dm10_13TeVDown", 1);
        
        wSF_OS_eleTau_TopPtUp = wjetsNorm((w_Control_base + cut_options_wCon)*cut_weights_TopPtUp, (w_Control_base + cut_options_wCon)*cut_weights_TopPtUp, (wSS_Control_base + cut_options_wCon)*cut_weights_TopPtUp, (wSS_Control_base_RelTauIso + cut_options_wCon)*cut_weights_TopPtUp, (w_Control_base + cut_options_wCon)*cut_weights_nom, (wSS_Control_base + cut_options_wCon)*cut_weights_nom, 0, "mt_tot", binnum_mt_tot, mt_totVarBinning , "eleTauWjetsOS", "et", "_CMS_xtt_ttbarShape_13TeVUp", 1);
        wSF_OS_eleTau_TopPtDown = wjetsNorm((w_Control_base + cut_options_wCon)*cut_weights_TopPtDown, (w_Control_base + cut_options_wCon)*cut_weights_TopPtDown, (wSS_Control_base + cut_options_wCon)*cut_weights_TopPtDown, (wSS_Control_base_RelTauIso + cut_options_wCon)*cut_weights_TopPtDown, (w_Control_base + cut_options_wCon)*cut_weights_nom, (wSS_Control_base + cut_options_wCon)*cut_weights_nom, 0, "mt_tot", binnum_mt_tot, mt_totVarBinning , "eleTauWjetsOS", "et", "_CMS_xtt_ttbarShape_13TeVDown", 1);
        
        wSF_OS_eleTau_ZZnloUp = wjetsNorm((w_Control_base + cut_options_wCon)*cut_weights_ZZUp, (w_Control_base + cut_options_wCon)*cut_weights_ZZUp, (wSS_Control_base + cut_options_wCon)*cut_weights_ZZUp, (wSS_Control_base_RelTauIso + cut_options_wCon)*cut_weights_ZZUp, (w_Control_base + cut_options_wCon)*cut_weights_nom, (wSS_Control_base + cut_options_wCon)*cut_weights_nom, 0, "mt_tot", binnum_mt_tot, mt_totVarBinning , "eleTauWjetsOS", "et", "_CMS_xtt_ZZNLOewk_13TeVUp", 1);
        wSF_OS_eleTau_ZZnloDown = wjetsNorm((w_Control_base + cut_options_wCon)*cut_weights_ZZDown, (w_Control_base + cut_options_wCon)*cut_weights_ZZDown, (wSS_Control_base + cut_options_wCon)*cut_weights_ZZDown, (wSS_Control_base_RelTauIso + cut_options_wCon)*cut_weights_ZZDown, (w_Control_base + cut_options_wCon)*cut_weights_nom, (wSS_Control_base + cut_options_wCon)*cut_weights_nom, 0, "mt_tot", binnum_mt_tot, mt_totVarBinning , "eleTauWjetsOS", "et", "_CMS_xtt_ZZNLOewk_13TeVDown", 1);
        
        wSF_OS_eleTau_WWnloUp = wjetsNorm((w_Control_base + cut_options_wCon)*cut_weights_WWUp, (w_Control_base + cut_options_wCon)*cut_weights_WWUp, (wSS_Control_base + cut_options_wCon)*cut_weights_WWUp, (wSS_Control_base_RelTauIso + cut_options_wCon)*cut_weights_WWUp, (w_Control_base + cut_options_wCon)*cut_weights_nom, (wSS_Control_base + cut_options_wCon)*cut_weights_nom, 0, "mt_tot", binnum_mt_tot, mt_totVarBinning , "eleTauWjetsOS", "et", "_CMS_xtt_WWNLOewk_13TeVUp", 1);
        wSF_OS_eleTau_WWnloDown = wjetsNorm((w_Control_base + cut_options_wCon)*cut_weights_WWDown, (w_Control_base + cut_options_wCon)*cut_weights_WWDown, (wSS_Control_base + cut_options_wCon)*cut_weights_WWDown, (wSS_Control_base_RelTauIso + cut_options_wCon)*cut_weights_WWDown, (w_Control_base + cut_options_wCon)*cut_weights_nom, (wSS_Control_base + cut_options_wCon)*cut_weights_nom, 0, "mt_tot", binnum_mt_tot, mt_totVarBinning , "eleTauWjetsOS", "et", "_CMS_xtt_WWNLOewk_13TeVDown", 1);
     }
    
    
    //double wSF_OS_eleTau = 0.90;
    
    if (sigChoice==1)
    {
        
        global_title = "QCD (Same Sign) Estimate Region";
        
        drawSignalRegion((signalCut_base + cut_options_nom)*cut_weights_nom, (rel_signalCut_base + cut_options_nom)*cut_weights_nom, (SScut_eleTau_base + cut_options_nom)*cut_weights_nom, (SScut_eleTau_base_RelTauIso + cut_options_nom)*cut_weights_nom, wSF_SS_eleTau, wSF_OS_eleTau, "met", metBinning, "sig region met", "et", "", 1, 1);
        reset_files();
    }
    else if (sigChoice==2)
    {

        global_title = "QCD (Same Sign) Estimate Region";
        
        drawSignalRegion((signalCut_base + TMVACutMap["ET_600"])*cut_weights_nom, (signalCut_base + TMVACutMap["ET_600"])*cut_weights_nom, (SScut_eleTau_base + TMVACutMap["ET_600"])*cut_weights_nom, (SScut_eleTau_base_RelTauIso + TMVACutMap["ET_600"])*cut_weights_nom, wSF_SS_eleTau, 1.0 , "mvaVar_et_MZP600A0400", tmvaBinning, "sig region Function TMVA variable", "et", "600", 1, 1);
        drawSignalRegion((signalCut_base + TMVACutMap["ET_800"])*cut_weights_nom, (signalCut_base + TMVACutMap["ET_800"])*cut_weights_nom, (SScut_eleTau_base + TMVACutMap["ET_800"])*cut_weights_nom, (SScut_eleTau_base_RelTauIso + TMVACutMap["ET_800"])*cut_weights_nom, wSF_SS_eleTau, 1.0 , "mvaVar_et_MZP800A0400", tmvaBinning, "sig region Function TMVA variable", "et", "800", 1, 1);
        drawSignalRegion((signalCut_base + TMVACutMap["ET_1000"])*cut_weights_nom, (signalCut_base + TMVACutMap["ET_1000"])*cut_weights_nom, (SScut_eleTau_base + TMVACutMap["ET_1000"])*cut_weights_nom, (SScut_eleTau_base_RelTauIso + TMVACutMap["ET_1000"])*cut_weights_nom, wSF_SS_eleTau, 1.0 , "mvaVar_et_MZP1000A0400", tmvaBinning, "sig region Function TMVA variable", "et", "1000", 1, 1);
        drawSignalRegion((signalCut_base + TMVACutMap["ET_1200"])*cut_weights_nom, (signalCut_base + TMVACutMap["ET_1200"])*cut_weights_nom, (SScut_eleTau_base + TMVACutMap["ET_1200"])*cut_weights_nom, (SScut_eleTau_base_RelTauIso + TMVACutMap["ET_1200"])*cut_weights_nom, wSF_SS_eleTau, 1.0 , "mvaVar_et_MZP1200A0400", tmvaBinning, "sig region Function TMVA variable", "et", "1200", 1, 1);
        
        if (doSyst)
        {
        
        }

        
    }
    else if (sigChoice==3)
    {
        
        global_title = "QCD (Same Sign) Estimate Region";

        drawSignalRegionVarBin((signalCut_base + cut_options_nom)*cut_weights_nom, (rel_signalCut_base + cut_options_nom)*cut_weights_nom, (SScut_eleTau_base + cut_options_nom)*cut_weights_nom, (SScut_eleTau_base_RelTauIso + cut_options_nom)*cut_weights_nom, wSF_SS_eleTau, (signalCut_base + cut_options_nom)*cut_weights_nom, (SScut_eleTau_base + cut_options_nom)*cut_weights_nom, 0, wSF_OS_eleTau, "mt_tot", binnum_mt_tot, mt_totVarBinning, "sig region Mt_total", "et", "", 1, 1);
        
        if (printSyst) countTotalsVarBin((signalCut_base + cut_options_nom)*cut_weights_nom, (rel_signalCut_base + cut_options_nom)*cut_weights_nom, (SScut_eleTau_base + cut_options_nom)*cut_weights_nom, (SScut_eleTau_base_RelTauIso + cut_options_nom)*cut_weights_nom, wSF_SS_eleTau, wSF_OS_eleTau, "mt_tot", binnum_mt_tot, mt_totVarBinning, "Mt_total", "et", "NOM");
        
        if (doSyst)
        {
            std::cout << "syst for et started" << std::endl;
            //syst
            
            drawSignalRegionVarBin((signalCut_base + cut_options_UESUp)*cut_weights_UESUp, (rel_signalCut_base + cut_options_UESUp)*cut_weights_UESUp, (SScut_eleTau_base + cut_options_UESUp)*cut_weights_UESUp, (SScut_eleTau_base_RelTauIso + cut_options_UESUp)*cut_weights_UESUp, wSF_SS_eleTau, (signalCut_base + cut_options_nom)*cut_weights_nom, (SScut_eleTau_base + cut_options_nom)*cut_weights_nom, 0, wSF_OS_eleTau_UESUp, "mt_tot_UESUp", binnum_mt_tot, mt_totVarBinning,  "", "et", "_CMS_scale_m_13TeVUp", 0, 1);
            drawSignalRegionVarBin((signalCut_base + cut_options_UESDown)*cut_weights_UESDown, (rel_signalCut_base + cut_options_UESDown)*cut_weights_UESDown, (SScut_eleTau_base + cut_options_UESDown)*cut_weights_UESDown, (SScut_eleTau_base_RelTauIso + cut_options_UESDown)*cut_weights_UESDown, wSF_SS_eleTau, (signalCut_base + cut_options_nom)*cut_weights_nom, (SScut_eleTau_base + cut_options_nom)*cut_weights_nom, 0, wSF_OS_eleTau_UESDown, "mt_tot_UESDown", binnum_mt_tot, mt_totVarBinning,  "", "et", "_CMS_scale_m_13TeVDown", 0, 1);
            drawSignalRegionVarBin((signalCut_base + cut_options_nom)*cut_weights_ZPtUp, (rel_signalCut_base + cut_options_nom)*cut_weights_ZPtUp, (SScut_eleTau_base + cut_options_nom)*cut_weights_ZPtUp, (SScut_eleTau_base_RelTauIso + cut_options_nom)*cut_weights_ZPtUp, wSF_SS_eleTau, (signalCut_base + cut_options_nom)*cut_weights_nom, (SScut_eleTau_base + cut_options_nom)*cut_weights_nom, 0, wSF_OS_eleTau_ZPtUp, "mt_tot", binnum_mt_tot, mt_totVarBinning,  "", "et", "_CMS_xtt_dyShape_13TeVUp", 0, 1);
            drawSignalRegionVarBin((signalCut_base + cut_options_nom)*cut_weights_ZPtDown, (rel_signalCut_base + cut_options_nom)*cut_weights_ZPtDown, (SScut_eleTau_base + cut_options_nom)*cut_weights_ZPtDown, (SScut_eleTau_base_RelTauIso + cut_options_nom)*cut_weights_ZPtDown, wSF_SS_eleTau, (signalCut_base + cut_options_nom)*cut_weights_nom, (SScut_eleTau_base + cut_options_nom)*cut_weights_nom, 0, wSF_OS_eleTau_ZPtDown, "mt_tot", binnum_mt_tot, mt_totVarBinning,  "", "et", "_CMS_xtt_dyShape_13TeVDown", 0, 1);
            drawSignalRegionVarBin((signalCut_base + cut_options_nom)*cut_weights_KPtUp, (rel_signalCut_base + cut_options_nom)*cut_weights_KPtUp, (SScut_eleTau_base + cut_options_nom)*cut_weights_KPtUp, (SScut_eleTau_base_RelTauIso + cut_options_nom)*cut_weights_KPtUp, wSF_SS_eleTau, (signalCut_base + cut_options_nom)*cut_weights_nom, (SScut_eleTau_base + cut_options_nom)*cut_weights_nom, 0, wSF_OS_eleTau_KPtUp, "mt_tot", binnum_mt_tot, mt_totVarBinning,  "", "et", "_CMS_xtt_wShape_13TeVUp", 0, 1);
            drawSignalRegionVarBin((signalCut_base + cut_options_nom)*cut_weights_KPtDown, (rel_signalCut_base + cut_options_nom)*cut_weights_KPtDown, (SScut_eleTau_base + cut_options_nom)*cut_weights_KPtDown, (SScut_eleTau_base_RelTauIso + cut_options_nom)*cut_weights_KPtDown, wSF_SS_eleTau, (signalCut_base + cut_options_nom)*cut_weights_nom, (SScut_eleTau_base + cut_options_nom)*cut_weights_nom, 0, wSF_OS_eleTau_KPtDown, "mt_tot", binnum_mt_tot, mt_totVarBinning,  "", "et", "_CMS_xtt_wShape_13TeVDown", 0, 1);
            drawSignalRegionVarBin((signalCut_base + cut_options_TauFakeUp)*cut_weights_TauFakeUp, (rel_signalCut_base + cut_options_TauFakeUp)*cut_weights_TauFakeUp, (SScut_eleTau_base + cut_options_TauFakeUp)*cut_weights_TauFakeUp, (SScut_eleTau_base_RelTauIso + cut_options_TauFakeUp)*cut_weights_TauFakeUp, wSF_SS_eleTau, (signalCut_base + cut_options_nom)*cut_weights_nom, (SScut_eleTau_base + cut_options_nom)*cut_weights_nom, 1, wSF_OS_eleTau_TauFakeUp, "mt_tot", binnum_mt_tot, mt_totVarBinning,  "", "et", "_CMS_xtt_jetToTauFake_13TeVUp", 0, 1);
            drawSignalRegionVarBin((signalCut_base + cut_options_TauFakeDown)*cut_weights_TauFakeDown, (rel_signalCut_base + cut_options_TauFakeDown)*cut_weights_TauFakeDown, (SScut_eleTau_base + cut_options_TauFakeDown)*cut_weights_TauFakeDown, (SScut_eleTau_base_RelTauIso + cut_options_TauFakeDown)*cut_weights_TauFakeDown, wSF_SS_eleTau, (signalCut_base + cut_options_nom)*cut_weights_nom, (SScut_eleTau_base + cut_options_nom)*cut_weights_nom, 1, wSF_OS_eleTau_TauFakeDown, "mt_tot", binnum_mt_tot, mt_totVarBinning,  "", "et", "_CMS_xtt_jetToTauFake_13TeVDown", 0, 1);
            drawSignalRegionVarBin((signalCut_base + cut_options_nom)*cut_weights_nom, (rel_signalCut_base + cut_options_nom)*cut_weights_nom, (SScut_eleTau_base + cut_options_nom)*cut_weights_nom, (SScut_eleTau_base_RelTauIso + cut_options_nom)*cut_weights_nom, wSF_SS_eleTau, (signalCut_base + cut_options_nom)*cut_weights_nom, (SScut_eleTau_base + cut_options_nom)*cut_weights_nom, 0, wSF_OS_eleTau_ZLUp, "mt_tot*0.98", binnum_mt_tot, mt_totVarBinning,  "", "et", "_CMS_xtt_ZLScale_etau_13TeVUp", 0, 1);
            drawSignalRegionVarBin((signalCut_base + cut_options_nom)*cut_weights_nom, (rel_signalCut_base + cut_options_nom)*cut_weights_nom, (SScut_eleTau_base + cut_options_nom)*cut_weights_nom, (SScut_eleTau_base_RelTauIso + cut_options_nom)*cut_weights_nom, wSF_SS_eleTau, (signalCut_base + cut_options_nom)*cut_weights_nom, (SScut_eleTau_base + cut_options_nom)*cut_weights_nom, 0, wSF_OS_eleTau_ZLDown, "mt_tot*1.02", binnum_mt_tot, mt_totVarBinning,  "", "et", "_CMS_xtt_ZLScale_etau_13TeVDown", 0, 1);
            drawSignalRegionVarBin((signalCut_base + cut_options_bTagUp)*cut_weights_bTagUp, (rel_signalCut_base + cut_options_bTagUp)*cut_weights_bTagUp, (SScut_eleTau_base + cut_options_bTagUp)*cut_weights_bTagUp, (SScut_eleTau_base_RelTauIso + cut_options_bTagUp)*cut_weights_bTagUp, wSF_SS_eleTau, (signalCut_base + cut_options_nom)*cut_weights_nom, (SScut_eleTau_base + cut_options_nom)*cut_weights_nom, 0, wSF_OS_eleTau_bTagUp, "mt_tot", binnum_mt_tot, mt_totVarBinning,  "", "et", "_CMS_xtt_bTag_13TeVUp", 0, 1);
            drawSignalRegionVarBin((signalCut_base + cut_options_bTagDown)*cut_weights_bTagDown, (rel_signalCut_base + cut_options_bTagDown)*cut_weights_bTagDown, (SScut_eleTau_base + cut_options_bTagDown)*cut_weights_bTagDown, (SScut_eleTau_base_RelTauIso + cut_options_bTagDown)*cut_weights_bTagDown, wSF_SS_eleTau, (signalCut_base + cut_options_nom)*cut_weights_nom, (SScut_eleTau_base + cut_options_nom)*cut_weights_nom, 0, wSF_OS_eleTau_bTagDown, "mt_tot", binnum_mt_tot, mt_totVarBinning,  "", "et", "_CMS_xtt_bTag_13TeVDown", 0, 1);
            drawSignalRegionVarBin((signalCut_base + cut_options_JECUp)*cut_weights_JECUp, (rel_signalCut_base + cut_options_JECUp)*cut_weights_JECUp, (SScut_eleTau_base + cut_options_JECUp)*cut_weights_JECUp, (SScut_eleTau_base_RelTauIso + cut_options_JECUp)*cut_weights_JECUp, wSF_SS_eleTau, (signalCut_base + cut_options_nom)*cut_weights_nom, (SScut_eleTau_base + cut_options_nom)*cut_weights_nom, 0, wSF_OS_eleTau_JECUp, "mt_tot_JEnUp", binnum_mt_tot, mt_totVarBinning,  "", "et", "_CMS_scale_j_13TeVUp", 0, 1);
            drawSignalRegionVarBin((signalCut_base + cut_options_JECDown)*cut_weights_JECDown, (rel_signalCut_base + cut_options_JECDown)*cut_weights_JECDown, (SScut_eleTau_base + cut_options_JECDown)*cut_weights_JECDown, (SScut_eleTau_base_RelTauIso + cut_options_JECDown)*cut_weights_JECDown, wSF_SS_eleTau, (signalCut_base + cut_options_nom)*cut_weights_nom, (SScut_eleTau_base + cut_options_nom)*cut_weights_nom, 0, wSF_OS_eleTau_JECDown, "mt_tot_JEnDown", binnum_mt_tot, mt_totVarBinning,  "", "et","_CMS_scale_j_13TeVDown", 0, 1);
            drawSignalRegionVarBin((signalCut_base + cut_options_HPTUp)*cut_weights_HPTUp, (rel_signalCut_base + cut_options_HPTUp)*cut_weights_HPTUp, (SScut_eleTau_base + cut_options_HPTUp)*cut_weights_HPTUp, (SScut_eleTau_base_RelTauIso + cut_options_HPTUp)*cut_weights_HPTUp, wSF_SS_eleTau, (signalCut_base + cut_options_nom)*cut_weights_nom, (SScut_eleTau_base + cut_options_nom)*cut_weights_nom, 0, wSF_OS_eleTau_HPTUp, "mt_tot", binnum_mt_tot, mt_totVarBinning,  "", "et", "_CMS_xtt_highTauEffi_13TeVUp", 0, 1);
            drawSignalRegionVarBin((signalCut_base + cut_options_HPTDown)*cut_weights_HPTDown, (rel_signalCut_base + cut_options_HPTDown)*cut_weights_HPTDown, (SScut_eleTau_base + cut_options_HPTDown)*cut_weights_HPTDown, (SScut_eleTau_base_RelTauIso + cut_options_HPTDown)*cut_weights_HPTDown, wSF_SS_eleTau, (signalCut_base + cut_options_nom)*cut_weights_nom, (SScut_eleTau_base + cut_options_nom)*cut_weights_nom, 0, wSF_OS_eleTau_HPTDown, "mt_tot", binnum_mt_tot, mt_totVarBinning,  "", "et", "_CMS_xtt_highTauEffi_13TeVDown", 0, 1);
            drawSignalRegionVarBin((signalCut_base + cut_options_TESUp)*cut_weights_nom, (rel_signalCut_base + cut_options_TESUp)*cut_weights_nom, (SScut_eleTau_base + cut_options_TESUp)*cut_weights_nom, (SScut_eleTau_base_RelTauIso + cut_options_TESUp)*cut_weights_nom, wSF_SS_eleTau, (signalCut_base + cut_options_nom)*cut_weights_nom, (SScut_eleTau_base + cut_options_nom)*cut_weights_nom, 0, wSF_OS_eleTau_TESUp, "mt_tot_TESUp", binnum_mt_tot, mt_totVarBinning,  "", "et", "_CMS_scale_t_13TeVUp", 0, 1);
            drawSignalRegionVarBin((signalCut_base + cut_options_TESDown)*cut_weights_nom, (rel_signalCut_base + cut_options_TESUp)*cut_weights_nom, (SScut_eleTau_base + cut_options_TESUp)*cut_weights_nom, (SScut_eleTau_base_RelTauIso + cut_options_TESUp)*cut_weights_nom, wSF_SS_eleTau, (signalCut_base + cut_options_nom)*cut_weights_nom, (SScut_eleTau_base + cut_options_nom)*cut_weights_nom, 0, wSF_OS_eleTau_TESDown, "mt_tot_TESDown", binnum_mt_tot, mt_totVarBinning,  "", "et","_CMS_scale_t_13TeVDown", 0, 1);
            drawSignalRegionVarBin((signalCut_base + cut_options_dm0TESUp)*cut_weights_nom, (rel_signalCut_base + cut_options_dm0TESUp)*cut_weights_nom, (SScut_eleTau_base + cut_options_dm0TESUp)*cut_weights_nom, (SScut_eleTau_base_RelTauIso + cut_options_dm0TESUp)*cut_weights_nom, wSF_SS_eleTau, (signalCut_base + cut_options_nom)*cut_weights_nom, (SScut_eleTau_base + cut_options_nom)*cut_weights_nom, 0, wSF_OS_eleTau_dm0TESUp, "mt_tot_dm0TESUp", binnum_mt_tot, mt_totVarBinning,  "", "et", "_CMS_scale_t_dm0_13TeVUp", 0, 1);
            drawSignalRegionVarBin((signalCut_base + cut_options_dm0TESDown)*cut_weights_nom, (rel_signalCut_base + cut_options_dm0TESUp)*cut_weights_nom, (SScut_eleTau_base + cut_options_dm0TESUp)*cut_weights_nom, (SScut_eleTau_base_RelTauIso + cut_options_dm0TESUp)*cut_weights_nom, wSF_SS_eleTau, (signalCut_base + cut_options_nom)*cut_weights_nom, (SScut_eleTau_base + cut_options_nom)*cut_weights_nom, 0, wSF_OS_eleTau_dm0TESDown, "mt_tot_dm0TESDown", binnum_mt_tot, mt_totVarBinning,  "", "et","_CMS_scale_t_dm0_13TeVDown", 0, 1);
            drawSignalRegionVarBin((signalCut_base + cut_options_dm1TESUp)*cut_weights_nom, (rel_signalCut_base + cut_options_dm1TESUp)*cut_weights_nom, (SScut_eleTau_base + cut_options_dm1TESUp)*cut_weights_nom, (SScut_eleTau_base_RelTauIso + cut_options_dm1TESUp)*cut_weights_nom, wSF_SS_eleTau, (signalCut_base + cut_options_nom)*cut_weights_nom, (SScut_eleTau_base + cut_options_nom)*cut_weights_nom, 0, wSF_OS_eleTau_dm1TESUp, "mt_tot_dm1TESUp", binnum_mt_tot, mt_totVarBinning,  "", "et", "_CMS_scale_t_dm1_13TeVUp", 0, 1);
            drawSignalRegionVarBin((signalCut_base + cut_options_dm1TESDown)*cut_weights_nom, (rel_signalCut_base + cut_options_dm1TESUp)*cut_weights_nom, (SScut_eleTau_base + cut_options_dm1TESUp)*cut_weights_nom, (SScut_eleTau_base_RelTauIso + cut_options_dm1TESUp)*cut_weights_nom, wSF_SS_eleTau, (signalCut_base + cut_options_nom)*cut_weights_nom, (SScut_eleTau_base + cut_options_nom)*cut_weights_nom, 0, wSF_OS_eleTau_dm1TESDown, "mt_tot_dm1TESDown", binnum_mt_tot, mt_totVarBinning,  "", "et","_CMS_scale_t_dm1_13TeVDown", 0, 1);
            drawSignalRegionVarBin((signalCut_base + cut_options_dm10TESUp)*cut_weights_nom, (rel_signalCut_base + cut_options_dm10TESUp)*cut_weights_nom, (SScut_eleTau_base + cut_options_dm10TESUp)*cut_weights_nom, (SScut_eleTau_base_RelTauIso + cut_options_dm10TESUp)*cut_weights_nom, wSF_SS_eleTau, (signalCut_base + cut_options_nom)*cut_weights_nom, (SScut_eleTau_base + cut_options_nom)*cut_weights_nom, 0, wSF_OS_eleTau_dm10TESUp, "mt_tot_dm10TESUp", binnum_mt_tot, mt_totVarBinning,  "", "et", "_CMS_scale_t_dm10_13TeVUp", 0, 1);
            drawSignalRegionVarBin((signalCut_base + cut_options_dm10TESDown)*cut_weights_nom, (rel_signalCut_base + cut_options_dm10TESUp)*cut_weights_nom, (SScut_eleTau_base + cut_options_dm10TESUp)*cut_weights_nom, (SScut_eleTau_base_RelTauIso + cut_options_dm10TESUp)*cut_weights_nom, wSF_SS_eleTau, (signalCut_base + cut_options_nom)*cut_weights_nom, (SScut_eleTau_base + cut_options_nom)*cut_weights_nom, 0, wSF_OS_eleTau_dm10TESDown, "mt_tot_dm10TESDown", binnum_mt_tot, mt_totVarBinning,  "", "et","_CMS_scale_t_dm10_13TeVDown", 0, 1);
            drawSignalRegionVarBin((signalCut_base + cut_options_nom)*cut_weights_TopPtUp, (rel_signalCut_base + cut_options_nom)*cut_weights_TopPtUp, (SScut_eleTau_base + cut_options_nom)*cut_weights_TopPtUp, (SScut_eleTau_base_RelTauIso + cut_options_nom)*cut_weights_TopPtUp, wSF_SS_eleTau, (signalCut_base + cut_options_nom)*cut_weights_nom, (SScut_eleTau_base + cut_options_nom)*cut_weights_nom, 0, wSF_OS_eleTau_TopPtUp, "mt_tot", binnum_mt_tot, mt_totVarBinning,  "", "et", "_CMS_xtt_ttbarShape_13TeVUp", 0, 1);
            drawSignalRegionVarBin((signalCut_base + cut_options_nom)*cut_weights_TopPtDown, (rel_signalCut_base + cut_options_nom)*cut_weights_TopPtDown, (SScut_eleTau_base + cut_options_nom)*cut_weights_TopPtDown, (SScut_eleTau_base_RelTauIso + cut_options_nom)*cut_weights_TopPtDown, wSF_SS_eleTau,(signalCut_base + cut_options_nom)*cut_weights_nom, (SScut_eleTau_base + cut_options_nom)*cut_weights_nom, 0, wSF_OS_eleTau_TopPtDown, "mt_tot", binnum_mt_tot, mt_totVarBinning,  "", "et","_CMS_xtt_ttbarShape_13TeVDown", 0, 1);
            drawSignalRegionVarBin((signalCut_base + cut_options_nom)*cut_weights_WWUp, (rel_signalCut_base + cut_options_nom)*cut_weights_WWUp, (SScut_eleTau_base + cut_options_nom)*cut_weights_WWUp, (SScut_eleTau_base_RelTauIso + cut_options_nom)*cut_weights_WWUp, wSF_SS_eleTau, (signalCut_base + cut_options_nom)*cut_weights_WWUp, (rel_signalCut_base + cut_options_nom)*cut_weights_WWUp, 0, wSF_OS_eleTau, "mt_tot", binnum_mt_tot, mt_totVarBinning,  "", "et", "_CMS_xtt_WWNLOewk_13TeVUp", 0, 1);
            drawSignalRegionVarBin((signalCut_base + cut_options_nom)*cut_weights_WWDown, (rel_signalCut_base + cut_options_nom)*cut_weights_WWDown, (SScut_eleTau_base + cut_options_nom)*cut_weights_WWDown, (SScut_eleTau_base_RelTauIso + cut_options_nom)*cut_weights_WWDown, wSF_SS_eleTau, (signalCut_base + cut_options_nom)*cut_weights_nom, (SScut_eleTau_base + cut_options_nom)*cut_weights_nom, 0, wSF_OS_eleTau, "mt_tot", binnum_mt_tot, mt_totVarBinning,  "", "et","_CMS_xtt_WWNLOewk_13TeVDown", 0, 1);
            drawSignalRegionVarBin((signalCut_base + cut_options_nom)*cut_weights_ZZUp, (rel_signalCut_base + cut_options_nom)*cut_weights_ZZUp, (SScut_eleTau_base + cut_options_nom)*cut_weights_ZZUp, (SScut_eleTau_base_RelTauIso + cut_options_nom)*cut_weights_ZZUp, wSF_SS_eleTau, (signalCut_base + cut_options_nom)*cut_weights_nom, (SScut_eleTau_base + cut_options_nom)*cut_weights_nom, 0, wSF_OS_eleTau, "mt_tot", binnum_mt_tot, mt_totVarBinning,  "", "et", "_CMS_xtt_ZZNLOewk_13TeVUp", 0, 1);
            drawSignalRegionVarBin((signalCut_base + cut_options_nom)*cut_weights_ZZDown, (rel_signalCut_base + cut_options_nom)*cut_weights_ZZDown, (SScut_eleTau_base + cut_options_nom)*cut_weights_ZZDown, (SScut_eleTau_base_RelTauIso + cut_options_nom)*cut_weights_ZZDown, wSF_SS_eleTau, (signalCut_base + cut_options_nom)*cut_weights_nom, (SScut_eleTau_base + cut_options_nom)*cut_weights_nom, 0, wSF_OS_eleTau, "mt_tot", binnum_mt_tot, mt_totVarBinning,  "", "et","_CMS_xtt_ZZNLOewk_13TeVDown", 0, 1);
            
            
            if(printSyst)
            {
            
                countTotalsVarBin((signalCut_base + cut_options_UESUp)*cut_weights_UESUp, (rel_signalCut_base + cut_options_UESUp)*cut_weights_UESUp, (SScut_eleTau_base + cut_options_UESUp)*cut_weights_UESUp, (SScut_eleTau_base_RelTauIso + cut_options_UESUp)*cut_weights_UESUp, wSF_SS_eleTau, wSF_OS_eleTau, "mt_tot_UESUp", binnum_mt_tot, mt_totVarBinning, "Mt_total", "et", "_CMS_scale_m_13TeVUp");
                countTotalsVarBin((signalCut_base + cut_options_UESDown)*cut_weights_UESDown, (rel_signalCut_base + cut_options_UESDown)*cut_weights_UESDown, (SScut_eleTau_base + cut_options_UESDown)*cut_weights_UESDown, (SScut_eleTau_base_RelTauIso + cut_options_UESDown)*cut_weights_UESDown, wSF_SS_eleTau, wSF_OS_eleTau, "mt_tot_UESDown", binnum_mt_tot, mt_totVarBinning, "Mt_total", "et", "_CMS_scale_m_13TeVDown");
                countTotalsVarBin((signalCut_base + cut_options_nom)*cut_weights_ZPtUp, (rel_signalCut_base + cut_options_nom)*cut_weights_ZPtUp, (SScut_eleTau_base + cut_options_nom)*cut_weights_ZPtUp, (SScut_eleTau_base_RelTauIso + cut_options_nom)*cut_weights_ZPtUp, wSF_SS_eleTau, wSF_OS_eleTau, "mt_tot", binnum_mt_tot, mt_totVarBinning, "Mt_total", "et", "_CMS_xtt_dyShape_13TeVUp");
                countTotalsVarBin((signalCut_base + cut_options_nom)*cut_weights_ZPtDown, (rel_signalCut_base + cut_options_nom)*cut_weights_ZPtDown, (SScut_eleTau_base + cut_options_nom)*cut_weights_ZPtDown, (SScut_eleTau_base_RelTauIso + cut_options_nom)*cut_weights_ZPtDown, wSF_SS_eleTau, wSF_OS_eleTau, "mt_tot", binnum_mt_tot, mt_totVarBinning, "Mt_total", "et", "_CMS_xtt_dyShape_13TeVDown");
                countTotalsVarBin((signalCut_base + cut_options_nom)*cut_weights_KPtUp, (rel_signalCut_base + cut_options_nom)*cut_weights_KPtUp, (SScut_eleTau_base + cut_options_nom)*cut_weights_KPtUp, (SScut_eleTau_base_RelTauIso + cut_options_nom)*cut_weights_KPtUp, wSF_SS_eleTau, wSF_OS_eleTau, "mt_tot", binnum_mt_tot, mt_totVarBinning, "Mt_total", "et", "_CMS_xtt_wShape_13TeVUp");
                countTotalsVarBin((signalCut_base + cut_options_nom)*cut_weights_KPtDown, (rel_signalCut_base + cut_options_nom)*cut_weights_KPtDown, (SScut_eleTau_base + cut_options_nom)*cut_weights_KPtDown, (SScut_eleTau_base_RelTauIso + cut_options_nom)*cut_weights_KPtDown, wSF_SS_eleTau, wSF_OS_eleTau, "mt_tot", binnum_mt_tot, mt_totVarBinning, "Mt_total", "et", "_CMS_xtt_wShape_13TeVDown");
                countTotalsVarBin((signalCut_base + cut_options_TauFakeUp)*cut_weights_TauFakeUp, (rel_signalCut_base + cut_options_TauFakeUp)*cut_weights_TauFakeUp, (SScut_eleTau_base + cut_options_TauFakeUp)*cut_weights_TauFakeUp, (SScut_eleTau_base_RelTauIso + cut_options_TauFakeUp)*cut_weights_TauFakeUp, wSF_SS_eleTau, wSF_OS_eleTau, "mt_tot", binnum_mt_tot, mt_totVarBinning, "Mt_total", "et", "_CMS_xtt_jetToTauFake_13TeVUp");
                countTotalsVarBin((signalCut_base + cut_options_TauFakeDown)*cut_weights_TauFakeDown, (rel_signalCut_base + cut_options_TauFakeDown)*cut_weights_TauFakeDown, (SScut_eleTau_base + cut_options_TauFakeDown)*cut_weights_TauFakeDown, (SScut_eleTau_base_RelTauIso + cut_options_TauFakeDown)*cut_weights_TauFakeDown, wSF_SS_eleTau, wSF_OS_eleTau, "mt_tot", binnum_mt_tot, mt_totVarBinning, "Mt_total", "et", "_CMS_xtt_jetToTauFake_13TeVDown");
                countTotalsVarBin((signalCut_base + cut_options_nom)*cut_weights_nom, (rel_signalCut_base + cut_options_nom)*cut_weights_nom, (SScut_eleTau_base + cut_options_nom)*cut_weights_nom, (SScut_eleTau_base_RelTauIso + cut_options_nom)*cut_weights_nom, wSF_SS_eleTau, wSF_OS_eleTau, "mt_tot*0.98", binnum_mt_tot, mt_totVarBinning, "Mt_total", "et", "_CMS_xtt_ZLScale_etau_13TeVUp");
                countTotalsVarBin((signalCut_base + cut_options_nom)*cut_weights_nom, (rel_signalCut_base + cut_options_nom)*cut_weights_nom, (SScut_eleTau_base + cut_options_nom)*cut_weights_nom, (SScut_eleTau_base_RelTauIso + cut_options_nom)*cut_weights_nom, wSF_SS_eleTau, wSF_OS_eleTau, "mt_tot*1.02", binnum_mt_tot, mt_totVarBinning, "Mt_total", "et", "_CMS_xtt_ZLScale_etau_13TeVDown");
                countTotalsVarBin((signalCut_base + cut_options_bTagUp)*cut_weights_bTagUp, (rel_signalCut_base + cut_options_bTagUp)*cut_weights_bTagUp, (SScut_eleTau_base + cut_options_bTagUp)*cut_weights_bTagUp, (SScut_eleTau_base_RelTauIso + cut_options_bTagUp)*cut_weights_bTagUp, wSF_SS_eleTau, wSF_OS_eleTau, "mt_tot", binnum_mt_tot, mt_totVarBinning, "Mt_total", "et", "_CMS_xtt_bTag_13TeVUp");
                countTotalsVarBin((signalCut_base + cut_options_bTagDown)*cut_weights_bTagDown, (rel_signalCut_base + cut_options_bTagDown)*cut_weights_bTagDown, (SScut_eleTau_base + cut_options_bTagDown)*cut_weights_bTagDown, (SScut_eleTau_base_RelTauIso + cut_options_bTagDown)*cut_weights_bTagDown, wSF_SS_eleTau, wSF_OS_eleTau, "mt_tot", binnum_mt_tot, mt_totVarBinning, "Mt_total", "et", "_CMS_xtt_bTag_13TeVDown");
                countTotalsVarBin((signalCut_base + cut_options_JECUp)*cut_weights_JECUp, (rel_signalCut_base + cut_options_JECUp)*cut_weights_JECUp, (SScut_eleTau_base + cut_options_JECUp)*cut_weights_JECUp, (SScut_eleTau_base_RelTauIso + cut_options_JECUp)*cut_weights_JECUp, wSF_SS_eleTau, wSF_OS_eleTau, "mt_tot_JEnUp", binnum_mt_tot, mt_totVarBinning, "Mt_total", "et", "_CMS_scale_j_13TeVUp");
                countTotalsVarBin((signalCut_base + cut_options_JECDown)*cut_weights_JECDown, (rel_signalCut_base + cut_options_JECDown)*cut_weights_JECDown, (SScut_eleTau_base + cut_options_JECDown)*cut_weights_JECDown, (SScut_eleTau_base_RelTauIso + cut_options_JECDown)*cut_weights_JECDown, wSF_SS_eleTau, wSF_OS_eleTau, "mt_tot_JEnDown", binnum_mt_tot, mt_totVarBinning, "Mt_total", "et", "_CMS_scale_j_13TeVDown");
                countTotalsVarBin((signalCut_base + cut_options_HPTUp)*cut_weights_HPTUp, (rel_signalCut_base + cut_options_HPTUp)*cut_weights_HPTUp, (SScut_eleTau_base + cut_options_HPTUp)*cut_weights_HPTUp, (SScut_eleTau_base_RelTauIso + cut_options_HPTUp)*cut_weights_HPTUp, wSF_SS_eleTau, wSF_OS_eleTau, "mt_tot", binnum_mt_tot, mt_totVarBinning, "Mt_total", "et", "_CMS_xtt_highTauEffi_13TeVUp");
                countTotalsVarBin((signalCut_base + cut_options_HPTDown)*cut_weights_HPTDown, (rel_signalCut_base + cut_options_HPTDown)*cut_weights_HPTDown, (SScut_eleTau_base + cut_options_HPTDown)*cut_weights_HPTDown, (SScut_eleTau_base_RelTauIso + cut_options_HPTDown)*cut_weights_HPTDown, wSF_SS_eleTau, wSF_OS_eleTau, "mt_tot", binnum_mt_tot, mt_totVarBinning, "Mt_total", "et", "_CMS_xtt_highTauEffi_13TeVDown");
                countTotalsVarBin((signalCut_base + cut_options_TESUp)*cut_weights_nom, (rel_signalCut_base + cut_options_TESUp)*cut_weights_nom, (SScut_eleTau_base + cut_options_nom)*cut_weights_nom, (SScut_eleTau_base_RelTauIso + cut_options_nom)*cut_weights_nom, wSF_SS_eleTau, wSF_OS_eleTau, "mt_tot", binnum_mt_tot, mt_totVarBinning, "Mt_total", "et", "_CMS_scale_t_13TeVUp");
                countTotalsVarBin((signalCut_base + cut_options_TESDown)*cut_weights_nom, (rel_signalCut_base + cut_options_TESDown)*cut_weights_nom, (SScut_eleTau_base + cut_options_nom)*cut_weights_nom, (SScut_eleTau_base_RelTauIso + cut_options_nom)*cut_weights_nom, wSF_SS_eleTau, wSF_OS_eleTau, "mt_tot", binnum_mt_tot, mt_totVarBinning, "Mt_total", "et", "_CMS_scale_t_13TeVDown");
                countTotalsVarBin((signalCut_base + cut_options_nom)*cut_weights_TopPtUp, (rel_signalCut_base + cut_options_nom)*cut_weights_TopPtUp, (SScut_eleTau_base + cut_options_nom)*cut_weights_TopPtUp, (SScut_eleTau_base_RelTauIso + cut_options_nom)*cut_weights_TopPtUp, wSF_SS_eleTau, wSF_OS_eleTau, "mt_tot", binnum_mt_tot, mt_totVarBinning, "Mt_total", "et", "_CMS_xtt_ttbarShape_13TeVUp");
                countTotalsVarBin((signalCut_base + cut_options_nom)*cut_weights_TopPtDown, (rel_signalCut_base + cut_options_nom)*cut_weights_TopPtDown, (SScut_eleTau_base + cut_options_nom)*cut_weights_TopPtDown, (SScut_eleTau_base_RelTauIso + cut_options_nom)*cut_weights_TopPtDown, wSF_SS_eleTau, wSF_OS_eleTau, "mt_tot", binnum_mt_tot, mt_totVarBinning, "Mt_total", "et", "_CMS_xtt_ttbarShape_13TeVDown");
                countTotalsVarBin((signalCut_base + cut_options_nom)*cut_weights_WWUp, (rel_signalCut_base + cut_options_nom)*cut_weights_WWUp, (SScut_eleTau_base + cut_options_nom)*cut_weights_WWUp, (SScut_eleTau_base_RelTauIso + cut_options_nom)*cut_weights_WWUp, wSF_SS_eleTau, wSF_OS_eleTau, "mt_tot", binnum_mt_tot, mt_totVarBinning, "Mt_total", "et", "_CMS_xtt_WWNLOewk_13TeVUp");
                countTotalsVarBin((signalCut_base + cut_options_nom)*cut_weights_WWDown, (rel_signalCut_base + cut_options_nom)*cut_weights_WWDown, (SScut_eleTau_base + cut_options_nom)*cut_weights_WWDown, (SScut_eleTau_base_RelTauIso + cut_options_nom)*cut_weights_WWDown, wSF_SS_eleTau, wSF_OS_eleTau, "mt_tot", binnum_mt_tot, mt_totVarBinning, "Mt_total", "et", "_CMS_xtt_WWNLOewk_13TeVDown");
                countTotalsVarBin((signalCut_base + cut_options_nom)*cut_weights_ZZUp, (rel_signalCut_base + cut_options_nom)*cut_weights_ZZUp, (SScut_eleTau_base + cut_options_nom)*cut_weights_ZZUp, (SScut_eleTau_base_RelTauIso + cut_options_nom)*cut_weights_ZZUp, wSF_SS_eleTau, wSF_OS_eleTau, "mt_tot", binnum_mt_tot, mt_totVarBinning, "Mt_total", "et", "_CMS_xtt_ZZNLOewk_13TeVUp");
                countTotalsVarBin((signalCut_base + cut_options_nom)*cut_weights_ZZDown, (rel_signalCut_base + cut_options_nom)*cut_weights_ZZDown, (SScut_eleTau_base + cut_options_nom)*cut_weights_ZZDown, (SScut_eleTau_base_RelTauIso + cut_options_nom)*cut_weights_ZZDown, wSF_SS_eleTau, wSF_OS_eleTau, "mt_tot", binnum_mt_tot, mt_totVarBinning, "Mt_total", "et", "_CMS_xtt_ZZNLOewk_13TeVDown");
                
            }
            
        }
        reset_files();
    }
    else if (choice2==2)
    {
        
        if (doControlTMVA)
        {
            
            global_title = "MET > " + control_met_cutoff + " GeV";
            drawSignalRegion((signalCut_base + cut_options_nom)*cut_weights_nom, (rel_signalCut_base + cut_options_nom)*cut_weights_nom, (SScut_eleTau_base + cut_options_nom)*cut_weights_nom, (SScut_eleTau_base_RelTauIso + cut_options_nom)*cut_weights_nom, wSF_SS_eleTau, 1.0 , "mvaVar_et_MZP600A0400", tmvaBinning, "sig region Function TMVA variable", "et", "", 1, 0);
            drawSignalRegion((signalCut_base + cut_options_nom)*cut_weights_nom, (rel_signalCut_base + cut_options_nom)*cut_weights_nom, (SScut_eleTau_base + cut_options_nom)*cut_weights_nom, (SScut_eleTau_base_RelTauIso + cut_options_nom)*cut_weights_nom, wSF_SS_eleTau, 1.0 , "mvaVar_et_MZP800A0400", tmvaBinning, "sig region Function TMVA variable", "et", "", 1, 0);
            drawSignalRegion((signalCut_base + cut_options_nom)*cut_weights_nom, (rel_signalCut_base + cut_options_nom)*cut_weights_nom, (SScut_eleTau_base + cut_options_nom)*cut_weights_nom, (SScut_eleTau_base_RelTauIso + cut_options_nom)*cut_weights_nom, wSF_SS_eleTau, 1.0 , "mvaVar_et_MZP1000A0400", tmvaBinning, "sig region Function TMVA variable", "et", "", 1, 0);
            drawSignalRegion((signalCut_base + cut_options_nom)*cut_weights_nom, (rel_signalCut_base + cut_options_nom)*cut_weights_nom, (SScut_eleTau_base + cut_options_nom)*cut_weights_nom, (SScut_eleTau_base_RelTauIso + cut_options_nom)*cut_weights_nom, wSF_SS_eleTau, 1.0 , "mvaVar_et_MZP1200A0400", tmvaBinning, "sig region Function TMVA variable", "et", "", 1, 0);

        }
        else
        {

            countTotalsVarBin((signalCut_base + cut_options_nom)*cut_weights_nom, (rel_signalCut_base + cut_options_nom)*cut_weights_nom, (SScut_eleTau_base + cut_options_nom)*cut_weights_nom, (SScut_eleTau_base_RelTauIso + cut_options_nom)*cut_weights_nom, wSF_SS_eleTau, wSF_OS_eleTau, "mt_tot", binnum_mt_tot, mt_totVarBinning, "control region Mt_total", "et", "");

            global_title = "MET > " + control_met_cutoff + " GeV";
            drawSignalRegion((signalCut_base + cut_options_nom)*cut_weights_nom, (rel_signalCut_base + cut_options_nom)*cut_weights_nom, (SScut_eleTau_base + cut_options_nom)*cut_weights_nom, (SScut_eleTau_base_RelTauIso + cut_options_nom)*cut_weights_nom, wSF_SS_eleTau, wSF_OS_eleTau, "mt_1", mtBinning, "Control mt", "et", "", 0, 1);
            drawSignalRegion((signalCut_base + cut_options_nom)*cut_weights_nom, (rel_signalCut_base + cut_options_nom)*cut_weights_nom, (SScut_eleTau_base + cut_options_nom)*cut_weights_nom, (SScut_eleTau_base_RelTauIso + cut_options_nom)*cut_weights_nom, wSF_SS_eleTau, wSF_OS_eleTau, "m_vis", m_visBinning, "Control mvis", "et", "", 0, 1);
            drawSignalRegion((signalCut_base + cut_options_nom)*cut_weights_nom, (rel_signalCut_base + cut_options_nom)*cut_weights_nom, (SScut_eleTau_base + cut_options_nom)*cut_weights_nom, (SScut_eleTau_base_RelTauIso + cut_options_nom)*cut_weights_nom, wSF_SS_eleTau, wSF_OS_eleTau, "nbtag", nbtagBinning, "Control nbtag", "et", "", 0, 1);
            drawSignalRegion((signalCut_base + cut_options_nom)*cut_weights_nom, (rel_signalCut_base + cut_options_nom)*cut_weights_nom, (SScut_eleTau_base + cut_options_nom)*cut_weights_nom, (SScut_eleTau_base_RelTauIso + cut_options_nom)*cut_weights_nom, wSF_SS_eleTau, wSF_OS_eleTau, "met", metBinning, "Control met", "et", "", 0, 1);
            drawSignalRegion((signalCut_base + cut_options_nom)*cut_weights_nom, (rel_signalCut_base + cut_options_nom)*cut_weights_nom, (SScut_eleTau_base + cut_options_nom)*cut_weights_nom, (SScut_eleTau_base_RelTauIso + cut_options_nom)*cut_weights_nom, wSF_SS_eleTau, wSF_OS_eleTau, "pt_1", ptBinning, "Control pt_1", "et", "", 0, 1);
            drawSignalRegion((signalCut_base + cut_options_nom)*cut_weights_nom, (rel_signalCut_base + cut_options_nom)*cut_weights_nom, (SScut_eleTau_base + cut_options_nom)*cut_weights_nom, (SScut_eleTau_base_RelTauIso + cut_options_nom)*cut_weights_nom, wSF_SS_eleTau, wSF_OS_eleTau, "pt_2", ptBinning, "Control pt_2", "et", "", 0, 1);
            drawSignalRegion((signalCut_base + cut_options_nom)*cut_weights_nom, (rel_signalCut_base + cut_options_nom)*cut_weights_nom, (SScut_eleTau_base + cut_options_nom)*cut_weights_nom, (SScut_eleTau_base_RelTauIso + cut_options_nom)*cut_weights_nom, wSF_SS_eleTau, wSF_OS_eleTau, "pt_tt", ptBinning, "Control pt_tt", "et", "", 0, 1);
            
            drawSignalRegion((signalCut_base + cut_options_nom)*cut_weights_nom, (rel_signalCut_base + cut_options_nom)*cut_weights_nom, (SScut_eleTau_base + cut_options_nom)*cut_weights_nom, (SScut_eleTau_base_RelTauIso + cut_options_nom)*cut_weights_nom, wSF_SS_eleTau, wSF_OS_eleTau, "DeltaR_leg1_leg2", drBinning, "Control DeltaR_leg1_leg2", "et", "", 0, 1);
            drawSignalRegion((signalCut_base + cut_options_nom)*cut_weights_nom, (rel_signalCut_base + cut_options_nom)*cut_weights_nom, (SScut_eleTau_base + cut_options_nom)*cut_weights_nom, (SScut_eleTau_base_RelTauIso + cut_options_nom)*cut_weights_nom, wSF_SS_eleTau, wSF_OS_eleTau, "m_1", mBinning, "Control m_1", "et", "", 0, 1);
            drawSignalRegion((signalCut_base + cut_options_nom)*cut_weights_nom, (rel_signalCut_base + cut_options_nom)*cut_weights_nom, (SScut_eleTau_base + cut_options_nom)*cut_weights_nom, (SScut_eleTau_base_RelTauIso + cut_options_nom)*cut_weights_nom, wSF_SS_eleTau, wSF_OS_eleTau, "m_2", mBinning, "Control m_2", "et", "", 0, 1);
            drawSignalRegion((signalCut_base + cut_options_nom)*cut_weights_nom, (rel_signalCut_base + cut_options_nom)*cut_weights_nom, (SScut_eleTau_base + cut_options_nom)*cut_weights_nom, (SScut_eleTau_base_RelTauIso + cut_options_nom)*cut_weights_nom, wSF_SS_eleTau, wSF_OS_eleTau, "mt_2", mtBinning, "Control mt_2", "et", "", 0, 1);
            drawSignalRegion((signalCut_base + cut_options_nom)*cut_weights_nom, (rel_signalCut_base + cut_options_nom)*cut_weights_nom, (SScut_eleTau_base + cut_options_nom)*cut_weights_nom, (SScut_eleTau_base_RelTauIso + cut_options_nom)*cut_weights_nom, wSF_SS_eleTau, wSF_OS_eleTau, "phi_1", phiBinning, "Control phi_1", "et", "", 0, 1);
            drawSignalRegion((signalCut_base + cut_options_nom)*cut_weights_nom, (rel_signalCut_base + cut_options_nom)*cut_weights_nom, (SScut_eleTau_base + cut_options_nom)*cut_weights_nom, (SScut_eleTau_base_RelTauIso + cut_options_nom)*cut_weights_nom, wSF_SS_eleTau, wSF_OS_eleTau, "phi_2", phiBinning, "Control phi_2", "et", "", 0, 1);
            drawSignalRegion((signalCut_base + cut_options_nom)*cut_weights_nom, (rel_signalCut_base + cut_options_nom)*cut_weights_nom, (SScut_eleTau_base + cut_options_nom)*cut_weights_nom, (SScut_eleTau_base_RelTauIso + cut_options_nom)*cut_weights_nom, wSF_SS_eleTau, wSF_OS_eleTau, "eta_1", etaBinning, "Control eta_1", "et", "", 0, 1);
            drawSignalRegion((signalCut_base + cut_options_nom)*cut_weights_nom, (rel_signalCut_base + cut_options_nom)*cut_weights_nom, (SScut_eleTau_base + cut_options_nom)*cut_weights_nom, (SScut_eleTau_base_RelTauIso + cut_options_nom)*cut_weights_nom, wSF_SS_eleTau, wSF_OS_eleTau, "eta_2", etaBinning, "Control eta_2", "et", "", 0, 1);
            drawSignalRegion((signalCut_base + cut_options_nom)*cut_weights_nom, (rel_signalCut_base + cut_options_nom)*cut_weights_nom, (SScut_eleTau_base + cut_options_nom)*cut_weights_nom, (SScut_eleTau_base_RelTauIso + cut_options_nom)*cut_weights_nom, wSF_SS_eleTau, wSF_OS_eleTau, "njetspt20", nbtagBinning, "Control njetspt20", "et", "", 0, 1);
            
        }
        reset_files();
    }
    
}

void doMuTau()
{
    setBkgErr("mt");

    float mt_totVarBinning[] = {0,20,40,60,90,120,150,180,210,235,260,285,325,400,2000};
    int binnum_mt_tot = sizeof(mt_totVarBinning)/sizeof(*mt_totVarBinning) - 1;
    
    /*
    TCut SScut_muTau_base_RelTauIso("isOsPair==0 && pairGoodForTrigger==1 && extramuon_veto==0  && extraelec_veto==0  && dilepton_veto==0 && iso_1 < 0.15 && againstElectronVLooseMVA6_2 > 0.5 && againstMuonTight3_2 > 0.5 && byLooseIsolationMVArun2v1DBoldDMwLT_2 > 0.5 && decayModeFinding_2 > 0.5");
    
    TCut rel_signalCut_base("isOsPair==1 && pairGoodForTrigger==1 && extramuon_veto==0  && extraelec_veto==0  && dilepton_veto==0 && iso_1 < 0.15 && againstElectronVLooseMVA6_2 > 0.5 && againstMuonTight3_2 > 0.5 && byLooseIsolationMVArun2v1DBoldDMwLT_2 > 0.5 && decayModeFinding_2 > 0.5");
    */
    
    TCut SScut_muTau_base_RelTauIso("isOsPair==0 && pairGoodForTrigger==1 && extramuon_veto==0  && extraelec_veto==0  && dilepton_veto==0 && iso_1 < 0.30 && againstElectronVLooseMVA6_2 > 0.5 && againstMuonTight3_2 > 0.5 && decayModeFinding_2 > 0.5");
    
    TCut rel_signalCut_base("isOsPair==1 && pairGoodForTrigger==1 && extramuon_veto==0  && extraelec_veto==0  && dilepton_veto==0 && iso_1 < 0.30 && againstElectronVLooseMVA6_2 > 0.5 && againstMuonTight3_2 > 0.5 && decayModeFinding_2 > 0.5");
    
    TCut wSS_Control_base_RelTauIso("isOsPair==0 && pairGoodForTrigger==1 && extramuon_veto==0  && extraelec_veto==0  && dilepton_veto==0 && iso_1 < 0.30 && againstElectronVLooseMVA6_2 > 0.5 && againstMuonTight3_2 > 0.5 && byTightIsolationMVArun2v1DBoldDMwLT_2 < 0.5 && decayModeFinding_2 > 0.5");
    
    TCut SScut_muTau_base("isOsPair==0 && pairGoodForTrigger==1 && extramuon_veto==0  && extraelec_veto==0  && dilepton_veto==0 && iso_1 < 0.15 && againstElectronVLooseMVA6_2 > 0.5 && againstMuonTight3_2 > 0.5 && byTightIsolationMVArun2v1DBoldDMwLT_2 > 0.5 && decayModeFinding_2 > 0.5");
    
    TCut signalCut_base("isOsPair==1 && pairGoodForTrigger==1 && extramuon_veto==0  && extraelec_veto==0  && dilepton_veto==0 && iso_1 < 0.15 && againstElectronVLooseMVA6_2 > 0.5 && againstMuonTight3_2 > 0.5 && byTightIsolationMVArun2v1DBoldDMwLT_2 > 0.5 && decayModeFinding_2 > 0.5");
    
    TCut w_Control_base("isOsPair==1 && pairGoodForTrigger==1 && extramuon_veto==0  && extraelec_veto==0  && dilepton_veto==0 && iso_1 < 0.15 && againstElectronVLooseMVA6_2 > 0.5 && againstMuonTight3_2 > 0.5 && byTightIsolationMVArun2v1DBoldDMwLT_2 < 0.5 && decayModeFinding_2 > 0.5");
    
    TCut wSS_Control_base("isOsPair==0 && pairGoodForTrigger==1 && extramuon_veto==0  && extraelec_veto==0  && dilepton_veto==0 && iso_1 < 0.15 && againstElectronVLooseMVA6_2 > 0.5 && againstMuonTight3_2 > 0.5 && byTightIsolationMVArun2v1DBoldDMwLT_2 < 0.5 && decayModeFinding_2 > 0.5");
    

    std::string metCutTmp;
    if (cRegion)
    {
        if(doControlTMVA){metCutTmp = "met > 0.";}
        else{metCutTmp = "met > " + control_met_cutoff;}
    }
    else
    {
       metCutTmp = met_options["NOM"];
    }
    TCut cut_options_nom(("pt_1 > 26. && pt_2 > 20. && " + drCutMap[drCut] + " && " + metCutTmp + " && " + mt_options["NOM"] + " && pt_tt > 65. && mt_tot <2000.").c_str());
    TCut cut_weights_nom = jetCutMap["NOM"];
    
    TCut mt_totBlind("mt_tot < 200");
    
    //Systematic Cuts
    TCut cut_options_UESUp(("pt_1 > 26. && pt_2 > 20. && " + drCutMap[drCut] + " && " + met_options["MET_UES_UP"] + " && " + mt_options["NOM"] + " && pt_tt > 65. && mt_tot <2000.").c_str());
    TCut cut_options_UESDown(("pt_1 > 26. && pt_2 > 20. && " + drCutMap[drCut] + " && " + met_options["MET_UES_DOWN"] + " && " + mt_options["NOM"] + " && pt_tt > 65. && mt_tot <2000.").c_str());
    
    TCut cut_options_TESUp(("pt_1_TESUp > 26. && pt_2_TESUp > 20. && " + drCutMap[drCut] + " && " + met_options["TES_UP"] + " && " + mt_options["TES_UP"] + " && pt_tt_TESUp > 65. && mt_tot <2000.").c_str());
    TCut cut_options_TESDown(("pt_1_TESDown > 26. && pt_2_TESDown > 20. && " + drCutMap[drCut] + " && " + met_options["TES_DOWN"] + " && " + mt_options["TES_DOWN"] + " && pt_tt_TESDown > 65. && mt_tot <2000.").c_str());
    
    TCut cut_options_dm0TESUp(("pt_1_dm0TESUp > 26. && pt_2_dm0TESUp > 20. && " + drCutMap[drCut] + " && " + met_options["TES_DM0UP"] + " && " + mt_options["TES_DM0UP"] + " && pt_tt_dm0TESUp > 65. && mt_tot <2000.").c_str());
    TCut cut_options_dm0TESDown(("pt_1_dm0TESDown > 26. && pt_2_dm0TESDown > 20. && " + drCutMap[drCut] + " && " + met_options["TES_DM0DOWN"] + " && " + mt_options["TES_DM0DOWN"] + " && pt_tt_dm0TESDown > 65. && mt_tot <2000.").c_str());
    
    TCut cut_options_dm1TESUp(("pt_1_dm1TESUp > 26. && pt_2_dm1TESUp > 20. && " + drCutMap[drCut] + " && " + met_options["TES_DM1UP"] + " && " + mt_options["TES_DM1UP"] + " && pt_tt_dm1TESUp > 65. && mt_tot <2000.").c_str());
    TCut cut_options_dm1TESDown(("pt_1_dm1TESDown > 26. && pt_2_dm1TESDown > 20. && " + drCutMap[drCut] + " && " + met_options["TES_DM1DOWN"] + " && " + mt_options["TES_DM1DOWN"] + " && pt_tt_dm1TESDown > 65. && mt_tot <2000.").c_str());
    
    TCut cut_options_dm10TESUp(("pt_1_dm10TESUp > 26. && pt_2_dm10TESUp > 20. && " + drCutMap[drCut] + " && " + met_options["TES_DM10UP"] + " && " + mt_options["TES_DM10UP"] + " && pt_tt_dm10TESUp > 65. && mt_tot <2000.").c_str());
    TCut cut_options_dm10TESDown(("pt_1_dm10TESDown > 26. && pt_2_dm10TESDown > 20. && " + drCutMap[drCut] + " && " + met_options["TES_DM10DOWN"] + " && " + mt_options["TES_DM10DOWN"] + " && pt_tt_dm10TESDown > 65. && mt_tot <2000.").c_str());
    
    TCut cut_options_TauFakeUp(("pt_1 > 26. && pt_2 > 20. && " + drCutMap[drCut] + " && " + met_options["NOM"] + " && " + mt_options["NOM"] + " && pt_tt > 65. && mt_tot <2000.").c_str());
    TCut cut_options_TauFakeDown(("pt_1 > 26. && pt_2 > 20. && " + drCutMap[drCut] + " && " + met_options["NOM"] + " && " + mt_options["NOM"] + " && pt_tt > 65. && mt_tot <2000.").c_str());
    
    TCut cut_options_JECUp(("pt_1 > 26. && pt_2 > 20. && " + drCutMap[drCut] + " && " + met_options["JECUP"] + " && " + mt_options["NOM"] + " && pt_tt > 65. && mt_tot <2000.").c_str());
    TCut cut_options_JECDown(("pt_1 > 26. && pt_2 > 20. && " + drCutMap[drCut] + " && " + met_options["JECDOWN"] + " && " + mt_options["NOM"] + " && pt_tt > 65. && mt_tot <2000.").c_str());
    
    TCut cut_options_HPTUp(("pt_1 > 26. && pt_2 > 20. && " + drCutMap[drCut] + " && " + met_options["NOM"] + " && " + mt_options["NOM"] + " && pt_tt > 65. && mt_tot <2000.").c_str());
    TCut cut_options_HPTDown(("pt_1 > 26. && pt_2 > 20. && " + drCutMap[drCut] + " && " + met_options["NOM"] + " && " + mt_options["NOM"] + " && pt_tt > 65. && mt_tot <2000.").c_str());
    
    TCut cut_options_bTagUp(("pt_1 > 26. && pt_2 > 20. && " + drCutMap[drCut] + " && " + met_options["NOM"] + " && " + mt_options["NOM"] + " && pt_tt > 65. && mt_tot <2000.").c_str());
    TCut cut_options_bTagDown(("pt_1 > 26. && pt_2 > 20. && " + drCutMap[drCut] + " && " + met_options["NOM"] + " && " + mt_options["NOM"] + " && pt_tt > 65. && mt_tot <2000.").c_str());
    
    //Systematic cuts for WCon

    TCut cut_options_wCon(("pt_1 > 26. && pt_2 > 20. && DeltaR_leg1_leg2 > 0.3 && met > 105. && mt_1 > 80. && pt_tt > 65. && mt_tot <2000."));
    
    TCut cut_options_wCon_UESUp(("pt_1 > 26. && pt_2 > 20. && DeltaR_leg1_leg2 > 0.3 && pfmet_type1_UnclusteredEnUp_Pt > 105. && pfmet_type1_UnclusteredEnUp_MT1 > 80. && pt_tt > 65. && mt_tot <2000."));
    TCut cut_options_wCon_UESDown(("pt_1 > 26. && pt_2 > 20. && DeltaR_leg1_leg2 > 0.3 && pfmet_type1_UnclusteredEnDown_Pt > 105. && pfmet_type1_UnclusteredEnDown_MT1 > 80. && pt_tt > 65. && mt_tot <2000."));
    
    TCut cut_options_wCon_TESUp(("pt_1_TESUp > 26. && pt_2_TESUp > 20. && DeltaR_leg1_leg2 > 0.3 && pfmet_type1_TESUp_Pt > 105. && mt_1_TESUp > 80. && pt_tt_TESUp > 65. && mt_tot <2000."));
    TCut cut_options_wCon_TESDown(("pt_1_TESDown > 26. && pt_2_TESDown > 20. && DeltaR_leg1_leg2 > 0.3 && pfmet_type1_TESDown_Pt > 105. && mt_1_TESDown > 80. && pt_tt_TESDown > 65. && mt_tot <2000."));
    
    TCut cut_options_wCon_dm0TESUp(("pt_1_dm0TESUp > 26. && pt_2_dm0TESUp > 20. && DeltaR_leg1_leg2 > 0.3 && pfmet_type1_dm0TESUp_Pt > 105. && mt_1_dm0TESUp > 80. && pt_tt_dm0TESUp > 65. && mt_tot <2000."));
    TCut cut_options_wCon_dm0TESDown(("pt_1_dm0TESDown > 26. && pt_2_dm0TESDown > 20. && DeltaR_leg1_leg2 > 0.3 && pfmet_type1_dm0TESDown_Pt > 105. && mt_1_dm0TESDown > 80. && pt_tt_dm0TESDown > 65. && mt_tot <2000."));
    
    TCut cut_options_wCon_dm1TESUp(("pt_1_dm1TESUp > 26. && pt_2_dm1TESUp > 20. && DeltaR_leg1_leg2 > 0.3 && pfmet_type1_dm1TESUp_Pt > 105. && mt_1_dm1TESUp > 80. && pt_tt_dm1TESUp > 65. && mt_tot <2000."));
    TCut cut_options_wCon_dm1TESDown(("pt_1_dm1TESDown > 26. && pt_2_dm1TESDown > 20. && DeltaR_leg1_leg2 > 0.3 && pfmet_type1_dm1TESDown_Pt > 105. && mt_1_dm1TESDown > 80. && pt_tt_dm1TESDown > 65. && mt_tot <2000."));
    
    TCut cut_options_wCon_dm10TESUp(("pt_1_dm10TESUp > 26. && pt_2_dm10TESUp > 20. && DeltaR_leg1_leg2 > 0.3 && pfmet_type1_dm10TESUp_Pt > 105. && mt_1_dm10TESUp > 80. && pt_tt_dm10TESUp > 65. && mt_tot <2000."));
    TCut cut_options_wCon_dm10TESDown(("pt_1_dm10TESDown > 26. && pt_2_dm10TESDown > 20. && DeltaR_leg1_leg2 > 0.3 && pfmet_type1_dm10TESDown_Pt > 105. && mt_1_dm10TESDown > 80. && pt_tt_dm10TESDown > 65. && mt_tot <2000."));
    
    TCut cut_options_wCon_JECUp(("pt_1 > 26. && pt_2 > 20. && DeltaR_leg1_leg2 > 0.3 && pfmet_type1_JetEnUp_Pt > 105. && pfmt_1_JEnUp > 80. && pt_tt > 65. && mt_tot <2000."));
    TCut cut_options_wCon_JECDown(("pt_1 > 26. && pt_2 > 20. && DeltaR_leg1_leg2 > 0.3 && pfmet_type1_JetEnDown_Pt > 105. && pfmt_1_JEnDown > 80. && pt_tt > 65. && mt_tot <2000."));
    
    //Weights
    TCut cut_weights_UESUp = jetCutMap["NOM"];
    TCut cut_weights_UESDown = jetCutMap["NOM"];
    
    TCut cut_weights_ZPtUp = jetCutMap["NOM"]; cut_weights_ZPtUp *= weightMap["ZUP"];
    TCut cut_weights_ZPtDown = jetCutMap["NOM"]; cut_weights_ZPtDown *= weightMap["ZDOWN"];
    
    TCut cut_weights_KPtUp = jetCutMap["NOM"]; cut_weights_KPtUp *= weightMap["KUP"];
    TCut cut_weights_KPtDown = jetCutMap["NOM"]; cut_weights_KPtDown *= weightMap["KDOWN"];
    
    TCut cut_weights_TauFakeUp = jetCutMap["NOM"]; cut_weights_TauFakeUp *= weightMap["JTFUP"];
    TCut cut_weights_TauFakeDown = jetCutMap["NOM"]; cut_weights_TauFakeDown *= weightMap["JTFDOWN"];
    
    TCut cut_weights_JECUp = jetCutMap["BSFJECUP"];
    TCut cut_weights_JECDown = jetCutMap["BSFJECDOWN"];
    
    TCut cut_weights_HPTUp = jetCutMap["NOM"]; cut_weights_HPTUp *= weightMap["HPTUP"];
    TCut cut_weights_HPTDown = jetCutMap["NOM"]; cut_weights_HPTDown *= weightMap["HPTDOWN"];
    
    TCut cut_weights_bTagUp = jetCutMap["BSFUP"];
    TCut cut_weights_bTagDown = jetCutMap["BSFDOWN"];
    
    TCut cut_weights_TopPtUp = jetCutMap["NOM"]; cut_weights_TopPtUp *= weightMap["TOPPTUP"];
    TCut cut_weights_TopPtDown = jetCutMap["NOM"];
    
    TCut cut_weights_WWUp = jetCutMap["NOM"]; cut_weights_WWUp *= weightMap["WWUP"];
    TCut cut_weights_WWDown = jetCutMap["NOM"]; cut_weights_WWDown *= weightMap["WWDOWN"];
    
    TCut cut_weights_ZZUp = jetCutMap["NOM"]; cut_weights_ZZUp *= weightMap["ZZUP"];
    TCut cut_weights_ZZDown = jetCutMap["NOM"]; cut_weights_ZZDown *= weightMap["ZZDown"];

	chan_label = "#mu+ #tau_{h}";

    /* declare the muTau qcd shapes for each binning */

    
    setup_files_muTau();
    
    double wSF_SS_muTau = 1.0; /* temp */

    double wSF_OS_muTau = wjetsNorm((w_Control_base + cut_options_wCon)*cut_weights_nom, (w_Control_base + cut_options_wCon)*cut_weights_nom, (wSS_Control_base + cut_options_wCon)*cut_weights_nom, (wSS_Control_base_RelTauIso + cut_options_wCon)*cut_weights_nom, (w_Control_base + cut_options_wCon)*cut_weights_nom, (wSS_Control_base + cut_options_wCon)*cut_weights_nom, 0, "mt_tot", binnum_mt_tot, mt_totVarBinning , "muTauWjetsOS", "mt", "", 1);
    
    std::cout << " WJets SF MuTau: " << wSF_OS_muTau << std::endl;
    
    double wSF_OS_muTau_UESUp = 1.0;
    double wSF_OS_muTau_UESDown = 1.0;
    double wSF_OS_muTau_ZPtUp = 1.0;
    double wSF_OS_muTau_ZPtDown = 1.0;
    double wSF_OS_muTau_KPtUp = 1.0;
    double wSF_OS_muTau_KPtDown = 1.0;
    double wSF_OS_muTau_TauFakeUp = 1.0;
    double wSF_OS_muTau_TauFakeDown = 1.0;
    double wSF_OS_muTau_ZLUp = 1.0;
    double wSF_OS_muTau_ZLDown = 1.0;
    double wSF_OS_muTau_bTagUp = 1.0;
    double wSF_OS_muTau_bTagDown = 1.0;
    double wSF_OS_muTau_JECUp = 1.0;
    double wSF_OS_muTau_JECDown = 1.0;
    double wSF_OS_muTau_HPTUp = 1.0;
    double wSF_OS_muTau_HPTDown = 1.0;
    double wSF_OS_muTau_TESUp = 1.0;
    double wSF_OS_muTau_TESDown = 1.0;
    double wSF_OS_muTau_dm0TESUp = 1.0;
    double wSF_OS_muTau_dm0TESDown = 1.0;
    double wSF_OS_muTau_dm1TESUp = 1.0;
    double wSF_OS_muTau_dm1TESDown = 1.0;
    double wSF_OS_muTau_dm10TESUp = 1.0;
    double wSF_OS_muTau_dm10TESDown = 1.0;
    double wSF_OS_muTau_TopPtUp = 1.0;
    double wSF_OS_muTau_TopPtDown = 1.0;
    double wSF_OS_muTau_ZZnloUp = 1.0;
    double wSF_OS_muTau_ZZnloDown = 1.0;
    double wSF_OS_muTau_WWnloUp = 1.0;
    double wSF_OS_muTau_WWnloDown = 1.0;
    
    if (doSyst)
    {
        wSF_OS_muTau_UESUp = wjetsNorm((w_Control_base + cut_options_wCon_UESUp)*cut_weights_nom, (w_Control_base + cut_options_wCon_UESUp)*cut_weights_nom, (wSS_Control_base + cut_options_wCon_UESUp)*cut_weights_nom, (wSS_Control_base_RelTauIso + cut_options_wCon_UESUp)*cut_weights_nom, (w_Control_base + cut_options_wCon)*cut_weights_nom, (wSS_Control_base + cut_options_wCon)*cut_weights_nom, 0, "mt_tot_UESUp", binnum_mt_tot, mt_totVarBinning , "muTauWjetsOS_UESUp", "mt", "_CMS_scale_m_13TeVUp", 1);
        wSF_OS_muTau_UESDown = wjetsNorm((w_Control_base + cut_options_wCon_UESDown)*cut_weights_nom, (w_Control_base + cut_options_wCon_UESDown)*cut_weights_nom, (wSS_Control_base + cut_options_wCon_UESDown)*cut_weights_nom, (wSS_Control_base_RelTauIso + cut_options_wCon_UESDown)*cut_weights_nom, (w_Control_base + cut_options_wCon)*cut_weights_nom, (wSS_Control_base + cut_options_wCon)*cut_weights_nom, 0, "mt_tot_UESDown", binnum_mt_tot, mt_totVarBinning , "muTauWjetsOS_UESDown", "mt", "_CMS_scale_m_13TeVDown", 1);
        
        wSF_OS_muTau_ZPtUp = wjetsNorm((w_Control_base + cut_options_wCon)*cut_weights_ZPtUp, (w_Control_base + cut_options_wCon)*cut_weights_ZPtUp, (wSS_Control_base + cut_options_wCon)*cut_weights_ZPtUp, (wSS_Control_base_RelTauIso + cut_options_wCon)*cut_weights_ZPtUp, (w_Control_base + cut_options_wCon)*cut_weights_nom, (wSS_Control_base + cut_options_wCon)*cut_weights_nom, 0, "mt_tot", binnum_mt_tot, mt_totVarBinning , "muTauWjetsOS", "mt", "_CMS_xtt_dyShape_13TeVUp", 1);
        wSF_OS_muTau_ZPtDown = wjetsNorm((w_Control_base + cut_options_wCon)*cut_weights_ZPtDown, (w_Control_base + cut_options_wCon)*cut_weights_ZPtDown, (wSS_Control_base + cut_options_wCon)*cut_weights_ZPtDown, (wSS_Control_base_RelTauIso + cut_options_wCon)*cut_weights_ZPtDown, (w_Control_base + cut_options_wCon)*cut_weights_nom, (wSS_Control_base + cut_options_wCon)*cut_weights_nom, 0, "mt_tot", binnum_mt_tot, mt_totVarBinning , "muTauWjetsOS", "mt", "_CMS_xtt_dyShape_13TeVDown", 1);
        
        wSF_OS_muTau_KPtUp = wjetsNorm((w_Control_base + cut_options_wCon)*cut_weights_KPtUp, (w_Control_base + cut_options_wCon)*cut_weights_KPtUp, (wSS_Control_base + cut_options_wCon)*cut_weights_ZPtUp, (wSS_Control_base_RelTauIso + cut_options_wCon)*cut_weights_ZPtUp, (w_Control_base + cut_options_wCon)*cut_weights_nom, (wSS_Control_base + cut_options_wCon)*cut_weights_nom, 0, "mt_tot", binnum_mt_tot, mt_totVarBinning , "muTauWjetsOS", "mt", "_CMS_xtt_wShape_13TeVUp", 1);
        wSF_OS_muTau_KPtDown = wjetsNorm((w_Control_base + cut_options_wCon)*cut_weights_KPtDown, (w_Control_base + cut_options_wCon)*cut_weights_KPtDown, (wSS_Control_base + cut_options_wCon)*cut_weights_KPtDown, (wSS_Control_base_RelTauIso + cut_options_wCon)*cut_weights_KPtDown, (w_Control_base + cut_options_wCon)*cut_weights_nom, (wSS_Control_base + cut_options_wCon)*cut_weights_nom, 0, "mt_tot", binnum_mt_tot, mt_totVarBinning , "muTauWjetsOS", "mt", "_CMS_xtt_wShape_13TeVDown", 1);
        
        wSF_OS_muTau_TauFakeUp = wjetsNorm((w_Control_base + cut_options_wCon)*cut_weights_TauFakeUp, (w_Control_base + cut_options_wCon)*cut_weights_TauFakeUp, (wSS_Control_base + cut_options_wCon)*cut_weights_ZPtUp, (wSS_Control_base_RelTauIso + cut_options_wCon)*cut_weights_ZPtUp, (w_Control_base + cut_options_wCon)*cut_weights_nom, (wSS_Control_base + cut_options_wCon)*cut_weights_nom, 1, "mt_tot", binnum_mt_tot, mt_totVarBinning , "muTauWjetsOS", "mt", "_CMS_xtt_jetToTauFake_13TeVUp", 1);
        wSF_OS_muTau_TauFakeDown = wjetsNorm((w_Control_base + cut_options_wCon)*cut_weights_TauFakeDown, (w_Control_base + cut_options_wCon)*cut_weights_TauFakeDown, (wSS_Control_base + cut_options_wCon)*cut_weights_TauFakeDown, (wSS_Control_base_RelTauIso + cut_options_wCon)*cut_weights_TauFakeDown, (w_Control_base + cut_options_wCon)*cut_weights_nom, (wSS_Control_base + cut_options_wCon)*cut_weights_nom, 1, "mt_tot", binnum_mt_tot, mt_totVarBinning , "muTauWjetsOS", "mt", "_CMS_xtt_jetToTauFake_13TeVDown", 1);
        
        wSF_OS_muTau_ZLUp = wjetsNorm((w_Control_base + cut_options_wCon)*cut_weights_nom, (w_Control_base + cut_options_wCon)*cut_weights_nom, (wSS_Control_base + cut_options_wCon)*cut_weights_nom, (wSS_Control_base_RelTauIso + cut_options_wCon)*cut_weights_nom, (w_Control_base + cut_options_wCon)*cut_weights_nom, (wSS_Control_base + cut_options_wCon)*cut_weights_nom, 0, "mt_tot*0.98", binnum_mt_tot, mt_totVarBinning , "muTauWjetsOS", "mt", "_CMS_xtt_ZLScale_mutau_13TeVUp", 1);
        wSF_OS_muTau_ZLDown = wjetsNorm((w_Control_base + cut_options_wCon)*cut_weights_nom, (w_Control_base + cut_options_wCon)*cut_weights_nom, (wSS_Control_base + cut_options_wCon)*cut_weights_nom, (wSS_Control_base_RelTauIso + cut_options_wCon)*cut_weights_nom, (w_Control_base + cut_options_wCon)*cut_weights_nom, (wSS_Control_base + cut_options_wCon)*cut_weights_nom, 0, "mt_tot*1.02", binnum_mt_tot, mt_totVarBinning , "muTauWjetsOS", "mt", "_CMS_xtt_ZLScale_mutau_13TeVDown", 1);
        
        wSF_OS_muTau_bTagUp = wjetsNorm((w_Control_base + cut_options_wCon)*cut_weights_bTagUp, (w_Control_base + cut_options_wCon)*cut_weights_bTagUp, (wSS_Control_base + cut_options_wCon)*cut_weights_bTagUp, (wSS_Control_base_RelTauIso + cut_options_wCon)*cut_weights_bTagUp, (w_Control_base + cut_options_wCon)*cut_weights_nom, (wSS_Control_base + cut_options_wCon)*cut_weights_nom, 0, "mt_tot", binnum_mt_tot, mt_totVarBinning , "muTauWjetsOS", "mt", "_CMS_xtt_bTag_13TeVUp", 1);
        wSF_OS_muTau_bTagDown = wjetsNorm((w_Control_base + cut_options_wCon)*cut_weights_bTagDown, (w_Control_base + cut_options_wCon)*cut_weights_bTagDown, (wSS_Control_base + cut_options_wCon)*cut_weights_bTagDown, (wSS_Control_base_RelTauIso + cut_options_wCon)*cut_weights_bTagDown, (w_Control_base + cut_options_wCon)*cut_weights_nom, (wSS_Control_base + cut_options_wCon)*cut_weights_nom, 0, "mt_tot", binnum_mt_tot, mt_totVarBinning , "muTauWjetsOS", "mt", "_CMS_xtt_bTag_13TeVDown", 1);
        
        wSF_OS_muTau_JECUp = wjetsNorm((w_Control_base + cut_options_wCon_JECUp)*cut_weights_JECUp, (w_Control_base + cut_options_wCon_JECUp)*cut_weights_JECUp, (wSS_Control_base + cut_options_wCon_JECUp)*cut_weights_JECUp, (wSS_Control_base_RelTauIso + cut_options_wCon_JECUp)*cut_weights_JECUp, (w_Control_base + cut_options_wCon)*cut_weights_nom, (wSS_Control_base + cut_options_wCon)*cut_weights_nom, 0, "mt_tot_JEnUp", binnum_mt_tot, mt_totVarBinning , "muTauWjetsOS", "mt", "_CMS_scale_j_13TeVUp", 1);
        wSF_OS_muTau_JECDown = wjetsNorm((w_Control_base + cut_options_wCon_JECDown)*cut_weights_JECDown, (w_Control_base + cut_options_wCon_JECDown)*cut_weights_JECDown, (wSS_Control_base + cut_options_wCon_JECDown)*cut_weights_JECDown, (wSS_Control_base_RelTauIso + cut_options_wCon_JECDown)*cut_weights_JECDown, (w_Control_base + cut_options_wCon)*cut_weights_nom, (wSS_Control_base + cut_options_wCon)*cut_weights_nom, 0, "mt_tot_JEnDown", binnum_mt_tot, mt_totVarBinning , "muTauWjetsOS", "mt", "_CMS_scale_j_13TeVDown", 1);
        
        wSF_OS_muTau_HPTUp = wjetsNorm((w_Control_base + cut_options_wCon)*cut_weights_HPTUp, (w_Control_base + cut_options_wCon)*cut_weights_HPTUp, (wSS_Control_base + cut_options_wCon)*cut_weights_HPTUp, (wSS_Control_base_RelTauIso + cut_options_wCon)*cut_weights_HPTUp, (w_Control_base + cut_options_wCon)*cut_weights_nom, (wSS_Control_base + cut_options_wCon)*cut_weights_nom, 0, "mt_tot", binnum_mt_tot, mt_totVarBinning , "muTauWjetsOS", "mt", "_CMS_xtt_highTauEffi_13TeVUp", 1);
        wSF_OS_muTau_HPTDown = wjetsNorm((w_Control_base + cut_options_wCon)*cut_weights_HPTDown, (w_Control_base + cut_options_wCon)*cut_weights_HPTDown, (wSS_Control_base + cut_options_wCon)*cut_weights_HPTDown, (wSS_Control_base_RelTauIso + cut_options_wCon)*cut_weights_HPTDown, (w_Control_base + cut_options_wCon)*cut_weights_nom, (wSS_Control_base + cut_options_wCon)*cut_weights_nom, 0, "mt_tot", binnum_mt_tot, mt_totVarBinning , "muTauWjetsOS", "mt", "_CMS_xtt_highTauEffi_13TeVDown", 1);
        
        wSF_OS_muTau_TESUp = wjetsNorm((w_Control_base + cut_options_wCon_TESUp)*cut_weights_nom, (w_Control_base + cut_options_wCon_TESUp)*cut_weights_nom, (wSS_Control_base + cut_options_wCon_TESUp)*cut_weights_nom, (wSS_Control_base_RelTauIso + cut_options_wCon_TESUp)*cut_weights_nom, (w_Control_base + cut_options_wCon)*cut_weights_nom, (wSS_Control_base + cut_options_wCon)*cut_weights_nom, 0, "mt_tot_TESUp", binnum_mt_tot, mt_totVarBinning , "muTauWjetsOS", "mt", "_CMS_scale_t_13TeVUp", 1);
        wSF_OS_muTau_TESDown = wjetsNorm((w_Control_base + cut_options_wCon_TESDown)*cut_weights_nom, (w_Control_base + cut_options_wCon_TESDown)*cut_weights_nom, (wSS_Control_base + cut_options_wCon_TESDown)*cut_weights_nom, (wSS_Control_base_RelTauIso + cut_options_wCon_TESDown)*cut_weights_nom, (w_Control_base + cut_options_wCon)*cut_weights_nom, (wSS_Control_base + cut_options_wCon)*cut_weights_nom, 0, "mt_tot_TESDown", binnum_mt_tot, mt_totVarBinning , "muTauWjetsOS", "mt", "_CMS_scale_t_13TeVDown", 1);
        
        wSF_OS_muTau_dm0TESUp = wjetsNorm((w_Control_base + cut_options_wCon_dm0TESUp)*cut_weights_nom, (w_Control_base + cut_options_wCon_dm0TESUp)*cut_weights_nom, (wSS_Control_base + cut_options_wCon_dm0TESUp)*cut_weights_nom, (wSS_Control_base_RelTauIso + cut_options_wCon_dm0TESUp)*cut_weights_nom, (w_Control_base + cut_options_wCon)*cut_weights_nom, (wSS_Control_base + cut_options_wCon)*cut_weights_nom, 0, "mt_tot_dm0TESUp", binnum_mt_tot, mt_totVarBinning , "muTauWjetsOS", "mt", "_CMS_scale_t_dm0_13TeVUp", 1);
        wSF_OS_muTau_dm0TESDown = wjetsNorm((w_Control_base + cut_options_wCon_dm0TESDown)*cut_weights_nom, (w_Control_base + cut_options_wCon_dm0TESDown)*cut_weights_nom, (wSS_Control_base + cut_options_wCon_dm0TESDown)*cut_weights_nom, (wSS_Control_base_RelTauIso + cut_options_wCon_dm0TESDown)*cut_weights_nom, (w_Control_base + cut_options_wCon)*cut_weights_nom, (wSS_Control_base + cut_options_wCon)*cut_weights_nom, 0, "mt_tot_dm0TESDown", binnum_mt_tot, mt_totVarBinning , "muTauWjetsOS", "mt", "_CMS_scale_t_dm0_13TeVDown", 1);
        
        wSF_OS_muTau_dm1TESUp = wjetsNorm((w_Control_base + cut_options_wCon_dm1TESUp)*cut_weights_nom, (w_Control_base + cut_options_wCon_dm1TESUp)*cut_weights_nom, (wSS_Control_base + cut_options_wCon_dm1TESUp)*cut_weights_nom, (wSS_Control_base_RelTauIso + cut_options_wCon_dm1TESUp)*cut_weights_nom, (w_Control_base + cut_options_wCon)*cut_weights_nom, (wSS_Control_base + cut_options_wCon)*cut_weights_nom, 0, "mt_tot_dm1TESUp", binnum_mt_tot, mt_totVarBinning , "muTauWjetsOS", "mt", "_CMS_scale_t_dm1_13TeVUp", 1);
        wSF_OS_muTau_dm1TESDown = wjetsNorm((w_Control_base + cut_options_wCon_dm1TESDown)*cut_weights_nom, (w_Control_base + cut_options_wCon_dm1TESDown)*cut_weights_nom, (wSS_Control_base + cut_options_wCon_dm1TESDown)*cut_weights_nom, (wSS_Control_base_RelTauIso + cut_options_wCon_dm1TESDown)*cut_weights_nom, (w_Control_base + cut_options_wCon)*cut_weights_nom, (wSS_Control_base + cut_options_wCon)*cut_weights_nom, 0, "mt_tot_dm1TESDown", binnum_mt_tot, mt_totVarBinning , "muTauWjetsOS", "mt", "_CMS_scale_t_dm1_13TeVDown", 1);
        
        wSF_OS_muTau_dm10TESUp = wjetsNorm((w_Control_base + cut_options_wCon_dm10TESUp)*cut_weights_nom, (w_Control_base + cut_options_wCon_dm10TESUp)*cut_weights_nom, (wSS_Control_base + cut_options_wCon_dm10TESUp)*cut_weights_nom, (wSS_Control_base_RelTauIso + cut_options_wCon_dm10TESUp)*cut_weights_nom, (w_Control_base + cut_options_wCon)*cut_weights_nom, (wSS_Control_base + cut_options_wCon)*cut_weights_nom, 0, "mt_tot_dm10TESUp", binnum_mt_tot, mt_totVarBinning , "muTauWjetsOS", "mt", "_CMS_scale_t_dm10_13TeVUp", 1);
        wSF_OS_muTau_dm10TESDown = wjetsNorm((w_Control_base + cut_options_wCon_dm10TESDown)*cut_weights_nom, (w_Control_base + cut_options_wCon_dm10TESDown)*cut_weights_nom, (wSS_Control_base + cut_options_wCon_dm10TESDown)*cut_weights_nom, (wSS_Control_base_RelTauIso + cut_options_wCon_dm10TESDown)*cut_weights_nom, (w_Control_base + cut_options_wCon)*cut_weights_nom, (wSS_Control_base + cut_options_wCon)*cut_weights_nom, 0, "mt_tot_dm10TESDown", binnum_mt_tot, mt_totVarBinning , "muTauWjetsOS", "mt", "_CMS_scale_t_dm10_13TeVDown", 1);
        
        wSF_OS_muTau_TopPtUp = wjetsNorm((w_Control_base + cut_options_wCon)*cut_weights_TopPtUp, (w_Control_base + cut_options_wCon)*cut_weights_TopPtUp, (wSS_Control_base + cut_options_wCon)*cut_weights_TopPtUp, (wSS_Control_base_RelTauIso + cut_options_wCon)*cut_weights_TopPtUp, (w_Control_base + cut_options_wCon)*cut_weights_nom, (wSS_Control_base + cut_options_wCon)*cut_weights_nom, 0, "mt_tot", binnum_mt_tot, mt_totVarBinning , "muTauWjetsOS", "mt", "_CMS_xtt_ttbarShape_13TeVUp", 1);
        wSF_OS_muTau_TopPtDown = wjetsNorm((w_Control_base + cut_options_wCon)*cut_weights_TopPtDown, (w_Control_base + cut_options_wCon)*cut_weights_TopPtDown, (wSS_Control_base + cut_options_wCon)*cut_weights_TopPtDown, (wSS_Control_base_RelTauIso + cut_options_wCon)*cut_weights_TopPtDown, (w_Control_base + cut_options_wCon)*cut_weights_nom, (wSS_Control_base + cut_options_wCon)*cut_weights_nom, 0, "mt_tot", binnum_mt_tot, mt_totVarBinning , "muTauWjetsOS", "mt", "_CMS_xtt_ttbarShape_13TeVDown", 1);
        
        wSF_OS_muTau_ZZnloUp = wjetsNorm((w_Control_base + cut_options_wCon)*cut_weights_ZZUp, (w_Control_base + cut_options_wCon)*cut_weights_ZZUp, (wSS_Control_base + cut_options_wCon)*cut_weights_ZZUp, (wSS_Control_base_RelTauIso + cut_options_wCon)*cut_weights_ZZUp, (w_Control_base + cut_options_wCon)*cut_weights_nom, (wSS_Control_base + cut_options_wCon)*cut_weights_nom, 0, "mt_tot", binnum_mt_tot, mt_totVarBinning , "muTauWjetsOS", "mt", "_CMS_xtt_ZZNLOewk_13TeVUp", 1);
        wSF_OS_muTau_ZZnloDown = wjetsNorm((w_Control_base + cut_options_wCon)*cut_weights_ZZDown, (w_Control_base + cut_options_wCon)*cut_weights_ZZDown, (wSS_Control_base + cut_options_wCon)*cut_weights_ZZDown, (wSS_Control_base_RelTauIso + cut_options_wCon)*cut_weights_ZZDown, (w_Control_base + cut_options_wCon)*cut_weights_nom, (wSS_Control_base + cut_options_wCon)*cut_weights_nom, 0, "mt_tot", binnum_mt_tot, mt_totVarBinning , "muTauWjetsOS", "mt", "_CMS_xtt_ZZNLOewk_13TeVDown", 1);
        
        wSF_OS_muTau_WWnloUp = wjetsNorm((w_Control_base + cut_options_wCon)*cut_weights_WWUp, (w_Control_base + cut_options_wCon)*cut_weights_WWUp, (wSS_Control_base + cut_options_wCon)*cut_weights_WWUp, (wSS_Control_base_RelTauIso + cut_options_wCon)*cut_weights_WWUp, (w_Control_base + cut_options_wCon)*cut_weights_nom, (wSS_Control_base + cut_options_wCon)*cut_weights_nom, 0, "mt_tot", binnum_mt_tot, mt_totVarBinning , "muTauWjetsOS", "mt", "_CMS_xtt_WWNLOewk_13TeVUp", 1);
        wSF_OS_muTau_WWnloDown = wjetsNorm((w_Control_base + cut_options_wCon)*cut_weights_WWDown, (w_Control_base + cut_options_wCon)*cut_weights_WWDown, (wSS_Control_base + cut_options_wCon)*cut_weights_WWDown, (wSS_Control_base_RelTauIso + cut_options_wCon)*cut_weights_WWDown, (w_Control_base + cut_options_wCon)*cut_weights_nom, (wSS_Control_base + cut_options_wCon)*cut_weights_nom, 0, "mt_tot", binnum_mt_tot, mt_totVarBinning , "muTauWjetsOS", "mt", "_CMS_xtt_WWNLOewk_13TeVDown", 1);
    }
    
    
    //double wSF_OS_muTau = 0.92;
    
    if (sigChoice==1)
    {
        
        global_title = "QCD (Same Sign) Estimate Region";
        drawSignalRegion((signalCut_base + cut_options_nom)*cut_weights_nom, (rel_signalCut_base + cut_options_nom)*cut_weights_nom, (SScut_muTau_base + cut_options_nom)*cut_weights_nom, (SScut_muTau_base_RelTauIso + cut_options_nom)*cut_weights_nom, wSF_SS_muTau, wSF_OS_muTau, "met", metBinning, "sig region met", "mt", "", 1, 1);
        
    }
    else if (sigChoice==2)
    {

        drawSignalRegion((signalCut_base + TMVACutMap["MT_600"])*cut_weights_nom, (signalCut_base + TMVACutMap["MT_600"])*cut_weights_nom, (SScut_muTau_base +  TMVACutMap["MT_600"])*cut_weights_nom, (SScut_muTau_base_RelTauIso + TMVACutMap["MT_600"])*cut_weights_nom, wSF_SS_muTau, 1.0 , "mvaVar_mt_MZP600A0400", tmvaBinning, "sig region Function TMVA variable", "mt", "600", 1, 1);
        drawSignalRegion((signalCut_base + TMVACutMap["MT_800"])*cut_weights_nom, (signalCut_base + TMVACutMap["MT_800"])*cut_weights_nom, (SScut_muTau_base + TMVACutMap["MT_800"])*cut_weights_nom, (SScut_muTau_base_RelTauIso + TMVACutMap["MT_800"])*cut_weights_nom , wSF_SS_muTau, 1.0 , "mvaVar_mt_MZP800A0400", tmvaBinning, "sig region Function TMVA variable", "mt", "800", 1, 1);
        drawSignalRegion((signalCut_base + TMVACutMap["MT_1000"])*cut_weights_nom, (signalCut_base + TMVACutMap["MT_1000"])*cut_weights_nom, (SScut_muTau_base + TMVACutMap["MT_1000"])*cut_weights_nom, (SScut_muTau_base_RelTauIso + TMVACutMap["MT_1000"])*cut_weights_nom, wSF_SS_muTau, 1.0 , "mvaVar_mt_MZP1000A0400", tmvaBinning, "sig region Function TMVA variable", "mt", "1000", 1, 1);
        drawSignalRegion((signalCut_base + TMVACutMap["MT_1200"])*cut_weights_nom, (signalCut_base + TMVACutMap["MT_1200"])*cut_weights_nom, (SScut_muTau_base + TMVACutMap["MT_1200"])*cut_weights_nom, (SScut_muTau_base_RelTauIso + TMVACutMap["MT_1200"])*cut_weights_nom, wSF_SS_muTau, 1.0 , "mvaVar_mt_MZP1200A0400", tmvaBinning, "sig region Function TMVA variable", "mt", "1200", 1, 1);
        
        if (doSyst)
        {
        }

    }
    else if (sigChoice==3)
    {
        //findBkgFractions(signalCut_base * cut_options_nom, 1.0, "mt_tot", mt_totVarBinning);
        global_title = "QCD (Same Sign) Estimate Region";
        drawSignalRegionVarBin((signalCut_base + cut_options_nom)*cut_weights_nom, (rel_signalCut_base + cut_options_nom)*cut_weights_nom, (SScut_muTau_base + cut_options_nom)*cut_weights_nom, (SScut_muTau_base_RelTauIso + cut_options_nom)*cut_weights_nom, wSF_SS_muTau, (signalCut_base + cut_options_nom)*cut_weights_nom, (SScut_muTau_base + cut_options_nom)*cut_weights_nom, 0, wSF_OS_muTau, "mt_tot", binnum_mt_tot, mt_totVarBinning, "sig region Mt_total", "mt", "", 1, 1);
        
        if (printSyst) countTotalsVarBin((signalCut_base + cut_options_nom)*cut_weights_nom, (rel_signalCut_base + cut_options_nom)*cut_weights_nom, (SScut_muTau_base + cut_options_nom)*cut_weights_nom, (SScut_muTau_base_RelTauIso + cut_options_nom)*cut_weights_nom, wSF_SS_muTau, wSF_OS_muTau, "mt_tot", binnum_mt_tot, mt_totVarBinning, "Mt_total", "mt", "NOM");
        
        if (doSyst)
        {
            std::cout << "syst for mt started" << std::endl;
            
            drawSignalRegionVarBin((signalCut_base + cut_options_UESUp)*cut_weights_UESUp, (rel_signalCut_base + cut_options_UESUp)*cut_weights_UESUp, (SScut_muTau_base + cut_options_UESUp)*cut_weights_UESUp, (SScut_muTau_base_RelTauIso + cut_options_UESUp)*cut_weights_UESUp, wSF_SS_muTau, (signalCut_base + cut_options_nom)*cut_weights_nom, (SScut_muTau_base + cut_options_nom)*cut_weights_nom, 0, wSF_OS_muTau_UESUp, "mt_tot_UESUp", binnum_mt_tot, mt_totVarBinning,  "", "mt", "_CMS_scale_m_13TeVUp", 0, 1);
            drawSignalRegionVarBin((signalCut_base + cut_options_UESDown)*cut_weights_UESDown, (rel_signalCut_base + cut_options_UESDown)*cut_weights_UESDown, (SScut_muTau_base + cut_options_UESDown)*cut_weights_UESDown, (SScut_muTau_base_RelTauIso + cut_options_UESDown)*cut_weights_UESDown, wSF_SS_muTau, (signalCut_base + cut_options_nom)*cut_weights_nom, (SScut_muTau_base + cut_options_nom)*cut_weights_nom, 0, wSF_OS_muTau_UESDown, "mt_tot_UESDown", binnum_mt_tot, mt_totVarBinning,  "", "mt", "_CMS_scale_m_13TeVDown", 0, 1);
            drawSignalRegionVarBin((signalCut_base + cut_options_nom)*cut_weights_ZPtUp, (rel_signalCut_base + cut_options_nom)*cut_weights_ZPtUp, (SScut_muTau_base + cut_options_nom)*cut_weights_ZPtUp, (SScut_muTau_base_RelTauIso + cut_options_nom)*cut_weights_ZPtUp, wSF_SS_muTau, (signalCut_base + cut_options_nom)*cut_weights_nom, (SScut_muTau_base + cut_options_nom)*cut_weights_nom, 0, wSF_OS_muTau_ZPtUp, "mt_tot", binnum_mt_tot, mt_totVarBinning,  "", "mt", "_CMS_xtt_dyShape_13TeVUp", 0, 1);
            drawSignalRegionVarBin((signalCut_base + cut_options_nom)*cut_weights_ZPtDown, (rel_signalCut_base + cut_options_nom)*cut_weights_ZPtDown, (SScut_muTau_base + cut_options_nom)*cut_weights_ZPtDown, (SScut_muTau_base_RelTauIso + cut_options_nom)*cut_weights_ZPtDown, wSF_SS_muTau, (signalCut_base + cut_options_nom)*cut_weights_nom, (SScut_muTau_base + cut_options_nom)*cut_weights_nom, 0, wSF_OS_muTau_ZPtDown, "mt_tot", binnum_mt_tot, mt_totVarBinning,  "", "mt", "_CMS_xtt_dyShape_13TeVDown", 0, 1);
            
            drawSignalRegionVarBin((signalCut_base + cut_options_nom)*cut_weights_KPtUp, (rel_signalCut_base + cut_options_nom)*cut_weights_KPtUp, (SScut_muTau_base + cut_options_nom)*cut_weights_KPtUp, (SScut_muTau_base_RelTauIso + cut_options_nom)*cut_weights_KPtUp, wSF_SS_muTau, (signalCut_base + cut_options_nom)*cut_weights_nom, (SScut_muTau_base + cut_options_nom)*cut_weights_nom, 0, wSF_OS_muTau_KPtUp, "mt_tot", binnum_mt_tot, mt_totVarBinning,  "", "mt", "_CMS_xtt_wShape_13TeVUp", 0, 1);
            drawSignalRegionVarBin((signalCut_base + cut_options_nom)*cut_weights_KPtDown, (rel_signalCut_base + cut_options_nom)*cut_weights_KPtDown, (SScut_muTau_base + cut_options_nom)*cut_weights_KPtDown, (SScut_muTau_base_RelTauIso + cut_options_nom)*cut_weights_KPtDown, wSF_SS_muTau, (signalCut_base + cut_options_nom)*cut_weights_nom, (SScut_muTau_base + cut_options_nom)*cut_weights_nom, 0, wSF_OS_muTau_KPtDown, "mt_tot", binnum_mt_tot, mt_totVarBinning,  "", "mt", "_CMS_xtt_wShape_13TeVDown", 0, 1);
            drawSignalRegionVarBin((signalCut_base + cut_options_TauFakeUp)*cut_weights_TauFakeUp, (rel_signalCut_base + cut_options_TauFakeUp)*cut_weights_TauFakeUp, (SScut_muTau_base + cut_options_TauFakeUp)*cut_weights_TauFakeUp, (SScut_muTau_base_RelTauIso + cut_options_TauFakeUp)*cut_weights_TauFakeUp, wSF_SS_muTau, (signalCut_base + cut_options_nom)*cut_weights_nom, (SScut_muTau_base + cut_options_nom)*cut_weights_nom, 1, wSF_OS_muTau_TauFakeUp, "mt_tot", binnum_mt_tot, mt_totVarBinning,  "", "mt", "_CMS_xtt_jetToTauFake_13TeVUp", 0, 1);
            drawSignalRegionVarBin((signalCut_base + cut_options_TauFakeDown)*cut_weights_TauFakeDown, (rel_signalCut_base + cut_options_TauFakeDown)*cut_weights_TauFakeDown, (SScut_muTau_base + cut_options_TauFakeDown)*cut_weights_TauFakeDown, (SScut_muTau_base_RelTauIso + cut_options_TauFakeDown)*cut_weights_TauFakeDown, wSF_SS_muTau, (signalCut_base + cut_options_nom)*cut_weights_nom, (SScut_muTau_base + cut_options_nom)*cut_weights_nom, 1, wSF_OS_muTau_TauFakeDown, "mt_tot", binnum_mt_tot, mt_totVarBinning,  "", "mt", "_CMS_xtt_jetToTauFake_13TeVDown", 0, 1);
            
            drawSignalRegionVarBin((signalCut_base + cut_options_nom)*cut_weights_nom, (rel_signalCut_base + cut_options_nom)*cut_weights_nom, (SScut_muTau_base + cut_options_nom)*cut_weights_nom, (SScut_muTau_base_RelTauIso + cut_options_nom)*cut_weights_nom, wSF_SS_muTau, (signalCut_base + cut_options_nom)*cut_weights_nom, (SScut_muTau_base + cut_options_nom)*cut_weights_nom, 0, wSF_OS_muTau_ZLUp , "mt_tot*0.98", binnum_mt_tot, mt_totVarBinning,  "", "mt", "_CMS_xtt_ZLScale_mutau_13TeVUp", 0, 1);
            drawSignalRegionVarBin((signalCut_base + cut_options_nom)*cut_weights_nom, (rel_signalCut_base + cut_options_nom)*cut_weights_nom, (SScut_muTau_base + cut_options_nom)*cut_weights_nom, (SScut_muTau_base_RelTauIso + cut_options_nom)*cut_weights_nom, wSF_SS_muTau, (signalCut_base + cut_options_nom)*cut_weights_nom, (SScut_muTau_base + cut_options_nom)*cut_weights_nom, 0, wSF_OS_muTau_ZLDown, "mt_tot*1.02", binnum_mt_tot, mt_totVarBinning,  "", "mt", "_CMS_xtt_ZLScale_mutau_13TeVDown", 0, 1);
            drawSignalRegionVarBin((signalCut_base + cut_options_bTagUp)*cut_weights_bTagUp, (rel_signalCut_base + cut_options_bTagUp)*cut_weights_bTagUp, (SScut_muTau_base + cut_options_bTagUp)*cut_weights_bTagUp, (SScut_muTau_base_RelTauIso + cut_options_bTagUp)*cut_weights_bTagUp, wSF_SS_muTau, (signalCut_base + cut_options_nom)*cut_weights_nom, (SScut_muTau_base + cut_options_nom)*cut_weights_nom, 0, wSF_OS_muTau_bTagUp, "mt_tot", binnum_mt_tot, mt_totVarBinning,  "", "mt", "_CMS_xtt_bTag_13TeVUp", 0, 1);
            drawSignalRegionVarBin((signalCut_base + cut_options_bTagDown)*cut_weights_bTagDown, (rel_signalCut_base + cut_options_bTagDown)*cut_weights_bTagDown, (SScut_muTau_base + cut_options_bTagDown)*cut_weights_bTagDown, (SScut_muTau_base_RelTauIso + cut_options_bTagDown)*cut_weights_bTagDown, wSF_SS_muTau, (signalCut_base + cut_options_nom)*cut_weights_nom, (SScut_muTau_base + cut_options_nom)*cut_weights_nom, 0, wSF_OS_muTau_bTagDown, "mt_tot", binnum_mt_tot, mt_totVarBinning,  "", "mt", "_CMS_xtt_bTag_13TeVDown", 0, 1);
            drawSignalRegionVarBin((signalCut_base + cut_options_JECUp)*cut_weights_JECUp, (rel_signalCut_base + cut_options_JECUp)*cut_weights_JECUp, (SScut_muTau_base + cut_options_JECUp)*cut_weights_JECUp, (SScut_muTau_base_RelTauIso + cut_options_JECUp)*cut_weights_JECUp, wSF_SS_muTau, (signalCut_base + cut_options_nom)*cut_weights_nom, (SScut_muTau_base + cut_options_nom)*cut_weights_nom, 0, wSF_OS_muTau_JECUp, "mt_tot_JEnUp", binnum_mt_tot, mt_totVarBinning,  "", "mt", "_CMS_scale_j_13TeVUp", 0, 1);
            drawSignalRegionVarBin((signalCut_base + cut_options_JECDown)*cut_weights_JECDown, (rel_signalCut_base + cut_options_JECDown)*cut_weights_JECDown, (SScut_muTau_base + cut_options_JECDown)*cut_weights_JECDown, (SScut_muTau_base_RelTauIso + cut_options_JECDown)*cut_weights_JECDown, wSF_SS_muTau, (signalCut_base + cut_options_nom)*cut_weights_nom, (SScut_muTau_base + cut_options_nom)*cut_weights_nom, 0, wSF_OS_muTau_JECDown, "mt_tot_JEnDown", binnum_mt_tot, mt_totVarBinning,  "", "mt","_CMS_scale_j_13TeVDown", 0, 1);
            drawSignalRegionVarBin((signalCut_base + cut_options_HPTUp)*cut_weights_HPTUp, (rel_signalCut_base + cut_options_HPTUp)*cut_weights_HPTUp, (SScut_muTau_base + cut_options_HPTUp)*cut_weights_HPTUp, (SScut_muTau_base_RelTauIso + cut_options_HPTUp)*cut_weights_HPTUp, wSF_SS_muTau, (signalCut_base + cut_options_nom)*cut_weights_nom, (SScut_muTau_base + cut_options_nom)*cut_weights_nom, 0, wSF_OS_muTau_HPTUp, "mt_tot", binnum_mt_tot, mt_totVarBinning,  "", "mt", "_CMS_xtt_highTauEffi_13TeVUp", 0, 1);
            drawSignalRegionVarBin((signalCut_base + cut_options_HPTDown)*cut_weights_HPTDown, (rel_signalCut_base + cut_options_HPTDown)*cut_weights_HPTDown, (SScut_muTau_base + cut_options_HPTDown)*cut_weights_HPTDown, (SScut_muTau_base_RelTauIso + cut_options_HPTDown)*cut_weights_HPTDown, wSF_SS_muTau, (signalCut_base + cut_options_nom)*cut_weights_nom, (SScut_muTau_base + cut_options_nom)*cut_weights_nom, 0, wSF_OS_muTau_HPTDown, "mt_tot", binnum_mt_tot, mt_totVarBinning,  "", "mt", "_CMS_xtt_highTauEffi_13TeVDown", 0, 1);
            drawSignalRegionVarBin((signalCut_base + cut_options_TESUp)*cut_weights_nom, (rel_signalCut_base + cut_options_TESUp)*cut_weights_nom, (SScut_muTau_base + cut_options_TESUp)*cut_weights_nom, (SScut_muTau_base_RelTauIso + cut_options_TESUp)*cut_weights_nom, wSF_SS_muTau, (signalCut_base + cut_options_nom)*cut_weights_nom, (SScut_muTau_base + cut_options_nom)*cut_weights_nom, 0, wSF_OS_muTau_TESUp, "mt_tot_TESUp", binnum_mt_tot, mt_totVarBinning,  "", "mt", "_CMS_scale_t_13TeVUp", 0, 1);
            drawSignalRegionVarBin((signalCut_base + cut_options_TESDown)*cut_weights_nom, (rel_signalCut_base + cut_options_TESUp)*cut_weights_nom, (SScut_muTau_base + cut_options_TESUp)*cut_weights_nom, (SScut_muTau_base_RelTauIso + cut_options_TESUp)*cut_weights_nom, wSF_SS_muTau, (signalCut_base + cut_options_nom)*cut_weights_nom, (SScut_muTau_base + cut_options_nom)*cut_weights_nom, 0, wSF_OS_muTau_TESDown, "mt_tot_TESDown", binnum_mt_tot, mt_totVarBinning,  "", "mt","_CMS_scale_t_13TeVDown", 0, 1);
            drawSignalRegionVarBin((signalCut_base + cut_options_dm0TESUp)*cut_weights_nom, (rel_signalCut_base + cut_options_dm0TESUp)*cut_weights_nom, (SScut_muTau_base + cut_options_dm0TESUp)*cut_weights_nom, (SScut_muTau_base_RelTauIso + cut_options_dm0TESUp)*cut_weights_nom, wSF_SS_muTau, (signalCut_base + cut_options_nom)*cut_weights_nom, (SScut_muTau_base + cut_options_nom)*cut_weights_nom, 0, wSF_OS_muTau_dm0TESUp, "mt_tot_dm0TESUp", binnum_mt_tot, mt_totVarBinning,  "", "mt", "_CMS_scale_t_dm0_13TeVUp", 0, 1);
            drawSignalRegionVarBin((signalCut_base + cut_options_dm0TESDown)*cut_weights_nom, (rel_signalCut_base + cut_options_dm0TESUp)*cut_weights_nom, (SScut_muTau_base + cut_options_dm0TESUp)*cut_weights_nom, (SScut_muTau_base_RelTauIso + cut_options_dm0TESUp)*cut_weights_nom, wSF_SS_muTau, (signalCut_base + cut_options_nom)*cut_weights_nom, (SScut_muTau_base + cut_options_nom)*cut_weights_nom, 0, wSF_OS_muTau_dm0TESDown, "mt_tot_dm0TESDown", binnum_mt_tot, mt_totVarBinning,  "", "mt","_CMS_scale_t_dm0_13TeVDown", 0, 1);
            drawSignalRegionVarBin((signalCut_base + cut_options_dm1TESUp)*cut_weights_nom, (rel_signalCut_base + cut_options_dm1TESUp)*cut_weights_nom, (SScut_muTau_base + cut_options_dm1TESUp)*cut_weights_nom, (SScut_muTau_base_RelTauIso + cut_options_dm1TESUp)*cut_weights_nom, wSF_SS_muTau, (signalCut_base + cut_options_nom)*cut_weights_nom, (SScut_muTau_base + cut_options_nom)*cut_weights_nom, 0, wSF_OS_muTau_dm1TESUp, "mt_tot_dm1TESUp", binnum_mt_tot, mt_totVarBinning,  "", "mt", "_CMS_scale_t_dm1_13TeVUp", 0, 1);
            drawSignalRegionVarBin((signalCut_base + cut_options_dm1TESDown)*cut_weights_nom, (rel_signalCut_base + cut_options_dm1TESUp)*cut_weights_nom, (SScut_muTau_base + cut_options_dm1TESUp)*cut_weights_nom, (SScut_muTau_base_RelTauIso + cut_options_dm1TESUp)*cut_weights_nom, wSF_SS_muTau, (signalCut_base + cut_options_nom)*cut_weights_nom, (SScut_muTau_base + cut_options_nom)*cut_weights_nom, 0, wSF_OS_muTau_dm1TESDown, "mt_tot_dm1TESDown", binnum_mt_tot, mt_totVarBinning,  "", "mt","_CMS_scale_t_dm1_13TeVDown", 0, 1);
            drawSignalRegionVarBin((signalCut_base + cut_options_dm10TESUp)*cut_weights_nom, (rel_signalCut_base + cut_options_dm10TESUp)*cut_weights_nom, (SScut_muTau_base + cut_options_dm10TESUp)*cut_weights_nom, (SScut_muTau_base_RelTauIso + cut_options_dm10TESUp)*cut_weights_nom, wSF_SS_muTau, (signalCut_base + cut_options_nom)*cut_weights_nom, (SScut_muTau_base + cut_options_nom)*cut_weights_nom, 0, wSF_OS_muTau_dm10TESUp, "mt_tot_dm10TESUp", binnum_mt_tot, mt_totVarBinning,  "", "mt", "_CMS_scale_t_dm10_13TeVUp", 0, 1);
            drawSignalRegionVarBin((signalCut_base + cut_options_dm10TESDown)*cut_weights_nom, (rel_signalCut_base + cut_options_dm10TESUp)*cut_weights_nom, (SScut_muTau_base + cut_options_dm10TESUp)*cut_weights_nom, (SScut_muTau_base_RelTauIso + cut_options_dm10TESUp)*cut_weights_nom, wSF_SS_muTau, (signalCut_base + cut_options_nom)*cut_weights_nom, (SScut_muTau_base + cut_options_nom)*cut_weights_nom, 0, wSF_OS_muTau_dm10TESDown, "mt_tot_dm10TESDown", binnum_mt_tot, mt_totVarBinning,  "", "mt","_CMS_scale_t_dm10_13TeVDown", 0, 1);
            drawSignalRegionVarBin((signalCut_base + cut_options_nom)*cut_weights_TopPtUp, (rel_signalCut_base + cut_options_nom)*cut_weights_TopPtUp, (SScut_muTau_base + cut_options_nom)*cut_weights_TopPtUp, (SScut_muTau_base_RelTauIso + cut_options_nom)*cut_weights_TopPtUp, wSF_SS_muTau, (signalCut_base + cut_options_nom)*cut_weights_nom, (SScut_muTau_base + cut_options_nom)*cut_weights_nom, 0, wSF_OS_muTau_TopPtUp, "mt_tot", binnum_mt_tot, mt_totVarBinning,  "", "mt", "_CMS_xtt_ttbarShape_13TeVUp", 0, 1);
            drawSignalRegionVarBin((signalCut_base + cut_options_nom)*cut_weights_TopPtDown, (rel_signalCut_base + cut_options_nom)*cut_weights_TopPtDown, (SScut_muTau_base + cut_options_nom)*cut_weights_TopPtDown, (SScut_muTau_base_RelTauIso + cut_options_nom)*cut_weights_TopPtDown, wSF_SS_muTau,(signalCut_base + cut_options_nom)*cut_weights_nom, (SScut_muTau_base + cut_options_nom)*cut_weights_nom, 0, wSF_OS_muTau_TopPtDown, "mt_tot", binnum_mt_tot, mt_totVarBinning,  "", "mt","_CMS_xtt_ttbarShape_13TeVDown", 0, 1);
            drawSignalRegionVarBin((signalCut_base + cut_options_nom)*cut_weights_WWUp, (rel_signalCut_base + cut_options_nom)*cut_weights_WWUp, (SScut_muTau_base + cut_options_nom)*cut_weights_WWUp, (SScut_muTau_base_RelTauIso + cut_options_nom)*cut_weights_WWUp, wSF_SS_muTau, (signalCut_base + cut_options_nom)*cut_weights_WWUp, (rel_signalCut_base + cut_options_nom)*cut_weights_WWUp, 0, wSF_OS_muTau, "mt_tot", binnum_mt_tot, mt_totVarBinning,  "", "mt", "_CMS_xtt_WWNLOewk_13TeVUp", 0, 1);
            drawSignalRegionVarBin((signalCut_base + cut_options_nom)*cut_weights_WWDown, (rel_signalCut_base + cut_options_nom)*cut_weights_WWDown, (SScut_muTau_base + cut_options_nom)*cut_weights_WWDown, (SScut_muTau_base_RelTauIso + cut_options_nom)*cut_weights_WWDown, wSF_SS_muTau, (signalCut_base + cut_options_nom)*cut_weights_nom, (SScut_muTau_base + cut_options_nom)*cut_weights_nom, 0, wSF_OS_muTau, "mt_tot", binnum_mt_tot, mt_totVarBinning,  "", "mt","_CMS_xtt_WWNLOewk_13TeVDown", 0, 1);
            drawSignalRegionVarBin((signalCut_base + cut_options_nom)*cut_weights_ZZUp, (rel_signalCut_base + cut_options_nom)*cut_weights_ZZUp, (SScut_muTau_base + cut_options_nom)*cut_weights_ZZUp, (SScut_muTau_base_RelTauIso + cut_options_nom)*cut_weights_ZZUp, wSF_SS_muTau, (signalCut_base + cut_options_nom)*cut_weights_nom, (SScut_muTau_base + cut_options_nom)*cut_weights_nom, 0, wSF_OS_muTau, "mt_tot", binnum_mt_tot, mt_totVarBinning,  "", "mt", "_CMS_xtt_ZZNLOewk_13TeVUp", 0, 1);
            drawSignalRegionVarBin((signalCut_base + cut_options_nom)*cut_weights_ZZDown, (rel_signalCut_base + cut_options_nom)*cut_weights_ZZDown, (SScut_muTau_base + cut_options_nom)*cut_weights_ZZDown, (SScut_muTau_base_RelTauIso + cut_options_nom)*cut_weights_ZZDown, wSF_SS_muTau, (signalCut_base + cut_options_nom)*cut_weights_nom, (SScut_muTau_base + cut_options_nom)*cut_weights_nom, 0, wSF_OS_muTau, "mt_tot", binnum_mt_tot, mt_totVarBinning,  "", "mt","_CMS_xtt_ZZNLOewk_13TeVDown", 0, 1);
            
            if(printSyst)
            {
            
                countTotalsVarBin((signalCut_base + cut_options_UESUp)*cut_weights_UESUp, (rel_signalCut_base + cut_options_UESUp)*cut_weights_UESUp, (SScut_muTau_base + cut_options_UESUp)*cut_weights_UESUp, (SScut_muTau_base_RelTauIso + cut_options_UESUp)*cut_weights_UESUp, wSF_SS_muTau, wSF_OS_muTau, "mt_tot_UESUp", binnum_mt_tot, mt_totVarBinning, "Mt_total", "mt", "_CMS_scale_m_13TeVUp");
                countTotalsVarBin((signalCut_base + cut_options_UESDown)*cut_weights_UESDown, (rel_signalCut_base + cut_options_UESDown)*cut_weights_UESDown, (SScut_muTau_base + cut_options_UESDown)*cut_weights_UESDown, (SScut_muTau_base_RelTauIso + cut_options_UESDown)*cut_weights_UESDown, wSF_SS_muTau, wSF_OS_muTau, "mt_tot_UESDown", binnum_mt_tot, mt_totVarBinning, "Mt_total", "mt", "_CMS_scale_m_13TeVDown");
                countTotalsVarBin((signalCut_base + cut_options_nom)*cut_weights_ZPtUp, (rel_signalCut_base + cut_options_nom)*cut_weights_ZPtUp, (SScut_muTau_base + cut_options_nom)*cut_weights_ZPtUp, (SScut_muTau_base_RelTauIso + cut_options_nom)*cut_weights_ZPtUp, wSF_SS_muTau, wSF_OS_muTau, "mt_tot", binnum_mt_tot, mt_totVarBinning, "Mt_total", "mt", "_CMS_xtt_dyShape_13TeVUp");
                countTotalsVarBin((signalCut_base + cut_options_nom)*cut_weights_ZPtDown, (rel_signalCut_base + cut_options_nom)*cut_weights_ZPtDown, (SScut_muTau_base + cut_options_nom)*cut_weights_ZPtDown, (SScut_muTau_base_RelTauIso + cut_options_nom)*cut_weights_ZPtDown, wSF_SS_muTau, wSF_OS_muTau, "mt_tot", binnum_mt_tot, mt_totVarBinning, "Mt_total", "mt", "_CMS_xtt_dyShape_13TeVDown");
                countTotalsVarBin((signalCut_base + cut_options_nom)*cut_weights_KPtUp, (rel_signalCut_base + cut_options_nom)*cut_weights_KPtUp, (SScut_muTau_base + cut_options_nom)*cut_weights_KPtUp, (SScut_muTau_base_RelTauIso + cut_options_nom)*cut_weights_KPtUp, wSF_SS_muTau, wSF_OS_muTau, "mt_tot", binnum_mt_tot, mt_totVarBinning, "Mt_total", "mt", "_CMS_xtt_wShape_13TeVUp");
                countTotalsVarBin((signalCut_base + cut_options_nom)*cut_weights_KPtDown, (rel_signalCut_base + cut_options_nom)*cut_weights_KPtDown, (SScut_muTau_base + cut_options_nom)*cut_weights_KPtDown, (SScut_muTau_base_RelTauIso + cut_options_nom)*cut_weights_KPtDown, wSF_SS_muTau, wSF_OS_muTau, "mt_tot", binnum_mt_tot, mt_totVarBinning, "Mt_total", "mt", "_CMS_xtt_wShape_13TeVDown");
                countTotalsVarBin((signalCut_base + cut_options_TauFakeUp)*cut_weights_TauFakeUp, (rel_signalCut_base + cut_options_TauFakeUp)*cut_weights_TauFakeUp, (SScut_muTau_base + cut_options_TauFakeUp)*cut_weights_TauFakeUp, (SScut_muTau_base_RelTauIso + cut_options_TauFakeUp)*cut_weights_TauFakeUp, wSF_SS_muTau, wSF_OS_muTau, "mt_tot", binnum_mt_tot, mt_totVarBinning, "Mt_total", "mt", "_CMS_xtt_jetToTauFake_13TeVUp");
                countTotalsVarBin((signalCut_base + cut_options_TauFakeDown)*cut_weights_TauFakeDown, (rel_signalCut_base + cut_options_TauFakeDown)*cut_weights_TauFakeDown, (SScut_muTau_base + cut_options_TauFakeDown)*cut_weights_TauFakeDown, (SScut_muTau_base_RelTauIso + cut_options_TauFakeDown)*cut_weights_TauFakeDown, wSF_SS_muTau, wSF_OS_muTau, "mt_tot", binnum_mt_tot, mt_totVarBinning, "Mt_total", "mt", "_CMS_xtt_jetToTauFake_13TeVDown");
                countTotalsVarBin((signalCut_base + cut_options_nom)*cut_weights_nom, (rel_signalCut_base + cut_options_nom)*cut_weights_nom, (SScut_muTau_base + cut_options_nom)*cut_weights_nom, (SScut_muTau_base_RelTauIso + cut_options_nom)*cut_weights_nom, wSF_SS_muTau, wSF_OS_muTau, "mt_tot*0.98", binnum_mt_tot, mt_totVarBinning, "Mt_total", "mt", "_CMS_xtt_ZLScale_mutau_13TeVUp");
                countTotalsVarBin((signalCut_base + cut_options_nom)*cut_weights_nom, (rel_signalCut_base + cut_options_nom)*cut_weights_nom, (SScut_muTau_base + cut_options_nom)*cut_weights_nom, (SScut_muTau_base_RelTauIso + cut_options_nom)*cut_weights_nom, wSF_SS_muTau, wSF_OS_muTau, "mt_tot*1.02", binnum_mt_tot, mt_totVarBinning, "Mt_total", "mt", "_CMS_xtt_ZLScale_mutau_13TeVDown");
                countTotalsVarBin((signalCut_base + cut_options_bTagUp)*cut_weights_bTagUp, (rel_signalCut_base + cut_options_bTagUp)*cut_weights_bTagUp, (SScut_muTau_base + cut_options_bTagUp)*cut_weights_bTagUp, (SScut_muTau_base_RelTauIso + cut_options_bTagUp)*cut_weights_bTagUp, wSF_SS_muTau, wSF_OS_muTau, "mt_tot", binnum_mt_tot, mt_totVarBinning, "Mt_total", "mt", "_CMS_xtt_bTag_13TeVUp");
                countTotalsVarBin((signalCut_base + cut_options_bTagDown)*cut_weights_bTagDown, (rel_signalCut_base + cut_options_bTagDown)*cut_weights_bTagDown, (SScut_muTau_base + cut_options_bTagDown)*cut_weights_bTagDown, (SScut_muTau_base_RelTauIso + cut_options_bTagDown)*cut_weights_bTagDown, wSF_SS_muTau, wSF_OS_muTau, "mt_tot", binnum_mt_tot, mt_totVarBinning, "Mt_total", "mt", "_CMS_xtt_bTag_13TeVDown");
                countTotalsVarBin((signalCut_base + cut_options_JECUp)*cut_weights_JECUp, (rel_signalCut_base + cut_options_JECUp)*cut_weights_JECUp, (SScut_muTau_base + cut_options_JECUp)*cut_weights_JECUp, (SScut_muTau_base_RelTauIso + cut_options_JECUp)*cut_weights_JECUp, wSF_SS_muTau, wSF_OS_muTau, "mt_tot_JEnUp", binnum_mt_tot, mt_totVarBinning, "Mt_total", "mt", "_CMS_scale_j_13TeVUp");
                countTotalsVarBin((signalCut_base + cut_options_JECDown)*cut_weights_JECDown, (rel_signalCut_base + cut_options_JECDown)*cut_weights_JECDown, (SScut_muTau_base + cut_options_JECDown)*cut_weights_JECDown, (SScut_muTau_base_RelTauIso + cut_options_JECDown)*cut_weights_JECDown, wSF_SS_muTau, wSF_OS_muTau, "mt_tot_JEnDown", binnum_mt_tot, mt_totVarBinning, "Mt_total", "mt", "_CMS_scale_j_13TeVDown");
                countTotalsVarBin((signalCut_base + cut_options_HPTUp)*cut_weights_HPTUp, (rel_signalCut_base + cut_options_HPTUp)*cut_weights_HPTUp, (SScut_muTau_base + cut_options_HPTUp)*cut_weights_HPTUp, (SScut_muTau_base_RelTauIso + cut_options_HPTUp)*cut_weights_HPTUp, wSF_SS_muTau, wSF_OS_muTau, "mt_tot", binnum_mt_tot, mt_totVarBinning, "Mt_total", "mt", "_CMS_xtt_highTauEffi_13TeVUp");
                countTotalsVarBin((signalCut_base + cut_options_HPTDown)*cut_weights_HPTDown, (rel_signalCut_base + cut_options_HPTDown)*cut_weights_HPTDown, (SScut_muTau_base + cut_options_HPTDown)*cut_weights_HPTDown, (SScut_muTau_base_RelTauIso + cut_options_HPTDown)*cut_weights_HPTDown, wSF_SS_muTau, wSF_OS_muTau, "mt_tot", binnum_mt_tot, mt_totVarBinning, "Mt_total", "mt", "_CMS_xtt_highTauEffi_13TeVDown");
                countTotalsVarBin((signalCut_base + cut_options_TESUp)*cut_weights_nom, (rel_signalCut_base + cut_options_TESUp)*cut_weights_nom, (SScut_muTau_base + cut_options_TESUp)*cut_weights_nom, (SScut_muTau_base_RelTauIso + cut_options_TESUp)*cut_weights_nom, wSF_SS_muTau, wSF_OS_muTau, "mt_tot", binnum_mt_tot, mt_totVarBinning, "Mt_total", "mt", "_CMS_scale_t_13TeVUp");
                countTotalsVarBin((signalCut_base + cut_options_TESDown)*cut_weights_nom, (rel_signalCut_base + cut_options_TESDown)*cut_weights_nom, (SScut_muTau_base + cut_options_nom)*cut_weights_nom, (SScut_muTau_base_RelTauIso + cut_options_nom)*cut_weights_nom, wSF_SS_muTau, wSF_OS_muTau, "mt_tot", binnum_mt_tot, mt_totVarBinning, "Mt_total", "mt", "_CMS_scale_t_13TeVDown");
                countTotalsVarBin((signalCut_base + cut_options_nom)*cut_weights_TopPtUp, (rel_signalCut_base + cut_options_nom)*cut_weights_TopPtUp, (SScut_muTau_base + cut_options_nom)*cut_weights_TopPtUp, (SScut_muTau_base_RelTauIso + cut_options_nom)*cut_weights_TopPtUp, wSF_SS_muTau, wSF_OS_muTau, "mt_tot", binnum_mt_tot, mt_totVarBinning, "Mt_total", "mt", "_CMS_xtt_ttbarShape_13TeVUp");
                countTotalsVarBin((signalCut_base + cut_options_nom)*cut_weights_TopPtDown, (rel_signalCut_base + cut_options_nom)*cut_weights_TopPtDown, (SScut_muTau_base + cut_options_nom)*cut_weights_TopPtDown, (SScut_muTau_base_RelTauIso + cut_options_nom)*cut_weights_TopPtDown, wSF_SS_muTau, wSF_OS_muTau, "mt_tot", binnum_mt_tot, mt_totVarBinning, "Mt_total", "mt", "_CMS_xtt_ttbarShape_13TeVDown");
                countTotalsVarBin((signalCut_base + cut_options_nom)*cut_weights_WWUp, (rel_signalCut_base + cut_options_nom)*cut_weights_WWUp, (SScut_muTau_base + cut_options_nom)*cut_weights_WWUp, (SScut_muTau_base_RelTauIso + cut_options_nom)*cut_weights_WWUp, wSF_SS_muTau, wSF_OS_muTau, "mt_tot", binnum_mt_tot, mt_totVarBinning, "Mt_total", "mt", "_CMS_xtt_WWNLOewk_13TeVUp");
                countTotalsVarBin((signalCut_base + cut_options_nom)*cut_weights_WWDown, (rel_signalCut_base + cut_options_nom)*cut_weights_WWDown, (SScut_muTau_base + cut_options_nom)*cut_weights_WWDown, (SScut_muTau_base_RelTauIso + cut_options_nom)*cut_weights_WWDown, wSF_SS_muTau, wSF_OS_muTau, "mt_tot", binnum_mt_tot, mt_totVarBinning, "Mt_total", "mt", "_CMS_xtt_WWNLOewk_13TeVDown");
                countTotalsVarBin((signalCut_base + cut_options_nom)*cut_weights_ZZUp, (rel_signalCut_base + cut_options_nom)*cut_weights_ZZUp, (SScut_muTau_base + cut_options_nom)*cut_weights_ZZUp, (SScut_muTau_base_RelTauIso + cut_options_nom)*cut_weights_ZZUp, wSF_SS_muTau, wSF_OS_muTau, "mt_tot", binnum_mt_tot, mt_totVarBinning, "Mt_total", "mt", "_CMS_xtt_ZZNLOewk_13TeVUp");
                countTotalsVarBin((signalCut_base + cut_options_nom)*cut_weights_ZZDown, (rel_signalCut_base + cut_options_nom)*cut_weights_ZZDown, (SScut_muTau_base + cut_options_nom)*cut_weights_ZZDown, (SScut_muTau_base_RelTauIso + cut_options_nom)*cut_weights_ZZDown, wSF_SS_muTau, wSF_OS_muTau, "mt_tot", binnum_mt_tot, mt_totVarBinning, "Mt_total", "mt", "_CMS_xtt_ZZNLOewk_13TeVDown");

            }

        }
        reset_files();

    }
    else if (choice2==2)
    {
        if (doControlTMVA)
        {
            global_title = "QCD (Same Sign) Estimate Region";
            
            global_title = "MET > " + control_met_cutoff + " GeV";
            drawSignalRegion((signalCut_base + cut_options_nom)*cut_weights_nom, (rel_signalCut_base + cut_options_nom)*cut_weights_nom, (SScut_muTau_base + cut_options_nom)*cut_weights_nom, (SScut_muTau_base_RelTauIso + cut_options_nom)*cut_weights_nom, wSF_SS_muTau, 1.0 , "mvaVar_mt_MZP600A0400", tmvaBinning, "sig region Function TMVA variable", "mt", "", 1, 0);
            drawSignalRegion((signalCut_base + cut_options_nom)*cut_weights_nom, (rel_signalCut_base + cut_options_nom)*cut_weights_nom, (SScut_muTau_base + cut_options_nom)*cut_weights_nom, (SScut_muTau_base_RelTauIso + cut_options_nom)*cut_weights_nom, wSF_SS_muTau, 1.0 , "mvaVar_mt_MZP800A0400", tmvaBinning, "sig region Function TMVA variable", "mt", "", 1, 0);
            drawSignalRegion((signalCut_base + cut_options_nom)*cut_weights_nom, (rel_signalCut_base + cut_options_nom)*cut_weights_nom, (SScut_muTau_base + cut_options_nom)*cut_weights_nom, (SScut_muTau_base_RelTauIso + cut_options_nom)*cut_weights_nom, wSF_SS_muTau, 1.0 , "mvaVar_mt_MZP1000A0400", tmvaBinning, "sig region Function TMVA variable", "mt", "", 1, 0);
            drawSignalRegion((signalCut_base + cut_options_nom)*cut_weights_nom, (rel_signalCut_base + cut_options_nom)*cut_weights_nom, (SScut_muTau_base + cut_options_nom)*cut_weights_nom, (SScut_muTau_base_RelTauIso + cut_options_nom)*cut_weights_nom, wSF_SS_muTau, 1.0 , "mvaVar_mt_MZP1200A0400", tmvaBinning, "sig region Function TMVA variable", "mt", "", 1, 0);

        }
        else
        {
            global_title = "QCD (Same Sign) Estimate Region";
            
            countTotalsVarBin((signalCut_base + cut_options_nom)*cut_weights_nom, (rel_signalCut_base + cut_options_nom)*cut_weights_nom, (SScut_muTau_base + cut_options_nom)*cut_weights_nom, (SScut_muTau_base_RelTauIso + cut_options_nom)*cut_weights_nom, wSF_SS_muTau, wSF_OS_muTau, "mt_tot", binnum_mt_tot, mt_totVarBinning, "control region Mt_total", "mt", "");
            
            global_title = "MET > " + control_met_cutoff + " GeV";
            drawSignalRegion((signalCut_base + cut_options_nom)*cut_weights_nom, (rel_signalCut_base + cut_options_nom)*cut_weights_nom, (SScut_muTau_base + cut_options_nom)*cut_weights_nom, (SScut_muTau_base_RelTauIso + cut_options_nom)*cut_weights_nom, wSF_SS_muTau, wSF_OS_muTau, "mt_1", mtBinning, "Control mt", "mt", "", 0, 1);
            drawSignalRegion((signalCut_base + cut_options_nom)*cut_weights_nom, (rel_signalCut_base + cut_options_nom)*cut_weights_nom, (SScut_muTau_base + cut_options_nom)*cut_weights_nom, (SScut_muTau_base_RelTauIso + cut_options_nom)*cut_weights_nom, wSF_SS_muTau, wSF_OS_muTau, "m_vis", m_visBinning, "Control mvis", "mt", "", 0, 1);
            drawSignalRegion((signalCut_base + cut_options_nom)*cut_weights_nom, (rel_signalCut_base + cut_options_nom)*cut_weights_nom, (SScut_muTau_base + cut_options_nom)*cut_weights_nom, (SScut_muTau_base_RelTauIso + cut_options_nom)*cut_weights_nom, wSF_SS_muTau, wSF_OS_muTau, "nbtag", nbtagBinning, "Control nbtag", "mt", "", 0, 1);
            drawSignalRegion((signalCut_base + cut_options_nom)*cut_weights_nom, (rel_signalCut_base + cut_options_nom)*cut_weights_nom, (SScut_muTau_base + cut_options_nom)*cut_weights_nom, (SScut_muTau_base_RelTauIso + cut_options_nom)*cut_weights_nom, wSF_SS_muTau, wSF_OS_muTau, "met", metBinning, "Control met", "mt", "", 0, 1);
            drawSignalRegion((signalCut_base + cut_options_nom)*cut_weights_nom, (rel_signalCut_base + cut_options_nom)*cut_weights_nom, (SScut_muTau_base + cut_options_nom)*cut_weights_nom, (SScut_muTau_base_RelTauIso + cut_options_nom)*cut_weights_nom, wSF_SS_muTau, wSF_OS_muTau, "pt_1", ptBinning, "Control pt_1", "mt", "", 0, 1);
            drawSignalRegion((signalCut_base + cut_options_nom)*cut_weights_nom, (rel_signalCut_base + cut_options_nom)*cut_weights_nom, (SScut_muTau_base + cut_options_nom)*cut_weights_nom, (SScut_muTau_base_RelTauIso + cut_options_nom)*cut_weights_nom, wSF_SS_muTau, wSF_OS_muTau, "pt_2", ptBinning, "Control pt_2", "mt", "", 0, 1);
            drawSignalRegion((signalCut_base + cut_options_nom)*cut_weights_nom, (rel_signalCut_base + cut_options_nom)*cut_weights_nom, (SScut_muTau_base + cut_options_nom)*cut_weights_nom, (SScut_muTau_base_RelTauIso + cut_options_nom)*cut_weights_nom, wSF_SS_muTau, wSF_OS_muTau, "pt_tt", ptBinning, "Control pt_tt", "mt", "", 0, 1);
            
            drawSignalRegion((signalCut_base + cut_options_nom)*cut_weights_nom, (rel_signalCut_base + cut_options_nom)*cut_weights_nom, (SScut_muTau_base + cut_options_nom)*cut_weights_nom, (SScut_muTau_base_RelTauIso + cut_options_nom)*cut_weights_nom, wSF_SS_muTau, wSF_SS_muTau, "DeltaR_leg1_leg2", drBinning, "Control DeltaR_leg1_leg2", "mt", "", 0, 1);
            drawSignalRegion((signalCut_base + cut_options_nom)*cut_weights_nom, (rel_signalCut_base + cut_options_nom)*cut_weights_nom, (SScut_muTau_base + cut_options_nom)*cut_weights_nom, (SScut_muTau_base_RelTauIso + cut_options_nom)*cut_weights_nom, wSF_SS_muTau, wSF_OS_muTau, "m_1", mBinning, "Control m_1", "mt", "", 0, 1);
            drawSignalRegion((signalCut_base + cut_options_nom)*cut_weights_nom, (rel_signalCut_base + cut_options_nom)*cut_weights_nom, (SScut_muTau_base + cut_options_nom)*cut_weights_nom, (SScut_muTau_base_RelTauIso + cut_options_nom)*cut_weights_nom, wSF_SS_muTau, wSF_OS_muTau, "m_2", mBinning, "Control m_2", "mt", "", 0, 1);
            drawSignalRegion((signalCut_base + cut_options_nom)*cut_weights_nom, (rel_signalCut_base + cut_options_nom)*cut_weights_nom, (SScut_muTau_base + cut_options_nom)*cut_weights_nom, (SScut_muTau_base_RelTauIso + cut_options_nom)*cut_weights_nom, wSF_SS_muTau, wSF_OS_muTau, "mt_2", mtBinning, "Control mt_2", "mt", "", 0, 1);
            drawSignalRegion((signalCut_base + cut_options_nom)*cut_weights_nom, (rel_signalCut_base + cut_options_nom)*cut_weights_nom, (SScut_muTau_base + cut_options_nom)*cut_weights_nom, (SScut_muTau_base_RelTauIso + cut_options_nom)*cut_weights_nom, wSF_SS_muTau, wSF_OS_muTau, "phi_1", phiBinning, "Control phi_1", "mt", "", 0, 1);
            drawSignalRegion((signalCut_base + cut_options_nom)*cut_weights_nom, (rel_signalCut_base + cut_options_nom)*cut_weights_nom, (SScut_muTau_base + cut_options_nom)*cut_weights_nom, (SScut_muTau_base_RelTauIso + cut_options_nom)*cut_weights_nom, wSF_SS_muTau, wSF_OS_muTau, "phi_2", phiBinning, "Control phi_2", "mt", "", 0, 1);
            drawSignalRegion((signalCut_base + cut_options_nom)*cut_weights_nom, (rel_signalCut_base + cut_options_nom)*cut_weights_nom, (SScut_muTau_base + cut_options_nom)*cut_weights_nom, (SScut_muTau_base_RelTauIso + cut_options_nom)*cut_weights_nom, wSF_SS_muTau, wSF_OS_muTau, "eta_1", etaBinning, "Control eta_1", "mt", "", 0, 1);
            drawSignalRegion((signalCut_base + cut_options_nom)*cut_weights_nom, (rel_signalCut_base + cut_options_nom)*cut_weights_nom, (SScut_muTau_base + cut_options_nom)*cut_weights_nom, (SScut_muTau_base_RelTauIso + cut_options_nom)*cut_weights_nom, wSF_SS_muTau, wSF_OS_muTau, "eta_2", etaBinning, "Control eta_2", "mt", "", 0, 1);
            drawSignalRegion((signalCut_base + cut_options_nom)*cut_weights_nom, (rel_signalCut_base + cut_options_nom)*cut_weights_nom, (SScut_muTau_base + cut_options_nom)*cut_weights_nom, (SScut_muTau_base_RelTauIso + cut_options_nom)*cut_weights_nom, wSF_SS_muTau, wSF_OS_muTau, "njetspt20", nbtagBinning, "Control njetspt20", "mt", "", 0, 1);
            
        }
        reset_files();
    }

}

void findBkgFractions(TCut cut, double sf, std::string parameter, float bin[3])
{

    TH1F * DYinc_ = new TH1F("DYinc_","DYinc_",bin[0],bin[1],bin[2]);
	DYinc_->Sumw2();

	TH1F * DYhigh_ = new TH1F("DYhigh_","DYhigh_",bin[0],bin[1],bin[2]);
	DYhigh_->Sumw2();
    
	TH1F * DYlow_ = new TH1F("DYlow_","DYlow_",bin[0],bin[1],bin[2]);
	DYlow_->Sumw2();
    
    TH1F * DY1Jets_ = new TH1F("DY1Jets_","DY1Jets_",bin[0],bin[1],bin[2]);
	DY1Jets_->Sumw2();
    
    TH1F * DY2Jets_ = new TH1F("DY2Jets_","DY2Jets_",bin[0],bin[1],bin[2]);
	DY2Jets_->Sumw2();
    
    TH1F * DY3Jets_ = new TH1F("DY3Jets_","DY3Jets_",bin[0],bin[1],bin[2]);
	DY3Jets_->Sumw2();
    
    TH1F * DY4Jets_ = new TH1F("DY4Jets_","DY4Jets_",bin[0],bin[1],bin[2]);
	DY4Jets_->Sumw2();
    
	TH1F * WJetsToLNu_ = new TH1F("WJetsToLNu_","WJetsToLNu_",bin[0],bin[1],bin[2]);
	WJetsToLNu_->Sumw2();
    
    TH1F * W1JetsToLNu_ = new TH1F("W1JetsToLNu_","W1JetsToLNu_",bin[0],bin[1],bin[2]);
	W1JetsToLNu_->Sumw2();
    
    TH1F * W2JetsToLNu_ = new TH1F("W2JetsToLNu_","W2JetsToLNu_",bin[0],bin[1],bin[2]);
	W2JetsToLNu_->Sumw2();
    
    TH1F * W3JetsToLNu_ = new TH1F("W3JetsToLNu_","W3JetsToLNu_",bin[0],bin[1],bin[2]);
	W3JetsToLNu_->Sumw2();
    
    TH1F * W4JetsToLNu_ = new TH1F("W4JetsToLNu_","W4JetsToLNu_",bin[0],bin[1],bin[2]);
	W4JetsToLNu_->Sumw2();
    
    TH1F * TT_ = new TH1F("TT_","TT_",bin[0],bin[1],bin[2]);
	TT_->Sumw2();

    TH1F * ST_t_channel_antitop_4f_leptonDecays_ = new TH1F("ST_t_channel_antitop_4f_leptonDecays_","ST_t_channel_antitop_4f_leptonDecays_",bin[0],bin[1],bin[2]);
	ST_t_channel_antitop_4f_leptonDecays_->Sumw2();
    
	TH1F * ST_t_channel_top_4f_leptonDecays_ = new TH1F("ST_t_channel_top_4f_leptonDecays_","ST_t_channel_top_4f_leptonDecays_",bin[0],bin[1],bin[2]);
	ST_t_channel_top_4f_leptonDecays_->Sumw2();
    
	TH1F * ST_tW_antitop_5f_inclusiveDecays_ = new TH1F("ST_tW_antitop_5f_inclusiveDecays_","ST_tW_antitop_5f_inclusiveDecays_",bin[0],bin[1],bin[2]);
	ST_tW_antitop_5f_inclusiveDecays_->Sumw2();
    
	TH1F * ST_tW_top_5f_inclusiveDecays_ = new TH1F("ST_tW_top_5f_inclusiveDecays_","ST_tW_top_5f_inclusiveDecays_",bin[0],bin[1],bin[2]);
	ST_tW_top_5f_inclusiveDecays_->Sumw2();
    
	TH1F * WWTo2L2Nu_ = new TH1F("WWTo2L2Nu_","WWTo2L2Nu_",bin[0],bin[1],bin[2]);
	WWTo2L2Nu_->Sumw2();
    
    TH1F * ZZTo2L2Nu_ = new TH1F("ZZTo2L2Nu_","ZZTo2L2Nu_",bin[0],bin[1],bin[2]);
	ZZTo2L2Nu_->Sumw2();
    
    TH1F * WWTo1L1Nu2Q_ = new TH1F("WWTo1L1Nu2Q_","WWTo1L1Nu2Q_",bin[0],bin[1],bin[2]);
	WWTo1L1Nu2Q_->Sumw2();

    TH1F * WZTo1L1Nu2Q_ = new TH1F("WZTo1L1Nu2Q_","WZTo1L1Nu2Q_",bin[0],bin[1],bin[2]);
	WZTo1L1Nu2Q_->Sumw2();
    
    TH1F * WZTo1L3Nu_ = new TH1F("WZTo1L3Nu_","WZTo1L3Nu_",bin[0],bin[1],bin[2]);
	WZTo1L3Nu_->Sumw2();
    
    TH1F * WZTo2L2Q_ = new TH1F("WZTo2L2Q_","WZTo2L2Q_",bin[0],bin[1],bin[2]);
	WZTo2L2Q_->Sumw2();
    
	TH1F * ZZTo2L2Q_ = new TH1F("ZZTo2L2Q_","ZZTo2L2Q_",bin[0],bin[1],bin[2]);
	ZZTo2L2Q_->Sumw2();
    
	TH1F * ZZTo2Q2Nu_ = new TH1F("ZZTo2Q2Nu_","ZZTo2Q2Nu_",bin[0],bin[1],bin[2]);
	ZZTo2Q2Nu_->Sumw2();
    
	TH1F * GluGluHTauTau_ = new TH1F("GluGluHTauTau_","GluGluHTauTau_",bin[0],bin[1],bin[2]);
	GluGluHTauTau_->Sumw2();
    
	TH1F * VBFHTauTau_ = new TH1F("VBFHTauTau_","VBFHTauTau_",bin[0],bin[1],bin[2]);
	VBFHTauTau_->Sumw2();
    
    TH1F * ZHTauTau_ = new TH1F("ZHTauTau_","ZHTauTau_",bin[0],bin[1],bin[2]);
	ZHTauTau_->Sumw2();
    
    TH1F * comb_ = new TH1F("comb_","comb_",bin[0],bin[1],bin[2]); 
	comb_->Sumw2();
    
    DYinc->Draw((parameter+">>DYinc_").c_str(),cut*weights["DY"].c_str());
    DY1Jets->Draw((parameter+">>DY1Jets_").c_str(),cut*weights["DY"].c_str());
    DY2Jets->Draw((parameter+">>DY2Jets_").c_str(),cut*weights["DY"].c_str());
    DY3Jets->Draw((parameter+">>DY3Jets_").c_str(),cut*weights["DY"].c_str());
    DY4Jets->Draw((parameter+">>DY4Jets_").c_str(),cut*weights["DY"].c_str());
    WJetsToLNu->Draw((parameter+">>WJetsToLNu_").c_str(),cut*weights["W"].c_str());
    W1JetsToLNu->Draw((parameter+">>W1JetsToLNu_").c_str(),cut*weights["W"].c_str());
    W2JetsToLNu->Draw((parameter+">>W2JetsToLNu_").c_str(),cut*weights["W"].c_str());
    W3JetsToLNu->Draw((parameter+">>W3JetsToLNu_").c_str(),cut*weights["W"].c_str());
    W4JetsToLNu->Draw((parameter+">>W4JetsToLNu_").c_str(),cut*weights["W"].c_str());
    TT->Draw((parameter+">>TT_").c_str(),cut*weights["TTT"].c_str());
    GluGluHTauTau->Draw((parameter+">>GluGluHTauTau_").c_str(),cut*weights["GluGluHTauTau"].c_str());
    VBFHTauTau->Draw((parameter+">>VBFHTauTau_").c_str(),cut*weights["VBFHTauTau"].c_str());
    ZHTauTau->Draw((parameter+">>ZHTauTau_").c_str(),cut*weights["ZHTauTau"].c_str());
    //ST_t_channel_antitop_4f_leptonDecays->Draw((parameter+">>ST_t_channel_antitop_4f_leptonDecays_").c_str(),cut*weights["VV"].c_str());
    //ST_t_channel_top_4f_leptonDecays->Draw((parameter+">>ST_t_channel_top_4f_leptonDecays_").c_str(),cut*weights["VV"].c_str());
    ST_tW_antitop_5f_inclusiveDecays->Draw((parameter+">>ST_tW_antitop_5f_inclusiveDecays_").c_str(),cut*weights["VV"].c_str());
    ST_tW_top_5f_inclusiveDecays->Draw((parameter+">>ST_tW_top_5f_inclusiveDecays_").c_str(),cut*weights["VV"].c_str());
    WWTo2L2Nu->Draw((parameter+">>WWTo2L2Nu_").c_str(),cut*weights["VV"].c_str());
    ZZTo2L2Nu->Draw((parameter+">>ZZTo2L2Nu_").c_str(),cut*weights["VV"].c_str());
    WWTo1L1Nu2Q->Draw((parameter+">>WWTo1L1Nu2Q_").c_str(),cut*weights["VV"].c_str());
    WZTo1L1Nu2Q->Draw((parameter+">>WZTo1L1Nu2Q_").c_str(),cut*weights["VV"].c_str());
    WZTo1L3Nu->Draw((parameter+">>WZTo1L3Nu_").c_str(),cut*weights["VV"].c_str());
    WZTo2L2Q->Draw((parameter+">>WZTo2L2Q_").c_str(),cut*weights["VV"].c_str());
    ZZTo2L2Q->Draw((parameter+">>ZZTo2L2Q_").c_str(),cut*weights["VV"].c_str());
    ZZTo2Q2Nu->Draw((parameter+">>ZZTo2Q2Nu_").c_str(),cut*weights["VV"].c_str());

    //not scaling by sf in datacards
    /*
    WJetsToLNu_->Scale(sf);
    W1JetsToLNu_->Scale(sf);
    W2JetsToLNu_->Scale(sf);
    W3JetsToLNu_->Scale(sf);
    W4JetsToLNu_->Scale(sf);
    */
    
    WJetsToLNu_->Scale(lumi_sf);
    W1JetsToLNu_->Scale(lumi_sf);
    W2JetsToLNu_->Scale(lumi_sf);
    W3JetsToLNu_->Scale(lumi_sf);
    W4JetsToLNu_->Scale(lumi_sf);
    DYhigh_->Scale(lumi_sf);
    DYlow_->Scale(lumi_sf);
    DYinc_->Scale(lumi_sf);
    DY1Jets_->Scale(lumi_sf);
    DY2Jets_->Scale(lumi_sf);
    DY3Jets_->Scale(lumi_sf);
    DY4Jets_->Scale(lumi_sf);
    TT_->Scale(lumi_sf);
    GluGluHTauTau_->Scale(lumi_sf);
    VBFHTauTau_->Scale(lumi_sf);
    ZHTauTau_->Scale(lumi_sf);
    ST_t_channel_antitop_4f_leptonDecays_->Scale(lumi_sf);
    ST_t_channel_top_4f_leptonDecays_->Scale(lumi_sf);
    ST_tW_antitop_5f_inclusiveDecays_->Scale(lumi_sf);
    ST_tW_top_5f_inclusiveDecays_->Scale(lumi_sf);
    WWTo2L2Nu_->Scale(lumi_sf);
    ZZTo2L2Nu_->Scale(lumi_sf);
    WWTo1L1Nu2Q_->Scale(lumi_sf);
    WZTo1L1Nu2Q_->Scale(lumi_sf);
    WZTo1L3Nu_->Scale(lumi_sf);
    WZTo2L2Q_->Scale(lumi_sf);
    ZZTo2L2Q_->Scale(lumi_sf);
    ZZTo2Q2Nu_->Scale(lumi_sf);
    
    comb_->Add(WJetsToLNu_);
    comb_->Add(W1JetsToLNu_);
    comb_->Add(W2JetsToLNu_);
    comb_->Add(W3JetsToLNu_);
    comb_->Add(W4JetsToLNu_);
    comb_->Add(DYhigh_);
    comb_->Add(DYlow_);
    comb_->Add(DYinc_);
    comb_->Add(DY1Jets_);
    comb_->Add(DY2Jets_);
    comb_->Add(DY3Jets_);
    comb_->Add(DY4Jets_);
    comb_->Add(TT_);
    comb_->Add(GluGluHTauTau_);
    comb_->Add(VBFHTauTau_);
    comb_->Add(ZHTauTau_);
    comb_->Add(ST_t_channel_antitop_4f_leptonDecays_);
    comb_->Add(ST_t_channel_top_4f_leptonDecays_);
    comb_->Add(ST_tW_antitop_5f_inclusiveDecays_);
    comb_->Add(ST_tW_top_5f_inclusiveDecays_);
    comb_->Add(WWTo2L2Nu_);
    comb_->Add(ZZTo2L2Nu_);
    comb_->Add(WWTo1L1Nu2Q_);
    comb_->Add(WZTo1L1Nu2Q_);
    comb_->Add(WZTo1L3Nu_);
    comb_->Add(WZTo2L2Q_);
    comb_->Add(ZZTo2L2Q_);
    comb_->Add(ZZTo2Q2Nu_);

	double error = 0;
    
    double totalBKGIntegral = comb_->Integral();
    
    std::cout << "\n";
    std::cout << "BKG Totals and Fractions" << "\n";
    std::cout << "\n";

	error = 0;
	WJetsToLNu_->IntegralAndError(-1,1e3, error);
	std::cout<<" WJetsToLNu................................"<<WJetsToLNu_->Integral()<<" +/- "<<error<<  "   Fraction: " <<  WJetsToLNu_->Integral()/totalBKGIntegral<< "\n";
	error = 0;
	W1JetsToLNu_->IntegralAndError(-1,1e3, error);
	std::cout<<" W1JetsToLNu..............................."<<W1JetsToLNu_->Integral()<<" +/- "<<error<<  "   Fraction: " <<  W1JetsToLNu_->Integral()/totalBKGIntegral<< "\n";
	error = 0;
	W2JetsToLNu_->IntegralAndError(-1,1e3, error);
	std::cout<<" W2JetsToLNu..............................."<<W2JetsToLNu_->Integral()<<" +/- "<<error<<  "   Fraction: " <<  W2JetsToLNu_->Integral()/totalBKGIntegral<< "\n";
	error = 0;
	W3JetsToLNu_->IntegralAndError(-1,1e3, error);
	std::cout<<" W3JetsToLNu..............................."<<W3JetsToLNu_->Integral()<<" +/- "<<error<<  "   Fraction: " <<  W3JetsToLNu_->Integral()/totalBKGIntegral<< "\n";
	error = 0;
	W4JetsToLNu_->IntegralAndError(-1,1e3, error);
	std::cout<<" W4JetsToLNu..............................."<<W4JetsToLNu_->Integral()<<" +/- "<<error<<  "   Fraction: " <<  W4JetsToLNu_->Integral()/totalBKGIntegral<< "\n";
    error = 0;
    DYhigh_->IntegralAndError(-1,1e3, error);
    std::cout<<" DYhigh...................................."<< DYhigh_->Integral()<<" +/- "<<error<<  "   Fraction: " <<  DYhigh_->Integral()/totalBKGIntegral<< "\n";
    error = 0;
    DYlow_->IntegralAndError(-1,1e3, error);
    std::cout<<" DYlow....................................."<< DYlow_->Integral()<<" +/- "<<error<< "   Fraction: " <<  DYlow_->Integral()/totalBKGIntegral<< "\n";
    error = 0;
    DYinc_->IntegralAndError(-1,1e3, error);
    std::cout<<" DYinc....................................."<< DYinc_->Integral()<<" +/- "<<error<< "   Fraction: " <<  DYinc_->Integral()/totalBKGIntegral<< "\n";
    error = 0;
    DY1Jets_->IntegralAndError(-1,1e3, error);
    std::cout<<" DY1Jets..................................."<< DY1Jets_->Integral()<<" +/- "<<error<< "   Fraction: " <<  DY1Jets_->Integral()/totalBKGIntegral<< "\n";
    error = 0;
    DY2Jets_->IntegralAndError(-1,1e3, error);
    std::cout<<" DY2Jets..................................."<< DY2Jets_->Integral()<<" +/- "<<error<< "   Fraction: " <<  DY2Jets_->Integral()/totalBKGIntegral<< "\n";
    error = 0;
    DY3Jets_->IntegralAndError(-1,1e3, error);
    std::cout<<" DY3Jets..................................."<< DY3Jets_->Integral()<<" +/- "<<error<< "   Fraction: " <<  DY3Jets_->Integral()/totalBKGIntegral<< "\n";
    error = 0;
    DY4Jets_->IntegralAndError(-1,1e3, error);
    std::cout<<" DY4Jets..................................."<< DY4Jets_->Integral()<<" +/- "<<error<< "   Fraction: " <<  DY4Jets_->Integral()/totalBKGIntegral<< "\n";
    error = 0;
    TT_->IntegralAndError(-1,1e3, error);
    std::cout<<" TT........................................"<< TT_->Integral()<<" +/- "<<error<< "   Fraction: " <<  TT_->Integral()/totalBKGIntegral<< "\n";
    error = 0;
    GluGluHTauTau_->IntegralAndError(-1,1e3, error);
    std::cout<<" GluGluHTauTau............................."<< GluGluHTauTau_->Integral()<<" +/- "<<error<< "   Fraction: " <<  GluGluHTauTau_->Integral()/totalBKGIntegral<< "\n";
    error = 0;
    VBFHTauTau_->IntegralAndError(-1,1e3, error);
    std::cout<<" VBFHTauTau................................"<< VBFHTauTau_->Integral()<<" +/- "<<error<< "   Fraction: " <<  VBFHTauTau_->Integral()/totalBKGIntegral<< "\n";
    error = 0;
    ZHTauTau_->IntegralAndError(-1,1e3, error);
    std::cout<<" ZHTauTau.................................."<< ZHTauTau_->Integral()<<" +/- "<<error<<  "   Fraction: " <<  ZHTauTau_->Integral()/totalBKGIntegral<< "\n";
    error = 0;
    ST_t_channel_antitop_4f_leptonDecays_->IntegralAndError(-1,1e3, error);
    std::cout<<" ST_t_channel_antitop_4f_leptonDecays......"<< ST_t_channel_antitop_4f_leptonDecays_->Integral()<<" +/- "<<error<<  "   Fraction: " <<  ST_t_channel_antitop_4f_leptonDecays_->Integral()/totalBKGIntegral<< "\n";
    error = 0;
    ST_t_channel_top_4f_leptonDecays_->IntegralAndError(-1,1e3, error);
    std::cout<<" ST_t_channel_top_4f_leptonDecays.........."<< ST_t_channel_top_4f_leptonDecays_->Integral()<<" +/- "<<error<<  "   Fraction: " <<  ST_t_channel_top_4f_leptonDecays_->Integral()/totalBKGIntegral<< "\n";
    error = 0;
    ST_tW_antitop_5f_inclusiveDecays_->IntegralAndError(-1,1e3, error);
    std::cout<<" ST_tW_antitop_5f_inclusiveDecays.........."<< ST_tW_antitop_5f_inclusiveDecays_->Integral()<<" +/- "<<error<<  "   Fraction: " <<  ST_tW_antitop_5f_inclusiveDecays_->Integral()/totalBKGIntegral<< "\n";
    error = 0;
    ST_tW_top_5f_inclusiveDecays_->IntegralAndError(-1,1e3, error);
    std::cout<<" ST_tW_top_5f_inclusiveDecays.............."<< ST_tW_top_5f_inclusiveDecays_->Integral()<<" +/- "<<error<<  "   Fraction: " <<  ST_tW_top_5f_inclusiveDecays_->Integral()/totalBKGIntegral<< "\n";
    error = 0;
    WWTo2L2Nu_->IntegralAndError(-1,1e3, error);
    std::cout<<" WWTo2L2Nu................................."<< WWTo2L2Nu_->Integral()<<" +/- "<<error<<  "   Fraction: " <<  WWTo2L2Nu_->Integral()/totalBKGIntegral<< "\n";
    error = 0;
    ZZTo2L2Nu_->IntegralAndError(-1,1e3, error);
    std::cout<<" ZZTo2L2Nu................................."<< ZZTo2L2Nu_->Integral()<<" +/- "<<error<<  "   Fraction: " <<  ZZTo2L2Nu_->Integral()/totalBKGIntegral<< "\n";
    error = 0;
    WWTo1L1Nu2Q_->IntegralAndError(-1,1e3, error);
    std::cout<<" WWTo1L1Nu2Q..............................."<< WWTo1L1Nu2Q_->Integral()<<" +/- "<<error<<  "   Fraction: " <<  WWTo1L1Nu2Q_->Integral()/totalBKGIntegral<< "\n";
    error = 0;
    WZTo1L1Nu2Q_->IntegralAndError(-1,1e3, error);
    std::cout<<" WZTo1L1Nu2Q..............................."<< WZTo1L1Nu2Q_->Integral()<<" +/- "<<error<<  "   Fraction: " <<  WZTo1L1Nu2Q_->Integral()/totalBKGIntegral<< "\n";
    error = 0;
    WZTo1L3Nu_->IntegralAndError(-1,1e3, error);
    std::cout<<" WZTo1L3Nu................................."<< WZTo1L3Nu_->Integral()<<" +/- "<<error<<  "   Fraction: " <<  WZTo1L3Nu_->Integral()/totalBKGIntegral<< "\n";
    error = 0;
    WZTo2L2Q_->IntegralAndError(-1,1e3, error);
    std::cout<<" WZTo2L2Q.................................."<< WZTo2L2Q_->Integral()<<" +/- "<<error<<  "   Fraction: " <<  WZTo2L2Q_->Integral()/totalBKGIntegral<< "\n";
    error = 0;
    ZZTo2L2Q_->IntegralAndError(-1,1e3, error);
    std::cout<<" ZZTo2L2Q.................................."<< ZZTo2L2Q_->Integral()<<" +/- "<<error<<  "   Fraction: " <<  ZZTo2L2Q_->Integral()/totalBKGIntegral<< "\n";
    error = 0;
    ZZTo2Q2Nu_->IntegralAndError(-1,1e3, error);
    std::cout<<" ZZTo2Q2Nu...................................."<< ZZTo2Q2Nu_->Integral()<<" +/- "<<error<< "   Fraction: " <<  ZZTo2Q2Nu_->Integral()/totalBKGIntegral<< "\n";
    
    error = 0;
    comb_->IntegralAndError(-1,1e3, error);
    std::cout<<" Combined Background Total................."<< comb_->Integral()<<" +/- "<<error<< "\n";
    
    delete WJetsToLNu_;
    delete W1JetsToLNu_;
    delete W2JetsToLNu_;
    delete W3JetsToLNu_;
    delete W4JetsToLNu_;
    delete DYhigh_;
    delete DYlow_;
    delete DYinc_;
    delete DY1Jets_;
    delete DY2Jets_;
    delete DY3Jets_;
    delete DY4Jets_;
    delete TT_;
    delete GluGluHTauTau_;
    delete VBFHTauTau_;
    delete ZHTauTau_;
    delete ST_t_channel_antitop_4f_leptonDecays_;
    delete ST_t_channel_top_4f_leptonDecays_;
    delete ST_tW_antitop_5f_inclusiveDecays_;
    delete ST_tW_top_5f_inclusiveDecays_;
    delete WWTo2L2Nu_;
    delete ZZTo2L2Nu_;
    delete WWTo1L1Nu2Q_;
    delete WZTo1L1Nu2Q_;
    delete WZTo1L3Nu_;
    delete WZTo2L2Q_;
    delete ZZTo2L2Q_;
    delete ZZTo2Q2Nu_;
    delete comb_;
}


void drawSignalRegion(TCut cut, TCut relCut, TCut SScut, TCut SScutRel, double sfQCD, double sf, std::string parameter, float bin[3], std::string can_name, std::string chan, std::string syst, bool plot, bool createOutputShapes)
{
    std::cout << "Drawing Parameter: " << parameter << std::endl;

    std::cout << "Yield CUT: " << cut.GetTitle() << std::endl;

	/* check if want blinded signal region met tail */
    if (plot)
    {
        if(choice2==1 && (parameter=="mt_tot"))
        {
            global_title = "Signal Region with Blinded Tail";
        }
        if(choice2==1 && (parameter=="met" || parameter=="mvamet"))
        {
            global_title = "Signal Region with Blind MET Tail";
        }
        if(choice2==1 && sigChoice==2 && (chan=="mt"||chan=="et" ||chan=="tt"))
        {
            global_title = "Signal Region in TMVA (blinded above 0.5)";
        }
    }
    
	/* create a TCanvas */
	TCanvas * Signal_region = new TCanvas(can_name.c_str(),can_name.c_str(),canDim[0],canDim[1]);
	Signal_region->Divide(1,2);

	/* create some histograms and format them */

	TH1F * data_ = new TH1F("data_","",bin[0],bin[1],bin[2]); 
	data_->Sumw2(); 
	data_->SetMarkerStyle(dataMarker);

	TH1F * qcd_ = new TH1F("qcd_","qcd_",bin[0],bin[1],bin[2]); 
	qcd_->Sumw2(); 
	qcd_->SetFillColor(colors["QCD"]);

	TH1F * ztt_ = new TH1F("ztt_","ztt_",bin[0],bin[1],bin[2]);
	ztt_->Sumw2();
	ztt_->SetFillColor(colors["ZTT"]);
    
    TH1F * zl_ = new TH1F("zl_","zl_",bin[0],bin[1],bin[2]);
	zl_->Sumw2();
	zl_->SetFillColor(colors["ZL"]);
    
    TH1F * zj_ = new TH1F("zj_","zj_",bin[0],bin[1],bin[2]);
	zj_->Sumw2();
	zj_->SetFillColor(colors["ZJ"]);
	
	TH1F * w_ = new TH1F("w_","w_",bin[0],bin[1],bin[2]); 
	w_->Sumw2(); 
	w_->SetFillColor(colors["W"]);
	
	TH1F * ttt_ = new TH1F("ttt_","ttt_",bin[0],bin[1],bin[2]);
	ttt_->Sumw2();
	ttt_->SetFillColor(colors["TTT"]);

    TH1F * ttj_ = new TH1F("ttj_","ttj_",bin[0],bin[1],bin[2]);
	ttj_->Sumw2();
	ttj_->SetFillColor(colors["TTJ"]);

	TH1F * vvt_ = new TH1F("vvt_","vvt_",bin[0],bin[1],bin[2]);
	vvt_->Sumw2();
	vvt_->SetFillColor(colors["VV"]);
    
    TH1F * vvj_ = new TH1F("vvj_","vvj_",bin[0],bin[1],bin[2]);
	vvj_->Sumw2();
	vvj_->SetFillColor(colors["VV"]);
    
    TH1F * ewk_ = new TH1F("ewk_","ewk_",bin[0],bin[1],bin[2]);
	ewk_->Sumw2();
	ewk_->SetFillColor(colors["EWK"]);
    
    TH1F * zvv_ = new TH1F("zvv_","zvv_",bin[0],bin[1],bin[2]);
	zvv_->Sumw2();
	zvv_->SetFillColor(colors["ZVV"]);
	
	TH1F * zhtautau_ = new TH1F("zhtautau_","zhtautau_",bin[0],bin[1],bin[2]);
	zhtautau_->Sumw2();
	zhtautau_->SetFillColor(colors["ZHTauTau"]);

	TH1F * glugluhtautau_ = new TH1F("glugluhtautau_","glugluhtautau_",bin[0],bin[1],bin[2]);
	glugluhtautau_->Sumw2();
	glugluhtautau_->SetFillColor(colors["GluGluHTauTau"]);
    
    TH1F * vbfhtautau_ = new TH1F("vbfhtautau_","vbfhtautau_",bin[0],bin[1],bin[2]);
	vbfhtautau_->Sumw2();
    vbfhtautau_->SetFillColor(colors["VBFHTauTau"]);
    
	TH1F * zp600a300_ = new TH1F("zp600a300_","zp600a300_",bin[0],bin[1],bin[2]);
	zp600a300_->Sumw2();
    TH1F * zp600a400_ = new TH1F("zp600a400_","zp600a400_",bin[0],bin[1],bin[2]);
	zp600a400_->Sumw2();

	TH1F * zp800a300_ = new TH1F("zp800a300_","zp800a300_",bin[0],bin[1],bin[2]);
	zp800a300_->Sumw2();
    TH1F * zp800a400_ = new TH1F("zp800a400_","zp800a400_",bin[0],bin[1],bin[2]);
	zp800a400_->Sumw2();
    TH1F * zp800a500_ = new TH1F("zp800a500_","zp800a500_",bin[0],bin[1],bin[2]);
	zp800a500_->Sumw2();
    TH1F * zp800a600_ = new TH1F("zp800a600_","zp800a600_",bin[0],bin[1],bin[2]);
	zp800a600_->Sumw2();

	TH1F * zp1000a300_ = new TH1F("zp1000a300_","zp1000a300_",bin[0],bin[1],bin[2]);
	zp1000a300_->Sumw2();
    TH1F * zp1000a400_ = new TH1F("zp1000a400_","zp1000a400_",bin[0],bin[1],bin[2]);
	zp1000a400_->Sumw2();
    TH1F * zp1000a500_ = new TH1F("zp1000a500_","zp1000a500_",bin[0],bin[1],bin[2]);
	zp1000a500_->Sumw2();
    TH1F * zp1000a600_ = new TH1F("zp1000a600_","zp1000a600_",bin[0],bin[1],bin[2]);
	zp1000a600_->Sumw2();
	TH1F * zp1000a700_ = new TH1F("zp1000a700_","zp1000a700_",bin[0],bin[1],bin[2]);
	zp1000a700_->Sumw2();
    TH1F * zp1000a800_ = new TH1F("zp1000a800_","zp1000a800_",bin[0],bin[1],bin[2]);
	zp1000a800_->Sumw2();
    
    TH1F * zp1200a300_ = new TH1F("zp1200a300_","zp1200a300_",bin[0],bin[1],bin[2]);
	zp1200a300_->Sumw2();
    TH1F * zp1200a400_ = new TH1F("zp1200a400_","zp1200a400_",bin[0],bin[1],bin[2]);
	zp1200a400_->Sumw2();
    TH1F * zp1200a500_ = new TH1F("zp1200a500_","zp1200a500_",bin[0],bin[1],bin[2]);
	zp1200a500_->Sumw2();
    TH1F * zp1200a600_ = new TH1F("zp1200a600_","zp1200a600_",bin[0],bin[1],bin[2]);
	zp1200a600_->Sumw2();
	TH1F * zp1200a700_ = new TH1F("zp1200a700_","zp1200a700_",bin[0],bin[1],bin[2]);
	zp1200a700_->Sumw2();
    TH1F * zp1200a800_ = new TH1F("zp1200a800_","zp1200a800_",bin[0],bin[1],bin[2]);
	zp1200a800_->Sumw2();
    
    TH1F * zp1400a300_ = new TH1F("zp1400a300_","zp1400a300_",bin[0],bin[1],bin[2]);
	zp1400a300_->Sumw2();
    TH1F * zp1400a400_ = new TH1F("zp1400a400_","zp1400a400_",bin[0],bin[1],bin[2]);
	zp1400a400_->Sumw2();
    TH1F * zp1400a500_ = new TH1F("zp1400a500_","zp1400a500_",bin[0],bin[1],bin[2]);
	zp1400a500_->Sumw2();
    TH1F * zp1400a600_ = new TH1F("zp1400a600_","zp1400a600_",bin[0],bin[1],bin[2]);
	zp1400a600_->Sumw2();
	TH1F * zp1400a700_ = new TH1F("zp1400a700_","zp1400a700_",bin[0],bin[1],bin[2]);
	zp1400a700_->Sumw2();
    TH1F * zp1400a800_ = new TH1F("zp1400a800_","zp1400a800_",bin[0],bin[1],bin[2]);
	zp1400a800_->Sumw2();
    
    TH1F * zp1700a300_ = new TH1F("zp1700a300_","zp1700a300_",bin[0],bin[1],bin[2]);
	zp1700a300_->Sumw2();
    TH1F * zp1700a400_ = new TH1F("zp1700a400_","zp1700a400_",bin[0],bin[1],bin[2]);
	zp1700a400_->Sumw2();
    TH1F * zp1700a500_ = new TH1F("zp1700a500_","zp1700a500_",bin[0],bin[1],bin[2]);
	zp1700a500_->Sumw2();
    TH1F * zp1700a600_ = new TH1F("zp1700a600_","zp1700a600_",bin[0],bin[1],bin[2]);
	zp1700a600_->Sumw2();
	TH1F * zp1700a700_ = new TH1F("zp1700a700_","zp1700a700_",bin[0],bin[1],bin[2]);
	zp1700a700_->Sumw2();
    TH1F * zp1700a800_ = new TH1F("zp1700a800_","zp1700a800_",bin[0],bin[1],bin[2]);
	zp1700a800_->Sumw2();
    
    TH1F * zp2000a300_ = new TH1F("zp2000a300_","zp2000a300_",bin[0],bin[1],bin[2]);
	zp2000a300_->Sumw2();
    TH1F * zp2000a400_ = new TH1F("zp2000a400_","zp2000a400_",bin[0],bin[1],bin[2]);
	zp2000a400_->Sumw2();
    TH1F * zp2000a500_ = new TH1F("zp2000a500_","zp2000a500_",bin[0],bin[1],bin[2]);
	zp2000a500_->Sumw2();
    TH1F * zp2000a600_ = new TH1F("zp2000a600_","zp2000a600_",bin[0],bin[1],bin[2]);
	zp2000a600_->Sumw2();
	TH1F * zp2000a700_ = new TH1F("zp2000a700_","zp2000a700_",bin[0],bin[1],bin[2]);
	zp2000a700_->Sumw2();
    TH1F * zp2000a800_ = new TH1F("zp2000a800_","zp2000a800_",bin[0],bin[1],bin[2]);
	zp2000a800_->Sumw2();
    
    TH1F * zp2500a300_ = new TH1F("zp2500a300_","zp2500a300_",bin[0],bin[1],bin[2]);
	zp2500a300_->Sumw2();
    TH1F * zp2500a400_ = new TH1F("zp2500a400_","zp2500a400_",bin[0],bin[1],bin[2]);
	zp2500a400_->Sumw2();
    TH1F * zp2500a500_ = new TH1F("zp2500a500_","zp2500a500_",bin[0],bin[1],bin[2]);
	zp2500a500_->Sumw2();
    TH1F * zp2500a600_ = new TH1F("zp2500a600_","zp2500a600_",bin[0],bin[1],bin[2]);
	zp2500a600_->Sumw2();
	TH1F * zp2500a700_ = new TH1F("zp2500a700_","zp2500a700_",bin[0],bin[1],bin[2]);
	zp2500a700_->Sumw2();
    TH1F * zp2500a800_ = new TH1F("zp2500a800_","zp2500a800_",bin[0],bin[1],bin[2]);
	zp2500a800_->Sumw2();
    
    TH1F * zp10x1_ = new TH1F("zp10x1_","zp10x1_",bin[0],bin[1],bin[2]);
	zp10x1_->Sumw2();
    TH1F * zp20x1_ = new TH1F("zp20x1_","zp20x1_",bin[0],bin[1],bin[2]);
	zp20x1_->Sumw2();
    TH1F * zp50x1_ = new TH1F("zp50x1_","zp50x1_",bin[0],bin[1],bin[2]);
	zp50x1_->Sumw2();
    TH1F * zp100x1_ = new TH1F("zp100x1_","zp100x1_",bin[0],bin[1],bin[2]);
	zp100x1_->Sumw2();
    TH1F * zp300x1_ = new TH1F("zp300x1_","zp300x1_",bin[0],bin[1],bin[2]);
	zp300x1_->Sumw2();
    TH1F * zp500x1_ = new TH1F("zp500x1_","zp500x1_",bin[0],bin[1],bin[2]);
	zp500x1_->Sumw2();
    TH1F * zp1000x1_ = new TH1F("zp1000x1_","zp1000x1_",bin[0],bin[1],bin[2]);
	zp1000x1_->Sumw2();
    TH1F * zp2000x1_ = new TH1F("zp2000x1_","zp2000x1_",bin[0],bin[1],bin[2]);
	zp2000x1_->Sumw2();
    TH1F * zp10000x1_ = new TH1F("zp10000x1_","zp10000x1_",bin[0],bin[1],bin[2]);
	zp10000x1_->Sumw2();
    
    TH1F * zp10x50_ = new TH1F("zp10x50_","zp10x50_",bin[0],bin[1],bin[2]);
	zp10x50_->Sumw2();
    TH1F * zp50x50_ = new TH1F("zp50x50_","zp50x50_",bin[0],bin[1],bin[2]);
	zp50x50_->Sumw2();
    TH1F * zp95x50_ = new TH1F("zp95x50_","zp95x50_",bin[0],bin[1],bin[2]);
	zp95x50_->Sumw2();
    TH1F * zp200x50_ = new TH1F("zp200x50_","zp200x50_",bin[0],bin[1],bin[2]);
	zp200x50_->Sumw2();
    TH1F * zp300x50_ = new TH1F("zp300x50_","zp300x50_",bin[0],bin[1],bin[2]);
	zp300x50_->Sumw2();
    TH1F * zp10000x50_ = new TH1F("zp10000x50_","zp10000x50_",bin[0],bin[1],bin[2]);
	zp10000x50_->Sumw2();

    TH1F * zp10x150_ = new TH1F("zp10x150_","zp10x150_",bin[0],bin[1],bin[2]);
	zp10x150_->Sumw2();
    TH1F * zp200x150_ = new TH1F("zp200x150_","zp200x150_",bin[0],bin[1],bin[2]);
	zp200x150_->Sumw2();
    TH1F * zp295x150_ = new TH1F("zp295x150_","zp295x150_",bin[0],bin[1],bin[2]);
	zp295x150_->Sumw2();
    TH1F * zp500x150_ = new TH1F("zp500x150_","zp500x150_",bin[0],bin[1],bin[2]);
	zp500x150_->Sumw2();
    TH1F * zp1000x150_ = new TH1F("zp1000x150_","zp1000x150_",bin[0],bin[1],bin[2]);
	zp1000x150_->Sumw2();
    TH1F * zp10000x150_ = new TH1F("zp10000x150_","zp10000x150_",bin[0],bin[1],bin[2]);
	zp10000x150_->Sumw2();
    
    TH1F * zp10x500_ = new TH1F("zp10x500_","zp10x500_",bin[0],bin[1],bin[2]);
	zp10x500_->Sumw2();
    TH1F * zp500x500_ = new TH1F("zp500x500_","zp500x500_",bin[0],bin[1],bin[2]);
	zp500x500_->Sumw2();
    TH1F * zp995x500_ = new TH1F("zp995x500_","zp995x500_",bin[0],bin[1],bin[2]);
	zp995x500_->Sumw2();
    TH1F * zp2000x500_ = new TH1F("zp2000x500_","zp2000x500_",bin[0],bin[1],bin[2]);
	zp2000x500_->Sumw2();
    TH1F * zp10000x500_ = new TH1F("zp10000x500_","zp10000x500_",bin[0],bin[1],bin[2]);
	zp10000x500_->Sumw2();
    
    TH1F * zp10x1000_ = new TH1F("zp10x1000_","zp10x1000_",bin[0],bin[1],bin[2]);
	zp10x1000_->Sumw2();
    TH1F * zp1000x1000_ = new TH1F("zp1000x1000_","zp1000x1000_",bin[0],bin[1],bin[2]);
	zp1000x1000_->Sumw2();
    TH1F * zp1995x1000_ = new TH1F("zp1995x1000_","zp1995x1000_",bin[0],bin[1],bin[2]);
	zp1995x1000_->Sumw2();
    TH1F * zp10000x1000_ = new TH1F("zp10000x1000_","zp10000x1000_",bin[0],bin[1],bin[2]);
	zp10000x1000_->Sumw2();

	/* sum of all bks */
	TH1F * comb_ = new TH1F("comb_","comb_",bin[0],bin[1],bin[2]); 
	comb_->Sumw2(); 

	TH1F * smh_ = new TH1F("smh_","smh_",bin[0],bin[1],bin[2]);
	smh_->Sumw2();
    smh_->SetFillColor(colors["SMH"]);
    
	/* fill the hists */
	if (plot) {Signal_region->cd(1);}

	DATA->Draw((parameter+">>data_").c_str(),cut);
    
	W->Draw((parameter+">>w_").c_str(),relCut*weights["W"].c_str());
	ZTT->Draw((parameter+">>ztt_").c_str(),relCut*weights["ZTT"].c_str());
    ZL->Draw((parameter+">>zl_").c_str(),relCut*weights["ZL"].c_str());
    ZJ->Draw((parameter+">>zj_").c_str(),relCut*weights["ZJ"].c_str());
	VVT->Draw((parameter+">>vvt_").c_str(),relCut*weights["VV"].c_str());
    VVJ->Draw((parameter+">>vvj_").c_str(),relCut*weights["VV"].c_str());
    ZVV->Draw((parameter+">>zvv_").c_str(),cut*weights["ZVV"].c_str());
    EWK->Draw((parameter+">>ewk_").c_str(),cut*weights["EWK"].c_str());
	TTT->Draw((parameter+">>ttt_").c_str(),relCut*weights["TTT"].c_str());
    TTJ->Draw((parameter+">>ttj_").c_str(),relCut*weights["TTJ"].c_str());
	ZHTauTau->Draw((parameter+">>zhtautau_").c_str(),cut*weights["ZHTauTau"].c_str());
    VBFHTauTau->Draw((parameter+">>vbfhtautau_").c_str(),cut*weights["VBFHTauTau"].c_str());
    GluGluHTauTau->Draw((parameter+">>glugluhtautau_").c_str(),cut*weights["GluGluHTauTau"].c_str());
    
    smh_->Add(zhtautau_);
    smh_->Add(glugluhtautau_);
    smh_->Add(vbfhtautau_);
    
	MZP600_MA0300->Draw((parameter+">>zp600a300_").c_str(),cut*weights["MONO"].c_str());
    MZP600_MA0400->Draw((parameter+">>zp600a400_").c_str(),cut*weights["MONO"].c_str());
    
	MZP800_MA0300->Draw((parameter+">>zp800a300_").c_str(),cut*weights["MONO"].c_str());
    MZP800_MA0400->Draw((parameter+">>zp800a400_").c_str(),cut*weights["MONO"].c_str());
    MZP800_MA0500->Draw((parameter+">>zp800a500_").c_str(),cut*weights["MONO"].c_str());
    MZP800_MA0600->Draw((parameter+">>zp800a600_").c_str(),cut*weights["MONO"].c_str());
    
	MZP1000_MA0300->Draw((parameter+">>zp1000a300_").c_str(),cut*weights["MONO"].c_str());
    MZP1000_MA0400->Draw((parameter+">>zp1000a400_").c_str(),cut*weights["MONO"].c_str());
    MZP1000_MA0500->Draw((parameter+">>zp1000a500_").c_str(),cut*weights["MONO"].c_str());
    MZP1000_MA0600->Draw((parameter+">>zp1000a600_").c_str(),cut*weights["MONO"].c_str());
    MZP1000_MA0700->Draw((parameter+">>zp1000a700_").c_str(),cut*weights["MONO"].c_str());
    MZP1000_MA0800->Draw((parameter+">>zp1000a800_").c_str(),cut*weights["MONO"].c_str());
    
    MZP1200_MA0300->Draw((parameter+">>zp1200a300_").c_str(),cut*weights["MONO"].c_str());
    MZP1200_MA0400->Draw((parameter+">>zp1200a400_").c_str(),cut*weights["MONO"].c_str());
    MZP1200_MA0500->Draw((parameter+">>zp1200a500_").c_str(),cut*weights["MONO"].c_str());
    MZP1200_MA0600->Draw((parameter+">>zp1200a600_").c_str(),cut*weights["MONO"].c_str());
    MZP1200_MA0700->Draw((parameter+">>zp1200a700_").c_str(),cut*weights["MONO"].c_str());
    MZP1200_MA0800->Draw((parameter+">>zp1200a800_").c_str(),cut*weights["MONO"].c_str());

	MZP1400_MA0300->Draw((parameter+">>zp1400a300_").c_str(),cut*weights["MONO"].c_str());
    MZP1400_MA0400->Draw((parameter+">>zp1400a400_").c_str(),cut*weights["MONO"].c_str());
    MZP1400_MA0500->Draw((parameter+">>zp1400a500_").c_str(),cut*weights["MONO"].c_str());
    MZP1400_MA0600->Draw((parameter+">>zp1400a600_").c_str(),cut*weights["MONO"].c_str());
    MZP1400_MA0700->Draw((parameter+">>zp1400a700_").c_str(),cut*weights["MONO"].c_str());
    MZP1400_MA0800->Draw((parameter+">>zp1400a800_").c_str(),cut*weights["MONO"].c_str());
    
	MZP1700_MA0300->Draw((parameter+">>zp1700a300_").c_str(),cut*weights["MONO"].c_str());
    MZP1700_MA0400->Draw((parameter+">>zp1700a400_").c_str(),cut*weights["MONO"].c_str());
    MZP1700_MA0500->Draw((parameter+">>zp1700a500_").c_str(),cut*weights["MONO"].c_str());
    MZP1700_MA0600->Draw((parameter+">>zp1700a600_").c_str(),cut*weights["MONO"].c_str());
    MZP1700_MA0700->Draw((parameter+">>zp1700a700_").c_str(),cut*weights["MONO"].c_str());
    MZP1700_MA0800->Draw((parameter+">>zp1700a800_").c_str(),cut*weights["MONO"].c_str());
    
	MZP2000_MA0300->Draw((parameter+">>zp2000a300_").c_str(),cut*weights["MONO"].c_str());
    MZP2000_MA0400->Draw((parameter+">>zp2000a400_").c_str(),cut*weights["MONO"].c_str());
    MZP2000_MA0500->Draw((parameter+">>zp2000a500_").c_str(),cut*weights["MONO"].c_str());
    MZP2000_MA0600->Draw((parameter+">>zp2000a600_").c_str(),cut*weights["MONO"].c_str());
    MZP2000_MA0700->Draw((parameter+">>zp2000a700_").c_str(),cut*weights["MONO"].c_str());
    MZP2000_MA0800->Draw((parameter+">>zp2000a800_").c_str(),cut*weights["MONO"].c_str());
    
	MZP2500_MA0300->Draw((parameter+">>zp2500a300_").c_str(),cut*weights["MONO"].c_str());
    MZP2500_MA0400->Draw((parameter+">>zp2500a400_").c_str(),cut*weights["MONO"].c_str());
    MZP2500_MA0500->Draw((parameter+">>zp2500a500_").c_str(),cut*weights["MONO"].c_str());
    MZP2500_MA0600->Draw((parameter+">>zp2500a600_").c_str(),cut*weights["MONO"].c_str());
    MZP2500_MA0700->Draw((parameter+">>zp2500a700_").c_str(),cut*weights["MONO"].c_str());
    MZP2500_MA0800->Draw((parameter+">>zp2500a800_").c_str(),cut*weights["MONO"].c_str());
    
    MZP10_MX1->Draw((parameter+">>zp10x1_").c_str(),cut*weights["MONO"].c_str());
    MZP20_MX1->Draw((parameter+">>zp20x1_").c_str(),cut*weights["MONO"].c_str());
    MZP50_MX1->Draw((parameter+">>zp50x1_").c_str(),cut*weights["MONO"].c_str());
    MZP100_MX1->Draw((parameter+">>zp100x1_").c_str(),cut*weights["MONO"].c_str());
    MZP300_MX1->Draw((parameter+">>zp300x1_").c_str(),cut*weights["MONO"].c_str());
    MZP500_MX1->Draw((parameter+">>zp500x1_").c_str(),cut*weights["MONO"].c_str());
    MZP1000_MX1->Draw((parameter+">>zp1000x1_").c_str(),cut*weights["MONO"].c_str());
    MZP2000_MX1->Draw((parameter+">>zp2000x1_").c_str(),cut*weights["MONO"].c_str());
    MZP10000_MX1->Draw((parameter+">>zp10000x1_").c_str(),cut*weights["MONO"].c_str());
    
    MZP10_MX50->Draw((parameter+">>zp10x50_").c_str(),cut*weights["MONO"].c_str());
    MZP50_MX50->Draw((parameter+">>zp50x50_").c_str(),cut*weights["MONO"].c_str());
    MZP95_MX50->Draw((parameter+">>zp95x50_").c_str(),cut*weights["MONO"].c_str());
    MZP200_MX50->Draw((parameter+">>zp200x50_").c_str(),cut*weights["MONO"].c_str());
    MZP300_MX50->Draw((parameter+">>zp300x50_").c_str(),cut*weights["MONO"].c_str());
    MZP10000_MX50->Draw((parameter+">>zp10000x50_").c_str(),cut*weights["MONO"].c_str());

    MZP10_MX150->Draw((parameter+">>zp10x150_").c_str(),cut*weights["MONO"].c_str());
    MZP200_MX150->Draw((parameter+">>zp200x150_").c_str(),cut*weights["MONO"].c_str());
    MZP295_MX150->Draw((parameter+">>zp295x150_").c_str(),cut*weights["MONO"].c_str());
    MZP500_MX150->Draw((parameter+">>zp500x150_").c_str(),cut*weights["MONO"].c_str());
    MZP1000_MX150->Draw((parameter+">>zp1000x150_").c_str(),cut*weights["MONO"].c_str());
    MZP10000_MX150->Draw((parameter+">>zp10000x150_").c_str(),cut*weights["MONO"].c_str());
    
    MZP10_MX500->Draw((parameter+">>zp10x500_").c_str(),cut*weights["MONO"].c_str());
    MZP500_MX500->Draw((parameter+">>zp500x500_").c_str(),cut*weights["MONO"].c_str());
    MZP995_MX500->Draw((parameter+">>zp995x500_").c_str(),cut*weights["MONO"].c_str());
    MZP2000_MX500->Draw((parameter+">>zp2000x500_").c_str(),cut*weights["MONO"].c_str());
    MZP10000_MX500->Draw((parameter+">>zp10000x500_").c_str(),cut*weights["MONO"].c_str());
    
    MZP10_MX1000->Draw((parameter+">>zp10x1000_").c_str(),cut*weights["MONO"].c_str());
    MZP1000_MX1000->Draw((parameter+">>zp1000x1000_").c_str(),cut*weights["MONO"].c_str());
    MZP1995_MX1000->Draw((parameter+">>zp1995x1000_").c_str(),cut*weights["MONO"].c_str());
    MZP10000_MX1000->Draw((parameter+">>zp10000x1000_").c_str(),cut*weights["MONO"].c_str());
    
    fillQCD_Shape(SScut, SScutRel, sfQCD, qcd_, parameter, bin, "QCD SS " + parameter + " sideband in " + chan + " channel", syst, chan, plot, createOutputShapes);
	qcd_->Scale(qcdOStoSS);
    
    std::cout << "Yields before rel scaling" << std::endl;
    std::cout << "data: " << data_->Integral() << std::endl;
    std::cout << "W: " << w_->Integral() << std::endl;
    std::cout << "ZTT: " << ztt_->Integral() << std::endl;
    std::cout << "ZL: " << zl_->Integral() << std::endl;
    std::cout << "ZJ: " << zj_->Integral() << std::endl;
    std::cout << "VVT: " << vvt_->Integral() << std::endl;
    std::cout << "VVJ: " << vvj_->Integral() << std::endl;
    std::cout << "ZVV: " << zvv_->Integral() << std::endl;
    std::cout << "EWK: " << ewk_->Integral() << std::endl;
    std::cout << "TTT: " << ttt_->Integral() << std::endl;
    std::cout << "TTJ: " << ttj_->Integral() << std::endl;
    std::cout << "ZH: " << zhtautau_->Integral() << std::endl;
    std::cout << "VBF H: " << glugluhtautau_->Integral() << std::endl;
    std::cout << "GG H: " << vbfhtautau_->Integral() << std::endl;
    
    // not scaling by sf in datacards
	//w_->Scale(sf);
    std::cout << "scaling W" << std::endl;
    w_->Scale(shapeScale(cut, relCut, W, weights["W"], parameter, bin));
    std::cout << "scaling VVT" << std::endl;
    vvt_->Scale(shapeScale(cut, relCut, VVT, weights["VV"], parameter, bin));
    std::cout << "scaling VVJ" << std::endl;
    vvj_->Scale(shapeScale(cut, relCut, VVJ, weights["VV"], parameter, bin));
    std::cout << "scaling TTT" << std::endl;
    ttt_->Scale(shapeScale(cut, relCut, TTT, weights["TTT"], parameter, bin));
    std::cout << "scaling TTJ" << std::endl;
    ttj_->Scale(shapeScale(cut, relCut, TTJ, weights["TTJ"], parameter, bin));
    std::cout << "scaling ZTT" << std::endl;
    ztt_->Scale(shapeScale(cut, relCut, ZTT, weights["ZTT"], parameter, bin));
    std::cout << "scaling ZL" << std::endl;
    zl_->Scale(shapeScale(cut, relCut, ZL, weights["ZL"], parameter, bin));
    std::cout << "scaling ZJ" << std::endl;
    zj_->Scale(shapeScale(cut, relCut, ZJ, weights["ZJ"], parameter, bin));


    
    std::cout << "Yields after rel scaling" << std::endl;
    std::cout << "data: " << data_->Integral() << std::endl;
    std::cout << "W: " << w_->Integral() << std::endl;
    std::cout << "ZTT: " << ztt_->Integral() << std::endl;
    std::cout << "ZL: " << zl_->Integral() << std::endl;
    std::cout << "ZJ: " << zj_->Integral() << std::endl;
    std::cout << "VVT: " << vvt_->Integral() << std::endl;
    std::cout << "VVJ: " << vvj_->Integral() << std::endl;
    std::cout << "ZVV: " << zvv_->Integral() << std::endl;
    std::cout << "EWK: " << ewk_->Integral() << std::endl;
    std::cout << "TTT: " << ttt_->Integral() << std::endl;
    std::cout << "TTJ: " << ttj_->Integral() << std::endl;
    std::cout << "ZH: " << zhtautau_->Integral() << std::endl;
    std::cout << "VBF H: " << glugluhtautau_->Integral() << std::endl;
    std::cout << "GG H: " << vbfhtautau_->Integral() << std::endl;
    
	// scale things to the projected lumi
	vvt_->Scale(lumi_sf);
    vvj_->Scale(lumi_sf);
    zvv_->Scale(lumi_sf);
    ewk_->Scale(lumi_sf);
	ttt_->Scale(lumi_sf);
    ttj_->Scale(lumi_sf);
	ztt_->Scale(lumi_sf);
    zl_->Scale(lumi_sf);
    zj_->Scale(lumi_sf);
	w_->Scale(lumi_sf);
	qcd_->Scale(lumi_sf);
    zhtautau_->Scale(lumi_sf);
    glugluhtautau_->Scale(lumi_sf);
    vbfhtautau_->Scale(lumi_sf);
    smh_->Scale(lumi_sf);
	zp600a300_->Scale(lumi_sf);
    zp600a400_->Scale(lumi_sf);
	zp800a300_->Scale(lumi_sf);
    zp800a400_->Scale(lumi_sf);
    zp800a500_->Scale(lumi_sf);
    zp800a600_->Scale(lumi_sf);
	zp1000a300_->Scale(lumi_sf);
    zp1000a400_->Scale(lumi_sf);
    zp1000a500_->Scale(lumi_sf);
    zp1000a600_->Scale(lumi_sf);
    zp1000a700_->Scale(lumi_sf);
    zp1000a800_->Scale(lumi_sf);
	zp1200a300_->Scale(lumi_sf);
    zp1200a400_->Scale(lumi_sf);
    zp1200a500_->Scale(lumi_sf);
    zp1200a600_->Scale(lumi_sf);
    zp1200a700_->Scale(lumi_sf);
    zp1200a800_->Scale(lumi_sf);
    zp1400a300_->Scale(lumi_sf);
    zp1400a400_->Scale(lumi_sf);
    zp1400a500_->Scale(lumi_sf);
    zp1400a600_->Scale(lumi_sf);
    zp1400a700_->Scale(lumi_sf);
    zp1400a800_->Scale(lumi_sf);
	zp1700a300_->Scale(lumi_sf);
    zp1700a400_->Scale(lumi_sf);
    zp1700a500_->Scale(lumi_sf);
    zp1700a600_->Scale(lumi_sf);
    zp1700a700_->Scale(lumi_sf);
    zp1700a800_->Scale(lumi_sf);
	zp2000a300_->Scale(lumi_sf);
    zp2000a400_->Scale(lumi_sf);
    zp2000a500_->Scale(lumi_sf);
    zp2000a600_->Scale(lumi_sf);
    zp2000a700_->Scale(lumi_sf);
    zp2000a800_->Scale(lumi_sf);
	zp2500a300_->Scale(lumi_sf);
    zp2500a400_->Scale(lumi_sf);
    zp2500a500_->Scale(lumi_sf);
    zp2500a600_->Scale(lumi_sf);
    zp2500a700_->Scale(lumi_sf);
    zp2500a800_->Scale(lumi_sf);
    
    zp10x1_->Scale(lumi_sf);
    zp20x1_->Scale(lumi_sf);
    zp50x1_->Scale(lumi_sf);
    zp100x1_->Scale(lumi_sf);
    zp300x1_->Scale(lumi_sf);
    zp500x1_->Scale(lumi_sf);
    zp1000x1_->Scale(lumi_sf);
    zp2000x1_->Scale(lumi_sf);
    zp10000x1_->Scale(lumi_sf);
    
    zp10x50_->Scale(lumi_sf);
    zp50x50_->Scale(lumi_sf);
    zp95x50_->Scale(lumi_sf);
    zp200x50_->Scale(lumi_sf);
    zp300x50_->Scale(lumi_sf);
    zp10000x50_->Scale(lumi_sf);

    zp10x150_->Scale(lumi_sf);
    zp200x150_->Scale(lumi_sf);
    zp295x150_->Scale(lumi_sf);
    zp500x150_->Scale(lumi_sf);
    zp1000x150_->Scale(lumi_sf);
    zp10000x150_->Scale(lumi_sf);
    
    zp10x500_->Scale(lumi_sf);
    zp500x500_->Scale(lumi_sf);
    zp995x500_->Scale(lumi_sf);
    zp2000x500_->Scale(lumi_sf);
    zp10000x500_->Scale(lumi_sf);
    
    zp10x1000_->Scale(lumi_sf);
    zp1000x1000_->Scale(lumi_sf);
    zp1995x1000_->Scale(lumi_sf);
    zp10000x1000_->Scale(lumi_sf);

    

    /* get a sum of bks */
	comb_->Add(vvt_);
    comb_->Add(vvj_);
    comb_->Add(zvv_);
    comb_->Add(ewk_);
	comb_->Add(ttt_);
    comb_->Add(ttj_);
	comb_->Add(w_);
	comb_->Add(ztt_);
    comb_->Add(zj_);
    comb_->Add(zl_);
    comb_->Add(smh_);
	comb_->Add(qcd_);
    
    std::cout << "Editing Shape File" << std::endl;
    
	if(createOutputShapes)
	{
		/* create a file to help with S/sqrt(S+B) cut opt */
        
        std::string paramName;
        if (parameter=="met" || parameter=="mvamet") {paramName = "met";}
        else if (parameter.substr(0,6)=="mvaVar") {paramName = "mva";}
        else if (parameter.substr(0,6) == "mt_tot") {paramName="mt";}
        else {paramName = parameter;}
        
        std::string catName = chan + "_inclusive";
		std::string optName = "xtt_" + chan + ".inputs-13TeV-" + paramName + ".root";
		TFile sbOp(optName.c_str(),"UPDATE");
        
        TDirectory * catdir = sbOp.GetDirectory(catName.c_str());
        if (catdir){catdir->cd();}
        else {sbOp.mkdir(catName.c_str()); TDirectory * catdir = sbOp.GetDirectory(catName.c_str()); catdir->cd();}
        
        std::string DATAname = "data_obs" + syst;
        const char *DATAsave = DATAname.c_str();
        std::string VVTname = "VVT" + syst;
        const char *VVTsave = VVTname.c_str();
        std::string VVJname = "VVJ" + syst;
        const char *VVJsave = VVJname.c_str();
        std::string ZVVname = "ZVV" + syst;
        const char *ZVVsave = ZVVname.c_str();
        std::string EWKname = "EWK" + syst;
        const char *EWKsave = EWKname.c_str();
        std::string TTTname = "TTT" + syst;
        const char *TTTsave = TTTname.c_str();
        std::string TTJname = "TTJ" + syst;
        const char *TTJsave = TTJname.c_str();
        std::string ZTTname = "ZTT" + syst;
        const char *ZTTsave = ZTTname.c_str();
        std::string ZLname = "ZL" + syst;
        const char *ZLsave = ZLname.c_str();
        std::string ZJname = "ZJ" + syst;
        const char *ZJsave = ZJname.c_str();
        std::string Wname = "W" + syst;
        const char *Wsave = Wname.c_str();
        std::string QCDname = "QCD" + syst;
        const char *QCDsave = QCDname.c_str();
        std::string SMHname = "SMH" + syst;
        const char *SMHsave = SMHname.c_str();
        std::string ZHname = "ZH" + syst;
        const char *ZHsave = ZHname.c_str();
        
        std::string combname = "comb" + syst;
        const char *combsave = combname.c_str();
        
        std::string ZprimeA400Zname = "ZprimeA400Z" + syst;
        const char *ZprimeA400Zsave = ZprimeA400Zname.c_str();
        
        std::string Zprime600A300name = "Zprime600A300" + syst;
        const char *Zprime600A300save = Zprime600A300name.c_str();
        std::string Zprime600A400name = "Zprime600A400" + syst;
        const char *Zprime600A400save = Zprime600A400name.c_str();
        
        std::string Zprime800A300name = "Zprime800A300" + syst;
        const char *Zprime800A300save = Zprime800A300name.c_str();
        std::string Zprime800A400name = "Zprime800A400" + syst;
        const char *Zprime800A400save = Zprime800A400name.c_str();
        std::string Zprime800A500name = "Zprime800A500" + syst;
        const char *Zprime800A500save = Zprime800A500name.c_str();
        std::string Zprime800A600name = "Zprime800A600" + syst;
        const char *Zprime800A600save = Zprime800A600name.c_str();
        
        std::string Zprime1000A300name = "Zprime1000A300" + syst;
        const char *Zprime1000A300save = Zprime1000A300name.c_str();
        std::string Zprime1000A400name = "Zprime1000A400" + syst;
        const char *Zprime1000A400save = Zprime1000A400name.c_str();
        std::string Zprime1000A500name = "Zprime1000A500" + syst;
        const char *Zprime1000A500save = Zprime1000A500name.c_str();
        std::string Zprime1000A600name = "Zprime1000A600" + syst;
        const char *Zprime1000A600save = Zprime1000A600name.c_str();
        std::string Zprime1000A700name = "Zprime1000A700" + syst;
        const char *Zprime1000A700save = Zprime1000A700name.c_str();
        std::string Zprime1000A800name = "Zprime1000A800" + syst;
        const char *Zprime1000A800save = Zprime1000A800name.c_str();
        
        std::string Zprime1200A300name = "Zprime1200A300" + syst;
        const char *Zprime1200A300save = Zprime1200A300name.c_str();
        std::string Zprime1200A400name = "Zprime1200A400" + syst;
        const char *Zprime1200A400save = Zprime1200A400name.c_str();
        std::string Zprime1200A500name = "Zprime1200A500" + syst;
        const char *Zprime1200A500save = Zprime1200A500name.c_str();
        std::string Zprime1200A600name = "Zprime1200A600" + syst;
        const char *Zprime1200A600save = Zprime1200A600name.c_str();
        std::string Zprime1200A700name = "Zprime1200A700" + syst;
        const char *Zprime1200A700save = Zprime1200A700name.c_str();
        std::string Zprime1200A800name = "Zprime1200A800" + syst;
        const char *Zprime1200A800save = Zprime1200A800name.c_str();
        
        std::string Zprime1400A300name = "Zprime1400A300" + syst;
        const char *Zprime1400A300save = Zprime1400A300name.c_str();
        std::string Zprime1400A400name = "Zprime1400A400" + syst;
        const char *Zprime1400A400save = Zprime1400A400name.c_str();
        std::string Zprime1400A500name = "Zprime1400A500" + syst;
        const char *Zprime1400A500save = Zprime1400A500name.c_str();
        std::string Zprime1400A600name = "Zprime1400A600" + syst;
        const char *Zprime1400A600save = Zprime1400A600name.c_str();
        std::string Zprime1400A700name = "Zprime1400A700" + syst;
        const char *Zprime1400A700save = Zprime1400A700name.c_str();
        std::string Zprime1400A800name = "Zprime1400A800" + syst;
        const char *Zprime1400A800save = Zprime1400A800name.c_str();
        
        std::string Zprime1700A300name = "Zprime1700A300" + syst;
        const char *Zprime1700A300save = Zprime1700A300name.c_str();
        std::string Zprime1700A400name = "Zprime1700A400" + syst;
        const char *Zprime1700A400save = Zprime1700A400name.c_str();
        std::string Zprime1700A500name = "Zprime1700A500" + syst;
        const char *Zprime1700A500save = Zprime1700A500name.c_str();
        std::string Zprime1700A600name = "Zprime1700A600" + syst;
        const char *Zprime1700A600save = Zprime1700A600name.c_str();
        std::string Zprime1700A700name = "Zprime1700A700" + syst;
        const char *Zprime1700A700save = Zprime1700A700name.c_str();
        std::string Zprime1700A800name = "Zprime1700A800" + syst;
        const char *Zprime1700A800save = Zprime1700A800name.c_str();
        
        std::string Zprime2000A300name = "Zprime2000A300" + syst;
        const char *Zprime2000A300save = Zprime2000A300name.c_str();
        std::string Zprime2000A400name = "Zprime2000A400" + syst;
        const char *Zprime2000A400save = Zprime2000A400name.c_str();
        std::string Zprime2000A500name = "Zprime2000A500" + syst;
        const char *Zprime2000A500save = Zprime2000A500name.c_str();
        std::string Zprime2000A600name = "Zprime2000A600" + syst;
        const char *Zprime2000A600save = Zprime2000A600name.c_str();
        std::string Zprime2000A700name = "Zprime2000A700" + syst;
        const char *Zprime2000A700save = Zprime2000A700name.c_str();
        std::string Zprime2000A800name = "Zprime2000A800" + syst;
        const char *Zprime2000A800save = Zprime2000A800name.c_str();
        
        std::string Zprime2500A300name = "Zprime2500A300" + syst;
        const char *Zprime2500A300save = Zprime2500A300name.c_str();
        std::string Zprime2500A400name = "Zprime2500A400" + syst;
        const char *Zprime2500A400save = Zprime2500A400name.c_str();
        std::string Zprime2500A500name = "Zprime2500A500" + syst;
        const char *Zprime2500A500save = Zprime2500A500name.c_str();
        std::string Zprime2500A600name = "Zprime2500A600" + syst;
        const char *Zprime2500A600save = Zprime2500A600name.c_str();
        std::string Zprime2500A700name = "Zprime2500A700" + syst;
        const char *Zprime2500A700save = Zprime2500A700name.c_str();
        std::string Zprime2500A800name = "Zprime2500A800" + syst;
        const char *Zprime2500A800save = Zprime2500A800name.c_str();
        
        std::string Zprime10X1name = "Zprime10X1" + syst;
        const char *Zprime10X1save = Zprime10X1name.c_str();
        std::string Zprime20X1name = "Zprime20X1" + syst;
        const char *Zprime20X1save = Zprime20X1name.c_str();
        std::string Zprime50X1name = "Zprime50X1" + syst;
        const char *Zprime50X1save = Zprime50X1name.c_str();
        std::string Zprime100X1name = "Zprime100X1" + syst;
        const char *Zprime100X1save = Zprime100X1name.c_str();
        std::string Zprime300X1name = "Zprime300X1" + syst;
        const char *Zprime300X1save = Zprime300X1name.c_str();
        std::string Zprime500X1name = "Zprime500X1" + syst;
        const char *Zprime500X1save = Zprime500X1name.c_str();
        std::string Zprime1000X1name = "Zprime1000X1" + syst;
        const char *Zprime1000X1save = Zprime1000X1name.c_str();
        std::string Zprime2000X1name = "Zprime2000X1" + syst;
        const char *Zprime2000X1save = Zprime2000X1name.c_str();
        std::string Zprime10000X1name = "Zprime10000X1" + syst;
        const char *Zprime10000X1save = Zprime10000X1name.c_str();
        
        std::string Zprime10X50name = "Zprime10X50" + syst;
        const char *Zprime10X50save = Zprime10X50name.c_str();
        std::string Zprime50X50name = "Zprime50X50" + syst;
        const char *Zprime50X50save = Zprime50X50name.c_str();
        std::string Zprime95X50name = "Zprime95X50" + syst;
        const char *Zprime95X50save = Zprime95X50name.c_str();
        std::string Zprime200X50name = "Zprime200X50" + syst;
        const char *Zprime200X50save = Zprime200X50name.c_str();
        std::string Zprime300X50name = "Zprime300X50" + syst;
        const char *Zprime300X50save = Zprime300X50name.c_str();
        std::string Zprime10000X50name = "Zprime10000X50" + syst;
        const char *Zprime10000X50save = Zprime10000X50name.c_str();
    
        std::string Zprime10X150name = "Zprime10X150" + syst;
        const char *Zprime10X150save = Zprime10X150name.c_str();
        std::string Zprime200X150name = "Zprime200X150" + syst;
        const char *Zprime200X150save = Zprime200X150name.c_str();
        std::string Zprime295X150name = "Zprime295X150" + syst;
        const char *Zprime295X150save = Zprime295X150name.c_str();
        std::string Zprime500X150name = "Zprime500X150" + syst;
        const char *Zprime500X150save = Zprime500X150name.c_str();
        std::string Zprime1000X150name = "Zprime1000X150" + syst;
        const char *Zprime1000X150save = Zprime1000X150name.c_str();
        std::string Zprime10000X150name = "Zprime10000X150" + syst;
        const char *Zprime10000X150save = Zprime10000X150name.c_str();
        
        std::string Zprime10X500name = "Zprime10X500" + syst;
        const char *Zprime10X500save = Zprime10X500name.c_str();
        std::string Zprime500X500name = "Zprime500X500" + syst;
        const char *Zprime500X500save = Zprime500X500name.c_str();
        std::string Zprime995X500name = "Zprime995X500" + syst;
        const char *Zprime995X500save = Zprime995X500name.c_str();
        std::string Zprime2000X500name = "Zprime2000X500" + syst;
        const char *Zprime2000X500save = Zprime2000X500name.c_str();
        std::string Zprime10000X500name = "Zprime10000X500" + syst;
        const char *Zprime10000X500save = Zprime10000X500name.c_str();
        
        std::string Zprime10X1000name = "Zprime10X1000" + syst;
        const char *Zprime10X1000save = Zprime10X1000name.c_str();
        std::string Zprime1000X1000name = "Zprime1000X1000" + syst;
        const char *Zprime1000X1000save = Zprime1000X1000name.c_str();
        std::string Zprime1995X1000name = "Zprime1995X1000" + syst;
        const char *Zprime1995X1000save = Zprime1995X1000name.c_str();
        std::string Zprime10000X1000name = "Zprime10000X1000" + syst;
        const char *Zprime10000X1000save = Zprime10000X1000name.c_str();

        data_->Write(DATAsave);
        
        vvt_->Write(VVTsave);
        vvj_->Write(VVJsave);
        zvv_->Write(ZVVsave);
        ewk_->Write(EWKsave);
        ttt_->Write(TTTsave);
        ttj_->Write(TTJsave);
        ztt_->Write(ZTTsave);
        zl_->Write(ZLsave);
        zj_->Write(ZJsave);
        w_->Write(Wsave);
        qcd_->Write(QCDsave);
        
        zhtautau_->Write(ZHsave);
        smh_->Write(SMHsave);
        
		comb_->Write(combsave);

        if (syst.substr(0,3) == "600"){zp600a400_->Write(ZprimeA400Zsave);}
        else if (syst.substr(0,3) == "800"){zp800a400_->Write(ZprimeA400Zsave);}
        else if (syst.substr(0,4) == "1000"){zp1000a400_->Write(ZprimeA400Zsave);}
        else if (syst.substr(0,4) == "1200"){zp1200a400_->Write(ZprimeA400Zsave);}
        else
        {
            zp600a300_->Write(Zprime600A300save);
            zp600a400_->Write(Zprime600A400save);
            
            zp800a300_->Write(Zprime800A300save);
            zp800a400_->Write(Zprime800A400save);
            zp800a500_->Write(Zprime800A500save);
            zp800a600_->Write(Zprime800A600save);
            
            zp1000a300_->Write(Zprime1000A300save);
            zp1000a400_->Write(Zprime1000A400save);
            zp1000a500_->Write(Zprime1000A500save);
            zp1000a600_->Write(Zprime1000A600save);
            zp1000a700_->Write(Zprime1000A700save);
            zp1000a800_->Write(Zprime1000A800save);
            
            zp1200a300_->Write(Zprime1200A300save);
            zp1200a400_->Write(Zprime1200A400save);
            zp1200a500_->Write(Zprime1200A500save);
            zp1200a600_->Write(Zprime1200A600save);
            zp1200a700_->Write(Zprime1200A700save);
            zp1200a800_->Write(Zprime1200A800save);
            
            zp1400a300_->Write(Zprime1400A300save);
            zp1400a400_->Write(Zprime1400A400save);
            zp1400a500_->Write(Zprime1400A500save);
            zp1400a600_->Write(Zprime1400A600save);
            zp1400a700_->Write(Zprime1400A700save);
            zp1400a800_->Write(Zprime1400A800save);
            
            zp1700a300_->Write(Zprime1700A300save);
            zp1700a400_->Write(Zprime1700A400save);
            zp1700a500_->Write(Zprime1700A500save);
            zp1700a600_->Write(Zprime1700A600save);
            zp1700a700_->Write(Zprime1700A700save);
            zp1700a800_->Write(Zprime1700A800save);
            
            zp2000a300_->Write(Zprime2000A300save);
            zp2000a400_->Write(Zprime2000A400save);
            zp2000a500_->Write(Zprime2000A500save);
            zp2000a600_->Write(Zprime2000A600save);
            zp2000a700_->Write(Zprime2000A700save);
            zp2000a800_->Write(Zprime2000A800save);
            
            zp2500a300_->Write(Zprime2500A300save);
            zp2500a400_->Write(Zprime2500A400save);
            zp2500a500_->Write(Zprime2500A500save);
            zp2500a600_->Write(Zprime2500A600save);
            zp2500a700_->Write(Zprime2500A700save);
            zp2500a800_->Write(Zprime2500A800save);
            
            zp10x1_->Write(Zprime10X1save);
            zp20x1_->Write(Zprime20X1save);
            zp50x1_->Write(Zprime50X1save);
            zp100x1_->Write(Zprime100X1save);
            zp300x1_->Write(Zprime300X1save);
            zp500x1_->Write(Zprime500X1save);
            zp1000x1_->Write(Zprime1000X1save);
            zp2000x1_->Write(Zprime2000X1save);
            zp10000x1_->Write(Zprime10000X1save);
            
            zp10x50_->Write(Zprime10X50save);
            zp50x50_->Write(Zprime50X50save);
            zp95x50_->Write(Zprime95X50save);
            zp200x50_->Write(Zprime200X50save);
            zp300x50_->Write(Zprime300X50save);
            zp10000x50_->Write(Zprime10000X50save);
        
            zp10x150_->Write(Zprime10X150save);
            zp200x150_->Write(Zprime200X150save);
            zp295x150_->Write(Zprime295X150save);
            zp500x150_->Write(Zprime500X150save);
            zp1000x150_->Write(Zprime1000X150save);
            zp10000x150_->Write(Zprime10000X150save);
            
            zp10x500_->Write(Zprime10X500save);
            zp500x500_->Write(Zprime500X500save);
            zp995x500_->Write(Zprime995X500save);
            zp2000x500_->Write(Zprime2000X500save);
            zp10000x500_->Write(Zprime10000X500save);
            
            zp10x1000_->Write(Zprime10X1000save);
            zp1000x1000_->Write(Zprime1000X1000save);
            zp1995x1000_->Write(Zprime1995X1000save);
            zp10000x1000_->Write(Zprime10000X1000save);

        }

        sbOp.cd();
		sbOp.Close();
	}
    std::cout << "FINISHED editing shape file" << std::endl;

    
    THStack * sigMTstack  = new THStack();
    
    if (plot)
    {
        if(parameter=="met")
        {
            for(int i=0; i<data_->GetNbinsX()+1; ++i)
            {
                if(data_->GetBinLowEdge(i)>= met_blind) 
                {
                    data_->SetBinContent(i,0);
                    data_->SetBinError(i,0);
                }
            }
        }
        else if(parameter.substr(0,6)=="mvaVar")
        {
            for(int i=0; i<data_->GetNbinsX()+1; ++i)
            {
                if(data_->GetBinLowEdge(i)>= 0.5)
                {
                    data_->SetBinContent(i,0);
                    data_->SetBinError(i,0);
                }
            }
        }
        else if(parameter=="mt_tot")
        {
            for(int i=0; i<data_->GetNbinsX()+1; ++i)
            {
                if(data_->GetBinLowEdge(i)>= 200)
                {
                    data_->SetBinContent(i,0);
                    data_->SetBinError(i,0);
                }
            }
        }
        
        sigMTstack->Add(ttt_);
        sigMTstack->Add(ttj_);
        sigMTstack->Add(zvv_);
        sigMTstack->Add(ewk_);
        sigMTstack->Add(qcd_);
        if(splitPlotSMH)
        {
            sigMTstack->Add(zhtautau_);
            sigMTstack->Add(vbfhtautau_);
            sigMTstack->Add(glugluhtautau_);
        }
        else
        {
            sigMTstack->Add(smh_);
        }
        sigMTstack->Add(vvt_);
        sigMTstack->Add(vvj_);
        sigMTstack->Add(w_);
        sigMTstack->Add(ztt_);
        sigMTstack->Add(zj_);
        sigMTstack->Add(zl_);

        Signal_region->cd(1);

        data_->SetMaximum(max(DATA_MAX*data_->GetMaximum(),30.0));
        data_->SetMinimum(DATA_MIN);

    //	data_->SetTitle((can_name+";XXX;Events/Bin").c_str());
        data_->SetTitle(";XXX;Events/Bin");
        data_->GetXaxis()->SetTitle(x_axisLabels[parameter].c_str());

        data_->DrawCopy("PE");
        sigMTstack->Draw("histsames");

        if(parameter.substr(0,6)=="mvaVar")
        {
            if (syst=="600"){zp600a400_->SetLineStyle(3); zp600a400_->SetLineWidth(3); zp600a400_->SetLineColor(50); std::cout << "DRAWING" << std::endl; zp600a400_->DrawCopy("hist sames");}
            else if (syst=="800"){zp800a400_->SetLineStyle(3); zp800a400_->SetLineWidth(3); zp800a400_->SetLineColor(60); std::cout << "DRAWING" << std::endl; zp800a400_->DrawCopy("hist sames");}
            else if (syst=="1000"){zp1000a400_->SetLineStyle(3); zp1000a400_->SetLineWidth(3); zp1000a400_->SetLineColor(70); std::cout << "DRAWING" << std::endl; zp1000a400_->DrawCopy("hist sames");}
            else if (syst=="1200"){zp1200a400_->SetLineStyle(3); zp1200a400_->SetLineWidth(3); zp1200a400_->SetLineColor(80); std::cout << "DRAWING" << std::endl; zp1200a400_->DrawCopy("hist sames");}
            else
            {
                zp600a400_->SetLineStyle(2); zp600a400_->SetLineWidth(3); zp600a400_->SetLineColor(50); std::cout << "DRAWING" << std::endl; zp600a400_->DrawCopy("hist sames");
                zp800a400_->SetLineStyle(2); zp800a400_->SetLineWidth(3); zp800a400_->SetLineColor(60); std::cout << "DRAWING" << std::endl; zp800a400_->DrawCopy("hist sames");
                zp1000a400_->SetLineStyle(2); zp1000a400_->SetLineWidth(3); zp1000a400_->SetLineColor(70); std::cout << "DRAWING" << std::endl; zp1000a400_->DrawCopy("hist sames");
                zp1200a400_->SetLineStyle(2); zp1200a400_->SetLineWidth(3); zp1200a400_->SetLineColor(80); std::cout << "DRAWING" << std::endl; zp1200a400_->DrawCopy("hist sames");
            }
        }
        else if (parameter == "mt_tot" || parameter == "met")
        {
            zp600a300_->SetLineStyle(2); zp600a300_->SetLineWidth(3); zp600a300_->SetLineColor(50); std::cout << "DRAWING" << std::endl; zp600a300_->DrawCopy("hist sames");
            zp800a300_->SetLineStyle(2); zp800a300_->SetLineWidth(3); zp800a300_->SetLineColor(60); std::cout << "DRAWING" << std::endl; zp800a300_->DrawCopy("hist sames");
            zp1000a300_->SetLineStyle(2); zp1000a300_->SetLineWidth(3); zp1000a300_->SetLineColor(70); std::cout << "DRAWING" << std::endl; zp1000a300_->DrawCopy("hist sames");
            zp1200a300_->SetLineStyle(2); zp1200a300_->SetLineWidth(3); zp1200a300_->SetLineColor(80); std::cout << "DRAWING" << std::endl; zp1200a300_->DrawCopy("hist sames");
            zp1400a300_->SetLineStyle(2); zp1400a300_->SetLineWidth(3); zp1400a300_->SetLineColor(50); std::cout << "DRAWING" << std::endl; zp1400a300_->DrawCopy("hist sames");
            zp1700a300_->SetLineStyle(2); zp1700a300_->SetLineWidth(3); zp1700a300_->SetLineColor(60); std::cout << "DRAWING" << std::endl; zp1700a300_->DrawCopy("hist sames");
            zp2000a300_->SetLineStyle(2); zp2000a300_->SetLineWidth(3); zp2000a300_->SetLineColor(70); std::cout << "DRAWING" << std::endl; zp2000a300_->DrawCopy("hist sames");
            zp2500a300_->SetLineStyle(2); zp2500a300_->SetLineWidth(3); zp2500a300_->SetLineColor(80); std::cout << "DRAWING" << std::endl; zp2500a300_->DrawCopy("hist sames");
        }
        data_->DrawCopy("PEsames");

        /* update canvas */
        Signal_region->Update();
        
        if(data_->GetSumOfWeights()==0)
        {
            //data_->SetMarkerColor(2);
            //std::cout<<" **** EMPTY BLINDED DATA ****, setting eq to bkg total ****\n";
            data_->SetMaximum(100.0);
            //data_->Add(comb_);
            //data_->DrawCopy("PE");
            sigMTstack->Draw("histsames");
            if(parameter=="met" || parameter=="mt_tot")
            {
            
                    zp600a300_->SetLineStyle(2); zp600a300_->SetLineWidth(3); zp600a300_->SetLineColor(50); std::cout << "DRAWING" << std::endl; zp600a300_->DrawCopy("hist sames");
                    zp800a300_->SetLineStyle(2); zp800a300_->SetLineWidth(3); zp800a300_->SetLineColor(60); std::cout << "DRAWING" << std::endl; zp800a300_->DrawCopy("hist sames");
                    zp1000a300_->SetLineStyle(2); zp1000a300_->SetLineWidth(3); zp1000a300_->SetLineColor(70); std::cout << "DRAWING" << std::endl; zp1000a300_->DrawCopy("hist sames");
                    zp1200a300_->SetLineStyle(2); zp1200a300_->SetLineWidth(3); zp1200a300_->SetLineColor(80); std::cout << "DRAWING" << std::endl; zp1200a300_->DrawCopy("hist sames");
                    zp1400a300_->SetLineStyle(2); zp1400a300_->SetLineWidth(3); zp1400a300_->SetLineColor(50); std::cout << "DRAWING" << std::endl; zp1400a300_->DrawCopy("hist sames");
                    zp1700a300_->SetLineStyle(2); zp1700a300_->SetLineWidth(3); zp1700a300_->SetLineColor(60); std::cout << "DRAWING" << std::endl; zp1700a300_->DrawCopy("hist sames");
                    zp2000a300_->SetLineStyle(2); zp2000a300_->SetLineWidth(3); zp2000a300_->SetLineColor(70); std::cout << "DRAWING" << std::endl; zp2000a300_->DrawCopy("hist sames");
                    zp2500a300_->SetLineStyle(2); zp2500a300_->SetLineWidth(3); zp2500a300_->SetLineColor(80); std::cout << "DRAWING" << std::endl; zp2500a300_->DrawCopy("hist sames");
            }
            else
            {
                if (syst=="600"){zp600a400_->SetLineStyle(3); zp600a400_->SetLineWidth(3); zp600a400_->SetLineColor(50); std::cout << "DRAWING" << std::endl; zp600a400_->DrawCopy("hist sames");}
                else if (syst=="800"){zp800a400_->SetLineStyle(3); zp800a400_->SetLineWidth(5); zp800a400_->SetLineColor(60); std::cout << "DRAWING" << std::endl; zp800a400_->DrawCopy("hist sames");}
                else if (syst=="1000"){zp1000a400_->SetLineStyle(3); zp1000a400_->SetLineWidth(5); zp1000a400_->SetLineColor(70); std::cout << "DRAWING" << std::endl; zp1000a400_->DrawCopy("hist sames");}
                else if (syst=="1200"){zp1200a400_->SetLineStyle(3); zp1200a400_->SetLineWidth(5); zp1200a400_->SetLineColor(80); std::cout << "DRAWING" << std::endl; zp1200a400_->DrawCopy("hist sames");}
                else if (parameter.substr(0,6)=="mvaVar")
                {
                    zp600a400_->SetLineStyle(3); zp600a400_->SetLineWidth(3); zp600a400_->SetLineColor(50); std::cout << "DRAWING" << std::endl; zp600a400_->DrawCopy("hist sames");
                    zp800a400_->SetLineStyle(3); zp800a400_->SetLineWidth(3); zp800a400_->SetLineColor(60); std::cout << "DRAWING" << std::endl; zp800a400_->DrawCopy("hist sames");
                    zp1000a400_->SetLineStyle(3); zp1000a400_->SetLineWidth(3); zp1000a400_->SetLineColor(70); std::cout << "DRAWING" << std::endl; zp1000a400_->DrawCopy("hist sames");
                    zp1200a400_->SetLineStyle(3); zp1200a400_->SetLineWidth(3); zp1200a400_->SetLineColor(80); std::cout << "DRAWING" << std::endl; zp1200a400_->DrawCopy("hist sames");
                }
            }
            data_->DrawCopy("PEsames");
        }

        drawTcanvasAndLegendAndRatio(Signal_region, data_, comb_, 1);

        Signal_region->DrawClone();

        std::string localSaveName = saveName + parameter + "_" + chan + "_" + drCutNameMap[drCut] + "Region.png";
        Signal_region->cd()->SaveAs(localSaveName.c_str());

        Signal_region->Close();
    }
    
	/* delete pointers */
	delete data_;
	delete vvt_;
    delete vvj_;
    delete zvv_;
    delete ewk_;
	delete ttt_;
    delete ttj_;
	delete ztt_;
    delete zl_;
    delete zj_;
	delete w_;
    delete zhtautau_;
    delete glugluhtautau_;
    delete vbfhtautau_;
	delete comb_;
    delete smh_;
	delete qcd_;
	delete zp600a300_;
    delete zp600a400_;
	delete zp800a300_;
    delete zp800a400_;
    delete zp800a500_;
    delete zp800a600_;
	delete zp1000a300_;
    delete zp1000a400_;
    delete zp1000a500_;
    delete zp1000a600_;
    delete zp1000a700_;
    delete zp1000a800_;
	delete zp1200a300_;
    delete zp1200a400_;
    delete zp1200a500_;
    delete zp1200a600_;
    delete zp1200a700_;
    delete zp1200a800_;
    delete zp1400a300_;
    delete zp1400a400_;
    delete zp1400a500_;
    delete zp1400a600_;
    delete zp1400a700_;
    delete zp1400a800_;
    delete zp1700a300_;
    delete zp1700a400_;
    delete zp1700a500_;
    delete zp1700a600_;
    delete zp1700a700_;
    delete zp1700a800_;
    delete zp2000a300_;
    delete zp2000a400_;
    delete zp2000a500_;
    delete zp2000a600_;
    delete zp2000a700_;
    delete zp2000a800_;
    delete zp2500a300_;
    delete zp2500a400_;
    delete zp2500a500_;
    delete zp2500a600_;
    delete zp2500a700_;
    delete zp2500a800_;
    
    delete zp10x1_;
    delete zp20x1_;
    delete zp50x1_;
    delete zp100x1_;
    delete zp300x1_;
    delete zp500x1_;
    delete zp1000x1_;
    delete zp2000x1_;
    delete zp10000x1_;
    
    delete zp10x50_;
    delete zp50x50_;
    delete zp95x50_;
    delete zp200x50_;
    delete zp300x50_;
    delete zp10000x50_;

    delete zp10x150_;
    delete zp200x150_;
    delete zp295x150_;
    delete zp500x150_;
    delete zp1000x150_;
    delete zp10000x150_;
    
    delete zp10x500_;
    delete zp500x500_;
    delete zp995x500_;
    delete zp2000x500_;
    delete zp10000x500_;
    
    delete zp10x1000_;
    delete zp1000x1000_;
    delete zp1995x1000_;
    delete zp10000x1000_;
    
	delete Signal_region;
	delete sigMTstack;

}

void drawSignalRegionVarBin(TCut cut, TCut relCut, TCut SScut, TCut SSrelCut, double sfQCD, TCut normCut, TCut SSnormCut, bool doNorm, double sf, std::string parameter, int binnum, float bin[], std::string can_name, std::string chan, std::string syst, bool plot, bool createOutputShapes)
{

    std::cout << "SIGNAL CUT: " << cut.GetTitle() << std::endl;

	/* check if want blinded signal region met tail */
    if (plot)
    {
        if(choice2==1 && (parameter=="mt_tot"))
        {
            global_title = "Signal Region with Blinded Tail";
        }
        if(choice2==1 && (parameter=="met" || parameter=="mvamet"))
        {
            global_title = "Signal Region with Blind MET Tail";
        }
        if(choice2==1 && sigChoice==2 && (chan=="mt"||chan=="et"))
        {
            global_title = "Signal Region in TMVA (blinded above 0.5)";
        }
        
    }
    
	/* create a TCanvas */
	TCanvas * Signal_region = new TCanvas(can_name.c_str(),can_name.c_str(),canDim[0],canDim[1]);
	Signal_region->Divide(1,2);
    
    if (plot)
    {
        Signal_region->cd(1);
    }

	/* create some histograms and format them */
    
	TH1F * data_ = new TH1F("data_","",binnum,bin);
	data_->Sumw2(); 
	data_->SetMarkerStyle(dataMarker);

	TH1F * qcd_ = new TH1F("qcd_","qcd_",binnum,bin); 
	qcd_->Sumw2(); 
	qcd_->SetFillColor(colors["QCD"]);

	TH1F * ztt_ = new TH1F("ztt_","ztt_",binnum,bin);
	ztt_->Sumw2();
	ztt_->SetFillColor(colors["ZTT"]);
    
    TH1F * zl_ = new TH1F("zl_","zl_",binnum,bin);
	zl_->Sumw2();
	zl_->SetFillColor(colors["ZL"]);
    
    TH1F * zj_ = new TH1F("zj_","zj_",binnum,bin);
	zj_->Sumw2();
	zj_->SetFillColor(colors["ZJ"]);
	
	TH1F * w_ = new TH1F("w_","w_",binnum,bin); 
	w_->Sumw2(); 
	w_->SetFillColor(colors["W"]);
	
	TH1F * ttj_ = new TH1F("ttj_","ttj_",binnum,bin);
	ttj_->Sumw2();
	ttj_->SetFillColor(colors["TTJ"]);
    
    TH1F * ttt_ = new TH1F("ttt_","ttt_",binnum,bin);
	ttt_->Sumw2();
	ttt_->SetFillColor(colors["TTT"]);
    
    TH1F * tt_ = new TH1F("tt_","tt_",binnum,bin);
	tt_->Sumw2();
	tt_->SetFillColor(colors["TT"]);
	
	TH1F * vvt_ = new TH1F("vvt_","vvt_",binnum,bin);
	vvt_->Sumw2();
	vvt_->SetFillColor(colors["VV"]);
    
    TH1F * vvj_ = new TH1F("vvj_","vvj_",binnum,bin);
	vvj_->Sumw2();
	vvj_->SetFillColor(colors["VV"]);
    
    TH1F * zvv_ = new TH1F("zvv_","zvv_",binnum,bin);
	zvv_->Sumw2();
	zvv_->SetFillColor(colors["ZVV"]);
    
    TH1F * ewk_ = new TH1F("ewk_","ewk_",binnum,bin);
	ewk_->Sumw2();
	ewk_->SetFillColor(colors["EWK"]);
	
	TH1F * zhtautau_ = new TH1F("zhtautau_","zhtautau_",binnum,bin);
	zhtautau_->Sumw2();
	zhtautau_->SetFillColor(colors["ZHTauTau"]);

	TH1F * glugluhtautau_ = new TH1F("glugluhtautau_","glugluhtautau_",binnum,bin);
	glugluhtautau_->Sumw2();
	glugluhtautau_->SetFillColor(colors["GluGluHTauTau"]);
    
    TH1F * vbfhtautau_ = new TH1F("vbfhtautau_","vbfhtautau_",binnum,bin);
	vbfhtautau_->Sumw2();
    vbfhtautau_->SetFillColor(colors["VBFHTauTau"]);
    
	TH1F * smh_ = new TH1F("smh_","smh_",binnum,bin);
	smh_->Sumw2();
    smh_->SetFillColor(colors["SMH"]);
    
	TH1F * comb_ = new TH1F("comb_","comb_",binnum,bin);
    
	TH1F * zp600a300_ = new TH1F("zp600a300_","zp600a300_",binnum,bin);
	zp600a300_->Sumw2();
    TH1F * zp600a400_ = new TH1F("zp600a400_","zp600a400_",binnum,bin);
	zp600a400_->Sumw2();

	TH1F * zp800a300_ = new TH1F("zp800a300_","zp800a300_",binnum,bin);
	zp800a300_->Sumw2();
    TH1F * zp800a400_ = new TH1F("zp800a400_","zp800a400_",binnum,bin);
	zp800a400_->Sumw2();
    TH1F * zp800a500_ = new TH1F("zp800a500_","zp800a500_",binnum,bin);
	zp800a500_->Sumw2();
    TH1F * zp800a600_ = new TH1F("zp800a600_","zp800a600_",binnum,bin);
	zp800a600_->Sumw2();

	TH1F * zp1000a300_ = new TH1F("zp1000a300_","zp1000a300_",binnum,bin);
	zp1000a300_->Sumw2();
    TH1F * zp1000a400_ = new TH1F("zp1000a400_","zp1000a400_",binnum,bin);
	zp1000a400_->Sumw2();
    TH1F * zp1000a500_ = new TH1F("zp1000a500_","zp1000a500_",binnum,bin);
	zp1000a500_->Sumw2();
    TH1F * zp1000a600_ = new TH1F("zp1000a600_","zp1000a600_",binnum,bin);
	zp1000a600_->Sumw2();
	TH1F * zp1000a700_ = new TH1F("zp1000a700_","zp1000a700_",binnum,bin);
	zp1000a700_->Sumw2();
    TH1F * zp1000a800_ = new TH1F("zp1000a800_","zp1000a800_",binnum,bin);
	zp1000a800_->Sumw2();
    
    TH1F * zp1200a300_ = new TH1F("zp1200a300_","zp1200a300_",binnum,bin);
	zp1200a300_->Sumw2();
    TH1F * zp1200a400_ = new TH1F("zp1200a400_","zp1200a400_",binnum,bin);
	zp1200a400_->Sumw2();
    TH1F * zp1200a500_ = new TH1F("zp1200a500_","zp1200a500_",binnum,bin);
	zp1200a500_->Sumw2();
    TH1F * zp1200a600_ = new TH1F("zp1200a600_","zp1200a600_",binnum,bin);
	zp1200a600_->Sumw2();
	TH1F * zp1200a700_ = new TH1F("zp1200a700_","zp1200a700_",binnum,bin);
	zp1200a700_->Sumw2();
    TH1F * zp1200a800_ = new TH1F("zp1200a800_","zp1200a800_",binnum,bin);
	zp1200a800_->Sumw2();
    
    TH1F * zp1400a300_ = new TH1F("zp1400a300_","zp1400a300_",binnum,bin);
	zp1400a300_->Sumw2();
    TH1F * zp1400a400_ = new TH1F("zp1400a400_","zp1400a400_",binnum,bin);
	zp1400a400_->Sumw2();
    TH1F * zp1400a500_ = new TH1F("zp1400a500_","zp1400a500_",binnum,bin);
	zp1400a500_->Sumw2();
    TH1F * zp1400a600_ = new TH1F("zp1400a600_","zp1400a600_",binnum,bin);
	zp1400a600_->Sumw2();
	TH1F * zp1400a700_ = new TH1F("zp1400a700_","zp1400a700_",binnum,bin);
	zp1400a700_->Sumw2();
    TH1F * zp1400a800_ = new TH1F("zp1400a800_","zp1400a800_",binnum,bin);
	zp1400a800_->Sumw2();
    
    TH1F * zp1700a300_ = new TH1F("zp1700a300_","zp1700a300_",binnum,bin);
	zp1700a300_->Sumw2();
    TH1F * zp1700a400_ = new TH1F("zp1700a400_","zp1700a400_",binnum,bin);
	zp1700a400_->Sumw2();
    TH1F * zp1700a500_ = new TH1F("zp1700a500_","zp1700a500_",binnum,bin);
	zp1700a500_->Sumw2();
    TH1F * zp1700a600_ = new TH1F("zp1700a600_","zp1700a600_",binnum,bin);
	zp1700a600_->Sumw2();
	TH1F * zp1700a700_ = new TH1F("zp1700a700_","zp1700a700_",binnum,bin);
	zp1700a700_->Sumw2();
    TH1F * zp1700a800_ = new TH1F("zp1700a800_","zp1700a800_",binnum,bin);
	zp1700a800_->Sumw2();
    
    TH1F * zp2000a300_ = new TH1F("zp2000a300_","zp2000a300_",binnum,bin);
	zp2000a300_->Sumw2();
    TH1F * zp2000a400_ = new TH1F("zp2000a400_","zp2000a400_",binnum,bin);
	zp2000a400_->Sumw2();
    TH1F * zp2000a500_ = new TH1F("zp2000a500_","zp2000a500_",binnum,bin);
	zp2000a500_->Sumw2();
    TH1F * zp2000a600_ = new TH1F("zp2000a600_","zp2000a600_",binnum,bin);
	zp2000a600_->Sumw2();
	TH1F * zp2000a700_ = new TH1F("zp2000a700_","zp2000a700_",binnum,bin);
	zp2000a700_->Sumw2();
    TH1F * zp2000a800_ = new TH1F("zp2000a800_","zp2000a800_",binnum,bin);
	zp2000a800_->Sumw2();
    
    TH1F * zp2500a300_ = new TH1F("zp2500a300_","zp2500a300_",binnum,bin);
	zp2500a300_->Sumw2();
    TH1F * zp2500a400_ = new TH1F("zp2500a400_","zp2500a400_",binnum,bin);
	zp2500a400_->Sumw2();
    TH1F * zp2500a500_ = new TH1F("zp2500a500_","zp2500a500_",binnum,bin);
	zp2500a500_->Sumw2();
    TH1F * zp2500a600_ = new TH1F("zp2500a600_","zp2500a600_",binnum,bin);
	zp2500a600_->Sumw2();
	TH1F * zp2500a700_ = new TH1F("zp2500a700_","zp2500a700_",binnum,bin);
	zp2500a700_->Sumw2();
    TH1F * zp2500a800_ = new TH1F("zp2500a800_","zp2500a800_",binnum,bin);
	zp2500a800_->Sumw2();
    
    TH1F * zp10x1_ = new TH1F("zp10x1_","zp10x1_",binnum,bin);
	zp10x1_->Sumw2();
    TH1F * zp20x1_ = new TH1F("zp20x1_","zp20x1_",binnum,bin);
	zp20x1_->Sumw2();
    TH1F * zp50x1_ = new TH1F("zp50x1_","zp50x1_",binnum,bin);
	zp50x1_->Sumw2();
    TH1F * zp100x1_ = new TH1F("zp100x1_","zp100x1_",binnum,bin);
	zp100x1_->Sumw2();
    TH1F * zp300x1_ = new TH1F("zp300x1_","zp300x1_",binnum,bin);
	zp300x1_->Sumw2();
    TH1F * zp500x1_ = new TH1F("zp500x1_","zp500x1_",binnum,bin);
	zp500x1_->Sumw2();
    TH1F * zp1000x1_ = new TH1F("zp1000x1_","zp1000x1_",binnum,bin);
	zp1000x1_->Sumw2();
    TH1F * zp2000x1_ = new TH1F("zp2000x1_","zp2000x1_",binnum,bin);
	zp2000x1_->Sumw2();
    TH1F * zp10000x1_ = new TH1F("zp10000x1_","zp10000x1_",binnum,bin);
	zp10000x1_->Sumw2();
    
    TH1F * zp10x50_ = new TH1F("zp10x50_","zp10x50_",binnum,bin);
	zp10x50_->Sumw2();
    TH1F * zp50x50_ = new TH1F("zp50x50_","zp50x50_",binnum,bin);
	zp50x50_->Sumw2();
    TH1F * zp95x50_ = new TH1F("zp95x50_","zp95x50_",binnum,bin);
	zp95x50_->Sumw2();
    TH1F * zp200x50_ = new TH1F("zp200x50_","zp200x50_",binnum,bin);
	zp200x50_->Sumw2();
    TH1F * zp300x50_ = new TH1F("zp300x50_","zp300x50_",binnum,bin);
	zp300x50_->Sumw2();
    TH1F * zp10000x50_ = new TH1F("zp10000x50_","zp10000x50_",binnum,bin);
	zp10000x50_->Sumw2();

    TH1F * zp10x150_ = new TH1F("zp10x150_","zp10x150_",binnum,bin);
	zp10x150_->Sumw2();
    TH1F * zp200x150_ = new TH1F("zp200x150_","zp200x150_",binnum,bin);
	zp200x150_->Sumw2();
    TH1F * zp295x150_ = new TH1F("zp295x150_","zp295x150_",binnum,bin);
	zp295x150_->Sumw2();
    TH1F * zp500x150_ = new TH1F("zp500x150_","zp500x150_",binnum,bin);
	zp500x150_->Sumw2();
    TH1F * zp1000x150_ = new TH1F("zp1000x150_","zp1000x150_",binnum,bin);
	zp1000x150_->Sumw2();
    TH1F * zp10000x150_ = new TH1F("zp10000x150_","zp10000x150_",binnum,bin);
	zp10000x150_->Sumw2();
    
    TH1F * zp10x500_ = new TH1F("zp10x500_","zp10x500_",binnum,bin);
	zp10x500_->Sumw2();
    TH1F * zp500x500_ = new TH1F("zp500x500_","zp500x500_",binnum,bin);
	zp500x500_->Sumw2();
    TH1F * zp995x500_ = new TH1F("zp995x500_","zp995x500_",binnum,bin);
	zp995x500_->Sumw2();
    TH1F * zp2000x500_ = new TH1F("zp2000x500_","zp2000x500_",binnum,bin);
	zp2000x500_->Sumw2();
    TH1F * zp10000x500_ = new TH1F("zp10000x500_","zp10000x500_",binnum,bin);
	zp10000x500_->Sumw2();
    
    TH1F * zp10x1000_ = new TH1F("zp10x1000_","zp10x1000_",binnum,bin);
	zp10x1000_->Sumw2();
    TH1F * zp1000x1000_ = new TH1F("zp1000x1000_","zp1000x1000_",binnum,bin);
	zp1000x1000_->Sumw2();
    TH1F * zp1995x1000_ = new TH1F("zp1995x1000_","zp1995x1000_",binnum,bin);
	zp1995x1000_->Sumw2();
    TH1F * zp10000x1000_ = new TH1F("zp10000x1000_","zp10000x1000_",binnum,bin);
	zp10000x1000_->Sumw2();

    
    /* fill the hists */

	DATA->Draw((parameter+">>data_").c_str(),cut);
    
	W->Draw((parameter+">>w_").c_str(),relCut*weights["W"].c_str());
	ZTT->Draw((parameter+">>ztt_").c_str(),relCut*weights["ZTT"].c_str());
    ZL->Draw((parameter+">>zl_").c_str(),relCut*weights["ZL"].c_str());
    ZJ->Draw((parameter+">>zj_").c_str(),relCut*weights["ZJ"].c_str());
	VVT->Draw((parameter+">>vvt_").c_str(),relCut*weights["VV"].c_str());
    VVJ->Draw((parameter+">>vvj_").c_str(),relCut*weights["VV"].c_str());
    ZVV->Draw((parameter+">>zvv_").c_str(),cut*weights["ZVV"].c_str());
    EWK->Draw((parameter+">>ewk_").c_str(),cut*weights["EWK"].c_str());
	TTT->Draw((parameter+">>ttt_").c_str(),relCut*weights["TTT"].c_str());
    TTJ->Draw((parameter+">>ttj_").c_str(),relCut*weights["TTJ"].c_str());
	ZHTauTau->Draw((parameter+">>zhtautau_").c_str(),cut*weights["ZHTauTau"].c_str());
    VBFHTauTau->Draw((parameter+">>vbfhtautau_").c_str(),cut*weights["VBFHTauTau"].c_str());
    GluGluHTauTau->Draw((parameter+">>glugluhtautau_").c_str(),cut*weights["GluGluHTauTau"].c_str());
    
	MZP600_MA0300->Draw((parameter+">>zp600a300_").c_str(),cut*weights["MONO"].c_str());
    MZP600_MA0400->Draw((parameter+">>zp600a400_").c_str(),cut*weights["MONO"].c_str());
    
	MZP800_MA0300->Draw((parameter+">>zp800a300_").c_str(),cut*weights["MONO"].c_str());
    MZP800_MA0400->Draw((parameter+">>zp800a400_").c_str(),cut*weights["MONO"].c_str());
    MZP800_MA0500->Draw((parameter+">>zp800a500_").c_str(),cut*weights["MONO"].c_str());
    MZP800_MA0600->Draw((parameter+">>zp800a600_").c_str(),cut*weights["MONO"].c_str());
    
	MZP1000_MA0300->Draw((parameter+">>zp1000a300_").c_str(),cut*weights["MONO"].c_str());
    MZP1000_MA0400->Draw((parameter+">>zp1000a400_").c_str(),cut*weights["MONO"].c_str());
    MZP1000_MA0500->Draw((parameter+">>zp1000a500_").c_str(),cut*weights["MONO"].c_str());
    MZP1000_MA0600->Draw((parameter+">>zp1000a600_").c_str(),cut*weights["MONO"].c_str());
    MZP1000_MA0700->Draw((parameter+">>zp1000a700_").c_str(),cut*weights["MONO"].c_str());
    MZP1000_MA0800->Draw((parameter+">>zp1000a800_").c_str(),cut*weights["MONO"].c_str());
    
    MZP1200_MA0300->Draw((parameter+">>zp1200a300_").c_str(),cut*weights["MONO"].c_str());
    MZP1200_MA0400->Draw((parameter+">>zp1200a400_").c_str(),cut*weights["MONO"].c_str());
    MZP1200_MA0500->Draw((parameter+">>zp1200a500_").c_str(),cut*weights["MONO"].c_str());
    MZP1200_MA0600->Draw((parameter+">>zp1200a600_").c_str(),cut*weights["MONO"].c_str());
    MZP1200_MA0700->Draw((parameter+">>zp1200a700_").c_str(),cut*weights["MONO"].c_str());
    MZP1200_MA0800->Draw((parameter+">>zp1200a800_").c_str(),cut*weights["MONO"].c_str());

	MZP1400_MA0300->Draw((parameter+">>zp1400a300_").c_str(),cut*weights["MONO"].c_str());
    MZP1400_MA0400->Draw((parameter+">>zp1400a400_").c_str(),cut*weights["MONO"].c_str());
    MZP1400_MA0500->Draw((parameter+">>zp1400a500_").c_str(),cut*weights["MONO"].c_str());
    MZP1400_MA0600->Draw((parameter+">>zp1400a600_").c_str(),cut*weights["MONO"].c_str());
    MZP1400_MA0700->Draw((parameter+">>zp1400a700_").c_str(),cut*weights["MONO"].c_str());
    MZP1400_MA0800->Draw((parameter+">>zp1400a800_").c_str(),cut*weights["MONO"].c_str());
    
	MZP1700_MA0300->Draw((parameter+">>zp1700a300_").c_str(),cut*weights["MONO"].c_str());
    MZP1700_MA0400->Draw((parameter+">>zp1700a400_").c_str(),cut*weights["MONO"].c_str());
    MZP1700_MA0500->Draw((parameter+">>zp1700a500_").c_str(),cut*weights["MONO"].c_str());
    MZP1700_MA0600->Draw((parameter+">>zp1700a600_").c_str(),cut*weights["MONO"].c_str());
    MZP1700_MA0700->Draw((parameter+">>zp1700a700_").c_str(),cut*weights["MONO"].c_str());
    MZP1700_MA0800->Draw((parameter+">>zp1700a800_").c_str(),cut*weights["MONO"].c_str());
    
	MZP2000_MA0300->Draw((parameter+">>zp2000a300_").c_str(),cut*weights["MONO"].c_str());
    MZP2000_MA0400->Draw((parameter+">>zp2000a400_").c_str(),cut*weights["MONO"].c_str());
    MZP2000_MA0500->Draw((parameter+">>zp2000a500_").c_str(),cut*weights["MONO"].c_str());
    MZP2000_MA0600->Draw((parameter+">>zp2000a600_").c_str(),cut*weights["MONO"].c_str());
    MZP2000_MA0700->Draw((parameter+">>zp2000a700_").c_str(),cut*weights["MONO"].c_str());
    MZP2000_MA0800->Draw((parameter+">>zp2000a800_").c_str(),cut*weights["MONO"].c_str());
    
	MZP2500_MA0300->Draw((parameter+">>zp2500a300_").c_str(),cut*weights["MONO"].c_str());
    MZP2500_MA0400->Draw((parameter+">>zp2500a400_").c_str(),cut*weights["MONO"].c_str());
    MZP2500_MA0500->Draw((parameter+">>zp2500a500_").c_str(),cut*weights["MONO"].c_str());
    MZP2500_MA0600->Draw((parameter+">>zp2500a600_").c_str(),cut*weights["MONO"].c_str());
    MZP2500_MA0700->Draw((parameter+">>zp2500a700_").c_str(),cut*weights["MONO"].c_str());
    MZP2500_MA0800->Draw((parameter+">>zp2500a800_").c_str(),cut*weights["MONO"].c_str());
    
    MZP10_MX1->Draw((parameter+">>zp10x1_").c_str(),cut*weights["MONO"].c_str());
    MZP20_MX1->Draw((parameter+">>zp20x1_").c_str(),cut*weights["MONO"].c_str());
    MZP50_MX1->Draw((parameter+">>zp50x1_").c_str(),cut*weights["MONO"].c_str());
    MZP100_MX1->Draw((parameter+">>zp100x1_").c_str(),cut*weights["MONO"].c_str());
    MZP300_MX1->Draw((parameter+">>zp300x1_").c_str(),cut*weights["MONO"].c_str());
    MZP500_MX1->Draw((parameter+">>zp500x1_").c_str(),cut*weights["MONO"].c_str());
    MZP1000_MX1->Draw((parameter+">>zp1000x1_").c_str(),cut*weights["MONO"].c_str());
    MZP2000_MX1->Draw((parameter+">>zp2000x1_").c_str(),cut*weights["MONO"].c_str());
    MZP10000_MX1->Draw((parameter+">>zp10000x1_").c_str(),cut*weights["MONO"].c_str());
    
    MZP10_MX50->Draw((parameter+">>zp10x50_").c_str(),cut*weights["MONO"].c_str());
    MZP50_MX50->Draw((parameter+">>zp50x50_").c_str(),cut*weights["MONO"].c_str());
    MZP95_MX50->Draw((parameter+">>zp95x50_").c_str(),cut*weights["MONO"].c_str());
    MZP200_MX50->Draw((parameter+">>zp200x50_").c_str(),cut*weights["MONO"].c_str());
    MZP300_MX50->Draw((parameter+">>zp300x50_").c_str(),cut*weights["MONO"].c_str());
    MZP10000_MX50->Draw((parameter+">>zp10000x50_").c_str(),cut*weights["MONO"].c_str());

    MZP10_MX150->Draw((parameter+">>zp10x150_").c_str(),cut*weights["MONO"].c_str());
    MZP200_MX150->Draw((parameter+">>zp200x150_").c_str(),cut*weights["MONO"].c_str());
    MZP295_MX150->Draw((parameter+">>zp295x150_").c_str(),cut*weights["MONO"].c_str());
    MZP500_MX150->Draw((parameter+">>zp500x150_").c_str(),cut*weights["MONO"].c_str());
    MZP1000_MX150->Draw((parameter+">>zp1000x150_").c_str(),cut*weights["MONO"].c_str());
    MZP10000_MX150->Draw((parameter+">>zp10000x150_").c_str(),cut*weights["MONO"].c_str());
    
    MZP10_MX500->Draw((parameter+">>zp10x500_").c_str(),cut*weights["MONO"].c_str());
    MZP500_MX500->Draw((parameter+">>zp500x500_").c_str(),cut*weights["MONO"].c_str());
    MZP995_MX500->Draw((parameter+">>zp995x500_").c_str(),cut*weights["MONO"].c_str());
    MZP2000_MX500->Draw((parameter+">>zp2000x500_").c_str(),cut*weights["MONO"].c_str());
    MZP10000_MX500->Draw((parameter+">>zp10000x500_").c_str(),cut*weights["MONO"].c_str());
    
    MZP10_MX1000->Draw((parameter+">>zp10x1000_").c_str(),cut*weights["MONO"].c_str());
    MZP1000_MX1000->Draw((parameter+">>zp1000x1000_").c_str(),cut*weights["MONO"].c_str());
    MZP1995_MX1000->Draw((parameter+">>zp1995x1000_").c_str(),cut*weights["MONO"].c_str());
    MZP10000_MX1000->Draw((parameter+">>zp10000x1000_").c_str(),cut*weights["MONO"].c_str());
    
    
	fillQCD_ShapeVarBin(SScut, SSrelCut, sfQCD, SSnormCut, doNorm, qcd_, parameter, binnum, bin, "QCD SS " + parameter + " sideband in " + chan + " channel", syst, chan, plot, createOutputShapes);
    
	qcd_->Scale(qcdOStoSS);
    
    if (plot)
    {
        Signal_region->cd(1);
    }
    
    //not scaling by SF in datacards
	//w_->Scale(sf);
    w_->Scale(shapeScaleVarBin(cut, relCut, W, weights["W"], binnum, bin));
    vvt_->Scale(shapeScaleVarBin(cut, relCut, VVT, weights["VV"], binnum, bin));
    vvj_->Scale(shapeScaleVarBin(cut, relCut, VVJ, weights["VV"], binnum, bin));
    ttt_->Scale(shapeScaleVarBin(cut, relCut, TTT, weights["TTT"], binnum, bin));
    ttj_->Scale(shapeScaleVarBin(cut, relCut, TTJ, weights["TTJ"], binnum, bin));
    ztt_->Scale(shapeScaleVarBin(cut, relCut, ZTT, weights["ZTT"], binnum, bin));
    zl_->Scale(shapeScaleVarBin(cut, relCut, ZL, weights["ZL"], binnum, bin));
    zj_->Scale(shapeScaleVarBin(cut, relCut, ZJ, weights["ZJ"], binnum, bin));

    
    // scale norms to nominal from variant if necessary
    if(doNorm)
    {
        
        vvt_->Scale(shapeScaleVarBin(normCut, cut, VVT, weights["VV"], binnum, bin));
        vvj_->Scale(shapeScaleVarBin(normCut, cut, VVJ, weights["VV"], binnum, bin));
        zvv_->Scale(shapeScaleVarBin(normCut, cut, ZVV, weights["ZVV"], binnum, bin));
        ewk_->Scale(shapeScaleVarBin(normCut, cut, EWK, weights["EWK"], binnum, bin));
        ttt_->Scale(shapeScaleVarBin(normCut, cut, TTT, weights["TTT"], binnum, bin));
        ttj_->Scale(shapeScaleVarBin(normCut, cut, TTJ, weights["TTJ"], binnum, bin));
        ztt_->Scale(shapeScaleVarBin(normCut, cut, ZTT, weights["ZTT"], binnum, bin));
        zl_->Scale(shapeScaleVarBin(normCut, cut, ZL, weights["ZL"], binnum, bin));
        zj_->Scale(shapeScaleVarBin(normCut, cut, ZJ, weights["ZJ"], binnum, bin));
        w_->Scale(shapeScaleVarBin(normCut, cut, W, weights["W"], binnum, bin));
        zhtautau_->Scale(shapeScaleVarBin(normCut, cut, ZHTauTau, weights["ZHTauTau"], binnum, bin));
        glugluhtautau_->Scale(shapeScaleVarBin(normCut, cut, GluGluHTauTau, weights["GluGluHTauTau"], binnum, bin));
        vbfhtautau_->Scale(shapeScaleVarBin(normCut, cut, VBFHTauTau, weights["VBFHTauTau"], binnum, bin));

        //QCD scaled in shape fillQCDShape MZP600_MA0400
        
        zp600a300_->Scale(shapeScaleVarBin(normCut, cut, MZP600_MA0300, weights["MONO"], binnum, bin));
        zp600a400_->Scale(shapeScaleVarBin(normCut, cut, MZP600_MA0400, weights["MONO"], binnum, bin));
        zp800a300_->Scale(shapeScaleVarBin(normCut, cut, MZP800_MA0300, weights["MONO"], binnum, bin));
        zp800a400_->Scale(shapeScaleVarBin(normCut, cut, MZP800_MA0400, weights["MONO"], binnum, bin));
        zp800a500_->Scale(shapeScaleVarBin(normCut, cut, MZP800_MA0500, weights["MONO"], binnum, bin));
        zp800a600_->Scale(shapeScaleVarBin(normCut, cut, MZP800_MA0600, weights["MONO"], binnum, bin));
        zp1000a300_->Scale(shapeScaleVarBin(normCut, cut, MZP1000_MA0300, weights["MONO"], binnum, bin));
        zp1000a400_->Scale(shapeScaleVarBin(normCut, cut, MZP1000_MA0400, weights["MONO"], binnum, bin));
        zp1000a500_->Scale(shapeScaleVarBin(normCut, cut, MZP1000_MA0500, weights["MONO"], binnum, bin));
        zp1000a600_->Scale(shapeScaleVarBin(normCut, cut, MZP1000_MA0600, weights["MONO"], binnum, bin));
        zp1000a700_->Scale(shapeScaleVarBin(normCut, cut, MZP1000_MA0700, weights["MONO"], binnum, bin));
        zp1000a800_->Scale(shapeScaleVarBin(normCut, cut, MZP1000_MA0800, weights["MONO"], binnum, bin));
        zp1200a300_->Scale(shapeScaleVarBin(normCut, cut, MZP1200_MA0300, weights["MONO"], binnum, bin));
        zp1200a400_->Scale(shapeScaleVarBin(normCut, cut, MZP1200_MA0400, weights["MONO"], binnum, bin));
        zp1200a500_->Scale(shapeScaleVarBin(normCut, cut, MZP1200_MA0500, weights["MONO"], binnum, bin));
        zp1200a600_->Scale(shapeScaleVarBin(normCut, cut, MZP1200_MA0600, weights["MONO"], binnum, bin));
        zp1200a700_->Scale(shapeScaleVarBin(normCut, cut, MZP1200_MA0700, weights["MONO"], binnum, bin));
        zp1200a800_->Scale(shapeScaleVarBin(normCut, cut, MZP1200_MA0800, weights["MONO"], binnum, bin));
        zp1400a300_->Scale(shapeScaleVarBin(normCut, cut, MZP1400_MA0300, weights["MONO"], binnum, bin));
        zp1400a400_->Scale(shapeScaleVarBin(normCut, cut, MZP1400_MA0400, weights["MONO"], binnum, bin));
        zp1400a500_->Scale(shapeScaleVarBin(normCut, cut, MZP1400_MA0500, weights["MONO"], binnum, bin));
        zp1400a600_->Scale(shapeScaleVarBin(normCut, cut, MZP1400_MA0600, weights["MONO"], binnum, bin));
        zp1400a700_->Scale(shapeScaleVarBin(normCut, cut, MZP1400_MA0700, weights["MONO"], binnum, bin));
        zp1400a800_->Scale(shapeScaleVarBin(normCut, cut, MZP1400_MA0800, weights["MONO"], binnum, bin));
        zp1700a300_->Scale(shapeScaleVarBin(normCut, cut, MZP1700_MA0300, weights["MONO"], binnum, bin));
        zp1700a400_->Scale(shapeScaleVarBin(normCut, cut, MZP1700_MA0400, weights["MONO"], binnum, bin));
        zp1700a500_->Scale(shapeScaleVarBin(normCut, cut, MZP1700_MA0500, weights["MONO"], binnum, bin));
        zp1700a600_->Scale(shapeScaleVarBin(normCut, cut, MZP1700_MA0600, weights["MONO"], binnum, bin));
        zp1700a700_->Scale(shapeScaleVarBin(normCut, cut, MZP1700_MA0700, weights["MONO"], binnum, bin));
        zp1700a800_->Scale(shapeScaleVarBin(normCut, cut, MZP1700_MA0800, weights["MONO"], binnum, bin));
        zp2000a300_->Scale(shapeScaleVarBin(normCut, cut, MZP2000_MA0300, weights["MONO"], binnum, bin));
        zp2000a400_->Scale(shapeScaleVarBin(normCut, cut, MZP2000_MA0400, weights["MONO"], binnum, bin));
        zp2000a500_->Scale(shapeScaleVarBin(normCut, cut, MZP2000_MA0500, weights["MONO"], binnum, bin));
        zp2000a600_->Scale(shapeScaleVarBin(normCut, cut, MZP2000_MA0600, weights["MONO"], binnum, bin));
        zp2000a700_->Scale(shapeScaleVarBin(normCut, cut, MZP2000_MA0700, weights["MONO"], binnum, bin));
        zp2000a800_->Scale(shapeScaleVarBin(normCut, cut, MZP2000_MA0800, weights["MONO"], binnum, bin));
        zp2500a300_->Scale(shapeScaleVarBin(normCut, cut, MZP2500_MA0300, weights["MONO"], binnum, bin));
        zp2500a400_->Scale(shapeScaleVarBin(normCut, cut, MZP2500_MA0400, weights["MONO"], binnum, bin));
        zp2500a500_->Scale(shapeScaleVarBin(normCut, cut, MZP2500_MA0500, weights["MONO"], binnum, bin));
        zp2500a600_->Scale(shapeScaleVarBin(normCut, cut, MZP2500_MA0600, weights["MONO"], binnum, bin));
        zp2500a700_->Scale(shapeScaleVarBin(normCut, cut, MZP2500_MA0700, weights["MONO"], binnum, bin));
        zp2500a800_->Scale(shapeScaleVarBin(normCut, cut, MZP2500_MA0800, weights["MONO"], binnum, bin));
        
        zp10x1_->Scale(shapeScaleVarBin(normCut, cut, MZP10_MX1, weights["MONO"], binnum, bin));
        zp20x1_->Scale(shapeScaleVarBin(normCut, cut, MZP20_MX1, weights["MONO"], binnum, bin));
        zp50x1_->Scale(shapeScaleVarBin(normCut, cut, MZP50_MX1, weights["MONO"], binnum, bin));
        zp100x1_->Scale(shapeScaleVarBin(normCut, cut, MZP100_MX1, weights["MONO"], binnum, bin));
        zp300x1_->Scale(shapeScaleVarBin(normCut, cut, MZP300_MX1, weights["MONO"], binnum, bin));
        zp500x1_->Scale(shapeScaleVarBin(normCut, cut, MZP500_MX1, weights["MONO"], binnum, bin));
        zp1000x1_->Scale(shapeScaleVarBin(normCut, cut, MZP1000_MX1, weights["MONO"], binnum, bin));
        zp2000x1_->Scale(shapeScaleVarBin(normCut, cut, MZP2000_MX1, weights["MONO"], binnum, bin));
        zp10000x1_->Scale(shapeScaleVarBin(normCut, cut, MZP10000_MX1, weights["MONO"], binnum, bin));
        
        zp10x50_->Scale(shapeScaleVarBin(normCut, cut, MZP10_MX50, weights["MONO"], binnum, bin));
        zp50x50_->Scale(shapeScaleVarBin(normCut, cut, MZP50_MX50, weights["MONO"], binnum, bin));
        zp95x50_->Scale(shapeScaleVarBin(normCut, cut, MZP95_MX50, weights["MONO"], binnum, bin));
        zp200x50_->Scale(shapeScaleVarBin(normCut, cut, MZP200_MX50, weights["MONO"], binnum, bin));
        zp300x50_->Scale(shapeScaleVarBin(normCut, cut, MZP300_MX50, weights["MONO"], binnum, bin));
        zp10000x50_->Scale(shapeScaleVarBin(normCut, cut, MZP10000_MX50, weights["MONO"], binnum, bin));

        zp10x150_->Scale(shapeScaleVarBin(normCut, cut, MZP10_MX150, weights["MONO"], binnum, bin));
        zp200x150_->Scale(shapeScaleVarBin(normCut, cut, MZP200_MX150, weights["MONO"], binnum, bin));
        zp295x150_->Scale(shapeScaleVarBin(normCut, cut, MZP295_MX150, weights["MONO"], binnum, bin));
        zp500x150_->Scale(shapeScaleVarBin(normCut, cut, MZP500_MX150, weights["MONO"], binnum, bin));
        zp1000x150_->Scale(shapeScaleVarBin(normCut, cut, MZP1000_MX150, weights["MONO"], binnum, bin));
        zp10000x150_->Scale(shapeScaleVarBin(normCut, cut, MZP10000_MX150, weights["MONO"], binnum, bin));
        
        zp10x500_->Scale(shapeScaleVarBin(normCut, cut, MZP10_MX500, weights["MONO"], binnum, bin));
        zp500x500_->Scale(shapeScaleVarBin(normCut, cut, MZP500_MX500, weights["MONO"], binnum, bin));
        zp995x500_->Scale(shapeScaleVarBin(normCut, cut, MZP995_MX500, weights["MONO"], binnum, bin));
        zp2000x500_->Scale(shapeScaleVarBin(normCut, cut, MZP2000_MX500, weights["MONO"], binnum, bin));
        zp10000x500_->Scale(shapeScaleVarBin(normCut, cut, MZP10000_MX500, weights["MONO"], binnum, bin));
        
        zp10x1000_->Scale(shapeScaleVarBin(normCut, cut, MZP10_MX1000, weights["MONO"], binnum, bin));
        zp1000x1000_->Scale(shapeScaleVarBin(normCut, cut, MZP1000_MX1000, weights["MONO"], binnum, bin));
        zp1995x1000_->Scale(shapeScaleVarBin(normCut, cut, MZP1995_MX1000, weights["MONO"], binnum, bin));
        zp10000x1000_->Scale(shapeScaleVarBin(normCut, cut, MZP10000_MX1000, weights["MONO"], binnum, bin));
        
        
    }
    
    // scale things to the projected lumi
	vvt_->Scale(lumi_sf);
    vvj_->Scale(lumi_sf);
    zvv_->Scale(lumi_sf);
    ewk_->Scale(lumi_sf);
	ttt_->Scale(lumi_sf);
    ttj_->Scale(lumi_sf);
	ztt_->Scale(lumi_sf);
    zl_->Scale(lumi_sf);
    zj_->Scale(lumi_sf);
	w_->Scale(lumi_sf);
    zhtautau_->Scale(lumi_sf);
    glugluhtautau_->Scale(lumi_sf);
    vbfhtautau_->Scale(lumi_sf);
    
    qcd_->Scale(lumi_sf);

	zp600a300_->Scale(lumi_sf);
    zp600a400_->Scale(lumi_sf);
	zp800a300_->Scale(lumi_sf);
    zp800a400_->Scale(lumi_sf);
    zp800a500_->Scale(lumi_sf);
    zp800a600_->Scale(lumi_sf);
	zp1000a300_->Scale(lumi_sf);
    zp1000a400_->Scale(lumi_sf);
    zp1000a500_->Scale(lumi_sf);
    zp1000a600_->Scale(lumi_sf);
    zp1000a700_->Scale(lumi_sf);
    zp1000a800_->Scale(lumi_sf);
	zp1200a300_->Scale(lumi_sf);
    zp1200a400_->Scale(lumi_sf);
    zp1200a500_->Scale(lumi_sf);
    zp1200a600_->Scale(lumi_sf);
    zp1200a700_->Scale(lumi_sf);
    zp1200a800_->Scale(lumi_sf);
    zp1400a300_->Scale(lumi_sf);
    zp1400a400_->Scale(lumi_sf);
    zp1400a500_->Scale(lumi_sf);
    zp1400a600_->Scale(lumi_sf);
    zp1400a700_->Scale(lumi_sf);
    zp1400a800_->Scale(lumi_sf);
	zp1700a300_->Scale(lumi_sf);
    zp1700a400_->Scale(lumi_sf);
    zp1700a500_->Scale(lumi_sf);
    zp1700a600_->Scale(lumi_sf);
    zp1700a700_->Scale(lumi_sf);
    zp1700a800_->Scale(lumi_sf);
	zp2000a300_->Scale(lumi_sf);
    zp2000a400_->Scale(lumi_sf);
    zp2000a500_->Scale(lumi_sf);
    zp2000a600_->Scale(lumi_sf);
    zp2000a700_->Scale(lumi_sf);
    zp2000a800_->Scale(lumi_sf);
	zp2500a300_->Scale(lumi_sf);
    zp2500a400_->Scale(lumi_sf);
    zp2500a500_->Scale(lumi_sf);
    zp2500a600_->Scale(lumi_sf);
    zp2500a700_->Scale(lumi_sf);
    zp2500a800_->Scale(lumi_sf);
    
    zp10x1_->Scale(lumi_sf);
    zp20x1_->Scale(lumi_sf);
    zp50x1_->Scale(lumi_sf);
    zp100x1_->Scale(lumi_sf);
    zp300x1_->Scale(lumi_sf);
    zp500x1_->Scale(lumi_sf);
    zp1000x1_->Scale(lumi_sf);
    zp2000x1_->Scale(lumi_sf);
    zp10000x1_->Scale(lumi_sf);
    
    zp10x50_->Scale(lumi_sf);
    zp50x50_->Scale(lumi_sf);
    zp95x50_->Scale(lumi_sf);
    zp200x50_->Scale(lumi_sf);
    zp300x50_->Scale(lumi_sf);
    zp10000x50_->Scale(lumi_sf);

    zp10x150_->Scale(lumi_sf);
    zp200x150_->Scale(lumi_sf);
    zp295x150_->Scale(lumi_sf);
    zp500x150_->Scale(lumi_sf);
    zp1000x150_->Scale(lumi_sf);
    zp10000x150_->Scale(lumi_sf);
    
    zp10x500_->Scale(lumi_sf);
    zp500x500_->Scale(lumi_sf);
    zp995x500_->Scale(lumi_sf);
    zp2000x500_->Scale(lumi_sf);
    zp10000x500_->Scale(lumi_sf);
    
    zp10x1000_->Scale(lumi_sf);
    zp1000x1000_->Scale(lumi_sf);
    zp1995x1000_->Scale(lumi_sf);
    zp10000x1000_->Scale(lumi_sf);

    // make combined histos and scale
    
    smh_->Add(zhtautau_);
    smh_->Add(glugluhtautau_);
    smh_->Add(vbfhtautau_);
    tt_->Add(ttt_);
    tt_->Add(ttj_);
    
    smh_->Scale(lumi_sf);
    tt_->Scale(lumi_sf);
    
    comb_->Add(vvt_);
    comb_->Add(vvj_);
    comb_->Add(zvv_);
    comb_->Add(ewk_);
    comb_->Add(ttt_);
    comb_->Add(ttj_);
    comb_->Add(w_);
    comb_->Add(ztt_);
    comb_->Add(zj_);
    comb_->Add(zl_);
    comb_->Add(smh_);
    comb_->Add(qcd_);

    std::cout << "Editing Shape File" << std::endl;
    
	if(createOutputShapes)
	{
		/* create a file to help with S/sqrt(S+B) cut opt */
        
        std::string paramName;
        if (parameter=="met" || parameter=="mvamet") {paramName = "met";}
        else if (parameter.substr(0,6)=="mvaVar") {paramName = "mva";}
        else if (parameter.substr(0,6) == "mt_tot") {paramName="mt";}
        else {paramName = parameter;}
        
        std::string catName = chan + "_inclusive";
		std::string optName = "xtt_" + chan + ".inputs-13TeV-" + paramName + ".root";
		TFile sbOp(optName.c_str(),"UPDATE");
        
        TDirectory * catdir = sbOp.GetDirectory(catName.c_str());
        if (catdir){catdir->cd();}
        else {sbOp.mkdir(catName.c_str()); TDirectory * catdir = sbOp.GetDirectory(catName.c_str()); catdir->cd();}
        
        std::string DATAname = "data_obs" + syst;
        const char *DATAsave = DATAname.c_str();
        std::string VVTname = "VVT" + syst;
        const char *VVTsave = VVTname.c_str();
        std::string VVJname = "VVJ" + syst;
        const char *VVJsave = VVJname.c_str();
        std::string ZVVname = "ZVV" + syst;
        const char *ZVVsave = ZVVname.c_str();
        std::string EWKname = "EWK" + syst;
        const char *EWKsave = EWKname.c_str();
        std::string TTTname = "TTT" + syst;
        const char *TTTsave = TTTname.c_str();
        std::string TTJname = "TTJ" + syst;
        const char *TTJsave = TTJname.c_str();
        std::string ZTTname = "ZTT" + syst;
        const char *ZTTsave = ZTTname.c_str();
        std::string ZLname = "ZL" + syst;
        const char *ZLsave = ZLname.c_str();
        std::string ZJname = "ZJ" + syst;
        const char *ZJsave = ZJname.c_str();
        std::string Wname = "W" + syst;
        const char *Wsave = Wname.c_str();
        std::string QCDname = "QCD" + syst;
        const char *QCDsave = QCDname.c_str();
        std::string SMHname = "SMH" + syst;
        const char *SMHsave = SMHname.c_str();
        std::string ZHname = "ZH" + syst;
        const char *ZHsave = ZHname.c_str();
        
        std::string combname = "comb" + syst;
        const char *combsave = combname.c_str();
        
        std::string ZprimeA300Zname = "Zprime" + syst + "A300";
        const char *ZprimeA300Zsave = ZprimeA300Zname.c_str();
        
        std::string Zprime600A300name = "Zprime600A300" + syst;
        const char *Zprime600A300save = Zprime600A300name.c_str();
        std::string Zprime600A400name = "Zprime600A400" + syst;
        const char *Zprime600A400save = Zprime600A400name.c_str();
        
        std::string Zprime800A300name = "Zprime800A300" + syst;
        const char *Zprime800A300save = Zprime800A300name.c_str();
        std::string Zprime800A400name = "Zprime800A400" + syst;
        const char *Zprime800A400save = Zprime800A400name.c_str();
        std::string Zprime800A500name = "Zprime800A500" + syst;
        const char *Zprime800A500save = Zprime800A500name.c_str();
        std::string Zprime800A600name = "Zprime800A600" + syst;
        const char *Zprime800A600save = Zprime800A600name.c_str();
        
        std::string Zprime1000A300name = "Zprime1000A300" + syst;
        const char *Zprime1000A300save = Zprime1000A300name.c_str();
        std::string Zprime1000A400name = "Zprime1000A400" + syst;
        const char *Zprime1000A400save = Zprime1000A400name.c_str();
        std::string Zprime1000A500name = "Zprime1000A500" + syst;
        const char *Zprime1000A500save = Zprime1000A500name.c_str();
        std::string Zprime1000A600name = "Zprime1000A600" + syst;
        const char *Zprime1000A600save = Zprime1000A600name.c_str();
        std::string Zprime1000A700name = "Zprime1000A700" + syst;
        const char *Zprime1000A700save = Zprime1000A700name.c_str();
        std::string Zprime1000A800name = "Zprime1000A800" + syst;
        const char *Zprime1000A800save = Zprime1000A800name.c_str();
        
        std::string Zprime1200A300name = "Zprime1200A300" + syst;
        const char *Zprime1200A300save = Zprime1200A300name.c_str();
        std::string Zprime1200A400name = "Zprime1200A400" + syst;
        const char *Zprime1200A400save = Zprime1200A400name.c_str();
        std::string Zprime1200A500name = "Zprime1200A500" + syst;
        const char *Zprime1200A500save = Zprime1200A500name.c_str();
        std::string Zprime1200A600name = "Zprime1200A600" + syst;
        const char *Zprime1200A600save = Zprime1200A600name.c_str();
        std::string Zprime1200A700name = "Zprime1200A700" + syst;
        const char *Zprime1200A700save = Zprime1200A700name.c_str();
        std::string Zprime1200A800name = "Zprime1200A800" + syst;
        const char *Zprime1200A800save = Zprime1200A800name.c_str();
        
        std::string Zprime1400A300name = "Zprime1400A300" + syst;
        const char *Zprime1400A300save = Zprime1400A300name.c_str();
        std::string Zprime1400A400name = "Zprime1400A400" + syst;
        const char *Zprime1400A400save = Zprime1400A400name.c_str();
        std::string Zprime1400A500name = "Zprime1400A500" + syst;
        const char *Zprime1400A500save = Zprime1400A500name.c_str();
        std::string Zprime1400A600name = "Zprime1400A600" + syst;
        const char *Zprime1400A600save = Zprime1400A600name.c_str();
        std::string Zprime1400A700name = "Zprime1400A700" + syst;
        const char *Zprime1400A700save = Zprime1400A700name.c_str();
        std::string Zprime1400A800name = "Zprime1400A800" + syst;
        const char *Zprime1400A800save = Zprime1400A800name.c_str();
        
        std::string Zprime1700A300name = "Zprime1700A300" + syst;
        const char *Zprime1700A300save = Zprime1700A300name.c_str();
        std::string Zprime1700A400name = "Zprime1700A400" + syst;
        const char *Zprime1700A400save = Zprime1700A400name.c_str();
        std::string Zprime1700A500name = "Zprime1700A500" + syst;
        const char *Zprime1700A500save = Zprime1700A500name.c_str();
        std::string Zprime1700A600name = "Zprime1700A600" + syst;
        const char *Zprime1700A600save = Zprime1700A600name.c_str();
        std::string Zprime1700A700name = "Zprime1700A700" + syst;
        const char *Zprime1700A700save = Zprime1700A700name.c_str();
        std::string Zprime1700A800name = "Zprime1700A800" + syst;
        const char *Zprime1700A800save = Zprime1700A800name.c_str();
        
        std::string Zprime2000A300name = "Zprime2000A300" + syst;
        const char *Zprime2000A300save = Zprime2000A300name.c_str();
        std::string Zprime2000A400name = "Zprime2000A400" + syst;
        const char *Zprime2000A400save = Zprime2000A400name.c_str();
        std::string Zprime2000A500name = "Zprime2000A500" + syst;
        const char *Zprime2000A500save = Zprime2000A500name.c_str();
        std::string Zprime2000A600name = "Zprime2000A600" + syst;
        const char *Zprime2000A600save = Zprime2000A600name.c_str();
        std::string Zprime2000A700name = "Zprime2000A700" + syst;
        const char *Zprime2000A700save = Zprime2000A700name.c_str();
        std::string Zprime2000A800name = "Zprime2000A800" + syst;
        const char *Zprime2000A800save = Zprime2000A800name.c_str();
        
        std::string Zprime2500A300name = "Zprime2500A300" + syst;
        const char *Zprime2500A300save = Zprime2500A300name.c_str();
        std::string Zprime2500A400name = "Zprime2500A400" + syst;
        const char *Zprime2500A400save = Zprime2500A400name.c_str();
        std::string Zprime2500A500name = "Zprime2500A500" + syst;
        const char *Zprime2500A500save = Zprime2500A500name.c_str();
        std::string Zprime2500A600name = "Zprime2500A600" + syst;
        const char *Zprime2500A600save = Zprime2500A600name.c_str();
        std::string Zprime2500A700name = "Zprime2500A700" + syst;
        const char *Zprime2500A700save = Zprime2500A700name.c_str();
        std::string Zprime2500A800name = "Zprime2500A800" + syst;
        const char *Zprime2500A800save = Zprime2500A800name.c_str();
        
        std::string Zprime10X1name = "Zprime10X1" + syst;
        const char *Zprime10X1save = Zprime10X1name.c_str();
        std::string Zprime20X1name = "Zprime20X1" + syst;
        const char *Zprime20X1save = Zprime20X1name.c_str();
        std::string Zprime50X1name = "Zprime50X1" + syst;
        const char *Zprime50X1save = Zprime50X1name.c_str();
        std::string Zprime100X1name = "Zprime100X1" + syst;
        const char *Zprime100X1save = Zprime100X1name.c_str();
        std::string Zprime300X1name = "Zprime300X1" + syst;
        const char *Zprime300X1save = Zprime300X1name.c_str();
        std::string Zprime500X1name = "Zprime500X1" + syst;
        const char *Zprime500X1save = Zprime500X1name.c_str();
        std::string Zprime1000X1name = "Zprime1000X1" + syst;
        const char *Zprime1000X1save = Zprime1000X1name.c_str();
        std::string Zprime2000X1name = "Zprime2000X1" + syst;
        const char *Zprime2000X1save = Zprime2000X1name.c_str();
        std::string Zprime10000X1name = "Zprime10000X1" + syst;
        const char *Zprime10000X1save = Zprime10000X1name.c_str();
        
        std::string Zprime10X50name = "Zprime10X50" + syst;
        const char *Zprime10X50save = Zprime10X50name.c_str();
        std::string Zprime50X50name = "Zprime50X50" + syst;
        const char *Zprime50X50save = Zprime50X50name.c_str();
        std::string Zprime95X50name = "Zprime95X50" + syst;
        const char *Zprime95X50save = Zprime95X50name.c_str();
        std::string Zprime200X50name = "Zprime200X50" + syst;
        const char *Zprime200X50save = Zprime200X50name.c_str();
        std::string Zprime300X50name = "Zprime300X50" + syst;
        const char *Zprime300X50save = Zprime300X50name.c_str();
        std::string Zprime10000X50name = "Zprime10000X50" + syst;
        const char *Zprime10000X50save = Zprime10000X50name.c_str();
    
        std::string Zprime10X150name = "Zprime10X150" + syst;
        const char *Zprime10X150save = Zprime10X150name.c_str();
        std::string Zprime200X150name = "Zprime200X150" + syst;
        const char *Zprime200X150save = Zprime200X150name.c_str();
        std::string Zprime295X150name = "Zprime295X150" + syst;
        const char *Zprime295X150save = Zprime295X150name.c_str();
        std::string Zprime500X150name = "Zprime500X150" + syst;
        const char *Zprime500X150save = Zprime500X150name.c_str();
        std::string Zprime1000X150name = "Zprime1000X150" + syst;
        const char *Zprime1000X150save = Zprime1000X150name.c_str();
        std::string Zprime10000X150name = "Zprime10000X150" + syst;
        const char *Zprime10000X150save = Zprime10000X150name.c_str();
        
        std::string Zprime10X500name = "Zprime10X500" + syst;
        const char *Zprime10X500save = Zprime10X500name.c_str();
        std::string Zprime500X500name = "Zprime500X500" + syst;
        const char *Zprime500X500save = Zprime500X500name.c_str();
        std::string Zprime995X500name = "Zprime995X500" + syst;
        const char *Zprime995X500save = Zprime995X500name.c_str();
        std::string Zprime2000X500name = "Zprime2000X500" + syst;
        const char *Zprime2000X500save = Zprime2000X500name.c_str();
        std::string Zprime10000X500name = "Zprime10000X500" + syst;
        const char *Zprime10000X500save = Zprime10000X500name.c_str();
        
        std::string Zprime10X1000name = "Zprime10X1000" + syst;
        const char *Zprime10X1000save = Zprime10X1000name.c_str();
        std::string Zprime1000X1000name = "Zprime1000X1000" + syst;
        const char *Zprime1000X1000save = Zprime1000X1000name.c_str();
        std::string Zprime1995X1000name = "Zprime1995X1000" + syst;
        const char *Zprime1995X1000save = Zprime1995X1000name.c_str();
        std::string Zprime10000X1000name = "Zprime10000X1000" + syst;
        const char *Zprime10000X1000save = Zprime10000X1000name.c_str();

        data_->Write(DATAsave);
        
        vvt_->Write(VVTsave);
        vvj_->Write(VVJsave);
        zvv_->Write(ZVVsave);
        ewk_->Write(EWKsave);
        ttt_->Write(TTTsave);
        ttj_->Write(TTJsave);
        ztt_->Write(ZTTsave);
        zl_->Write(ZLsave);
        zj_->Write(ZJsave);
        w_->Write(Wsave);
        qcd_->Write(QCDsave);
        
        zhtautau_->Write(ZHsave);
        smh_->Write(SMHsave);
        
		comb_->Write(combsave);

        if (syst == "600"){zp600a300_->Write(ZprimeA300Zsave);}
        else if (syst == "800"){zp800a300_->Write(ZprimeA300Zsave);}
        else if (syst == "1000"){zp1000a300_->Write(ZprimeA300Zsave);}
        else if (syst == "1200"){zp1200a300_->Write(ZprimeA300Zsave);}
        else
        {
            zp600a300_->Write(Zprime600A300save);
            zp600a400_->Write(Zprime600A400save);
            
            zp800a300_->Write(Zprime800A300save);
            zp800a400_->Write(Zprime800A400save);
            zp800a500_->Write(Zprime800A500save);
            zp800a600_->Write(Zprime800A600save);
            
            zp1000a300_->Write(Zprime1000A300save);
            zp1000a400_->Write(Zprime1000A400save);
            zp1000a500_->Write(Zprime1000A500save);
            zp1000a600_->Write(Zprime1000A600save);
            zp1000a700_->Write(Zprime1000A700save);
            zp1000a800_->Write(Zprime1000A800save);
            
            zp1200a300_->Write(Zprime1200A300save);
            zp1200a400_->Write(Zprime1200A400save);
            zp1200a500_->Write(Zprime1200A500save);
            zp1200a600_->Write(Zprime1200A600save);
            zp1200a700_->Write(Zprime1200A700save);
            zp1200a800_->Write(Zprime1200A800save);
            
            zp1400a300_->Write(Zprime1400A300save);
            zp1400a400_->Write(Zprime1400A400save);
            zp1400a500_->Write(Zprime1400A500save);
            zp1400a600_->Write(Zprime1400A600save);
            zp1400a700_->Write(Zprime1400A700save);
            zp1400a800_->Write(Zprime1400A800save);
            
            zp1700a300_->Write(Zprime1700A300save);
            zp1700a400_->Write(Zprime1700A400save);
            zp1700a500_->Write(Zprime1700A500save);
            zp1700a600_->Write(Zprime1700A600save);
            zp1700a700_->Write(Zprime1700A700save);
            zp1700a800_->Write(Zprime1700A800save);
            
            zp2000a300_->Write(Zprime2000A300save);
            zp2000a400_->Write(Zprime2000A400save);
            zp2000a500_->Write(Zprime2000A500save);
            zp2000a600_->Write(Zprime2000A600save);
            zp2000a700_->Write(Zprime2000A700save);
            zp2000a800_->Write(Zprime2000A800save);
            
            zp2500a300_->Write(Zprime2500A300save);
            zp2500a400_->Write(Zprime2500A400save);
            zp2500a500_->Write(Zprime2500A500save);
            zp2500a600_->Write(Zprime2500A600save);
            zp2500a700_->Write(Zprime2500A700save);
            zp2500a800_->Write(Zprime2500A800save);
            
            zp10x1_->Write(Zprime10X1save);
            zp20x1_->Write(Zprime20X1save);
            zp50x1_->Write(Zprime50X1save);
            zp100x1_->Write(Zprime100X1save);
            zp300x1_->Write(Zprime300X1save);
            zp500x1_->Write(Zprime500X1save);
            zp1000x1_->Write(Zprime1000X1save);
            zp2000x1_->Write(Zprime2000X1save);
            zp10000x1_->Write(Zprime10000X1save);
            
            zp10x50_->Write(Zprime10X50save);
            zp50x50_->Write(Zprime50X50save);
            zp95x50_->Write(Zprime95X50save);
            zp200x50_->Write(Zprime200X50save);
            zp300x50_->Write(Zprime300X50save);
            zp10000x50_->Write(Zprime10000X50save);
        
            zp10x150_->Write(Zprime10X150save);
            zp200x150_->Write(Zprime200X150save);
            zp295x150_->Write(Zprime295X150save);
            zp500x150_->Write(Zprime500X150save);
            zp1000x150_->Write(Zprime1000X150save);
            zp10000x150_->Write(Zprime10000X150save);
            
            zp10x500_->Write(Zprime10X500save);
            zp500x500_->Write(Zprime500X500save);
            zp995x500_->Write(Zprime995X500save);
            zp2000x500_->Write(Zprime2000X500save);
            zp10000x500_->Write(Zprime10000X500save);
            
            zp10x1000_->Write(Zprime10X1000save);
            zp1000x1000_->Write(Zprime1000X1000save);
            zp1995x1000_->Write(Zprime1995X1000save);
            zp10000x1000_->Write(Zprime10000X1000save);

        }

        sbOp.cd();
		sbOp.Close();
	}
    std::cout << "FINISHED editing shape file" << std::endl;
    
    THStack * sigMTstack  = new THStack();
    
    if (plot)
    {
        if(parameter=="met")
        {
            for(int i=0; i<data_->GetNbinsX()+1; ++i)
            {
                if(data_->GetBinLowEdge(i)>= met_blind) 
                {
                    data_->SetBinContent(i,0);
                    data_->SetBinError(i,0);
                }
            }
        }
        else if(parameter.substr(0,6)=="mvaVar")
        {
            for(int i=0; i<data_->GetNbinsX()+1; ++i)
            {
                if(data_->GetBinLowEdge(i)>= 0.5)
                {
                    data_->SetBinContent(i,0);
                    data_->SetBinError(i,0);
                }
            }
        }
        if(parameter=="mt_tot")
        {
            for(int i=0; i<data_->GetNbinsX()+1; ++i)
            {
                if(data_->GetBinLowEdge(i) >= 200)
                {
                    data_->SetBinContent(i,0);
                    data_->SetBinError(i,0);
                }
            }
        }
    
        sigMTstack->Add(tt_);
        sigMTstack->Add(zvv_);
        sigMTstack->Add(ewk_);
        sigMTstack->Add(qcd_);
        if(splitPlotSMH)
        {
            sigMTstack->Add(zhtautau_);
            sigMTstack->Add(vbfhtautau_);
            sigMTstack->Add(glugluhtautau_);
        }
        else
        {
            sigMTstack->Add(smh_);
        }
        sigMTstack->Add(vvt_);
        sigMTstack->Add(vvj_);
        sigMTstack->Add(w_);
        sigMTstack->Add(ztt_);
        sigMTstack->Add(zj_);
        sigMTstack->Add(zl_);
        
        Signal_region->cd(1);

        data_->SetMaximum(max(DATA_MAX*data_->GetMaximum(),30.0));
        data_->SetMinimum(DATA_MIN);

    //	data_->SetTitle((can_name+";XXX;Events/Bin").c_str());
        data_->SetTitle(";XXX;Events/Bin");
        data_->GetXaxis()->SetTitle(x_axisLabels[parameter].c_str());

        data_->DrawCopy("PE");
        sigMTstack->Draw("histsames");

        if(parameter.substr(0,6)=="mvaVar")
        {
            if (syst=="600"){zp600a300_->SetLineStyle(3); zp600a300_->SetLineWidth(3); zp600a300_->SetLineColor(50); std::cout << "DRAWING" << std::endl; zp600a300_->DrawCopy("hist sames");}
            else if (syst=="800"){zp800a300_->SetLineStyle(3); zp800a300_->SetLineWidth(3); zp800a300_->SetLineColor(60); std::cout << "DRAWING" << std::endl; zp800a300_->DrawCopy("hist sames");}
            else if (syst=="1000"){zp1000a300_->SetLineStyle(3); zp1000a300_->SetLineWidth(3); zp1000a300_->SetLineColor(70); std::cout << "DRAWING" << std::endl; zp1000a300_->DrawCopy("hist sames");}
            else if (syst=="1200"){zp1200a300_->SetLineStyle(3); zp1200a300_->SetLineWidth(3); zp1200a300_->SetLineColor(80); std::cout << "DRAWING" << std::endl; zp1200a300_->DrawCopy("hist sames");}
            else
            {
                zp600a300_->SetLineStyle(2); zp600a300_->SetLineWidth(3); zp600a300_->SetLineColor(50); std::cout << "DRAWING" << std::endl; zp600a300_->DrawCopy("hist sames");
                zp800a300_->SetLineStyle(2); zp800a300_->SetLineWidth(3); zp800a300_->SetLineColor(60); std::cout << "DRAWING" << std::endl; zp800a300_->DrawCopy("hist sames");
                zp1000a300_->SetLineStyle(2); zp1000a300_->SetLineWidth(3); zp1000a300_->SetLineColor(70); std::cout << "DRAWING" << std::endl; zp1000a300_->DrawCopy("hist sames");
                zp1200a300_->SetLineStyle(2); zp1200a300_->SetLineWidth(3); zp1200a300_->SetLineColor(80); std::cout << "DRAWING" << std::endl; zp1200a300_->DrawCopy("hist sames");
            }
        }
        else if (parameter == "mt_tot" || parameter == "met")
        {
            zp600a300_->SetLineStyle(2); zp600a300_->SetLineWidth(3); zp600a300_->SetLineColor(50); std::cout << "DRAWING" << std::endl; zp600a300_->DrawCopy("hist sames");
            zp800a300_->SetLineStyle(2); zp800a300_->SetLineWidth(3); zp800a300_->SetLineColor(60); std::cout << "DRAWING" << std::endl; zp800a300_->DrawCopy("hist sames");
            zp1000a300_->SetLineStyle(2); zp1000a300_->SetLineWidth(3); zp1000a300_->SetLineColor(70); std::cout << "DRAWING" << std::endl; zp1000a300_->DrawCopy("hist sames");
            zp1200a300_->SetLineStyle(2); zp1200a300_->SetLineWidth(3); zp1200a300_->SetLineColor(80); std::cout << "DRAWING" << std::endl; zp1200a300_->DrawCopy("hist sames");
        }
        data_->DrawCopy("PEsames");

        /* update canvas */
        Signal_region->Update();

        /* small fix for when blinded data is empty */
        if(data_->GetSumOfWeights()==0)
        {
            data_->SetMaximum(100.0);
            sigMTstack->Draw("histsames");
            if(parameter=="met" || parameter=="mt_tot")
            {
            
                    zp600a300_->SetLineStyle(2); zp600a300_->SetLineWidth(3); zp600a300_->SetLineColor(50); std::cout << "DRAWING" << std::endl; zp600a300_->DrawCopy("hist sames");
                    zp800a300_->SetLineStyle(2); zp800a300_->SetLineWidth(3); zp800a300_->SetLineColor(60); std::cout << "DRAWING" << std::endl; zp800a300_->DrawCopy("hist sames");
                    zp1000a300_->SetLineStyle(2); zp1000a300_->SetLineWidth(3); zp1000a300_->SetLineColor(70); std::cout << "DRAWING" << std::endl; zp1000a300_->DrawCopy("hist sames");
                    zp1200a300_->SetLineStyle(2); zp1200a300_->SetLineWidth(3); zp1200a300_->SetLineColor(80); std::cout << "DRAWING" << std::endl; zp1200a300_->DrawCopy("hist sames");

            }
            else
            {
                if (syst=="600"){zp600a300_->SetLineStyle(3); zp600a300_->SetLineWidth(3); zp600a300_->SetLineColor(50); std::cout << "DRAWING" << std::endl; zp600a300_->DrawCopy("hist sames");}
                else if (syst=="800"){zp800a300_->SetLineStyle(3); zp800a300_->SetLineWidth(5); zp800a300_->SetLineColor(60); std::cout << "DRAWING" << std::endl; zp800a300_->DrawCopy("hist sames");}
                else if (syst=="1000"){zp1000a300_->SetLineStyle(3); zp1000a300_->SetLineWidth(5); zp1000a300_->SetLineColor(70); std::cout << "DRAWING" << std::endl; zp1000a300_->DrawCopy("hist sames");}
                else if (syst=="1200"){zp1200a300_->SetLineStyle(3); zp1200a300_->SetLineWidth(5); zp1200a300_->SetLineColor(80); std::cout << "DRAWING" << std::endl; zp1200a300_->DrawCopy("hist sames");}
                else if (parameter.substr(0,6)=="mvaVar")
                {
                    zp600a300_->SetLineStyle(3); zp600a300_->SetLineWidth(3); zp600a300_->SetLineColor(50); std::cout << "DRAWING" << std::endl; zp600a300_->DrawCopy("hist sames");
                    zp800a300_->SetLineStyle(3); zp800a300_->SetLineWidth(3); zp800a300_->SetLineColor(60); std::cout << "DRAWING" << std::endl; zp800a300_->DrawCopy("hist sames");
                    zp1000a300_->SetLineStyle(3); zp1000a300_->SetLineWidth(3); zp1000a300_->SetLineColor(70); std::cout << "DRAWING" << std::endl; zp1000a300_->DrawCopy("hist sames");
                    zp1200a300_->SetLineStyle(3); zp1200a300_->SetLineWidth(3); zp1200a300_->SetLineColor(80); std::cout << "DRAWING" << std::endl; zp1200a300_->DrawCopy("hist sames");
                }
            }
            data_->DrawCopy("PEsames");
        }

        drawTcanvasAndLegendAndRatio(Signal_region, data_, comb_, 1);

        Signal_region->DrawClone();

        std::string localSaveName = saveName + parameter + "_" + chan + "_" + drCutNameMap[drCut] + "Region.png";
        Signal_region->cd()->SaveAs(localSaveName.c_str());

        Signal_region->Close();
    }
    
	/* delete pointers */
	delete data_;
	delete vvt_;
    delete vvj_;
    delete zvv_;
    delete ewk_;
	delete ttt_;
    delete ttj_;
    delete tt_;
	delete ztt_;
    delete zl_;
    delete zj_;
	delete w_;
    delete zhtautau_;
    delete glugluhtautau_;
    delete vbfhtautau_;
	delete comb_;
    delete smh_;
	delete qcd_;
	delete zp600a300_;
    delete zp600a400_;
	delete zp800a300_;
    delete zp800a400_;
    delete zp800a500_;
    delete zp800a600_;
	delete zp1000a300_;
    delete zp1000a400_;
    delete zp1000a500_;
    delete zp1000a600_;
    delete zp1000a700_;
    delete zp1000a800_;
	delete zp1200a300_;
    delete zp1200a400_;
    delete zp1200a500_;
    delete zp1200a600_;
    delete zp1200a700_;
    delete zp1200a800_;
    delete zp1400a300_;
    delete zp1400a400_;
    delete zp1400a500_;
    delete zp1400a600_;
    delete zp1400a700_;
    delete zp1400a800_;
    delete zp1700a300_;
    delete zp1700a400_;
    delete zp1700a500_;
    delete zp1700a600_;
    delete zp1700a700_;
    delete zp1700a800_;
    delete zp2000a300_;
    delete zp2000a400_;
    delete zp2000a500_;
    delete zp2000a600_;
    delete zp2000a700_;
    delete zp2000a800_;
    delete zp2500a300_;
    delete zp2500a400_;
    delete zp2500a500_;
    delete zp2500a600_;
    delete zp2500a700_;
    delete zp2500a800_;
    
    delete zp10x1_;
    delete zp20x1_;
    delete zp50x1_;
    delete zp100x1_;
    delete zp300x1_;
    delete zp500x1_;
    delete zp1000x1_;
    delete zp2000x1_;
    delete zp10000x1_;
    
    delete zp10x50_;
    delete zp50x50_;
    delete zp95x50_;
    delete zp200x50_;
    delete zp300x50_;
    delete zp10000x50_;

    delete zp10x150_;
    delete zp200x150_;
    delete zp295x150_;
    delete zp500x150_;
    delete zp1000x150_;
    delete zp10000x150_;
    
    delete zp10x500_;
    delete zp500x500_;
    delete zp995x500_;
    delete zp2000x500_;
    delete zp10000x500_;
    
    delete zp10x1000_;
    delete zp1000x1000_;
    delete zp1995x1000_;
    delete zp10000x1000_;

	delete Signal_region;
	delete sigMTstack;
}


void fillQCD_Shape(TCut SScut, TCut SScutRel, double WnormForQCD, TH1F * QCDshape, std::string parameter, float bin[3], std::string can_name, std::string syst, std::string chan, bool plots, bool createOutputShapes)
{

    /* create a TCanvas */
    TCanvas * QCD_SS = new TCanvas(can_name.c_str(),can_name.c_str(),canDim[0],canDim[1]);
    QCD_SS->Divide(1,2);
    	/* fill the hists */
	if (plots) {QCD_SS->cd(1);}
    
    //First create the shapes from SS + relaxed shape
    
	/* create some histograms and format them */
	TH1F * dataShape_ = new TH1F("dataShape_","dataShape_",bin[0],bin[1],bin[2]); 
	dataShape_->Sumw2(); 
	dataShape_->SetMarkerStyle(dataMarker);

	TH1F * qcdShape_ = new TH1F("qcdShape_","qcdShape_",bin[0],bin[1],bin[2]); 
	qcdShape_->Sumw2(); 
	qcdShape_->SetFillColor(colors["QCD"]);

	TH1F * zttShape_ = new TH1F("zttShape_","zttShape_",bin[0],bin[1],bin[2]);
	zttShape_->Sumw2();
	zttShape_->SetFillColor(colors["ZTT"]);
    
    TH1F * zlShape_ = new TH1F("zlShape_","zlShape_",bin[0],bin[1],bin[2]);
	zlShape_->Sumw2();
	zlShape_->SetFillColor(colors["ZL"]);
    
    TH1F * zjShape_ = new TH1F("zjShape_","zjShape_",bin[0],bin[1],bin[2]);
	zjShape_->Sumw2();
	zjShape_->SetFillColor(colors["ZJ"]);
	
	TH1F * wShape_ = new TH1F("wShape_","wShape_",bin[0],bin[1],bin[2]); 
	wShape_->Sumw2(); 
	wShape_->SetFillColor(colors["W"]);
	
	TH1F * tttShape_ = new TH1F("tttShape_","tttShape_",bin[0],bin[1],bin[2]);
	tttShape_->Sumw2();
	tttShape_->SetFillColor(colors["TTT"]);
    
    TH1F * ttjShape_ = new TH1F("ttjShape_","ttjShape_",bin[0],bin[1],bin[2]);
	ttjShape_->Sumw2();
	ttjShape_->SetFillColor(colors["TTJ"]);
	
	TH1F * vvtShape_ = new TH1F("vvtShape_","vvtShape_",bin[0],bin[1],bin[2]);
	vvtShape_->Sumw2();
	vvtShape_->SetFillColor(colors["VV"]);
    
    TH1F * vvjShape_ = new TH1F("vvjShape_","vvjShape_",bin[0],bin[1],bin[2]);
	vvjShape_->Sumw2();
	vvjShape_->SetFillColor(colors["VV"]);
    
    TH1F * zvvShape_ = new TH1F("zvvShape_","zvvShape_",bin[0],bin[1],bin[2]);
	zvvShape_->Sumw2();
	zvvShape_->SetFillColor(colors["ZVV"]);
    
    TH1F * ewkShape_ = new TH1F("ewkShape_","ewkShape_",bin[0],bin[1],bin[2]);
	ewkShape_->Sumw2();
	ewkShape_->SetFillColor(colors["EWK"]);
	
	TH1F * zhtautauShape_ = new TH1F("zhtautauShape_","zhtautauShape_",bin[0],bin[1],bin[2]);
	zhtautauShape_->Sumw2();
	zhtautauShape_->SetFillColor(colors["ZHTauTau"]);

	TH1F * glugluhtautauShape_ = new TH1F("glugluhtautauShape_","glugluhtautauShape_",bin[0],bin[1],bin[2]);
	glugluhtautauShape_->Sumw2();
	glugluhtautauShape_->SetFillColor(colors["GluGluHTauTau"]);
    
    TH1F * vbfhtautauShape_ = new TH1F("vbfhtautauShape_","vbfhtautauShape_",bin[0],bin[1],bin[2]);
	vbfhtautauShape_->Sumw2();
    vbfhtautauShape_->SetFillColor(colors["VBFHTauTau"]);
    
    TH1F * smhShape_ = new TH1F("smhShape_","smhShape_",bin[0],bin[1],bin[2]);
	smhShape_->Sumw2();
    smhShape_->SetFillColor(colors["SMH"]);

	/* fill the hists */
	//if (plots) {QCD_SS->cd(1);}

	DATA->Draw((parameter+">>dataShape_").c_str(),SScutRel);
	W->Draw((parameter+">>wShape_").c_str(),SScutRel*weights["W"].c_str());
	ZTT->Draw((parameter+">>zttShape_").c_str(),SScutRel*weights["ZTT"].c_str());
    ZL->Draw((parameter+">>zlShape_").c_str(),SScutRel*weights["ZL"].c_str());
    ZJ->Draw((parameter+">>zjShape_").c_str(),SScutRel*weights["ZJ"].c_str());
	VVT->Draw((parameter+">>vvtShape_").c_str(),SScutRel*weights["VV"].c_str());
    VVJ->Draw((parameter+">>vvjShape_").c_str(),SScutRel*weights["VV"].c_str());
    ZVV->Draw((parameter+">>zvvShape_").c_str(),SScutRel*weights["ZVV"].c_str());
    EWK->Draw((parameter+">>ewkShape_").c_str(),SScutRel*weights["EWK"].c_str());
	TTT->Draw((parameter+">>tttShape_").c_str(),SScutRel*weights["TTT"].c_str());
    TTJ->Draw((parameter+">>ttjShape_").c_str(),SScutRel*weights["TTJ"].c_str());
	ZHTauTau->Draw((parameter+">>zhtautauShape_").c_str(),SScutRel*weights["ZHTauTau"].c_str());
    GluGluHTauTau->Draw((parameter+">>glugluhtautauShape_").c_str(),SScutRel*weights["GluGluHTauTau"].c_str());
    VBFHTauTau->Draw((parameter+">>vbfhtautauShape_").c_str(),SScutRel*weights["VBFHTauTau"].c_str());

    smhShape_->Add(zhtautauShape_);
    smhShape_->Add(glugluhtautauShape_);
    smhShape_->Add(vbfhtautauShape_);
    
	/* QCD is data - all mc bks */
	wShape_->Scale(WnormForQCD);
	qcdShape_->Add(dataShape_);
    for(int i = 0; i<qcdShape_->GetNbinsX()+1; i++)
	{
		if(qcdShape_->GetBinContent(i)<0.) {qcdShape_->SetBinContent(i,0.);}
	}
	qcdShape_->Add(wShape_,-1);
    for(int i = 0; i<qcdShape_->GetNbinsX()+1; i++)
	{
		if(qcdShape_->GetBinContent(i)<0.) {qcdShape_->SetBinContent(i,0.);}
	}
	qcdShape_->Add(zttShape_,-1);
    for(int i = 0; i<qcdShape_->GetNbinsX()+1; i++)
	{
		if(qcdShape_->GetBinContent(i)<0.) {qcdShape_->SetBinContent(i,0.);}
	}
    qcdShape_->Add(zlShape_,-1);
    for(int i = 0; i<qcdShape_->GetNbinsX()+1; i++)
	{
		if(qcdShape_->GetBinContent(i)<0.) {qcdShape_->SetBinContent(i,0.);}
	}
    qcdShape_->Add(zjShape_,-1);
    for(int i = 0; i<qcdShape_->GetNbinsX()+1; i++)
	{
		if(qcdShape_->GetBinContent(i)<0.) {qcdShape_->SetBinContent(i,0.);}
	}
	qcdShape_->Add(tttShape_,-1);
    for(int i = 0; i<qcdShape_->GetNbinsX()+1; i++)
	{
		if(qcdShape_->GetBinContent(i)<0.) {qcdShape_->SetBinContent(i,0.);}
	}
    qcdShape_->Add(ttjShape_,-1);
    for(int i = 0; i<qcdShape_->GetNbinsX()+1; i++)
	{
		if(qcdShape_->GetBinContent(i)<0.) {qcdShape_->SetBinContent(i,0.);}
	}
	qcdShape_->Add(vvtShape_,-1);
    for(int i = 0; i<qcdShape_->GetNbinsX()+1; i++)
	{
		if(qcdShape_->GetBinContent(i)<0.) {qcdShape_->SetBinContent(i,0.);}
	}
    qcdShape_->Add(vvjShape_,-1);
    for(int i = 0; i<qcdShape_->GetNbinsX()+1; i++)
	{
		if(qcdShape_->GetBinContent(i)<0.) {qcdShape_->SetBinContent(i,0.);}
	}
    qcdShape_->Add(zvvShape_,-1);
    for(int i = 0; i<qcdShape_->GetNbinsX()+1; i++)
	{
		if(qcdShape_->GetBinContent(i)<0.) {qcdShape_->SetBinContent(i,0.);}
	}
    qcdShape_->Add(ewkShape_,-1);
    for(int i = 0; i<qcdShape_->GetNbinsX()+1; i++)
	{
		if(qcdShape_->GetBinContent(i)<0.) {qcdShape_->SetBinContent(i,0.);}
	}
	qcdShape_->Add(smhShape_,-1);
    for(int i = 0; i<qcdShape_->GetNbinsX()+1; i++)
	{
		if(qcdShape_->GetBinContent(i)<0.) {qcdShape_->SetBinContent(i,0.);}
	}
    
    delete dataShape_;
	delete vvtShape_;
    delete vvjShape_;
    delete zvvShape_;
    delete ewkShape_;
	delete tttShape_;
    delete ttjShape_;
	delete zttShape_;
    delete zlShape_;
    delete zjShape_;
	delete wShape_;
    delete zhtautauShape_;
    delete glugluhtautauShape_;
    delete vbfhtautauShape_;
    delete smhShape_;
    
    //now create the normalization for regular selection
	/* create some histograms and format them */
	TH1F * dataY_ = new TH1F("dataY_","dataY_",bin[0],bin[1],bin[2]); 
	dataY_->Sumw2(); 
	dataY_->SetMarkerStyle(dataMarker);

	TH1F * zttY_ = new TH1F("zttY_","zttY_",bin[0],bin[1],bin[2]);
	zttY_->Sumw2();
	zttY_->SetFillColor(colors["ZTT"]);
    
    TH1F * zlY_ = new TH1F("zlY_","zlY_",bin[0],bin[1],bin[2]);
	zlY_->Sumw2();
	zlY_->SetFillColor(colors["ZL"]);
    
    TH1F * zjY_ = new TH1F("zjY_","zjY_",bin[0],bin[1],bin[2]);
	zjY_->Sumw2();
	zjY_->SetFillColor(colors["ZJ"]);
	
	TH1F * wY_ = new TH1F("wY_","wY_",bin[0],bin[1],bin[2]); 
	wY_->Sumw2(); 
	wY_->SetFillColor(colors["W"]);
	
	TH1F * tttY_ = new TH1F("tttY_","tttY_",bin[0],bin[1],bin[2]);
	tttY_->Sumw2();
	tttY_->SetFillColor(colors["TTT"]);
    
    TH1F * ttjY_ = new TH1F("ttjY_","ttjY_",bin[0],bin[1],bin[2]);
	ttjY_->Sumw2();
	ttjY_->SetFillColor(colors["TTJ"]);
	
	TH1F * vvtY_ = new TH1F("vvtY_","vvtY_",bin[0],bin[1],bin[2]);
	vvtY_->Sumw2();
	vvtY_->SetFillColor(colors["VV"]);
    
    TH1F * vvjY_ = new TH1F("vvjY_","vvjY_",bin[0],bin[1],bin[2]);
	vvjY_->Sumw2();
	vvjY_->SetFillColor(colors["VV"]);
    
    TH1F * zvvY_ = new TH1F("zvvY_","zvvY_",bin[0],bin[1],bin[2]);
	zvvY_->Sumw2();
	zvvY_->SetFillColor(colors["ZVV"]);
    
    TH1F * ewkY_ = new TH1F("ewkY_","ewkY_",bin[0],bin[1],bin[2]);
	ewkY_->Sumw2();
	ewkY_->SetFillColor(colors["EWK"]);
	
	TH1F * zhtautauY_ = new TH1F("zhtautauY_","zhtautauY_",bin[0],bin[1],bin[2]);
	zhtautauY_->Sumw2();
	zhtautauY_->SetFillColor(colors["ZHTauTau"]);

	TH1F * glugluhtautauY_ = new TH1F("glugluhtautauY_","glugluhtautauY_",bin[0],bin[1],bin[2]);
	glugluhtautauY_->Sumw2();
	glugluhtautauY_->SetFillColor(colors["GluGluHTauTau"]);
    
    TH1F * vbfhtautauY_ = new TH1F("vbfhtautauY_","vbfhtautauY_",bin[0],bin[1],bin[2]);
	vbfhtautauY_->Sumw2();
    vbfhtautauY_->SetFillColor(colors["VBFHTauTau"]);
    
    TH1F * smhY_ = new TH1F("smhY_","smhY_",bin[0],bin[1],bin[2]);
	smhY_->Sumw2();
    smhY_->SetFillColor(colors["SMH"]);

	TH1F * combY_ = new TH1F("combY_","combY_",bin[0],bin[1],bin[2]); 
	combY_->Sumw2();

	DATA->Draw((parameter+">>dataY_").c_str(),SScut);
	W->Draw((parameter+">>wY_").c_str(),SScutRel*weights["W"].c_str());
	ZTT->Draw((parameter+">>zttY_").c_str(),SScutRel*weights["ZTT"].c_str());
    ZL->Draw((parameter+">>zlY_").c_str(),SScutRel*weights["ZL"].c_str());
    ZJ->Draw((parameter+">>zjY_").c_str(),SScutRel*weights["ZJ"].c_str());
	VVT->Draw((parameter+">>vvtY_").c_str(),SScutRel*weights["VV"].c_str());
    VVJ->Draw((parameter+">>vvjY_").c_str(),SScutRel*weights["VV"].c_str());
    ZVV->Draw((parameter+">>zvvY_").c_str(),SScut*weights["ZVV"].c_str());
    EWK->Draw((parameter+">>ewkY_").c_str(),SScut*weights["EWK"].c_str());
	TTT->Draw((parameter+">>tttY_").c_str(),SScutRel*weights["TTT"].c_str());
    TTJ->Draw((parameter+">>ttjY_").c_str(),SScutRel*weights["TTJ"].c_str());
	ZHTauTau->Draw((parameter+">>zhtautauY_").c_str(),SScut*weights["ZHTauTau"].c_str());
    GluGluHTauTau->Draw((parameter+">>glugluhtautauY_").c_str(),SScut*weights["GluGluHTauTau"].c_str());
    VBFHTauTau->Draw((parameter+">>vbfhtautauY_").c_str(),SScut*weights["VBFHTauTau"].c_str());
    
    wY_->Scale(shapeScale(SScut, SScutRel, W, weights["W"], parameter, bin));
    vvtY_->Scale(shapeScale(SScut, SScutRel, VVT, weights["VV"], parameter, bin));
    vvjY_->Scale(shapeScale(SScut, SScutRel, VVJ, weights["VV"], parameter, bin));
    tttY_->Scale(shapeScale(SScut, SScutRel, TTT, weights["TTT"], parameter, bin));
    ttjY_->Scale(shapeScale(SScut, SScutRel, TTJ, weights["TTJ"], parameter, bin));
    zttY_->Scale(shapeScale(SScut, SScutRel, ZTT, weights["ZTT"], parameter, bin));
    zlY_->Scale(shapeScale(SScut, SScutRel, ZL, weights["ZL"], parameter, bin));
    zjY_->Scale(shapeScale(SScut, SScutRel, ZJ, weights["ZJ"], parameter, bin));
    
    smhY_->Add(zhtautauY_);
    smhY_->Add(glugluhtautauY_);
    smhY_->Add(vbfhtautauY_);

	/* QCD is data - all mc bks */
    
    double qcdYield;
    qcdYield = TMath::Nint(dataY_->Integral() - wY_->Integral() - zttY_->Integral() - zlY_->Integral() - zjY_->Integral() - tttY_->Integral() - ttjY_->Integral() - vvtY_->Integral() - vvjY_->Integral() - zvvY_->Integral() - ewkY_->Integral() - smhY_->Integral());
    
    std::cout << "QCD: " << qcdYield << std::endl;
    std::cout << "QCD Shape: " << qcdShape_->Integral() << std::endl;
    
    double corrFactor;
    if (qcdShape_->Integral() == 0.) {corrFactor = 1.0;}
    else {corrFactor = qcdYield/qcdShape_->Integral();}
    std::cout << "correction factor for QCD relaxed to full yield: " << corrFactor << std::endl;
    qcdShape_->Scale(corrFactor);
    std::cout << "QCD after scaling: " << qcdShape_->Integral() << std::endl;
    
    /* eliminate negative bins */
    for(int i = 0; i<qcdShape_->GetNbinsX()+1; i++)
	{
		if(qcdShape_->GetBinContent(i)<0.) {qcdShape_->SetBinContent(i,0.);}
	}
    std::cout << "QCD after rm negative bins: " << qcdShape_->Integral() << std::endl;
    
    if(createOutputShapes)
	{
		/* create a file to help with S/sqrt(S+B) cut opt */
        
        std::string paramName;
        if (parameter=="met" || parameter=="mvamet") {paramName = "met";}
        else if (parameter.substr(0,6)=="mvaVar") {paramName = "mva";}
        else if (parameter.substr(0,6) == "mt_tot") {paramName="mt";}
        else {paramName = parameter;}
        
        std::string catName = chan + "_QCD_inclusive_cr";
		std::string optName = "xtt_" + chan + ".inputs-13TeV-" + paramName + ".root";
		TFile sbOp(optName.c_str(),"UPDATE");
        
        TDirectory * catdir = sbOp.GetDirectory(catName.c_str());
        if (catdir){catdir->cd();}
        else {sbOp.mkdir(catName.c_str()); TDirectory * catdir = sbOp.GetDirectory(catName.c_str()); catdir->cd();}
        
        std::string DATAname = "data_obs" + syst;
        const char *DATAsave = DATAname.c_str();
        std::string VVTname = "VVT" + syst;
        const char *VVTsave = VVTname.c_str();
        std::string VVJname = "VVJ" + syst;
        const char *VVJsave = VVJname.c_str();
        std::string ZVVname = "ZVV" + syst;
        const char *ZVVsave = ZVVname.c_str();
        std::string EWKname = "EWK" + syst;
        const char *EWKsave = EWKname.c_str();
        std::string TTTname = "TTT" + syst;
        const char *TTTsave = TTTname.c_str();
        std::string TTJname = "TTJ" + syst;
        const char *TTJsave = TTJname.c_str();
        std::string ZTTname = "ZTT" + syst;
        const char *ZTTsave = ZTTname.c_str();
        std::string ZLname = "ZL" + syst;
        const char *ZLsave = ZLname.c_str();
        std::string ZJname = "ZJ" + syst;
        const char *ZJsave = ZJname.c_str();
        std::string Wname = "W" + syst;
        const char *Wsave = Wname.c_str();
        std::string QCDname = "QCD" + syst;
        const char *QCDsave = QCDname.c_str();
        std::string SMHname = "SMH" + syst;
        const char *SMHsave = SMHname.c_str();
        std::string ZHname = "ZH" + syst;
        const char *ZHsave = ZHname.c_str();
        
        std::string combname = "comb" + syst;
        const char *combsave = combname.c_str();

        dataY_->Write(DATAsave);
        
        vvtY_->Write(VVTsave);
        vvjY_->Write(VVJsave);
        zvvY_->Write(ZVVsave);
        ewkY_->Write(EWKsave);
        tttY_->Write(TTTsave);
        ttjY_->Write(TTJsave);
        zttY_->Write(ZTTsave);
        zlY_->Write(ZLsave);
        zjY_->Write(ZJsave);
        wY_->Write(Wsave);
        qcdShape_->Write(QCDsave);
        
        zhtautauY_->Write(ZHsave);
        smhY_->Write(SMHsave);
        
		combY_->Write(combsave);

        sbOp.cd();
		sbOp.Close();
	}
    std::cout << "FINISHED editing shape file" << std::endl;

    THStack * QCDStack  = new THStack();

    if (plots)
    {

        QCDStack->Add(tttY_);
        QCDStack->Add(ttjY_);
        QCDStack->Add(zvvY_);
        QCDStack->Add(ewkY_);
        QCDStack->Add(qcdShape_);
        QCDStack->Add(zhtautauY_);
        QCDStack->Add(vbfhtautauY_);
        QCDStack->Add(glugluhtautauY_);
        QCDStack->Add(vvtY_);
        QCDStack->Add(vvjY_);
        QCDStack->Add(wY_);
        QCDStack->Add(zttY_);
        QCDStack->Add(zlY_);
        QCDStack->Add(zjY_);

        QCD_SS->cd(1);dataY_->SetMaximum(DATA_MAX*dataY_->GetMaximum());dataY_->SetTitle(";XXX;Events/Bin");dataY_->GetXaxis()->SetTitle(x_axisLabels[parameter].c_str());dataY_->DrawCopy("PE");QCDStack->Draw("histsames");dataY_->DrawCopy("PEsames");QCD_SS->Update();

        combY_->Add(vvtY_);
        combY_->Add(vvjY_);
        combY_->Add(zvvY_);
        combY_->Add(ewkY_);
        combY_->Add(tttY_);
        combY_->Add(ttjY_);
        combY_->Add(wY_);
        combY_->Add(zttY_);
        combY_->Add(zlY_);
        combY_->Add(zjY_);
        combY_->Add(smhY_);
        combY_->Add(qcdShape_);

        drawTcanvasAndLegendAndRatio(QCD_SS, dataY_, combY_, 0);
        QCD_SS->DrawClone();std::string localSaveName = "QCD_"+saveName + parameter + "_" + chan + "_" + drCutNameMap[drCut] + "Region.png";QCD_SS->cd()->SaveAs(localSaveName.c_str());QCD_SS->Close();
    }
    
	//std::cout<<" ***** "<<qcd_->Integral()<<"\n";

	QCDshape->Add(qcdShape_);

	/* delete pointers */
	delete dataY_;
	delete vvtY_;
    delete vvjY_;
    delete zvvY_;
    delete ewkY_;
	delete tttY_;
    delete ttjY_;
	delete zttY_;
    delete zlY_;
    delete zjY_;
	delete wY_;
    delete zhtautauY_;
    delete glugluhtautauY_;
    delete vbfhtautauY_;
    delete smhY_;
	delete combY_;
    delete qcdShape_;
	delete QCDStack;
	delete QCD_SS;

}

void fillQCD_ShapeVarBin(TCut SScut, TCut SScutRel, double WnormForQCD, TCut SSnormCut, bool doNorm, TH1F * QCDshape, std::string parameter, int binnum, float bin[], std::string can_name,  std::string syst, std::string chan, bool plots, bool createOutputShapes)
{
    
    std::cout << "SS Full CUT: " << SScut.GetTitle() << std::endl;
    
    std::cout << "SS Relaxed CUT: " << SScutRel.GetTitle() << std::endl;

    /* create a TCanvas */
    TCanvas * QCD_SS = new TCanvas(can_name.c_str(),can_name.c_str(),canDim[0],canDim[1]);
    QCD_SS->Divide(1,2);
    
    /* fill the hists */
	if (plots) {QCD_SS->cd(1);}

    //First create the shapes from SS + relaxed shape
    
	/* create some histograms and format them */
	
    TH1F * dataShape_ = new TH1F("dataShape_","dataShape_",binnum,bin);
	dataShape_->Sumw2(); 
	dataShape_->SetMarkerStyle(dataMarker);

	TH1F * qcdShape_ = new TH1F("qcdShape_","qcdShape_",binnum,bin); 
	qcdShape_->Sumw2(); 
	qcdShape_->SetFillColor(colors["QCD"]);

	TH1F * zttShape_ = new TH1F("zttShape_","zttShape_",binnum,bin);
	zttShape_->Sumw2();
	zttShape_->SetFillColor(colors["ZTT"]);
    
    TH1F * zlShape_ = new TH1F("zlShape_","zlShape_",binnum,bin);
	zlShape_->Sumw2();
	zlShape_->SetFillColor(colors["ZL"]);
    
    TH1F * zjShape_ = new TH1F("zjShape_","zjShape_",binnum,bin);
	zjShape_->Sumw2();
	zjShape_->SetFillColor(colors["ZJ"]);
	
	TH1F * wShape_ = new TH1F("wShape_","wShape_",binnum,bin); 
	wShape_->Sumw2(); 
	wShape_->SetFillColor(colors["W"]);
	
	TH1F * tttShape_ = new TH1F("tttShape_","tttShape_",binnum,bin);
	tttShape_->Sumw2();
	tttShape_->SetFillColor(colors["TTT"]);
    
    TH1F * ttjShape_ = new TH1F("ttjShape_","ttjShape_",binnum,bin);
	ttjShape_->Sumw2();
	ttjShape_->SetFillColor(colors["TTJ"]);
	
	TH1F * vvtShape_ = new TH1F("vvtShape_","vvtShape_",binnum,bin);
	vvtShape_->Sumw2();
	vvtShape_->SetFillColor(colors["VV"]);
    
    TH1F * vvjShape_ = new TH1F("vvjShape_","vvjShape_",binnum,bin);
	vvjShape_->Sumw2();
	vvjShape_->SetFillColor(colors["VV"]);
    
    TH1F * zvvShape_ = new TH1F("zvvShape_","zvvShape_",binnum,bin);
	zvvShape_->Sumw2();
	zvvShape_->SetFillColor(colors["ZVV"]);
    
    TH1F * ewkShape_ = new TH1F("ewkShape_","ewkShape_",binnum,bin);
	ewkShape_->Sumw2();
	ewkShape_->SetFillColor(colors["EWK"]);
	
	TH1F * zhtautauShape_ = new TH1F("zhtautauShape_","zhtautauShape_",binnum,bin);
	zhtautauShape_->Sumw2();
	zhtautauShape_->SetFillColor(colors["ZHTauTau"]);

	TH1F * glugluhtautauShape_ = new TH1F("glugluhtautauShape_","glugluhtautauShape_",binnum,bin);
	glugluhtautauShape_->Sumw2();
	glugluhtautauShape_->SetFillColor(colors["GluGluHTauTau"]);
    
    TH1F * vbfhtautauShape_ = new TH1F("vbfhtautauShape_","vbfhtautauShape_",binnum,bin);
	vbfhtautauShape_->Sumw2();
    vbfhtautauShape_->SetFillColor(colors["VBFHTauTau"]);
    
    TH1F * preQCDShape_ = new TH1F("preQCDShape_","preQCDShape_",binnum,bin);
    
	/* fill the hists */
	//if (plots) {QCD_SS->cd(1);}
	DATA->Draw((parameter+">>dataShape_").c_str(),SScutRel);
	W->Draw((parameter+">>wShape_").c_str(),SScutRel*weights["W"].c_str());
	ZTT->Draw((parameter+">>zttShape_").c_str(),SScutRel*weights["ZTT"].c_str());
    ZL->Draw((parameter+">>zlShape_").c_str(),SScutRel*weights["ZL"].c_str());
    ZJ->Draw((parameter+">>zjShape_").c_str(),SScutRel*weights["ZJ"].c_str());
	VVT->Draw((parameter+">>vvtShape_").c_str(),SScutRel*weights["VV"].c_str());
    VVJ->Draw((parameter+">>vvjShape_").c_str(),SScutRel*weights["VV"].c_str());
    ZVV->Draw((parameter+">>zvvShape_").c_str(),SScutRel*weights["ZVV"].c_str());
    EWK->Draw((parameter+">>ewkShape_").c_str(),SScutRel*weights["EWK"].c_str());
	TTT->Draw((parameter+">>tttShape_").c_str(),SScutRel*weights["TTT"].c_str());
    TTJ->Draw((parameter+">>ttjShape_").c_str(),SScutRel*weights["TTJ"].c_str());
	ZHTauTau->Draw((parameter+">>zhtautauShape_").c_str(),SScutRel*weights["ZHTauTau"].c_str());
    GluGluHTauTau->Draw((parameter+">>glugluhtautauShape_").c_str(),SScutRel*weights["GluGluHTauTau"].c_str());
    VBFHTauTau->Draw((parameter+">>vbfhtautauShape_").c_str(),SScutRel*weights["VBFHTauTau"].c_str());

    std::cout << "SS relaxed" << std::endl;
    std::cout << "data: " << dataShape_->Integral() << std::endl;
    std::cout << "W: " << wShape_->Integral() << std::endl;
    std::cout << "ZTT: " << zttShape_->Integral() << std::endl;
    std::cout << "ZL: " << zlShape_->Integral() << std::endl;
    std::cout << "ZJ: " << zjShape_->Integral() << std::endl;
    std::cout << "VVT: " << vvtShape_->Integral() << std::endl;
    std::cout << "VVJ: " << vvjShape_->Integral() << std::endl;
    std::cout << "ZVV: " << zvvShape_->Integral() << std::endl;
    std::cout << "EWK: " << ewkShape_->Integral() << std::endl;
    std::cout << "TTT: " << tttShape_->Integral() << std::endl;
    std::cout << "TTJ: " << ttjShape_->Integral() << std::endl;
    std::cout << "ZH: " << zhtautauShape_->Integral() << std::endl;
    std::cout << "VBF H: " << glugluhtautauShape_->Integral() << std::endl;
    std::cout << "GG H: " << vbfhtautauShape_->Integral() << std::endl;

	/* QCD is data - all mc bks */

	wShape_->Scale(WnormForQCD);
	qcdShape_->Add(dataShape_);
    for(int i = 0; i<qcdShape_->GetNbinsX()+1; i++)
	{
		if(qcdShape_->GetBinContent(i)<0.) {qcdShape_->SetBinContent(i,0.);}
	}
	qcdShape_->Add(wShape_,-1);
    for(int i = 0; i<qcdShape_->GetNbinsX()+1; i++)
	{
		if(qcdShape_->GetBinContent(i)<0.) {qcdShape_->SetBinContent(i,0.);}
	}
	qcdShape_->Add(zttShape_,-1);
    for(int i = 0; i<qcdShape_->GetNbinsX()+1; i++)
	{
		if(qcdShape_->GetBinContent(i)<0.) {qcdShape_->SetBinContent(i,0.);}
	}
    qcdShape_->Add(zlShape_,-1);
    for(int i = 0; i<qcdShape_->GetNbinsX()+1; i++)
	{
		if(qcdShape_->GetBinContent(i)<0.) {qcdShape_->SetBinContent(i,0.);}
	}
    qcdShape_->Add(zjShape_,-1);
    for(int i = 0; i<qcdShape_->GetNbinsX()+1; i++)
	{
		if(qcdShape_->GetBinContent(i)<0.) {qcdShape_->SetBinContent(i,0.);}
	}
	qcdShape_->Add(tttShape_,-1);
    for(int i = 0; i<qcdShape_->GetNbinsX()+1; i++)
	{
		if(qcdShape_->GetBinContent(i)<0.) {qcdShape_->SetBinContent(i,0.);}
	}
    qcdShape_->Add(ttjShape_,-1);
    for(int i = 0; i<qcdShape_->GetNbinsX()+1; i++)
	{
		if(qcdShape_->GetBinContent(i)<0.) {qcdShape_->SetBinContent(i,0.);}
	}
	qcdShape_->Add(vvtShape_,-1);
    for(int i = 0; i<qcdShape_->GetNbinsX()+1; i++)
	{
		if(qcdShape_->GetBinContent(i)<0.) {qcdShape_->SetBinContent(i,0.);}
	}
    qcdShape_->Add(vvjShape_,-1);
    for(int i = 0; i<qcdShape_->GetNbinsX()+1; i++)
	{
		if(qcdShape_->GetBinContent(i)<0.) {qcdShape_->SetBinContent(i,0.);}
	}
    qcdShape_->Add(zvvShape_,-1);
    for(int i = 0; i<qcdShape_->GetNbinsX()+1; i++)
	{
		if(qcdShape_->GetBinContent(i)<0.) {qcdShape_->SetBinContent(i,0.);}
	}
    qcdShape_->Add(ewkShape_,-1);
    for(int i = 0; i<qcdShape_->GetNbinsX()+1; i++)
	{
		if(qcdShape_->GetBinContent(i)<0.) {qcdShape_->SetBinContent(i,0.);}
	}
	qcdShape_->Add(zhtautauShape_,-1);
    for(int i = 0; i<qcdShape_->GetNbinsX()+1; i++)
	{
		if(qcdShape_->GetBinContent(i)<0.) {qcdShape_->SetBinContent(i,0.);}
	}
    qcdShape_->Add(glugluhtautauShape_,-1);
    for(int i = 0; i<qcdShape_->GetNbinsX()+1; i++)
	{
		if(qcdShape_->GetBinContent(i)<0.) {qcdShape_->SetBinContent(i,0.);}
	}
    qcdShape_->Add(vbfhtautauShape_,-1);
    for(int i = 0; i<qcdShape_->GetNbinsX()+1; i++)
	{
		if(qcdShape_->GetBinContent(i)<0.) {qcdShape_->SetBinContent(i,0.);}
	}
    
    preQCDShape_->Add(qcdShape_);

    std::cout << "QCD Shape: " << qcdShape_->Integral() << std::endl;
    
	delete dataShape_;
	delete vvtShape_;
    delete vvjShape_;
    delete zvvShape_;
    delete ewkShape_;
	delete tttShape_;
    delete ttjShape_;
	delete zttShape_;
    delete zlShape_;
    delete zjShape_;
	delete wShape_;
    delete zhtautauShape_;
    delete glugluhtautauShape_;
    delete vbfhtautauShape_;

    //now create the normalization for regular selection
	/* create some histograms and format them */
    
	TH1F * dataY_ = new TH1F("dataY_","dataY_",binnum,bin); 
	dataY_->Sumw2(); 
	dataY_->SetMarkerStyle(dataMarker);

	TH1F * zttY_ = new TH1F("zttY_","zttY_",binnum,bin);
	zttY_->Sumw2();
	zttY_->SetFillColor(colors["ZTT"]);
    
    TH1F * zlY_ = new TH1F("zlY_","zlY_",binnum,bin);
	zlY_->Sumw2();
	zlY_->SetFillColor(colors["ZL"]);
    
    TH1F * zjY_ = new TH1F("zjY_","zjY_",binnum,bin);
	zjY_->Sumw2();
	zjY_->SetFillColor(colors["ZJ"]);
	
	TH1F * wY_ = new TH1F("wY_","wY_",binnum,bin); 
	wY_->Sumw2(); 
	wY_->SetFillColor(colors["W"]);
	
	TH1F * tttY_ = new TH1F("tttY_","tttY_",binnum,bin);
	tttY_->Sumw2();
	tttY_->SetFillColor(colors["TTT"]);
    
    TH1F * ttjY_ = new TH1F("ttjY_","ttjY_",binnum,bin);
	ttjY_->Sumw2();
	ttjY_->SetFillColor(colors["TTJ"]);
	
	TH1F * vvtY_ = new TH1F("vvtY_","vvtY_",binnum,bin);
	vvtY_->Sumw2();
	vvtY_->SetFillColor(colors["VV"]);
    
    TH1F * vvjY_ = new TH1F("vvjY_","vvjY_",binnum,bin);
	vvjY_->Sumw2();
	vvjY_->SetFillColor(colors["VV"]);
    
    TH1F * zvvY_ = new TH1F("zvvY_","zvvY_",binnum,bin);
	zvvY_->Sumw2();
	zvvY_->SetFillColor(colors["ZVV"]);
    
    TH1F * ewkY_ = new TH1F("ewkY_","ewkY_",binnum,bin);
	ewkY_->Sumw2();
	ewkY_->SetFillColor(colors["EWK"]);
	
	TH1F * zhtautauY_ = new TH1F("zhtautauY_","zhtautauY_",binnum,bin);
	zhtautauY_->Sumw2();
	zhtautauY_->SetFillColor(colors["ZHTauTau"]);

	TH1F * glugluhtautauY_ = new TH1F("glugluhtautauY_","glugluhtautauY_",binnum,bin);
	glugluhtautauY_->Sumw2();
	glugluhtautauY_->SetFillColor(colors["GluGluHTauTau"]);
    
    TH1F * vbfhtautauY_ = new TH1F("vbfhtautauY_","vbfhtautauY_",binnum,bin);
	vbfhtautauY_->Sumw2();
    vbfhtautauY_->SetFillColor(colors["VBFHTauTau"]);

    TH1F * smhY_ = new TH1F("smhY_","smhY_",binnum,bin);
	smhY_->Sumw2();
    smhY_->SetFillColor(colors["SMH"]);
    
    TH1F * combY_ = new TH1F("combY_","combY_",binnum,bin);
	combY_->Sumw2();

	DATA->Draw((parameter+">>dataY_").c_str(),SScut);
	W->Draw((parameter+">>wY_").c_str(),SScutRel*weights["W"].c_str());
	ZTT->Draw((parameter+">>zttY_").c_str(),SScutRel*weights["ZTT"].c_str());
    ZL->Draw((parameter+">>zlY_").c_str(),SScutRel*weights["ZL"].c_str());
    ZJ->Draw((parameter+">>zjY_").c_str(),SScutRel*weights["ZJ"].c_str());
	VVT->Draw((parameter+">>vvtY_").c_str(),SScutRel*weights["VV"].c_str());
    VVJ->Draw((parameter+">>vvjY_").c_str(),SScutRel*weights["VV"].c_str());
    ZVV->Draw((parameter+">>zvvY_").c_str(),SScut*weights["ZVV"].c_str());
    EWK->Draw((parameter+">>ewkY_").c_str(),SScut*weights["EWK"].c_str());
	TTT->Draw((parameter+">>tttY_").c_str(),SScutRel*weights["TTT"].c_str());
    TTJ->Draw((parameter+">>ttjY_").c_str(),SScutRel*weights["TTJ"].c_str());
	ZHTauTau->Draw((parameter+">>zhtautauY_").c_str(),SScut*weights["ZHTauTau"].c_str());
    GluGluHTauTau->Draw((parameter+">>glugluhtautauY_").c_str(),SScut*weights["GluGluHTauTau"].c_str());
    VBFHTauTau->Draw((parameter+">>vbfhtautauY_").c_str(),SScut*weights["VBFHTauTau"].c_str());
    
    wY_->Scale(shapeScaleVarBin(SScut, SScutRel, W, weights["W"], binnum, bin));
    vvtY_->Scale(shapeScaleVarBin(SScut, SScutRel, VVT, weights["VV"], binnum, bin));
    vvjY_->Scale(shapeScaleVarBin(SScut, SScutRel, VVJ, weights["VV"], binnum, bin));
    tttY_->Scale(shapeScaleVarBin(SScut, SScutRel, TTT, weights["TTT"], binnum, bin));
    ttjY_->Scale(shapeScaleVarBin(SScut, SScutRel, TTJ, weights["TTJ"], binnum, bin));
    zttY_->Scale(shapeScaleVarBin(SScut, SScutRel, ZTT, weights["ZTT"], binnum, bin));
    zlY_->Scale(shapeScaleVarBin(SScut, SScutRel, ZL, weights["ZL"], binnum, bin));
    zjY_->Scale(shapeScaleVarBin(SScut, SScutRel, ZJ, weights["ZJ"], binnum, bin));
    
    // scale norms to nominal from variant if necessary
    if(doNorm)
    {
        dataY_->Scale(shapeScaleVarBin(SSnormCut, SScut, DATA, "1.0", binnum, bin));
        vvtY_->Scale(shapeScaleVarBin(SSnormCut, SScut, VVT, weights["VV"], binnum, bin));
        vvjY_->Scale(shapeScaleVarBin(SSnormCut, SScut, VVJ, weights["VV"], binnum, bin));
        zvvY_->Scale(shapeScaleVarBin(SSnormCut, SScut, ZVV, weights["ZVV"], binnum, bin));
        ewkY_->Scale(shapeScaleVarBin(SSnormCut, SScut, EWK, weights["EWK"], binnum, bin));
        tttY_->Scale(shapeScaleVarBin(SSnormCut, SScut, TTT, weights["TTT"], binnum, bin));
        ttjY_->Scale(shapeScaleVarBin(SSnormCut, SScut, TTJ, weights["TTJ"], binnum, bin));
        zttY_->Scale(shapeScaleVarBin(SSnormCut, SScut, ZTT, weights["ZTT"], binnum, bin));
        zlY_->Scale(shapeScaleVarBin(SSnormCut, SScut, ZL, weights["ZL"], binnum, bin));
        zjY_->Scale(shapeScaleVarBin(SSnormCut, SScut, ZJ, weights["ZJ"], binnum, bin));
        wY_->Scale(shapeScaleVarBin(SSnormCut, SScut, W, weights["W"], binnum, bin));
        zhtautauY_->Scale(shapeScaleVarBin(SSnormCut, SScut, ZHTauTau, weights["ZHTauTau"], binnum, bin));
        glugluhtautauY_->Scale(shapeScaleVarBin(SSnormCut, SScut, GluGluHTauTau, weights["GluGluHTauTau"], binnum, bin));
        vbfhtautauY_->Scale(shapeScaleVarBin(SSnormCut, SScut, VBFHTauTau, weights["VBFHTauTau"], binnum, bin));
    }
    
    smhY_->Add(zhtautauY_);
    smhY_->Add(glugluhtautauY_);
    smhY_->Add(vbfhtautauY_);

    std::cout << "SS Normal Cut" << std::endl;
    std::cout << "data: " << dataY_->Integral() << std::endl;
    std::cout << "W: " << wY_->Integral() << std::endl;
    std::cout << "ZTT: " << zttY_->Integral() << std::endl;
    std::cout << "ZL: " << zlY_->Integral() << std::endl;
    std::cout << "ZJ: " << zjY_->Integral() << std::endl;
    std::cout << "VVT: " << vvtY_->Integral() << std::endl;
    std::cout << "VVJ: " << vvjY_->Integral() << std::endl;
    std::cout << "ZVV: " << zvvY_->Integral() << std::endl;
    std::cout << "EWK: " << ewkY_->Integral() << std::endl;
    std::cout << "TTT: " << tttY_->Integral() << std::endl;
    std::cout << "TTJ: " << ttjY_->Integral() << std::endl;
    std::cout << "ZH: " << zhtautauY_->Integral() << std::endl;
    std::cout << "VBF H: " << glugluhtautauY_->Integral() << std::endl;
    std::cout << "GG H: " << vbfhtautauY_->Integral() << std::endl;

	/* QCD is data - all mc bks */
    double qcdYield;
    qcdYield = TMath::Nint(dataY_->Integral() - wY_->Integral() - zttY_->Integral() - zlY_->Integral() - zjY_->Integral() - tttY_->Integral() - ttjY_->Integral() - vvtY_->Integral() - vvjY_->Integral() - zvvY_->Integral() - ewkY_->Integral() - smhY_->Integral());
    
    std::cout << "QCD: " << qcdYield << std::endl;
    std::cout << "QCD Shape: " << qcdShape_->Integral() << std::endl;
    
    double corrFactor;
    if (qcdShape_->Integral() == 0.) {corrFactor = 1.0;}
    else {corrFactor = qcdYield/qcdShape_->Integral();}
    std::cout << "correction factor for QCD relaxed to full yield: " << corrFactor << std::endl;
    qcdShape_->Scale(corrFactor);
    std::cout << "QCD after scaling: " << qcdShape_->Integral() << std::endl;
    
    /* eliminate negative bins */
    for(int i = 0; i<qcdShape_->GetNbinsX()+1; i++)
	{
		if(qcdShape_->GetBinContent(i)<0.) {qcdShape_->SetBinContent(i,0.);}
	}
    
    std::cout << "QCD after rm negative bins and final scaling: " << qcdShape_->Integral() << std::endl;
    
    std::cout << "Editing Shape File" << std::endl;
    
	if(createOutputShapes)
	{
        
        std::string paramName;
        if (parameter=="met" || parameter=="mvamet") {paramName = "met";}
        else if (parameter.substr(0,6)=="mvaVar") {paramName = "mva";}
        else if (parameter.substr(0,6) == "mt_tot") {paramName="mt";}
        else {paramName = parameter;}
        
        std::string catName = chan + "_QCD_inclusive_cr";
		std::string optName = "xtt_" + chan + ".inputs-13TeV-" + paramName + ".root";
		TFile sbOp(optName.c_str(),"UPDATE");
        
        TDirectory * catdir = sbOp.GetDirectory(catName.c_str());
        if (catdir){catdir->cd();}
        else {sbOp.mkdir(catName.c_str()); TDirectory * catdir = sbOp.GetDirectory(catName.c_str()); catdir->cd();}
        
        std::string DATAname = "data_obs" + syst;
        const char *DATAsave = DATAname.c_str();
        std::string VVTname = "VVT" + syst;
        const char *VVTsave = VVTname.c_str();
        std::string VVJname = "VVJ" + syst;
        const char *VVJsave = VVJname.c_str();
        std::string ZVVname = "ZVV" + syst;
        const char *ZVVsave = ZVVname.c_str();
        std::string EWKname = "EWK" + syst;
        const char *EWKsave = EWKname.c_str();
        std::string TTTname = "TTT" + syst;
        const char *TTTsave = TTTname.c_str();
        std::string TTJname = "TTJ" + syst;
        const char *TTJsave = TTJname.c_str();
        std::string ZTTname = "ZTT" + syst;
        const char *ZTTsave = ZTTname.c_str();
        std::string ZLname = "ZL" + syst;
        const char *ZLsave = ZLname.c_str();
        std::string ZJname = "ZJ" + syst;
        const char *ZJsave = ZJname.c_str();
        std::string Wname = "W" + syst;
        const char *Wsave = Wname.c_str();
        std::string QCDname = "QCD" + syst;
        const char *QCDsave = QCDname.c_str();
        std::string preQCDname = "preQCD" + syst;
        const char *preQCDsave = preQCDname.c_str();

        std::string SMHname = "SMH" + syst;
        const char *SMHsave = SMHname.c_str();
        std::string ZHname = "ZH" + syst;
        const char *ZHsave = ZHname.c_str();

        dataY_->Write(DATAsave);
        
        vvtY_->Write(VVTsave);
        vvjY_->Write(VVJsave);
        zvvY_->Write(ZVVsave);
        ewkY_->Write(EWKsave);
        tttY_->Write(TTTsave);
        ttjY_->Write(TTJsave);
        zttY_->Write(ZTTsave);
        zlY_->Write(ZLsave);
        zjY_->Write(ZJsave);
        wY_->Write(Wsave);
        qcdShape_->Write(QCDsave);
        preQCDShape_->Write(preQCDsave);
        
        zhtautauY_->Write(ZHsave);
        smhY_->Write(SMHsave);

        sbOp.cd();
		sbOp.Close();
	}
    std::cout << "FINISHED editing shape file" << std::endl;
    
    THStack * QCDStack  = new THStack();

    if (plots)
    {

        QCDStack->Add(tttY_);
        QCDStack->Add(ttjY_);
        QCDStack->Add(zvvY_);
        QCDStack->Add(ewkY_);
        QCDStack->Add(qcdShape_);
        QCDStack->Add(smhY_);
        QCDStack->Add(vvtY_);
        QCDStack->Add(vvjY_);
        QCDStack->Add(wY_);
        QCDStack->Add(zttY_);
        QCDStack->Add(zlY_);
        QCDStack->Add(zjY_);

        QCD_SS->cd(1);dataY_->SetMaximum(DATA_MAX*dataY_->GetMaximum());dataY_->SetTitle(";XXX;Events/Bin");dataY_->GetXaxis()->SetTitle(x_axisLabels[parameter].c_str());dataY_->DrawCopy("PE");QCDStack->Draw("histsames");dataY_->DrawCopy("PEsames");QCD_SS->Update();

        combY_->Add(vvtY_);
        combY_->Add(vvjY_);
        combY_->Add(zvvY_);
        combY_->Add(ewkY_);
        combY_->Add(tttY_);
        combY_->Add(ttjY_);
        combY_->Add(wY_);
        combY_->Add(zttY_);
        combY_->Add(zlY_);
        combY_->Add(zjY_);
        combY_->Add(smhY_);
        combY_->Add(qcdShape_);

        drawTcanvasAndLegendAndRatio(QCD_SS, dataY_, combY_, 0);QCD_SS->DrawClone();std::string localSaveName = "QCD_"+saveName + parameter + "_" + chan + "_" + drCutNameMap[drCut] + "Region.png";QCD_SS->cd()->SaveAs(localSaveName.c_str());QCD_SS->Close();
    }
    
	//std::cout<<" ***** "<<qcd_->Integral()<<"\n";

    //finally returns shape used for QCD in SR plots
	QCDshape->Add(qcdShape_);

	/* delete pointers */
	delete dataY_;
	delete vvtY_;
    delete vvjY_;
    delete zvvY_;
    delete ewkY_;
	delete tttY_;
    delete ttjY_;
	delete zttY_;
    delete zlY_;
    delete zjY_;
	delete wY_;
    delete zhtautauY_;
    delete glugluhtautauY_;
    delete vbfhtautauY_;
    delete smhY_;
    delete qcdShape_;
    delete preQCDShape_;
	delete combY_;
	delete QCDStack;
	delete QCD_SS;
}

double wjetsNorm(TCut wCut, TCut wCutRel, TCut wSSCut, TCut wSSCutRel, TCut wNormCut, TCut wSSNormCut, bool doNorm, std::string parameter, int binnum, float bin[], std::string can_name, std::string chan, std::string syst, bool createOutputShapes) //not ready for 80X yet, may be used
{

    //not currenlty using relaxed shape in anti-iso, statistics are high as is

	/* create a TCanvas */
    TCanvas * wjetNorm_can = new TCanvas(can_name.c_str(),can_name.c_str(),canDim[0],canDim[1]);
    wjetNorm_can->Divide(1,2);

	/* create some histograms and format them */

	TH1F * data_ = new TH1F("data_","",binnum,bin);
	data_->Sumw2(); 
	data_->SetMarkerStyle(dataMarker);

	TH1F * ztt_ = new TH1F("ztt_","ztt_",binnum,bin); 
	ztt_->Sumw2();
	ztt_->SetFillColor(colors["ZTT"]);
    
    TH1F * zl_ = new TH1F("zl_","zl_",binnum,bin); 
	zl_->Sumw2();
	zl_->SetFillColor(colors["ZL"]);
    
    TH1F * zj_ = new TH1F("zj_","zj_",binnum,bin); 
	zj_->Sumw2();
	zj_->SetFillColor(colors["ZJ"]);
	
	TH1F * w_ = new TH1F("w_","w_",binnum,bin);  
	w_->Sumw2(); 
	w_->SetFillColor(colors["W"]);
	
	TH1F * ttj_ = new TH1F("ttj_","ttj_",binnum,bin); 
	ttj_->Sumw2();
	ttj_->SetFillColor(colors["TTJ"]);
    
    TH1F * ttt_ = new TH1F("ttt_","ttt_",binnum,bin); 
	ttt_->Sumw2();
	ttt_->SetFillColor(colors["TTT"]);
	
	TH1F * vvt_ = new TH1F("vvt_","vvt_",binnum,bin); 
	vvt_->Sumw2();
	vvt_->SetFillColor(colors["VV"]);
    
    TH1F * vvj_ = new TH1F("vvj_","vvj_",binnum,bin); 
	vvj_->Sumw2();
	vvj_->SetFillColor(colors["VV"]);
    
    TH1F * zvv_ = new TH1F("zvv_","zvv_",binnum,bin); 
	zvv_->Sumw2();
	zvv_->SetFillColor(colors["ZVV"]);
    
    TH1F * ewk_ = new TH1F("ewk_","ewk_",binnum,bin); 
	ewk_->Sumw2();
	ewk_->SetFillColor(colors["EWK"]);
    
    TH1F * qcd_ = new TH1F("qcd_","qcd_",binnum,bin); 
	qcd_->Sumw2();
	qcd_->SetFillColor(colors["QCD"]);
	
	TH1F * zhtautau_ = new TH1F("zhtautau_","zhtautau_",binnum,bin); 
	zhtautau_->Sumw2();
	zhtautau_->SetFillColor(colors["ZHTauTau"]);

	TH1F * glugluhtautau_ = new TH1F("glugluhtautau_","glugluhtautau_",binnum,bin); 
	glugluhtautau_->Sumw2();
	glugluhtautau_->SetFillColor(colors["GluGluHTauTau"]);
    
    TH1F * vbfhtautau_ = new TH1F("vbfhtautau_","vbfhtautau_",binnum,bin); 
	vbfhtautau_->Sumw2();
    vbfhtautau_->SetFillColor(colors["VBFHTauTau"]);
    
    TH1F * smh_ = new TH1F("smh_","smh_",binnum,bin); 
	smh_->Sumw2();
    smh_->SetFillColor(colors["SMH"]);

	TH1F * comb_ = new TH1F("comb_","comb_",binnum,bin);  
	comb_->Sumw2(); 

	TH1F * combSUM_ = new TH1F("combSUM_","combSUM_",binnum,bin);
	combSUM_->Sumw2();

	/* fill the hists */
	wjetNorm_can->cd(1);

	DATA->Draw((parameter+">>data_").c_str(),wCut);
	W->Draw((parameter+">>w_").c_str(),wCut*weights["W"].c_str());
	ZTT->Draw((parameter+">>ztt_").c_str(),wCut*weights["ZTT"].c_str());
    ZJ->Draw((parameter+">>zj_").c_str(),wCut*weights["ZJ"].c_str());
    ZL->Draw((parameter+">>zl_").c_str(),wCut*weights["ZL"].c_str());
	VVT->Draw((parameter+">>vvt_").c_str(),wCut*weights["VV"].c_str());
    VVJ->Draw((parameter+">>vvj_").c_str(),wCut*weights["VV"].c_str());
    ZVV->Draw((parameter+">>zvv_").c_str(),wCut*weights["ZVV"].c_str());
    EWK->Draw((parameter+">>ewk_").c_str(),wCut*weights["EWK"].c_str());
	TTT->Draw((parameter+">>ttt_").c_str(),wCut*weights["TTT"].c_str());
    TTJ->Draw((parameter+">>ttj_").c_str(),wCut*weights["TTJ"].c_str());
	ZHTauTau->Draw((parameter+">>zhtautau_").c_str(),wCut*weights["ZHTauTau"].c_str());
    GluGluHTauTau->Draw((parameter+">>glugluhtautau_").c_str(),wCut*weights["GluGluHTauTau"].c_str());
    VBFHTauTau->Draw((parameter+">>vbfhtautau_").c_str(),wCut*weights["VBFHTauTau"].c_str());
    
    fillQCD_ShapeVarBin(wSSCut, wSSCutRel, 1.0, wSSNormCut, doNorm, qcd_, parameter, binnum, bin, "QCD AI SS sideband in " + chan + " channel", syst, chan, 0, 0);
    
    //qcd_->Scale(qcdOStoSS);
    
    // scale norms to nominal from variant if necessary
    if(doNorm)
    {
        vvt_->Scale(shapeScaleVarBin(wNormCut, wCut, VVT, weights["VV"], binnum, bin));
        vvj_->Scale(shapeScaleVarBin(wNormCut, wCut, VVJ, weights["VV"], binnum, bin));
        zvv_->Scale(shapeScaleVarBin(wNormCut, wCut, ZVV, weights["ZVV"], binnum, bin));
        ewk_->Scale(shapeScaleVarBin(wNormCut, wCut, EWK, weights["EWK"], binnum, bin));
        ttt_->Scale(shapeScaleVarBin(wNormCut, wCut, TTT, weights["TTT"], binnum, bin));
        ttj_->Scale(shapeScaleVarBin(wNormCut, wCut, TTJ, weights["TTJ"], binnum, bin));
        ztt_->Scale(shapeScaleVarBin(wNormCut, wCut, ZTT, weights["ZTT"], binnum, bin));
        zl_->Scale(shapeScaleVarBin(wNormCut, wCut, ZL, weights["ZL"], binnum, bin));
        zj_->Scale(shapeScaleVarBin(wNormCut, wCut, ZJ, weights["ZJ"], binnum, bin));
        w_->Scale(shapeScaleVarBin(wNormCut, wCut, W, weights["W"], binnum, bin));
        zhtautau_->Scale(shapeScaleVarBin(wNormCut, wCut, ZHTauTau, weights["ZHTauTau"], binnum, bin));
        glugluhtautau_->Scale(shapeScaleVarBin(wNormCut, wCut, GluGluHTauTau, weights["GluGluHTauTau"], binnum, bin));
        vbfhtautau_->Scale(shapeScaleVarBin(wNormCut, wCut, VBFHTauTau, weights["VBFHTauTau"], binnum, bin));
    }
    
    smh_->Add(zhtautau_);
    smh_->Add(glugluhtautau_);
    smh_->Add(vbfhtautau_);
    
	THStack * wnormStack  = new THStack();

	/* data minus all non-W */
	comb_->Add(data_);
    comb_->Add(vvt_,-1);
    comb_->Add(vvj_,-1);
    comb_->Add(zvv_,-1);
    comb_->Add(ewk_,-1);
    comb_->Add(qcd_,-1);
	comb_->Add(ttt_,-1);
    comb_->Add(ttj_,-1);
	comb_->Add(ztt_,-1);
    comb_->Add(zj_,-1);
    comb_->Add(zl_,-1);
    comb_->Add(smh_,-1);

    /* eliminate negative bins */
	std::cout<<" for w norm,  data-(non-W) = "<<comb_->Integral()<<"\n";
	std::cout<<" while W mc raw norm = "<<w_->Integral()<<"\n";
	double sf = comb_->Integral()/w_->Integral();

	std::cout<<" SF (not currently used in datacards) is "<<sf<<" for systematic: " << syst << std::endl;
    if (sf < 0)
    {
        std::cout<<" SF is negative, setting to 1.0" << std::endl;
        sf = 1.0;
    }
	//w_->Scale(sf);

    if(createOutputShapes)
	{
		/* create a file to help with S/sqrt(S+B) cut opt */
        
        std::string paramName;
        paramName = "mt";
        
        std::string catName = chan + "_W_inclusive_cr";
		std::string optName = "xtt_" + chan + ".inputs-13TeV-mt.root";
		TFile sbOp(optName.c_str(),"UPDATE");
        
        TDirectory * catdir = sbOp.GetDirectory(catName.c_str());
        if (catdir){catdir->cd();}
        else {sbOp.mkdir(catName.c_str()); TDirectory * catdir = sbOp.GetDirectory(catName.c_str()); catdir->cd();}
        
        std::string DATAname = "data_obs" + syst;
        const char *DATAsave = DATAname.c_str();
        std::string VVTname = "VVT" + syst;
        const char *VVTsave = VVTname.c_str();
        std::string VVJname = "VVJ" + syst;
        const char *VVJsave = VVJname.c_str();
        std::string ZVVname = "ZVV" + syst;
        const char *ZVVsave = ZVVname.c_str();
        std::string EWKname = "EWK" + syst;
        const char *EWKsave = EWKname.c_str();
        std::string TTTname = "TTT" + syst;
        const char *TTTsave = TTTname.c_str();
        std::string TTJname = "TTJ" + syst;
        const char *TTJsave = TTJname.c_str();
        std::string ZTTname = "ZTT" + syst;
        const char *ZTTsave = ZTTname.c_str();
        std::string ZLname = "ZL" + syst;
        const char *ZLsave = ZLname.c_str();
        std::string ZJname = "ZJ" + syst;
        const char *ZJsave = ZJname.c_str();
        std::string Wname = "W" + syst;
        const char *Wsave = Wname.c_str();
        std::string QCDname = "QCD" + syst;
        const char *QCDsave = QCDname.c_str();
        std::string SMHname = "SMH" + syst;
        const char *SMHsave = SMHname.c_str();
        std::string ZHname = "ZH" + syst;
        const char *ZHsave = ZHname.c_str();
        
        std::string combname = "comb" + syst;
        const char *combsave = combname.c_str();

        data_->Write(DATAsave);
        
        vvt_->Write(VVTsave);
        vvj_->Write(VVJsave);
        zvv_->Write(ZVVsave);
        ewk_->Write(EWKsave);
        ttt_->Write(TTTsave);
        ttj_->Write(TTJsave);
        ztt_->Write(ZTTsave);
        zl_->Write(ZLsave);
        zj_->Write(ZJsave);
        w_->Write(Wsave);
        qcd_->Write(QCDsave);
        
        zhtautau_->Write(ZHsave);
        smh_->Write(SMHsave);

        sbOp.cd();
		sbOp.Close();
	}
    std::cout << "FINISHED editing AI shape file" << std::endl;


	std::cout<<" post sf W norm = "<<w_->Integral()<<"\n";

	wnormStack->Add(vvt_);
    wnormStack->Add(vvj_);
    wnormStack->Add(zvv_);
    wnormStack->Add(ewk_);
    wnormStack->Add(qcd_);
	wnormStack->Add(ttt_);
    wnormStack->Add(ttj_);
	wnormStack->Add(w_);
	wnormStack->Add(ztt_);
    wnormStack->Add(zj_);
    wnormStack->Add(zl_);
    wnormStack->Add(smh_);

    wjetNorm_can->cd(1); data_->SetMaximum(DATA_MAX*data_->GetMaximum());/*data_->SetTitle((can_name+";XXX;Events/Bin").c_str());*/data_->SetTitle(";XXX;Events/Bin");data_->GetXaxis()->SetTitle(x_axisLabels["m_vis"].c_str());data_->DrawCopy("PE");wnormStack->Draw("histsames");data_->DrawCopy("PEsames");wjetNorm_can->Update();

	combSUM_->Add(vvt_);
    combSUM_->Add(vvj_);
    combSUM_->Add(zvv_);
    combSUM_->Add(ewk_);
    combSUM_->Add(qcd_);
	combSUM_->Add(ttt_);
    combSUM_->Add(ttj_);
	combSUM_->Add(w_);
	combSUM_->Add(ztt_);
    combSUM_->Add(zl_);
    combSUM_->Add(zj_);
    combSUM_->Add(smh_);

    drawTcanvasAndLegendAndRatio(wjetNorm_can, data_, combSUM_,0);wjetNorm_can->DrawClone();std::string localSaveName = can_name+"_"+saveName + "MT.png";wjetNorm_can->cd()->SaveAs(localSaveName.c_str());wjetNorm_can->Close();

	/* delete pointers */
	delete data_;
    delete comb_;
	delete combSUM_;
	delete vvt_;
    delete vvj_;
    delete zvv_;
    delete ewk_;
    delete qcd_;
	delete ttt_;
    delete ttj_;
	delete ztt_;
    delete zl_;
    delete zj_;
	delete w_;
    delete zhtautau_;
    delete glugluhtautau_;
    delete vbfhtautau_;
    delete smh_;

	delete wnormStack;
	delete wjetNorm_can;

	return sf;
}

double shapeScaleVarBin(TCut NomCut,TCut RelCut, TChain * CHAIN, std::string weight, int binnum, float bin[])
{

	/* create some histograms and format them */
	TH1F * T_ = new TH1F("T_","T_",binnum,bin);
	T_->Sumw2();
    
	TH1F * Rel_ = new TH1F("Rel_","Rel_",binnum,bin);
	Rel_->Sumw2();
    
	CHAIN->Draw("mt_tot>>T_",NomCut*weight.c_str());
    CHAIN->Draw("mt_tot>>Rel_",RelCut*weight.c_str());

    double sf = T_->Integral()/Rel_->Integral();
    if (Rel_->Integral() == 0.0) {sf = 1.0;}
    
    /*std::cout << "Nominal Integral: " << T_->Integral() << std::endl;
    std::cout << "Relaxed Integral: " << Rel_->Integral() << std::endl;
    std::cout << "Scale Factor: " << sf << std::endl;*/
    
    delete T_;
    delete Rel_;
    
	return sf;
}

double shapeScale(TCut NomCut,TCut RelCut, TChain * CHAIN, std::string weight, std::string parameter, float bin[3])
{

	/* create some histograms and format them */
	TH1F * T_ = new TH1F("T_","T_", bin[0], bin[1], bin[2]);
	T_->Sumw2();
    
	TH1F * Rel_ = new TH1F("Rel_","Rel_", bin[0], bin[1], bin[2]);
	Rel_->Sumw2();
    
	CHAIN->Draw((parameter+">>T_").c_str(),NomCut*weight.c_str());
    CHAIN->Draw((parameter+">>Rel_").c_str(),RelCut*weight.c_str());

    double sf = T_->Integral()/Rel_->Integral();
    if (Rel_->Integral() == 0.0) {sf = 1.0;}
    
    std::cout << "Nominal Integral: " << T_->Integral() << std::endl;
    std::cout << "Relaxed Integral: " << Rel_->Integral() << std::endl;
    std::cout << "Scale Factor: " << sf << std::endl;
    
    delete T_;
    delete Rel_;
    
	return sf;
}

void countTotalsVarBin(TCut cut, TCut relCut, TCut SScut, TCut SSrelCut, double sfQCD, double sf, std::string parameter, int binnum, float bin[], std::string can_name, std::string chan, std::string syst)

{

    std::cout << "CUT: " << cut.GetTitle() << std::endl;
    
	/* create a TCanvas */
	TCanvas * Signal_region = new TCanvas(can_name.c_str(),can_name.c_str(),canDim[0],canDim[1]);
	Signal_region->Divide(1,2);

	/* create some histograms and format them */

	TH1F * data_ = new TH1F("data_","",binnum,bin);
	data_->Sumw2(); 
	data_->SetMarkerStyle(dataMarker);

	TH1F * qcd_ = new TH1F("qcd_","qcd_",binnum,bin); 
	qcd_->Sumw2(); 
	qcd_->SetFillColor(colors["QCD"]);

	TH1F * ztt_ = new TH1F("ztt_","ztt_",binnum,bin);
	ztt_->Sumw2();
	ztt_->SetFillColor(colors["ZTT"]);
    
    TH1F * zl_ = new TH1F("zl_","zl_",binnum,bin);
	zl_->Sumw2();
	zl_->SetFillColor(colors["ZL"]);
    
    TH1F * zj_ = new TH1F("zj_","zj_",binnum,bin);
	zj_->Sumw2();
	zj_->SetFillColor(colors["ZJ"]);
	
	TH1F * w_ = new TH1F("w_","w_",binnum,bin); 
	w_->Sumw2(); 
	w_->SetFillColor(colors["W"]);
	
	TH1F * ttj_ = new TH1F("ttj_","ttj_",binnum,bin);
	ttj_->Sumw2();
	ttj_->SetFillColor(colors["TTJ"]);
    
    TH1F * ttt_ = new TH1F("ttt_","ttt_",binnum,bin);
	ttt_->Sumw2();
	ttt_->SetFillColor(colors["TTT"]);
    
    TH1F * tt_ = new TH1F("tt_","tt_",binnum,bin);
	tt_->Sumw2();
	tt_->SetFillColor(colors["TT"]);
	
	TH1F * vvt_ = new TH1F("vvt_","vvt_",binnum,bin);
	vvt_->Sumw2();
	vvt_->SetFillColor(colors["VV"]);
    
    TH1F * vvj_ = new TH1F("vvj_","vvj_",binnum,bin);
	vvj_->Sumw2();
	vvj_->SetFillColor(colors["VV"]);
    
    TH1F * zvv_ = new TH1F("zvv_","zvv_",binnum,bin);
	zvv_->Sumw2();
	zvv_->SetFillColor(colors["ZVV"]);
    
    TH1F * ewk_ = new TH1F("ewk_","ewk_",binnum,bin);
	ewk_->Sumw2();
	ewk_->SetFillColor(colors["EWK"]);
	
	TH1F * zhtautau_ = new TH1F("zhtautau_","zhtautau_",binnum,bin);
	zhtautau_->Sumw2();
	zhtautau_->SetFillColor(colors["ZHTauTau"]);

	TH1F * glugluhtautau_ = new TH1F("glugluhtautau_","glugluhtautau_",binnum,bin);
	glugluhtautau_->Sumw2();
	glugluhtautau_->SetFillColor(colors["GluGluHTauTau"]);
    
    TH1F * vbfhtautau_ = new TH1F("vbfhtautau_","vbfhtautau_",binnum,bin);
	vbfhtautau_->Sumw2();
    vbfhtautau_->SetFillColor(colors["VBFHTauTau"]);
    
	TH1F * zp600a300_ = new TH1F("zp600a300_","zp600a300_",binnum,bin);
	zp600a300_->Sumw2();
    TH1F * zp600a400_ = new TH1F("zp600a400_","zp600a400_",binnum,bin);
	zp600a400_->Sumw2();

	TH1F * zp800a300_ = new TH1F("zp800a300_","zp800a300_",binnum,bin);
	zp800a300_->Sumw2();
    TH1F * zp800a400_ = new TH1F("zp800a400_","zp800a400_",binnum,bin);
	zp800a400_->Sumw2();
    TH1F * zp800a500_ = new TH1F("zp800a500_","zp800a500_",binnum,bin);
	zp800a500_->Sumw2();
    TH1F * zp800a600_ = new TH1F("zp800a600_","zp800a600_",binnum,bin);
	zp800a600_->Sumw2();

	TH1F * zp1000a300_ = new TH1F("zp1000a300_","zp1000a300_",binnum,bin);
	zp1000a300_->Sumw2();
    TH1F * zp1000a400_ = new TH1F("zp1000a400_","zp1000a400_",binnum,bin);
	zp1000a400_->Sumw2();
    TH1F * zp1000a500_ = new TH1F("zp1000a500_","zp1000a500_",binnum,bin);
	zp1000a500_->Sumw2();
    TH1F * zp1000a600_ = new TH1F("zp1000a600_","zp1000a600_",binnum,bin);
	zp1000a600_->Sumw2();
	TH1F * zp1000a700_ = new TH1F("zp1000a700_","zp1000a700_",binnum,bin);
	zp1000a700_->Sumw2();
    TH1F * zp1000a800_ = new TH1F("zp1000a800_","zp1000a800_",binnum,bin);
	zp1000a800_->Sumw2();
    
    TH1F * zp1200a300_ = new TH1F("zp1200a300_","zp1200a300_",binnum,bin);
	zp1200a300_->Sumw2();
    TH1F * zp1200a400_ = new TH1F("zp1200a400_","zp1200a400_",binnum,bin);
	zp1200a400_->Sumw2();
    TH1F * zp1200a500_ = new TH1F("zp1200a500_","zp1200a500_",binnum,bin);
	zp1200a500_->Sumw2();
    TH1F * zp1200a600_ = new TH1F("zp1200a600_","zp1200a600_",binnum,bin);
	zp1200a600_->Sumw2();
	TH1F * zp1200a700_ = new TH1F("zp1200a700_","zp1200a700_",binnum,bin);
	zp1200a700_->Sumw2();
    TH1F * zp1200a800_ = new TH1F("zp1200a800_","zp1200a800_",binnum,bin);
	zp1200a800_->Sumw2();
    
    TH1F * zp1400a300_ = new TH1F("zp1400a300_","zp1400a300_",binnum,bin);
	zp1400a300_->Sumw2();
    TH1F * zp1400a400_ = new TH1F("zp1400a400_","zp1400a400_",binnum,bin);
	zp1400a400_->Sumw2();
    TH1F * zp1400a500_ = new TH1F("zp1400a500_","zp1400a500_",binnum,bin);
	zp1400a500_->Sumw2();
    TH1F * zp1400a600_ = new TH1F("zp1400a600_","zp1400a600_",binnum,bin);
	zp1400a600_->Sumw2();
	TH1F * zp1400a700_ = new TH1F("zp1400a700_","zp1400a700_",binnum,bin);
	zp1400a700_->Sumw2();
    TH1F * zp1400a800_ = new TH1F("zp1400a800_","zp1400a800_",binnum,bin);
	zp1400a800_->Sumw2();
    
    TH1F * zp1700a300_ = new TH1F("zp1700a300_","zp1700a300_",binnum,bin);
	zp1700a300_->Sumw2();
    TH1F * zp1700a400_ = new TH1F("zp1700a400_","zp1700a400_",binnum,bin);
	zp1700a400_->Sumw2();
    TH1F * zp1700a500_ = new TH1F("zp1700a500_","zp1700a500_",binnum,bin);
	zp1700a500_->Sumw2();
    TH1F * zp1700a600_ = new TH1F("zp1700a600_","zp1700a600_",binnum,bin);
	zp1700a600_->Sumw2();
	TH1F * zp1700a700_ = new TH1F("zp1700a700_","zp1700a700_",binnum,bin);
	zp1700a700_->Sumw2();
    TH1F * zp1700a800_ = new TH1F("zp1700a800_","zp1700a800_",binnum,bin);
	zp1700a800_->Sumw2();
    
    TH1F * zp2000a300_ = new TH1F("zp2000a300_","zp2000a300_",binnum,bin);
	zp2000a300_->Sumw2();
    TH1F * zp2000a400_ = new TH1F("zp2000a400_","zp2000a400_",binnum,bin);
	zp2000a400_->Sumw2();
    TH1F * zp2000a500_ = new TH1F("zp2000a500_","zp2000a500_",binnum,bin);
	zp2000a500_->Sumw2();
    TH1F * zp2000a600_ = new TH1F("zp2000a600_","zp2000a600_",binnum,bin);
	zp2000a600_->Sumw2();
	TH1F * zp2000a700_ = new TH1F("zp2000a700_","zp2000a700_",binnum,bin);
	zp2000a700_->Sumw2();
    TH1F * zp2000a800_ = new TH1F("zp2000a800_","zp2000a800_",binnum,bin);
	zp2000a800_->Sumw2();
    
    TH1F * zp2500a300_ = new TH1F("zp2500a300_","zp2500a300_",binnum,bin);
	zp2500a300_->Sumw2();
    TH1F * zp2500a400_ = new TH1F("zp2500a400_","zp2500a400_",binnum,bin);
	zp2500a400_->Sumw2();
    TH1F * zp2500a500_ = new TH1F("zp2500a500_","zp2500a500_",binnum,bin);
	zp2500a500_->Sumw2();
    TH1F * zp2500a600_ = new TH1F("zp2500a600_","zp2500a600_",binnum,bin);
	zp2500a600_->Sumw2();
	TH1F * zp2500a700_ = new TH1F("zp2500a700_","zp2500a700_",binnum,bin);
	zp2500a700_->Sumw2();
    TH1F * zp2500a800_ = new TH1F("zp2500a800_","zp2500a800_",binnum,bin);
	zp2500a800_->Sumw2();

	/* sum of all bks */
	TH1F * comb_ = new TH1F("comb_","comb_",binnum,bin); 
	comb_->Sumw2(); 

	TH1F * smh_ = new TH1F("smh_","smh_",binnum,bin);
	smh_->Sumw2();
    smh_->SetFillColor(colors["SMH"]);

    DATA->Draw((parameter+">>data_").c_str(),cut);
    
	W->Draw((parameter+">>w_").c_str(),relCut*weights["W"].c_str());
	ZTT->Draw((parameter+">>ztt_").c_str(),relCut*weights["ZTT"].c_str());
    ZL->Draw((parameter+">>zl_").c_str(),relCut*weights["ZL"].c_str());
    ZJ->Draw((parameter+">>zj_").c_str(),relCut*weights["ZJ"].c_str());
	VVT->Draw((parameter+">>vvt_").c_str(),relCut*weights["VV"].c_str());
    VVJ->Draw((parameter+">>vvj_").c_str(),relCut*weights["VV"].c_str());
    ZVV->Draw((parameter+">>zvv_").c_str(),cut*weights["ZVV"].c_str());
    EWK->Draw((parameter+">>ewk_").c_str(),cut*weights["EWK"].c_str());
	TTT->Draw((parameter+">>ttt_").c_str(),relCut*weights["TTT"].c_str());
    TTJ->Draw((parameter+">>ttj_").c_str(),relCut*weights["TTJ"].c_str());
	ZHTauTau->Draw((parameter+">>zhtautau_").c_str(),cut*weights["ZHTauTau"].c_str());
    VBFHTauTau->Draw((parameter+">>vbfhtautau_").c_str(),cut*weights["VBFHTauTau"].c_str());
    GluGluHTauTau->Draw((parameter+">>glugluhtautau_").c_str(),cut*weights["GluGluHTauTau"].c_str());
    
    smh_->Add(zhtautau_);
    smh_->Add(glugluhtautau_);
    smh_->Add(vbfhtautau_);
    
    tt_->Add(ttt_);
    tt_->Add(ttj_);
    
	MZP600_MA0300->Draw((parameter+">>zp600a300_").c_str(),cut*weights["MONO"].c_str());
    MZP600_MA0400->Draw((parameter+">>zp600a400_").c_str(),cut*weights["MONO"].c_str());
    
	MZP800_MA0300->Draw((parameter+">>zp800a300_").c_str(),cut*weights["MONO"].c_str());
    MZP800_MA0400->Draw((parameter+">>zp800a400_").c_str(),cut*weights["MONO"].c_str());
    MZP800_MA0500->Draw((parameter+">>zp800a500_").c_str(),cut*weights["MONO"].c_str());
    MZP800_MA0600->Draw((parameter+">>zp800a600_").c_str(),cut*weights["MONO"].c_str());
    
	MZP1000_MA0300->Draw((parameter+">>zp1000a300_").c_str(),cut*weights["MONO"].c_str());
    MZP1000_MA0400->Draw((parameter+">>zp1000a400_").c_str(),cut*weights["MONO"].c_str());
    MZP1000_MA0500->Draw((parameter+">>zp1000a500_").c_str(),cut*weights["MONO"].c_str());
    MZP1000_MA0600->Draw((parameter+">>zp1000a600_").c_str(),cut*weights["MONO"].c_str());
    MZP1000_MA0700->Draw((parameter+">>zp1000a700_").c_str(),cut*weights["MONO"].c_str());
    MZP1000_MA0800->Draw((parameter+">>zp1000a800_").c_str(),cut*weights["MONO"].c_str());
    
    MZP1200_MA0300->Draw((parameter+">>zp1200a300_").c_str(),cut*weights["MONO"].c_str());
    MZP1200_MA0400->Draw((parameter+">>zp1200a400_").c_str(),cut*weights["MONO"].c_str());
    MZP1200_MA0500->Draw((parameter+">>zp1200a500_").c_str(),cut*weights["MONO"].c_str());
    MZP1200_MA0600->Draw((parameter+">>zp1200a600_").c_str(),cut*weights["MONO"].c_str());
    MZP1200_MA0700->Draw((parameter+">>zp1200a700_").c_str(),cut*weights["MONO"].c_str());
    MZP1200_MA0800->Draw((parameter+">>zp1200a800_").c_str(),cut*weights["MONO"].c_str());

	MZP1400_MA0300->Draw((parameter+">>zp1400a300_").c_str(),cut*weights["MONO"].c_str());
    MZP1400_MA0400->Draw((parameter+">>zp1400a400_").c_str(),cut*weights["MONO"].c_str());
    MZP1400_MA0500->Draw((parameter+">>zp1400a500_").c_str(),cut*weights["MONO"].c_str());
    MZP1400_MA0600->Draw((parameter+">>zp1400a600_").c_str(),cut*weights["MONO"].c_str());
    MZP1400_MA0700->Draw((parameter+">>zp1400a700_").c_str(),cut*weights["MONO"].c_str());
    MZP1400_MA0800->Draw((parameter+">>zp1400a800_").c_str(),cut*weights["MONO"].c_str());
    
	MZP1700_MA0300->Draw((parameter+">>zp1700a300_").c_str(),cut*weights["MONO"].c_str());
    MZP1700_MA0400->Draw((parameter+">>zp1700a400_").c_str(),cut*weights["MONO"].c_str());
    MZP1700_MA0500->Draw((parameter+">>zp1700a500_").c_str(),cut*weights["MONO"].c_str());
    MZP1700_MA0600->Draw((parameter+">>zp1700a600_").c_str(),cut*weights["MONO"].c_str());
    MZP1700_MA0700->Draw((parameter+">>zp1700a700_").c_str(),cut*weights["MONO"].c_str());
    MZP1700_MA0800->Draw((parameter+">>zp1700a800_").c_str(),cut*weights["MONO"].c_str());
    
	MZP2000_MA0300->Draw((parameter+">>zp2000a300_").c_str(),cut*weights["MONO"].c_str());
    MZP2000_MA0400->Draw((parameter+">>zp2000a400_").c_str(),cut*weights["MONO"].c_str());
    MZP2000_MA0500->Draw((parameter+">>zp2000a500_").c_str(),cut*weights["MONO"].c_str());
    MZP2000_MA0600->Draw((parameter+">>zp2000a600_").c_str(),cut*weights["MONO"].c_str());
    MZP2000_MA0700->Draw((parameter+">>zp2000a700_").c_str(),cut*weights["MONO"].c_str());
    MZP2000_MA0800->Draw((parameter+">>zp2000a800_").c_str(),cut*weights["MONO"].c_str());
    
	MZP2500_MA0300->Draw((parameter+">>zp2500a300_").c_str(),cut*weights["MONO"].c_str());
    MZP2500_MA0400->Draw((parameter+">>zp2500a400_").c_str(),cut*weights["MONO"].c_str());
    MZP2500_MA0500->Draw((parameter+">>zp2500a500_").c_str(),cut*weights["MONO"].c_str());
    MZP2500_MA0600->Draw((parameter+">>zp2500a600_").c_str(),cut*weights["MONO"].c_str());
    MZP2500_MA0700->Draw((parameter+">>zp2500a700_").c_str(),cut*weights["MONO"].c_str());
    MZP2500_MA0800->Draw((parameter+">>zp2500a800_").c_str(),cut*weights["MONO"].c_str());
    
	fillQCD_ShapeVarBin(SScut, SSrelCut, sfQCD, SScut, 0.0, qcd_, parameter, binnum, bin, "QCD SS " + parameter + " sideband in " + chan + " channel", syst, chan, 0, 0);
	qcd_->Scale(qcdOStoSS);
    
    //not W scaling in datacards
	//w_->Scale(sf);
    w_->Scale(shapeScaleVarBin(cut, relCut, W, weights["W"], binnum, bin));
    vvt_->Scale(shapeScaleVarBin(cut, relCut, VVT, weights["VV"], binnum, bin));
    vvj_->Scale(shapeScaleVarBin(cut, relCut, VVJ, weights["VV"], binnum, bin));
    ttt_->Scale(shapeScaleVarBin(cut, relCut, TTT, weights["TTT"], binnum, bin));
    ttj_->Scale(shapeScaleVarBin(cut, relCut, TTJ, weights["TTJ"], binnum, bin));
    ztt_->Scale(shapeScaleVarBin(cut, relCut, ZTT, weights["ZTT"], binnum, bin));
    zl_->Scale(shapeScaleVarBin(cut, relCut, ZL, weights["ZL"], binnum, bin));
    zj_->Scale(shapeScaleVarBin(cut, relCut, ZJ, weights["ZJ"], binnum, bin));

    comb_->Add(vvt_);
    comb_->Add(vvj_);
    comb_->Add(zvv_);
    comb_->Add(ewk_);
	comb_->Add(ttt_);
    comb_->Add(ttj_);
	comb_->Add(w_);
	comb_->Add(ztt_);
    comb_->Add(zj_);
    comb_->Add(zl_);
    comb_->Add(zhtautau_);
    comb_->Add(glugluhtautau_);
    comb_->Add(vbfhtautau_);
	comb_->Add(qcd_);
    
	// scale things to the projected lumi
	vvt_->Scale(lumi_sf);
    vvj_->Scale(lumi_sf);
    zvv_->Scale(lumi_sf);
    ewk_->Scale(lumi_sf);
	ttt_->Scale(lumi_sf);
    ttj_->Scale(lumi_sf);
	ztt_->Scale(lumi_sf);
    zl_->Scale(lumi_sf);
    zj_->Scale(lumi_sf);
	w_->Scale(lumi_sf);
	qcd_->Scale(lumi_sf);
    zhtautau_->Scale(lumi_sf);
    glugluhtautau_->Scale(lumi_sf);
    vbfhtautau_->Scale(lumi_sf);
    smh_->Scale(lumi_sf);
    
	zp600a300_->Scale(lumi_sf);
    zp600a400_->Scale(lumi_sf);
	zp800a300_->Scale(lumi_sf);
    zp800a400_->Scale(lumi_sf);
    zp800a500_->Scale(lumi_sf);
    zp800a600_->Scale(lumi_sf);
	zp1000a300_->Scale(lumi_sf);
    zp1000a400_->Scale(lumi_sf);
    zp1000a500_->Scale(lumi_sf);
    zp1000a600_->Scale(lumi_sf);
    zp1000a700_->Scale(lumi_sf);
    zp1000a800_->Scale(lumi_sf);
	zp1200a300_->Scale(lumi_sf);
    zp1200a400_->Scale(lumi_sf);
    zp1200a500_->Scale(lumi_sf);
    zp1200a600_->Scale(lumi_sf);
    zp1200a700_->Scale(lumi_sf);
    zp1200a800_->Scale(lumi_sf);
    zp1400a300_->Scale(lumi_sf);
    zp1400a400_->Scale(lumi_sf);
    zp1400a500_->Scale(lumi_sf);
    zp1400a600_->Scale(lumi_sf);
    zp1400a700_->Scale(lumi_sf);
    zp1400a800_->Scale(lumi_sf);
	zp1700a300_->Scale(lumi_sf);
    zp1700a400_->Scale(lumi_sf);
    zp1700a500_->Scale(lumi_sf);
    zp1700a600_->Scale(lumi_sf);
    zp1700a700_->Scale(lumi_sf);
    zp1700a800_->Scale(lumi_sf);
	zp2000a300_->Scale(lumi_sf);
    zp2000a400_->Scale(lumi_sf);
    zp2000a500_->Scale(lumi_sf);
    zp2000a600_->Scale(lumi_sf);
    zp2000a700_->Scale(lumi_sf);
    zp2000a800_->Scale(lumi_sf);
	zp2500a300_->Scale(lumi_sf);
    zp2500a400_->Scale(lumi_sf);
    zp2500a500_->Scale(lumi_sf);
    zp2500a600_->Scale(lumi_sf);
    zp2500a700_->Scale(lumi_sf);
    zp2500a800_->Scale(lumi_sf);
    
	comb_->Scale(lumi_sf);
    
    // print signal region bkg yields (under whatever cut is provided)

    std::cout<<" bkg totals under "<<can_name<<"\n";

    double error = 0;

    error = 0;
    vvt_->IntegralAndError(-1,1e3, error);
    std::cout<<" VVT ... "<<vvt_->Integral()<<" +/- "<<error<<"\n";
    
    error = 0;
    vvj_->IntegralAndError(-1,1e3, error);
    std::cout<<" VVJ ... "<<vvj_->Integral()<<" +/- "<<error<<"\n";
    
    error = 0;
    zvv_->IntegralAndError(-1,1e3, error);
    std::cout<<" ZVV ... "<< zvv_->Integral()<<" +/- "<<error<<"\n";
    
    error = 0;
    ewk_->IntegralAndError(-1,1e3, error);
    std::cout<<" EWK ... "<< ewk_->Integral()<<" +/- "<<error<<"\n";

    error = 0;
    ttt_->IntegralAndError(-1,1e3, error);
    std::cout<<" TTT ... "<<ttt_->Integral()<<" +/- "<<error<<"\n";
    
    error = 0;
    ttj_->IntegralAndError(-1,1e3, error);
    std::cout<<" TTJ ... "<<ttj_->Integral()<<" +/- "<<error<<"\n";

    error = 0;
    ztt_->IntegralAndError(-1,1e3, error);
    std::cout<<" ZTT ... "<<ztt_->Integral()<<" +/- "<<error<<"\n";

    error = 0;
    zl_->IntegralAndError(-1,1e3, error);
    std::cout<<" ZL ... "<<zl_->Integral()<<" +/- "<<error<<"\n";
    
    error = 0;
    zj_->IntegralAndError(-1,1e3, error);
    std::cout<<" ZJ ... "<<zj_->Integral()<<" +/- "<<error<<"\n";

    error = 0;
    w_->IntegralAndError(-1,1e3, error);
    std::cout<<" W ... "<<w_->Integral()<<" +/- "<<error<<"\n";
    
    error = 0;
    zhtautau_->IntegralAndError(-1,1e3, error);
    std::cout<<" ZHTauTau ... "<<zhtautau_->Integral()<<" +/- "<<error<<"\n";

    error = 0;
    glugluhtautau_->IntegralAndError(-1,1e3, error);
    std::cout<<" GluGluHTauTau ... "<<glugluhtautau_->Integral()<<" +/- "<<error<<"\n";
    
    error = 0;
    vbfhtautau_->IntegralAndError(-1,1e3, error);
    std::cout<<" VBFHTauTau ... "<<vbfhtautau_->Integral()<<" +/- "<<error<<"\n";

    error = 0;
    qcd_->IntegralAndError(-1,1e3, error);
    error = sqrt(qcd_->GetEntries());
    std::cout<<" QCD ... "<<qcd_->Integral()<<" +/- "<<error<<"\n";
    
    error = 0;
    smh_->IntegralAndError(-1,1e3, error);
    std::cout<<" SUM SMH ... "<<smh_->Integral()<<" +/- "<<error<<"\n";

    error = 0;
    zp600a300_->IntegralAndError(-1,1e3, error);
    std::cout<<" ZP600 ... "<<zp600a300_->Integral()<<" +/- "<<error<<"\n";

    error = 0;
    zp800a300_->IntegralAndError(-1,1e3, error);
    std::cout<<" ZP800 ... "<<zp800a300_->Integral()<<" +/- "<<error<<"\n";

    error = 0;
    zp1000a300_->IntegralAndError(-1,1e3, error);
    std::cout<<" ZP1000 ... "<<zp1000a300_->Integral()<<" +/- "<<error<<"\n";

    error = 0;
    zp1200a300_->IntegralAndError(-1,1e3, error);
    std::cout<<" ZP1200 ... "<<zp1200a300_->Integral()<<" +/- "<<error<<"\n";

    error = 0;
    zp1400a300_->IntegralAndError(-1,1e3, error);
    std::cout<<" ZP1400 ... "<<zp1400a300_->Integral()<<" +/- "<<error<<"\n";
    
    error = 0;
    zp1700a300_->IntegralAndError(-1,1e3, error);
    std::cout<<" ZP1700 ... "<<zp1700a300_->Integral()<<" +/- "<<error<<"\n";
    
    error = 0;
    zp2000a300_->IntegralAndError(-1,1e3, error);
    std::cout<<" ZP2000 ... "<<zp2000a300_->Integral()<<" +/- "<<error<<"\n";
    
    error = 0;
    zp2500a300_->IntegralAndError(-1,1e3, error);
    std::cout<<" ZP2500 ... "<<zp2500a300_->Integral()<<" +/- "<<error<<"\n";

    error = 0;
    comb_->IntegralAndError(-1,1e3, error);
    std::cout<<" totalBkg ... "<<comb_->Integral()<<" +/- "<<error<<"\n";
    
    error = 0;
    data_->IntegralAndError(-1,1e3, error);
    std::cout<<" totalObs ... "<<data_->Integral()<<" +/- "<<error<<"\n";

    // test

    double mvaChoiceI;
    double mvaChoiceE;

    double wI, wE;
    double vvtI, vvtE;
    double vvjI, vvjE;
    double zvvI, zvvE;
    double ewkI, ewkE;
    double zttI, zttE;
    double zlI, zlE;
    double zjI, zjE;
    double qcdI, qcdE;
    double tttI, tttE;
    double ttjI, ttjE;
    double zhtautauI, zhtautauE;
    double vbfhtautauI, vbfhtautauE;
    double glugluhtautauI, glugluhtautauE;
    double smhI, smhE;
    double zp600I, zp600E;
    double zp800I, zp800E;
    double zp1000I, zp1000E;
    double zp1200I, zp1200E;
    double zp1400I, zp1400E;
    double zp1700I, zp1700E;
    double zp2000I, zp2000E;
    double zp2500I, zp2500E;

    zp600I = zp600a300_->Integral(-1, zp600a300_->GetNbinsX()+1);
    zp600a300_->IntegralAndError(-1, zp600a300_->GetNbinsX()+1,zp600E);

    zp800I = zp800a300_->Integral(-1, zp800a300_->GetNbinsX()+1);
    zp800a300_->IntegralAndError(-1, zp800a300_->GetNbinsX()+1,zp800E);

    zp1000I = zp1000a300_->Integral(-1, zp1000a300_->GetNbinsX()+1);
    zp1000a300_->IntegralAndError(-1, zp1000a300_->GetNbinsX()+1,zp1000E);

    zp1200I = zp1200a300_->Integral(-1, zp1200a300_->GetNbinsX()+1);
    zp1200a300_->IntegralAndError(-1, zp1200a300_->GetNbinsX()+1,zp1200E);

    zp1400I = zp1400a300_->Integral(-1, zp1400a300_->GetNbinsX()+1);
    zp1400a300_->IntegralAndError(-1, zp1400a300_->GetNbinsX()+1,zp1400E);
    
    zp1700I = zp1700a300_->Integral(-1, zp1700a300_->GetNbinsX()+1);
    zp1700a300_->IntegralAndError(-1, zp1700a300_->GetNbinsX()+1,zp1700E);
    
    zp2000I = zp2000a300_->Integral(-1, zp2000a300_->GetNbinsX()+1);
    zp2000a300_->IntegralAndError(-1, zp2000a300_->GetNbinsX()+1,zp2000E);
    
    zp2500I = zp2500a300_->Integral(-1, zp2500a300_->GetNbinsX()+1);
    zp2500a300_->IntegralAndError(-1, zp2500a300_->GetNbinsX()+1,zp2500E);

    wI = w_->Integral(-1, w_->GetNbinsX()+1);
    w_->IntegralAndError(-1, w_->GetNbinsX()+1,wE);

    vvtI = vvt_->Integral(-1, vvt_->GetNbinsX()+1);
    vvt_->IntegralAndError(-1, vvt_->GetNbinsX()+1,vvtE);
    
    vvjI = vvj_->Integral(-1, vvj_->GetNbinsX()+1);
    vvj_->IntegralAndError(-1, vvj_->GetNbinsX()+1,vvjE);
    
    zvvI = zvv_->Integral(-1, zvv_->GetNbinsX()+1);
    zvv_->IntegralAndError(-1, zvv_->GetNbinsX()+1,zvvE);
    
    ewkI = ewk_->Integral(-1, ewk_->GetNbinsX()+1);
    ewk_->IntegralAndError(-1, ewk_->GetNbinsX()+1,ewkE);

    zttI = ztt_->Integral(-1, ztt_->GetNbinsX()+1);
    ztt_->IntegralAndError(-1, ztt_->GetNbinsX()+1,zttE);
    
    zlI = zl_->Integral(-1, zl_->GetNbinsX()+1);
    zl_->IntegralAndError(-1, zl_->GetNbinsX()+1,zlE);
    
    zjI = zj_->Integral(-1, zj_->GetNbinsX()+1);
    zj_->IntegralAndError(-1, zj_->GetNbinsX()+1,zjE);

    qcdI = qcd_->Integral(-1, qcd_->GetNbinsX()+1);
    qcd_->IntegralAndError(-1, qcd_->GetNbinsX()+1,qcdE);

    tttI = ttt_->Integral(-1, ttt_->GetNbinsX()+1);
    ttt_->IntegralAndError(-1, ttt_->GetNbinsX()+1,tttE);
    
    ttjI = ttj_->Integral(-1, ttj_->GetNbinsX()+1);
    ttj_->IntegralAndError(-1, ttj_->GetNbinsX()+1,ttjE);

    zhtautauI = zhtautau_->Integral(-1, zhtautau_->GetNbinsX()+1);
    zhtautau_->IntegralAndError(-1, zhtautau_->GetNbinsX()+1,zhtautauE);
    
    glugluhtautauI = glugluhtautau_->Integral(-1, glugluhtautau_->GetNbinsX()+1);
    glugluhtautau_->IntegralAndError(-1, glugluhtautau_->GetNbinsX()+1,glugluhtautauE);
    
    vbfhtautauI = vbfhtautau_->Integral(-1, vbfhtautau_->GetNbinsX()+1);
    vbfhtautau_->IntegralAndError(-1, vbfhtautau_->GetNbinsX()+1,vbfhtautauE);
    
    smhI = smh_->Integral(-1, smh_->GetNbinsX()+1);
    smh_->IntegralAndError(-1, smh_->GetNbinsX()+1,smhE);

    if(wI>0.0){ wE = 1 + wE/wI; } else {wI = 0.0; wE = 1.0;}
    if(vvtI>0.0){ vvtE = 1 + vvtE/vvtI; } else {vvtI = 0.0; vvtE = 1.0;}
    if(vvjI>0.0){ vvjE = 1 + vvjE/vvjI; } else {vvjI = 0.0; vvjE = 1.0;}
    if(zvvI>0.0){ zvvE = 1 + zvvE/zvvI; } else {zvvI = 0.0; zvvE = 1.0;}
    if(ewkI>0.0){ ewkE = 1 + ewkE/ewkI; } else {ewkI = 0.0; ewkE = 1.0;}
    if(zttI>0.0){ zttE = 1 + zttE/zttI; } else {zttI = 0.0; zttE = 1.0;}
    if(zlI>0.0){ zlE = 1 + zlE/zlI; } else {zlI = 0.0; zlE = 1.0;}
    if(zjI>0.0){ zjE = 1 + zjE/zjI; } else {zjI = 0.0; zjE = 1.0;}
    if(qcdI>0.0){ qcdE = 1 + qcdE/qcdI; } else {qcdI = 0.0; qcdE = 1.0;}
    if(tttI>0.0){ tttE = 1 + tttE/tttI; } else {tttI = 0.0; tttE = 1.0;}
    if(ttjI>0.0){ ttjE = 1 + ttjE/ttjI; } else {ttjI = 0.0; ttjE = 1.0;}
    if(zhtautauI>0.0){ zhtautauE = 1 + zhtautauE/zhtautauI; } else {zhtautauI = 0.0; zhtautauE = 1.0;}
    if(vbfhtautauI>0.0){ vbfhtautauE = 1 + vbfhtautauE/vbfhtautauI; } else {vbfhtautauI = 0.0; vbfhtautauE = 1.0;}
    if(glugluhtautauI>0.0){ glugluhtautauE = 1 + glugluhtautauE/glugluhtautauI; } else {glugluhtautauI = 0.0; glugluhtautauE = 1.0;}
    if(smhI>0.0){ smhE = 1 + smhE/smhI; } else {smhI = 0.0; smhE = 1.0;}
    if(zp600I>0.0){ zp600E = 1 + zp600E/zp600I; } else {zp600I = 0.0; zp600E = 1.0;}
    if(zp800I>0.0){ zp800E = 1 + zp800E/zp800I; } else {zp800I = 0.0; zp800E = 1.0;}
    if(zp1000I>0.0){ zp1000E = 1 + zp1000E/zp1000I; } else {zp1000I = 0.0; zp1000E = 1.0;}
    if(zp1200I>0.0){ zp1200E = 1 + zp1200E/zp1200I; } else {zp1200I = 0.0; zp1200E = 1.0;}
    if(zp1400I>0.0){ zp1400E = 1 + zp1400E/zp1400I; } else {zp1400I = 0.0; zp1400E = 1.0;}
    if(zp1700I>0.0){ zp1700E = 1 + zp1700E/zp1700I; } else {zp1700I = 0.0; zp1700E = 1.0;}
    if(zp2000I>0.0){ zp2000E = 1 + zp2000E/zp2000I; } else {zp2000I = 0.0; zp2000E = 1.0;}
    if(zp2500I>0.0){ zp2500E = 1 + zp2500E/zp2500I; } else {zp2500I = 0.0; zp2500E = 1.0;}
    if(comb_->Integral()!=0) std::cout<<" @ ZP = 600 s/rt(b) = "<<zp600I/sqrt(comb_->Integral())<<"\n";
    if(comb_->Integral()!=0) std::cout<<" @ ZP = 800 s/rt(b) = "<<zp800I/sqrt(comb_->Integral())<<"\n";
    if(comb_->Integral()!=0) std::cout<<" @ ZP = 1000 s/rt(b) = "<<zp1000I/sqrt(comb_->Integral())<<"\n";
    if(comb_->Integral()!=0) std::cout<<" @ ZP = 1200 s/rt(b) = "<<zp1200I/sqrt(comb_->Integral())<<"\n";
    if(comb_->Integral()!=0) std::cout<<" @ ZP = 1400 s/rt(b) = "<<zp1400I/sqrt(comb_->Integral())<<"\n";
    if(comb_->Integral()!=0) std::cout<<" @ ZP = 1700 s/rt(b) = "<<zp1700I/sqrt(comb_->Integral())<<"\n";
    if(comb_->Integral()!=0) std::cout<<" @ ZP = 2000 s/rt(b) = "<<zp2000I/sqrt(comb_->Integral())<<"\n";
    if(comb_->Integral()!=0) std::cout<<" @ ZP = 2500 s/rt(b) = "<<zp2500I/sqrt(comb_->Integral())<<"\n";

    if (choice2==1 && sigChoice==2 && (chan=="mt" || chan=="et"))
    {
        if (parameter.substr(15,4)=="600A"){mvaChoiceI = zp600I; mvaChoiceE = zp600E; }
        if (parameter.substr(15,4)=="800A"){mvaChoiceI = zp800I; mvaChoiceE = zp800E; }
        if (parameter.substr(15,4)=="1000"){mvaChoiceI = zp1000I; mvaChoiceE = zp1000E; }
        if (parameter.substr(15,4)=="1200"){mvaChoiceI = zp1200I; mvaChoiceE = zp1200E; }
        
        out_data << "QUANTITY,SIGNAL,W,VVT,VVJ,ZVV,EWK,ZTT,ZJ,ZL,QCD,TTT,TTJ,ZHTauTau,GluGluHTauTau,VBFHTauTau" << std::endl;
        
        out_data << drCutNameMap[drCut] + chan + parameter.substr(10,16) + " rate," << mvaChoiceI << "," << wI << "," << vvtI << "," << vvjI << "," << zvvI << "," << ewkI << "," << zttI << "," << zjI << "," << zlI << "," << qcdI << "," << tttI << "," << ttjI << "," << zhtautauI << "," << glugluhtautauI << "," << vbfhtautauI << std::endl;
        
        out_data << drCutNameMap[drCut] + chan + parameter.substr(10,16) + " stat.err," << mvaChoiceE << "," << wE << "," << vvtE << "," << vvjE << "," << zvvI << "," << ewkI << "," << zttE << "," << zjE << "," << zlE << "," << qcdE << "," << ttjE << "," << ttjE << "," << zhtautauE << "," << glugluhtautauE << "," << vbfhtautauE << std::endl;
    }
    else
    {
        out_data << "QUANTITY,W,VVT,VVJ,ZVV,EWK,ZTT,ZJ,ZL,QCD,TTT,TTJ,SMH" << std::endl;
        
        out_data << syst << "_" <<  chan << "," << wI << "," << vvtI << "," << vvjI << "," << zvvI << "," << ewkI << "," << zttI << "," << zjI << "," << zlI << "," << qcdI << "," << tttI << "," << ttjI << "," << smhI << std::endl;
        
        /*out_data << "QUANTITY,SIGNAL,W,VV,ZVV,EWK,ZTT,ZJ,ZL,QCD,TTT,TTJ,ZHTauTau,GluGluHTauTau,VBFHTauTau" << std::endl;
        out_data << syst + chan + "ZP600rate," << zp600I << "," << wI << "," << vvI << "," << zvvI << "," << ewkI << "," << zttI << "," << zjI << "," << zlI << "," << qcdI << "," << tttI << "," << ttjI << "," << zhtautauI << "," << glugluhtautauI << "," << vbfhtautauI << std::endl;
        out_data << drCutNameMap[drCut] + chan + "ZP600stat.err," << zp600E << "," << wE << "," << vvE << "," << zvvE << "," << ewkE << "," << zttE << "," << zjE << "," << zlE << "," << qcdE << "," << tttE << "," << ttjE << "," << zhtautauE << "," << glugluhtautauE << "," << vbfhtautauE << std::endl;
        
        out_data << drCutNameMap[drCut] + chan + "ZP800rate," << zp800I << "," << wI << "," << vvI << "," << zvvI << "," << ewkI << "," << zttI << "," << zjI << "," << zlI << "," << qcdI << "," << tttI << "," << ttjI << "," << zhtautauI << "," << glugluhtautauI << "," << vbfhtautauI << std::endl;
        out_data << drCutNameMap[drCut] + chan + "ZP800stat.err," << zp800E << "," << wE << "," << vvE << "," << zvvE << "," << ewkE << "," << zttE << "," << zjE << "," << zlE << "," << qcdE << "," << tttE << "," << ttjE << "," << zhtautauE << "," << glugluhtautauE << "," << vbfhtautauE << std::endl;
        
        out_data << drCutNameMap[drCut] + chan + "ZP1000rate," << zp1000I << "," << wI << "," << vvI << "," << zvvI << "," << ewkI << "," << zttI << "," << zjI << "," << zlI << "," << qcdI << "," << tttI << "," << ttjI << "," << zhtautauI << "," << glugluhtautauI << "," << vbfhtautauI << std::endl;
        out_data << drCutNameMap[drCut] + chan + "ZP1000stat.err," << zp1000E << "," << wE << "," << vvE << "," << zvvE << "," << ewkE << "," << zttE << "," << zjE << "," << zlE << "," << qcdE << "," << tttE << "," << ttjE << "," << zhtautauE << "," << glugluhtautauE << "," << vbfhtautauE << std::endl;
        
        out_data << drCutNameMap[drCut] + chan + "ZP1200rate," << zp1200I << "," << wI << "," << vvI << "," << zvvI << "," << ewkI << "," << zttI << "," << zjI << "," << zlI << "," << qcdI << "," << tttI << "," << ttjI << "," << zhtautauI << "," << glugluhtautauI << "," << vbfhtautauI << std::endl;
        out_data << drCutNameMap[drCut] + chan + "ZP1200stat.err," << zp1200E << "," << wE << "," << vvE << "," << zvvE << "," << ewkE << "," << zttE << "," << zjE << "," << zlE << "," << qcdE << "," << tttE << "," << ttjE << "," << zhtautauE << "," << glugluhtautauE << "," << vbfhtautauE << std::endl;*/
    }
    
	delete data_;
	delete vvt_;
    delete vvj_;
    delete zvv_;
    delete ewk_;
	delete ttt_;
    delete ttj_;
    delete tt_;
	delete ztt_;
    delete zl_;
    delete zj_;
	delete w_;
    delete zhtautau_;
    delete glugluhtautau_;
    delete vbfhtautau_;
	delete comb_;
    delete smh_;
	delete qcd_;
	delete zp600a300_;
    delete zp600a400_;
	delete zp800a300_;
    delete zp800a400_;
    delete zp800a500_;
    delete zp800a600_;
	delete zp1000a300_;
    delete zp1000a400_;
    delete zp1000a500_;
    delete zp1000a600_;
    delete zp1000a700_;
    delete zp1000a800_;
	delete zp1200a300_;
    delete zp1200a400_;
    delete zp1200a500_;
    delete zp1200a600_;
    delete zp1200a700_;
    delete zp1200a800_;
    delete zp1400a300_;
    delete zp1400a400_;
    delete zp1400a500_;
    delete zp1400a600_;
    delete zp1400a700_;
    delete zp1400a800_;
    delete zp1700a300_;
    delete zp1700a400_;
    delete zp1700a500_;
    delete zp1700a600_;
    delete zp1700a700_;
    delete zp1700a800_;
    delete zp2000a300_;
    delete zp2000a400_;
    delete zp2000a500_;
    delete zp2000a600_;
    delete zp2000a700_;
    delete zp2000a800_;
    delete zp2500a300_;
    delete zp2500a400_;
    delete zp2500a500_;
    delete zp2500a600_;
    delete zp2500a700_;
    delete zp2500a800_;

}

void countTotals(TCut cut, TCut relCut, TCut SScut, TCut SScutRel, double sfQCD, double sf, TH1F * QCDshape, std::string parameter, float bin[3], std::string can_name, std::string chan, std::string syst)
{
    std::cout << "CUT: " << cut.GetTitle() << std::endl;

	/* create some histograms and format them */

	TH1F * data_ = new TH1F("data_","",bin[0],bin[1],bin[2]);
	data_->Sumw2(); 
	data_->SetMarkerStyle(dataMarker);

	TH1F * qcd_ = new TH1F("qcd_","qcd_",bin[0],bin[1],bin[2]); 
	qcd_->Sumw2(); 
	qcd_->SetFillColor(colors["QCD"]);

	TH1F * ztt_ = new TH1F("ztt_","ztt_",bin[0],bin[1],bin[2]);
	ztt_->Sumw2();
	ztt_->SetFillColor(colors["ZTT"]);
    
    TH1F * zl_ = new TH1F("zl_","zl_",bin[0],bin[1],bin[2]);
	zl_->Sumw2();
	zl_->SetFillColor(colors["ZL"]);
    
    TH1F * zj_ = new TH1F("zj_","zj_",bin[0],bin[1],bin[2]);
	zj_->Sumw2();
	zj_->SetFillColor(colors["ZJ"]);
	
	TH1F * w_ = new TH1F("w_","w_",bin[0],bin[1],bin[2]); 
	w_->Sumw2(); 
	w_->SetFillColor(colors["W"]);
	
	TH1F * ttj_ = new TH1F("ttj_","ttj_",bin[0],bin[1],bin[2]);
	ttj_->Sumw2();
	ttj_->SetFillColor(colors["TTJ"]);
    
    TH1F * ttt_ = new TH1F("ttt_","ttt_",bin[0],bin[1],bin[2]);
	ttt_->Sumw2();
	ttt_->SetFillColor(colors["TTT"]);
    
    TH1F * tt_ = new TH1F("tt_","tt_",bin[0],bin[1],bin[2]);
	tt_->Sumw2();
	tt_->SetFillColor(colors["TT"]);
	
	TH1F * vvj_ = new TH1F("vvj_","vvj_",bin[0],bin[1],bin[2]);
	vvj_->Sumw2();
	vvj_->SetFillColor(colors["VV"]);
    
    TH1F * vvt_ = new TH1F("vvt_","vvt_",bin[0],bin[1],bin[2]);
	vvt_->Sumw2();
	vvt_->SetFillColor(colors["VV"]);
    
    TH1F * zvv_ = new TH1F("zvv_","zvv_",bin[0],bin[1],bin[2]);
	zvv_->Sumw2();
	zvv_->SetFillColor(colors["ZVV"]);
    
    TH1F * ewk_ = new TH1F("ewk_","ewk_",bin[0],bin[1],bin[2]);
	ewk_->Sumw2();
	ewk_->SetFillColor(colors["EWK"]);
	
	TH1F * zhtautau_ = new TH1F("zhtautau_","zhtautau_",bin[0],bin[1],bin[2]);
	zhtautau_->Sumw2();
	zhtautau_->SetFillColor(colors["ZHTauTau"]);

	TH1F * glugluhtautau_ = new TH1F("glugluhtautau_","glugluhtautau_",bin[0],bin[1],bin[2]);
	glugluhtautau_->Sumw2();
	glugluhtautau_->SetFillColor(colors["GluGluHTauTau"]);
    
    TH1F * vbfhtautau_ = new TH1F("vbfhtautau_","vbfhtautau_",bin[0],bin[1],bin[2]);
	vbfhtautau_->Sumw2();
    vbfhtautau_->SetFillColor(colors["VBFHTauTau"]);
    
	TH1F * zp600a300_ = new TH1F("zp600a300_","zp600a300_",bin[0],bin[1],bin[2]);
	zp600a300_->Sumw2();
    TH1F * zp600a400_ = new TH1F("zp600a400_","zp600a400_",bin[0],bin[1],bin[2]);
	zp600a400_->Sumw2();

	TH1F * zp800a300_ = new TH1F("zp800a300_","zp800a300_",bin[0],bin[1],bin[2]);
	zp800a300_->Sumw2();
    TH1F * zp800a400_ = new TH1F("zp800a400_","zp800a400_",bin[0],bin[1],bin[2]);
	zp800a400_->Sumw2();
    TH1F * zp800a500_ = new TH1F("zp800a500_","zp800a500_",bin[0],bin[1],bin[2]);
	zp800a500_->Sumw2();
    TH1F * zp800a600_ = new TH1F("zp800a600_","zp800a600_",bin[0],bin[1],bin[2]);
	zp800a600_->Sumw2();

	TH1F * zp1000a300_ = new TH1F("zp1000a300_","zp1000a300_",bin[0],bin[1],bin[2]);
	zp1000a300_->Sumw2();
    TH1F * zp1000a400_ = new TH1F("zp1000a400_","zp1000a400_",bin[0],bin[1],bin[2]);
	zp1000a400_->Sumw2();
    TH1F * zp1000a500_ = new TH1F("zp1000a500_","zp1000a500_",bin[0],bin[1],bin[2]);
	zp1000a500_->Sumw2();
    TH1F * zp1000a600_ = new TH1F("zp1000a600_","zp1000a600_",bin[0],bin[1],bin[2]);
	zp1000a600_->Sumw2();
	TH1F * zp1000a700_ = new TH1F("zp1000a700_","zp1000a700_",bin[0],bin[1],bin[2]);
	zp1000a700_->Sumw2();
    TH1F * zp1000a800_ = new TH1F("zp1000a800_","zp1000a800_",bin[0],bin[1],bin[2]);
	zp1000a800_->Sumw2();
    
    TH1F * zp1200a300_ = new TH1F("zp1200a300_","zp1200a300_",bin[0],bin[1],bin[2]);
	zp1200a300_->Sumw2();
    TH1F * zp1200a400_ = new TH1F("zp1200a400_","zp1200a400_",bin[0],bin[1],bin[2]);
	zp1200a400_->Sumw2();
    TH1F * zp1200a500_ = new TH1F("zp1200a500_","zp1200a500_",bin[0],bin[1],bin[2]);
	zp1200a500_->Sumw2();
    TH1F * zp1200a600_ = new TH1F("zp1200a600_","zp1200a600_",bin[0],bin[1],bin[2]);
	zp1200a600_->Sumw2();
	TH1F * zp1200a700_ = new TH1F("zp1200a700_","zp1200a700_",bin[0],bin[1],bin[2]);
	zp1200a700_->Sumw2();
    TH1F * zp1200a800_ = new TH1F("zp1200a800_","zp1200a800_",bin[0],bin[1],bin[2]);
	zp1200a800_->Sumw2();
    
    TH1F * zp1400a300_ = new TH1F("zp1400a300_","zp1400a300_",bin[0],bin[1],bin[2]);
	zp1400a300_->Sumw2();
    TH1F * zp1400a400_ = new TH1F("zp1400a400_","zp1400a400_",bin[0],bin[1],bin[2]);
	zp1400a400_->Sumw2();
    TH1F * zp1400a500_ = new TH1F("zp1400a500_","zp1400a500_",bin[0],bin[1],bin[2]);
	zp1400a500_->Sumw2();
    TH1F * zp1400a600_ = new TH1F("zp1400a600_","zp1400a600_",bin[0],bin[1],bin[2]);
	zp1400a600_->Sumw2();
	TH1F * zp1400a700_ = new TH1F("zp1400a700_","zp1400a700_",bin[0],bin[1],bin[2]);
	zp1400a700_->Sumw2();
    TH1F * zp1400a800_ = new TH1F("zp1400a800_","zp1400a800_",bin[0],bin[1],bin[2]);
	zp1400a800_->Sumw2();
    
    TH1F * zp1700a300_ = new TH1F("zp1700a300_","zp1700a300_",bin[0],bin[1],bin[2]);
	zp1700a300_->Sumw2();
    TH1F * zp1700a400_ = new TH1F("zp1700a400_","zp1700a400_",bin[0],bin[1],bin[2]);
	zp1700a400_->Sumw2();
    TH1F * zp1700a500_ = new TH1F("zp1700a500_","zp1700a500_",bin[0],bin[1],bin[2]);
	zp1700a500_->Sumw2();
    TH1F * zp1700a600_ = new TH1F("zp1700a600_","zp1700a600_",bin[0],bin[1],bin[2]);
	zp1700a600_->Sumw2();
	TH1F * zp1700a700_ = new TH1F("zp1700a700_","zp1700a700_",bin[0],bin[1],bin[2]);
	zp1700a700_->Sumw2();
    TH1F * zp1700a800_ = new TH1F("zp1700a800_","zp1700a800_",bin[0],bin[1],bin[2]);
	zp1700a800_->Sumw2();
    
    TH1F * zp2000a300_ = new TH1F("zp2000a300_","zp2000a300_",bin[0],bin[1],bin[2]);
	zp2000a300_->Sumw2();
    TH1F * zp2000a400_ = new TH1F("zp2000a400_","zp2000a400_",bin[0],bin[1],bin[2]);
	zp2000a400_->Sumw2();
    TH1F * zp2000a500_ = new TH1F("zp2000a500_","zp2000a500_",bin[0],bin[1],bin[2]);
	zp2000a500_->Sumw2();
    TH1F * zp2000a600_ = new TH1F("zp2000a600_","zp2000a600_",bin[0],bin[1],bin[2]);
	zp2000a600_->Sumw2();
	TH1F * zp2000a700_ = new TH1F("zp2000a700_","zp2000a700_",bin[0],bin[1],bin[2]);
	zp2000a700_->Sumw2();
    TH1F * zp2000a800_ = new TH1F("zp2000a800_","zp2000a800_",bin[0],bin[1],bin[2]);
	zp2000a800_->Sumw2();
    
    TH1F * zp2500a300_ = new TH1F("zp2500a300_","zp2500a300_",bin[0],bin[1],bin[2]);
	zp2500a300_->Sumw2();
    TH1F * zp2500a400_ = new TH1F("zp2500a400_","zp2500a400_",bin[0],bin[1],bin[2]);
	zp2500a400_->Sumw2();
    TH1F * zp2500a500_ = new TH1F("zp2500a500_","zp2500a500_",bin[0],bin[1],bin[2]);
	zp2500a500_->Sumw2();
    TH1F * zp2500a600_ = new TH1F("zp2500a600_","zp2500a600_",bin[0],bin[1],bin[2]);
	zp2500a600_->Sumw2();
	TH1F * zp2500a700_ = new TH1F("zp2500a700_","zp2500a700_",bin[0],bin[1],bin[2]);
	zp2500a700_->Sumw2();
    TH1F * zp2500a800_ = new TH1F("zp2500a800_","zp2500a800_",bin[0],bin[1],bin[2]);
	zp2500a800_->Sumw2();

	/* sum of all bks */
	TH1F * comb_ = new TH1F("comb_","comb_",bin[0],bin[1],bin[2]); 
	comb_->Sumw2(); 

	TH1F * smh_ = new TH1F("smh_","smh_",bin[0],bin[1],bin[2]);
	smh_->Sumw2();
    smh_->SetFillColor(colors["SMH"]);

	W->Draw((parameter+">>w_").c_str(),relCut*weights["W"].c_str());
	ZTT->Draw((parameter+">>ztt_").c_str(),relCut*weights["ZTT"].c_str());
    ZL->Draw((parameter+">>zl_").c_str(),relCut*weights["ZL"].c_str());
    ZJ->Draw((parameter+">>zj_").c_str(),relCut*weights["ZJ"].c_str());
	VVT->Draw((parameter+">>vvt_").c_str(),relCut*weights["VV"].c_str());
    VVJ->Draw((parameter+">>vvj_").c_str(),relCut*weights["VV"].c_str());
    ZVV->Draw((parameter+">>zvv_").c_str(),cut*weights["ZVV"].c_str());
    EWK->Draw((parameter+">>ewk_").c_str(),cut*weights["EWK"].c_str());
	TTT->Draw((parameter+">>ttt_").c_str(),relCut*weights["TTT"].c_str());
    TTJ->Draw((parameter+">>ttj_").c_str(),relCut*weights["TTJ"].c_str());
	ZHTauTau->Draw((parameter+">>zhtautau_").c_str(),cut*weights["ZHTauTau"].c_str());
    VBFHTauTau->Draw((parameter+">>vbfhtautau_").c_str(),cut*weights["VBFHTauTau"].c_str());
    GluGluHTauTau->Draw((parameter+">>glugluhtautau_").c_str(),cut*weights["GluGluHTauTau"].c_str());
    
    smh_->Add(zhtautau_);
    smh_->Add(glugluhtautau_);
    smh_->Add(vbfhtautau_);
    
    tt_->Add(ttt_);
    tt_->Add(ttj_);
    
	MZP600_MA0300->Draw((parameter+">>zp600a300_").c_str(),cut*weights["MONO"].c_str());
    MZP600_MA0400->Draw((parameter+">>zp600a400_").c_str(),cut*weights["MONO"].c_str());
    
	MZP800_MA0300->Draw((parameter+">>zp800a300_").c_str(),cut*weights["MONO"].c_str());
    MZP800_MA0400->Draw((parameter+">>zp800a400_").c_str(),cut*weights["MONO"].c_str());
    MZP800_MA0500->Draw((parameter+">>zp800a500_").c_str(),cut*weights["MONO"].c_str());
    MZP800_MA0600->Draw((parameter+">>zp800a600_").c_str(),cut*weights["MONO"].c_str());
    
	MZP1000_MA0300->Draw((parameter+">>zp1000a300_").c_str(),cut*weights["MONO"].c_str());
    MZP1000_MA0400->Draw((parameter+">>zp1000a400_").c_str(),cut*weights["MONO"].c_str());
    MZP1000_MA0500->Draw((parameter+">>zp1000a500_").c_str(),cut*weights["MONO"].c_str());
    MZP1000_MA0600->Draw((parameter+">>zp1000a600_").c_str(),cut*weights["MONO"].c_str());
    MZP1000_MA0700->Draw((parameter+">>zp1000a700_").c_str(),cut*weights["MONO"].c_str());
    MZP1000_MA0800->Draw((parameter+">>zp1000a800_").c_str(),cut*weights["MONO"].c_str());
    
    MZP1200_MA0300->Draw((parameter+">>zp1200a300_").c_str(),cut*weights["MONO"].c_str());
    MZP1200_MA0400->Draw((parameter+">>zp1200a400_").c_str(),cut*weights["MONO"].c_str());
    MZP1200_MA0500->Draw((parameter+">>zp1200a500_").c_str(),cut*weights["MONO"].c_str());
    MZP1200_MA0600->Draw((parameter+">>zp1200a600_").c_str(),cut*weights["MONO"].c_str());
    MZP1200_MA0700->Draw((parameter+">>zp1200a700_").c_str(),cut*weights["MONO"].c_str());
    MZP1200_MA0800->Draw((parameter+">>zp1200a800_").c_str(),cut*weights["MONO"].c_str());

	MZP1400_MA0300->Draw((parameter+">>zp1400a300_").c_str(),cut*weights["MONO"].c_str());
    MZP1400_MA0400->Draw((parameter+">>zp1400a400_").c_str(),cut*weights["MONO"].c_str());
    MZP1400_MA0500->Draw((parameter+">>zp1400a500_").c_str(),cut*weights["MONO"].c_str());
    MZP1400_MA0600->Draw((parameter+">>zp1400a600_").c_str(),cut*weights["MONO"].c_str());
    MZP1400_MA0700->Draw((parameter+">>zp1400a700_").c_str(),cut*weights["MONO"].c_str());
    MZP1400_MA0800->Draw((parameter+">>zp1400a800_").c_str(),cut*weights["MONO"].c_str());
    
	MZP1700_MA0300->Draw((parameter+">>zp1700a300_").c_str(),cut*weights["MONO"].c_str());
    MZP1700_MA0400->Draw((parameter+">>zp1700a400_").c_str(),cut*weights["MONO"].c_str());
    MZP1700_MA0500->Draw((parameter+">>zp1700a500_").c_str(),cut*weights["MONO"].c_str());
    MZP1700_MA0600->Draw((parameter+">>zp1700a600_").c_str(),cut*weights["MONO"].c_str());
    MZP1700_MA0700->Draw((parameter+">>zp1700a700_").c_str(),cut*weights["MONO"].c_str());
    MZP1700_MA0800->Draw((parameter+">>zp1700a800_").c_str(),cut*weights["MONO"].c_str());
    
	MZP2000_MA0300->Draw((parameter+">>zp2000a300_").c_str(),cut*weights["MONO"].c_str());
    MZP2000_MA0400->Draw((parameter+">>zp2000a400_").c_str(),cut*weights["MONO"].c_str());
    MZP2000_MA0500->Draw((parameter+">>zp2000a500_").c_str(),cut*weights["MONO"].c_str());
    MZP2000_MA0600->Draw((parameter+">>zp2000a600_").c_str(),cut*weights["MONO"].c_str());
    MZP2000_MA0700->Draw((parameter+">>zp2000a700_").c_str(),cut*weights["MONO"].c_str());
    MZP2000_MA0800->Draw((parameter+">>zp2000a800_").c_str(),cut*weights["MONO"].c_str());
    
	MZP2500_MA0300->Draw((parameter+">>zp2500a300_").c_str(),cut*weights["MONO"].c_str());
    MZP2500_MA0400->Draw((parameter+">>zp2500a400_").c_str(),cut*weights["MONO"].c_str());
    MZP2500_MA0500->Draw((parameter+">>zp2500a500_").c_str(),cut*weights["MONO"].c_str());
    MZP2500_MA0600->Draw((parameter+">>zp2500a600_").c_str(),cut*weights["MONO"].c_str());
    MZP2500_MA0700->Draw((parameter+">>zp2500a700_").c_str(),cut*weights["MONO"].c_str());
    MZP2500_MA0800->Draw((parameter+">>zp2500a800_").c_str(),cut*weights["MONO"].c_str());
    
	fillQCD_Shape(SScut, SScutRel, sfQCD, qcd_, parameter, bin, "QCD SS " + parameter + " sideband in " + chan + " channel", syst, chan, 0, 0);
	qcd_->Scale(qcdOStoSS);
    
    //not scaling in datacards
	//w_->Scale(sf);
    w_->Scale(shapeScale(cut, relCut, W, weights["W"], parameter, bin));
    vvt_->Scale(shapeScale(cut, relCut, VVT, weights["VV"], parameter, bin));
    vvj_->Scale(shapeScale(cut, relCut, VVJ, weights["VV"], parameter, bin));
    ttt_->Scale(shapeScale(cut, relCut, TTT, weights["TTT"], parameter, bin));
    ttj_->Scale(shapeScale(cut, relCut, TTJ, weights["TTJ"], parameter, bin));
    ztt_->Scale(shapeScale(cut, relCut, ZTT, weights["ZTT"], parameter, bin));
    zl_->Scale(shapeScale(cut, relCut, ZL, weights["ZL"], parameter, bin));
    zj_->Scale(shapeScale(cut, relCut, ZJ, weights["ZJ"], parameter, bin));

    comb_->Add(vvt_);
    comb_->Add(vvj_);
    comb_->Add(zvv_);
    comb_->Add(ewk_);
	comb_->Add(ttt_);
    comb_->Add(ttj_);
	comb_->Add(w_);
	comb_->Add(ztt_);
    comb_->Add(zj_);
    comb_->Add(zl_);
    comb_->Add(zhtautau_);
    comb_->Add(glugluhtautau_);
    comb_->Add(vbfhtautau_);
	comb_->Add(qcd_);
    
	// scale things to the projected lumi
	vvt_->Scale(lumi_sf);
    vvj_->Scale(lumi_sf);
    zvv_->Scale(lumi_sf);
    ewk_->Scale(lumi_sf);
	ttt_->Scale(lumi_sf);
    ttj_->Scale(lumi_sf);
	ztt_->Scale(lumi_sf);
    zl_->Scale(lumi_sf);
    zj_->Scale(lumi_sf);
	w_->Scale(lumi_sf);
	qcd_->Scale(lumi_sf);
    zhtautau_->Scale(lumi_sf);
    glugluhtautau_->Scale(lumi_sf);
    vbfhtautau_->Scale(lumi_sf);
    smh_->Scale(lumi_sf);
    
	zp600a300_->Scale(lumi_sf);
    zp600a400_->Scale(lumi_sf);
	zp800a300_->Scale(lumi_sf);
    zp800a400_->Scale(lumi_sf);
    zp800a500_->Scale(lumi_sf);
    zp800a600_->Scale(lumi_sf);
	zp1000a300_->Scale(lumi_sf);
    zp1000a400_->Scale(lumi_sf);
    zp1000a500_->Scale(lumi_sf);
    zp1000a600_->Scale(lumi_sf);
    zp1000a700_->Scale(lumi_sf);
    zp1000a800_->Scale(lumi_sf);
	zp1200a300_->Scale(lumi_sf);
    zp1200a400_->Scale(lumi_sf);
    zp1200a500_->Scale(lumi_sf);
    zp1200a600_->Scale(lumi_sf);
    zp1200a700_->Scale(lumi_sf);
    zp1200a800_->Scale(lumi_sf);
    zp1400a300_->Scale(lumi_sf);
    zp1400a400_->Scale(lumi_sf);
    zp1400a500_->Scale(lumi_sf);
    zp1400a600_->Scale(lumi_sf);
    zp1400a700_->Scale(lumi_sf);
    zp1400a800_->Scale(lumi_sf);
	zp1700a300_->Scale(lumi_sf);
    zp1700a400_->Scale(lumi_sf);
    zp1700a500_->Scale(lumi_sf);
    zp1700a600_->Scale(lumi_sf);
    zp1700a700_->Scale(lumi_sf);
    zp1700a800_->Scale(lumi_sf);
	zp2000a300_->Scale(lumi_sf);
    zp2000a400_->Scale(lumi_sf);
    zp2000a500_->Scale(lumi_sf);
    zp2000a600_->Scale(lumi_sf);
    zp2000a700_->Scale(lumi_sf);
    zp2000a800_->Scale(lumi_sf);
	zp2500a300_->Scale(lumi_sf);
    zp2500a400_->Scale(lumi_sf);
    zp2500a500_->Scale(lumi_sf);
    zp2500a600_->Scale(lumi_sf);
    zp2500a700_->Scale(lumi_sf);
    zp2500a800_->Scale(lumi_sf);
    
	comb_->Scale(lumi_sf);
    
    // print signal region bkg yields (under whatever cut is provided)

    std::cout<<" bkg totals under "<<can_name<<"\n";

    double error = 0;

    error = 0;
    vvt_->IntegralAndError(-1,1e3, error);
    std::cout<<" VVT ... "<<vvt_->Integral()<<" +/- "<<error<<"\n";
    
    error = 0;
    vvj_->IntegralAndError(-1,1e3, error);
    std::cout<<" VVJ ... "<<vvj_->Integral()<<" +/- "<<error<<"\n";
    
    error = 0;
    zvv_->IntegralAndError(-1,1e3, error);
    std::cout<<" ZVV ... "<< zvv_->Integral()<<" +/- "<<error<<"\n";
    
    error = 0;
    ewk_->IntegralAndError(-1,1e3, error);
    std::cout<<" EWK ... "<< ewk_->Integral()<<" +/- "<<error<<"\n";

    error = 0;
    ttt_->IntegralAndError(-1,1e3, error);
    std::cout<<" TTT ... "<<ttt_->Integral()<<" +/- "<<error<<"\n";
    
    error = 0;
    ttj_->IntegralAndError(-1,1e3, error);
    std::cout<<" TTJ ... "<<ttj_->Integral()<<" +/- "<<error<<"\n";

    error = 0;
    ztt_->IntegralAndError(-1,1e3, error);
    std::cout<<" ZTT ... "<<ztt_->Integral()<<" +/- "<<error<<"\n";

    error = 0;
    zl_->IntegralAndError(-1,1e3, error);
    std::cout<<" ZL ... "<<zl_->Integral()<<" +/- "<<error<<"\n";
    
    error = 0;
    zj_->IntegralAndError(-1,1e3, error);
    std::cout<<" ZJ ... "<<zj_->Integral()<<" +/- "<<error<<"\n";

    error = 0;
    w_->IntegralAndError(-1,1e3, error);
    std::cout<<" W ... "<<w_->Integral()<<" +/- "<<error<<"\n";
    
    error = 0;
    zhtautau_->IntegralAndError(-1,1e3, error);
    std::cout<<" ZHTauTau ... "<<zhtautau_->Integral()<<" +/- "<<error<<"\n";

    error = 0;
    glugluhtautau_->IntegralAndError(-1,1e3, error);
    std::cout<<" GluGluHTauTau ... "<<glugluhtautau_->Integral()<<" +/- "<<error<<"\n";
    
    error = 0;
    vbfhtautau_->IntegralAndError(-1,1e3, error);
    std::cout<<" VBFHTauTau ... "<<vbfhtautau_->Integral()<<" +/- "<<error<<"\n";

    error = 0;
    qcd_->IntegralAndError(-1,1e3, error);
    error = sqrt(qcd_->GetEntries());
    std::cout<<" QCD ... "<<qcd_->Integral()<<" +/- "<<error<<"\n";

    error = 0;
    zp600a300_->IntegralAndError(-1,1e3, error);
    std::cout<<" ZP600 ... "<<zp600a300_->Integral()<<" +/- "<<error<<"\n";

    error = 0;
    zp800a300_->IntegralAndError(-1,1e3, error);
    std::cout<<" ZP800 ... "<<zp800a300_->Integral()<<" +/- "<<error<<"\n";

    error = 0;
    zp1000a300_->IntegralAndError(-1,1e3, error);
    std::cout<<" ZP1000 ... "<<zp1000a300_->Integral()<<" +/- "<<error<<"\n";

    error = 0;
    zp1200a300_->IntegralAndError(-1,1e3, error);
    std::cout<<" ZP1200 ... "<<zp1200a300_->Integral()<<" +/- "<<error<<"\n";

    error = 0;
    zp1400a300_->IntegralAndError(-1,1e3, error);
    std::cout<<" ZP1400 ... "<<zp1400a300_->Integral()<<" +/- "<<error<<"\n";
    
    error = 0;
    zp1700a300_->IntegralAndError(-1,1e3, error);
    std::cout<<" ZP1700 ... "<<zp1700a300_->Integral()<<" +/- "<<error<<"\n";
    
    error = 0;
    zp2000a300_->IntegralAndError(-1,1e3, error);
    std::cout<<" ZP2000 ... "<<zp2000a300_->Integral()<<" +/- "<<error<<"\n";
    
    error = 0;
    zp2500a300_->IntegralAndError(-1,1e3, error);
    std::cout<<" ZP2500 ... "<<zp2500a300_->Integral()<<" +/- "<<error<<"\n";

    error = 0;
    comb_->IntegralAndError(-1,1e3, error);
    std::cout<<" totalBkg ... "<<comb_->Integral()<<" +/- "<<error<<"\n";

    // test

    double mvaChoiceI;
    double mvaChoiceE;

    double wI, wE;
    double vvtI, vvtE;
    double vvjI, vvjE;
    double zvvI, zvvE;
    double ewkI, ewkE;
    double zttI, zttE;
    double zlI, zlE;
    double zjI, zjE;
    double qcdI, qcdE;
    double tttI, tttE;
    double ttjI, ttjE;
    double zhtautauI, zhtautauE;
    double vbfhtautauI, vbfhtautauE;
    double glugluhtautauI, glugluhtautauE;

    double zp600I, zp600E;
    double zp800I, zp800E;
    double zp1000I, zp1000E;
    double zp1200I, zp1200E;
    double zp1400I, zp1400E;
    double zp1700I, zp1700E;
    double zp2000I, zp2000E;
    double zp2500I, zp2500E;

    zp600I = zp600a300_->Integral(-1, zp600a300_->GetNbinsX()+1);
    zp600a300_->IntegralAndError(-1, zp600a300_->GetNbinsX()+1,zp600E);

    zp800I = zp800a300_->Integral(-1, zp800a300_->GetNbinsX()+1);
    zp800a300_->IntegralAndError(-1, zp800a300_->GetNbinsX()+1,zp800E);

    zp1000I = zp1000a300_->Integral(-1, zp1000a300_->GetNbinsX()+1);
    zp1000a300_->IntegralAndError(-1, zp1000a300_->GetNbinsX()+1,zp1000E);

    zp1200I = zp1200a300_->Integral(-1, zp1200a300_->GetNbinsX()+1);
    zp1200a300_->IntegralAndError(-1, zp1200a300_->GetNbinsX()+1,zp1200E);

    zp1400I = zp1400a300_->Integral(-1, zp1400a300_->GetNbinsX()+1);
    zp1400a300_->IntegralAndError(-1, zp1400a300_->GetNbinsX()+1,zp1400E);
    
    zp1700I = zp1700a300_->Integral(-1, zp1700a300_->GetNbinsX()+1);
    zp1700a300_->IntegralAndError(-1, zp1700a300_->GetNbinsX()+1,zp1700E);
    
    zp2000I = zp2000a300_->Integral(-1, zp2000a300_->GetNbinsX()+1);
    zp2000a300_->IntegralAndError(-1, zp2000a300_->GetNbinsX()+1,zp2000E);
    
    zp2500I = zp2500a300_->Integral(-1, zp2500a300_->GetNbinsX()+1);
    zp2500a300_->IntegralAndError(-1, zp2500a300_->GetNbinsX()+1,zp2500E);

    wI = w_->Integral(-1, w_->GetNbinsX()+1);
    w_->IntegralAndError(-1, w_->GetNbinsX()+1,wE);

    vvtI = vvt_->Integral(-1, vvt_->GetNbinsX()+1);
    vvt_->IntegralAndError(-1, vvt_->GetNbinsX()+1,vvtE);
    
    vvjI = vvj_->Integral(-1, vvj_->GetNbinsX()+1);
    vvj_->IntegralAndError(-1, vvj_->GetNbinsX()+1,vvjE);
    
    zvvI = zvv_->Integral(-1, zvv_->GetNbinsX()+1);
    zvv_->IntegralAndError(-1, zvv_->GetNbinsX()+1,zvvE);
    
    ewkI = ewk_->Integral(-1, ewk_->GetNbinsX()+1);
    ewk_->IntegralAndError(-1, ewk_->GetNbinsX()+1,ewkE);

    zttI = ztt_->Integral(-1, ztt_->GetNbinsX()+1);
    ztt_->IntegralAndError(-1, ztt_->GetNbinsX()+1,zttE);
    
    zlI = zl_->Integral(-1, zl_->GetNbinsX()+1);
    zl_->IntegralAndError(-1, zl_->GetNbinsX()+1,zlE);
    
    zjI = zj_->Integral(-1, zj_->GetNbinsX()+1);
    zj_->IntegralAndError(-1, zj_->GetNbinsX()+1,zjE);

    qcdI = qcd_->Integral(-1, qcd_->GetNbinsX()+1);
    qcd_->IntegralAndError(-1, qcd_->GetNbinsX()+1,qcdE);

    tttI = ttt_->Integral(-1, ttt_->GetNbinsX()+1);
    ttt_->IntegralAndError(-1, ttt_->GetNbinsX()+1,tttE);
    
    ttjI = ttj_->Integral(-1, ttj_->GetNbinsX()+1);
    ttj_->IntegralAndError(-1, ttj_->GetNbinsX()+1,ttjE);

    zhtautauI = zhtautau_->Integral(-1, zhtautau_->GetNbinsX()+1);
    zhtautau_->IntegralAndError(-1, zhtautau_->GetNbinsX()+1,zhtautauE);
    
    glugluhtautauI = glugluhtautau_->Integral(-1, glugluhtautau_->GetNbinsX()+1);
    glugluhtautau_->IntegralAndError(-1, glugluhtautau_->GetNbinsX()+1,glugluhtautauE);
    
    vbfhtautauI = vbfhtautau_->Integral(-1, vbfhtautau_->GetNbinsX()+1);
    vbfhtautau_->IntegralAndError(-1, vbfhtautau_->GetNbinsX()+1,vbfhtautauE);

    if(wI>0.0){ wE = 1 + wE/wI; } else {wI = 0.0; wE = 1.0;}
    if(vvtI>0.0){ vvtE = 1 + vvtE/vvtI; } else {vvtI = 0.0; vvtE = 1.0;}
    if(vvjI>0.0){ vvjE = 1 + vvjE/vvjI; } else {vvjI = 0.0; vvjE = 1.0;}
    if(zvvI>0.0){ zvvE = 1 + zvvE/zvvI; } else {zvvI = 0.0; zvvE = 1.0;}
    if(ewkI>0.0){ ewkE = 1 + ewkE/ewkI; } else {ewkI = 0.0; ewkE = 1.0;}
    if(zttI>0.0){ zttE = 1 + zttE/zttI; } else {zttI = 0.0; zttE = 1.0;}
    if(zlI>0.0){ zlE = 1 + zlE/zlI; } else {zlI = 0.0; zlE = 1.0;}
    if(zjI>0.0){ zjE = 1 + zjE/zjI; } else {zjI = 0.0; zjE = 1.0;}
    if(qcdI>0.0){ qcdE = 1 + qcdE/qcdI; } else {qcdI = 0.0; qcdE = 1.0;}
    if(tttI>0.0){ tttE = 1 + tttE/tttI; } else {tttI = 0.0; tttE = 1.0;}
    if(ttjI>0.0){ ttjE = 1 + ttjE/ttjI; } else {ttjI = 0.0; ttjE = 1.0;}
    if(zhtautauI>0.0){ zhtautauE = 1 + zhtautauE/zhtautauI; } else {zhtautauI = 0.0; zhtautauE = 1.0;}
    if(vbfhtautauI>0.0){ vbfhtautauE = 1 + vbfhtautauE/vbfhtautauI; } else {vbfhtautauI = 0.0; vbfhtautauE = 1.0;}
    if(glugluhtautauI>0.0){ glugluhtautauE = 1 + glugluhtautauE/glugluhtautauI; } else {glugluhtautauI = 0.0; glugluhtautauE = 1.0;}
    if(zp600I>0.0){ zp600E = 1 + zp600E/zp600I; } else {zp600I = 0.0; zp600E = 1.0;}
    if(zp800I>0.0){ zp800E = 1 + zp800E/zp800I; } else {zp800I = 0.0; zp800E = 1.0;}
    if(zp1000I>0.0){ zp1000E = 1 + zp1000E/zp1000I; } else {zp1000I = 0.0; zp1000E = 1.0;}
    if(zp1200I>0.0){ zp1200E = 1 + zp1200E/zp1200I; } else {zp1200I = 0.0; zp1200E = 1.0;}
    if(zp1400I>0.0){ zp1400E = 1 + zp1400E/zp1400I; } else {zp1400I = 0.0; zp1400E = 1.0;}
    if(zp1700I>0.0){ zp1700E = 1 + zp1700E/zp1700I; } else {zp1700I = 0.0; zp1700E = 1.0;}
    if(zp2000I>0.0){ zp2000E = 1 + zp2000E/zp2000I; } else {zp2000I = 0.0; zp2000E = 1.0;}
    if(zp2500I>0.0){ zp2500E = 1 + zp2500E/zp2500I; } else {zp2500I = 0.0; zp2500E = 1.0;}
    if(comb_->Integral()!=0) std::cout<<" @ ZP = 600 s/rt(b) = "<<zp600I/sqrt(comb_->Integral())<<"\n";
    if(comb_->Integral()!=0) std::cout<<" @ ZP = 800 s/rt(b) = "<<zp800I/sqrt(comb_->Integral())<<"\n";
    if(comb_->Integral()!=0) std::cout<<" @ ZP = 1000 s/rt(b) = "<<zp1000I/sqrt(comb_->Integral())<<"\n";
    if(comb_->Integral()!=0) std::cout<<" @ ZP = 1200 s/rt(b) = "<<zp1200I/sqrt(comb_->Integral())<<"\n";
    if(comb_->Integral()!=0) std::cout<<" @ ZP = 1400 s/rt(b) = "<<zp1400I/sqrt(comb_->Integral())<<"\n";
    if(comb_->Integral()!=0) std::cout<<" @ ZP = 1700 s/rt(b) = "<<zp1700I/sqrt(comb_->Integral())<<"\n";
    if(comb_->Integral()!=0) std::cout<<" @ ZP = 2000 s/rt(b) = "<<zp2000I/sqrt(comb_->Integral())<<"\n";
    if(comb_->Integral()!=0) std::cout<<" @ ZP = 2500 s/rt(b) = "<<zp2500I/sqrt(comb_->Integral())<<"\n";

    //out_data
    if (choice2==1 && sigChoice==2 && (chan=="mt" || chan=="et"))
    {
        if (parameter.substr(15,4)=="600A"){mvaChoiceI = zp600I; mvaChoiceE = zp600E; }
        if (parameter.substr(15,4)=="800A"){mvaChoiceI = zp800I; mvaChoiceE = zp800E; }
        if (parameter.substr(15,4)=="1000"){mvaChoiceI = zp1000I; mvaChoiceE = zp1000E; }
        if (parameter.substr(15,4)=="1200"){mvaChoiceI = zp1200I; mvaChoiceE = zp1200E; }
        
        out_data << "QUANTITY       SIGNAL     W     VVT     VVJ     ZVV     EWK     ZTT     ZJ     ZL     QCD     TTT     TTJ     ZHTauTau     GluGluHTauTau     VBFHTauTau" << std::endl;
        
        out_data << drCutNameMap[drCut] + chan + parameter.substr(10,16) + " rate," << mvaChoiceI << "," << wI << "," << vvtI << "," << vvjI << "," << zvvI << "," << ewkI << "," << zttI << "," << zjI << "," << zlI << "," << qcdI << "," << tttI << "," << ttjI << "," << zhtautauI << "," << glugluhtautauI << "," << vbfhtautauI << std::endl;
        
        out_data << drCutNameMap[drCut] + chan + parameter.substr(10,16) + " stat.err    " << mvaChoiceE << "," << wE << "," << vvtE << "," << vvjE << "," << zvvI << "," << ewkI << "," << zttE << "," << zjE << "," << zlE << "," << qcdE << "," << tttE << "," << ttjE << "," << zhtautauE << "," << glugluhtautauE << "," << vbfhtautauE << std::endl;
    }
    else
    {
    
        out_data << "QUANTITY       SIGNAL     W     VVT     VVJ     ZVV     EWK     ZTT     ZJ     ZL     QCD     TTT     TTJ     ZHTauTau     GluGluHTauTau     VBFHTauTau" << std::endl;
        
        out_data << drCutNameMap[drCut] + chan + "ZP600 rate," << zp600I << "," << wI << "," << vvtI << "," << vvjI << "," << zvvI << "," << ewkI << "," << zttI << "," << zjI << "," << zlI << "," << qcdI << "," << tttI << "," << ttjI << "," << zhtautauI << "," << glugluhtautauI << "," << vbfhtautauI << std::endl;
        out_data << drCutNameMap[drCut] + chan + "ZP600 stat.err    " << zp600E << "," << wE << "," << vvtE << "," << vvjE << "," << zvvE << "," << ewkE << "," << zttE << "," << zjE << "," << zlE << "," << qcdE << "," << tttE << "," << ttjE << "," << zhtautauE << "," << glugluhtautauE << "," << vbfhtautauE << std::endl;
        
        out_data << drCutNameMap[drCut] + chan + "ZP800 rate," << zp800I << "," << wI << "," << vvtI << "," << vvjI << "," << zvvI << "," << ewkI << "," << zttI << "," << zjI << "," << zlI << "," << qcdI << "," << tttI << "," << ttjI << "," << zhtautauI << "," << glugluhtautauI << "," << vbfhtautauI << std::endl;
        out_data << drCutNameMap[drCut] + chan + "ZP800 stat.err    " << zp800E << "," << wE << "," << vvtE << "," << vvjE << "," << zvvE << "," << ewkE << "," << zttE << "," << zjE << "," << zlE << "," << qcdE << "," << tttE << "," << ttjE << "," << zhtautauE << "," << glugluhtautauE << "," << vbfhtautauE << std::endl;
        
        out_data << drCutNameMap[drCut] + chan + "ZP1000 rate," << zp1000I << "," << wI << "," << vvtI << "," << vvjI << "," << zvvI << "," << ewkI << "," << zttI << "," << zjI << "," << zlI << "," << qcdI << "," << tttI << "," << ttjI << "," << zhtautauI << "," << glugluhtautauI << "," << vbfhtautauI << std::endl;
        out_data << drCutNameMap[drCut] + chan + "ZP1000 stat.err    " << zp1000E << "," << wE << "," << vvtE << "," << vvjE << "," << zvvE << "," << ewkE << "," << zttE << "," << zjE << "," << zlE << "," << qcdE << "," << tttE << "," << ttjE << "," << zhtautauE << "," << glugluhtautauE << "," << vbfhtautauE << std::endl;
        
        out_data << drCutNameMap[drCut] + chan + "ZP1200 rate," << zp1200I << "," << wI << "," << vvtI << "," << vvjI << "," << zvvI << "," << ewkI << "," << zttI << "," << zjI << "," << zlI << "," << qcdI << "," << tttI << "," << ttjI << "," << zhtautauI << "," << glugluhtautauI << "," << vbfhtautauI << std::endl;
        out_data << drCutNameMap[drCut] + chan + "ZP1200 stat.err    " << zp1200E << "," << wE << "," << vvtE << "," << vvjE << "," << zvvE << "," << ewkE << "," << zttE << "," << zjE << "," << zlE << "," << qcdE << "," << tttE << "," << ttjE << "," << zhtautauE << "," << glugluhtautauE << "," << vbfhtautauE << std::endl;
    }
}

void setBkgErr(std::string chan)
{
    if(chan == "et")
    {
        vvtError = 0.10;
        vvjError = 0.10;
        wError = 0.32;
        zttError = 0.12;
        zlError = 0.18;
        zjError = 0.16;
        zvvError = 0.065;
        ewkError = 0.094;
        tttError = 0.16;
        ttjError = 0.17;
        smhError = 0.086;
        qcdError = 0.31;
    }
    else if (chan == "mt")
    {
        vvtError = 0.097;
        vvjError = 0.097;
        wError = 0.32;
        zttError = 0.012;
        zlError = 0.27;
        zjError = 0.16;
        zvvError = 0.065;
        ewkError = 0.12;
        tttError = 0.17;
        ttjError = 0.17;
        smhError = 0.094;
        qcdError = 0.31;
    }
    else if (chan == "tt")
    {
        vvtError = 0.13;
        vvjError = 0.13;
        wError = 0.062;
        zttError = 0.17;
        zlError = 0.10;
        zjError = 0.18;
        zvvError = 0.062;
        ewkError = 0.17;
        tttError = 0.21;
        ttjError = 0.20;
        smhError = 0.12;
        qcdError = 0.31;
    }
}

void drawTcanvasAndLegendAndRatio(TCanvas * C, TH1F * num__, TH1F* den__, bool drawSignals)
{

	C->cd(1);

	gPad->SetGridx(0);
	gPad->SetGridy(0);
	gPad->SetLogy(setLog);

	TLegend leg(.58,.50,0.90,.90);
    leg.SetTextSize(.04);
    leg.SetTextFont(22);
    leg.SetBorderSize(0);
    leg.SetFillColor(0);

	TH1F * data__ = new TH1F("data__","data__",mtBinning[0],mtBinning[1],mtBinning[2]); 
	data__->SetMarkerStyle(dataMarker);

	TH1F * ztt__ = new TH1F("ztt__","ztt__",mtBinning[0],mtBinning[1],mtBinning[2]);
	ztt__->SetFillColor(colors["ZTT"]);
    
    TH1F * zl__ = new TH1F("zl__","zl__",mtBinning[0],mtBinning[1],mtBinning[2]);
	zl__->SetFillColor(colors["ZL"]);
    
    TH1F * zj__ = new TH1F("zj__","zj__",mtBinning[0],mtBinning[1],mtBinning[2]);
	zj__->SetFillColor(colors["ZJ"]);
	
	TH1F * w__ = new TH1F("w__","w__",mtBinning[0],mtBinning[1],mtBinning[2]); 
	w__->SetFillColor(colors["W"]);
	
	TH1F * tt__ = new TH1F("tt__","tt__",mtBinning[0],mtBinning[1],mtBinning[2]); 
	tt__->SetFillColor(colors["TTT"]);
	
	TH1F * vv__ = new TH1F("vv__","vv__",mtBinning[0],mtBinning[1],mtBinning[2]); 
	vv__->SetFillColor(colors["VV"]);
    
    TH1F * zvv__ = new TH1F("zvv__","zvv__",mtBinning[0],mtBinning[1],mtBinning[2]);
	zvv__->SetFillColor(colors["ZVV"]);
    
    TH1F * ewk__ = new TH1F("ewk__","ewk__",mtBinning[0],mtBinning[1],mtBinning[2]);
	ewk__->SetFillColor(colors["EWK"]);
	
	TH1F * zhtautau__ = new TH1F("zhtautau__","zhtautau__",mtBinning[0],mtBinning[1],mtBinning[2]);
	zhtautau__->SetFillColor(colors["ZHTauTau"]);

    TH1F * glugluhtautau__ = new TH1F("glugluhtautau__","glugluhtautau__",mtBinning[0],mtBinning[1],mtBinning[2]);
	glugluhtautau__->SetFillColor(colors["GluGluHTauTau"]);
    
    TH1F * vbfhtautau__ = new TH1F("vbfhtautau__","vbfhtautau__",mtBinning[0],mtBinning[1],mtBinning[2]);
	vbfhtautau__->SetFillColor(colors["VBFHTauTau"]);

    TH1F * smh__ = new TH1F("smh__","smh__",mtBinning[0],mtBinning[1],mtBinning[2]);
    smh__->SetFillColor(colors["SMH"]);
    
	TH1F * qcd__ = new TH1F("qcd__","qcd__",mtBinning[0],mtBinning[1],mtBinning[2]); 
	qcd__->SetFillColor(colors["QCD"]);

	TH1F * zp600__ = new TH1F("zp600__","zp600__",mtBinning[0],mtBinning[1],mtBinning[2]);
	TH1F * zp1200__ = new TH1F("zp1200__","zp1200__",mtBinning[0],mtBinning[1],mtBinning[2]);
    TH1F * zp1700__ = new TH1F("zp1700__","zp1700__",mtBinning[0],mtBinning[1],mtBinning[2]);
    TH1F * zp2000__ = new TH1F("zp2000__","zp2000__",mtBinning[0],mtBinning[1],mtBinning[2]);
    TH1F * zp2500__ = new TH1F("zp2500__","zp2500__",mtBinning[0],mtBinning[1],mtBinning[2]);
    
    zp600__->SetLineStyle(2);
	zp1200__->SetLineStyle(2);
	zp1700__->SetLineStyle(2);
	zp2000__->SetLineStyle(2);
	zp2500__->SetLineStyle(2);
    
    zp600__->SetLineWidth(3);
	zp1200__->SetLineWidth(3);
	zp1700__->SetLineWidth(3);
	zp2000__->SetLineWidth(3);
	zp2500__->SetLineWidth(3);

	zp600__->SetLineColor(50);
	zp1200__->SetLineColor(60);
	zp1700__->SetLineColor(70);
    zp2000__->SetLineColor(80);
	zp2500__->SetLineColor(6);

    
    leg.AddEntry(data__, "Observed","PE");
    leg.AddEntry(ztt__, "ZTT","F");
    leg.AddEntry(zl__, "ZL","F");
    leg.AddEntry(zj__, "ZJ","F");
    leg.AddEntry(w__, "W+jets","F");
    leg.AddEntry(vv__, "VV, SingleTop","F");
    if (splitPlotSMH)
    {
        leg.AddEntry(zhtautau__, "ZHTauTau(125 GeV)","F");
        leg.AddEntry(glugluhtautau__, "GluGluHTauTau(125 GeV)","F");
        leg.AddEntry(vbfhtautau__, "VBFHTauTau(125 GeV)","F");
    }
    else {leg.AddEntry(smh__, "SM Higgs: ggH, vbfH, ZH","F");}
    leg.AddEntry(qcd__, "QCD","F");
    leg.AddEntry(tt__, "t#bar{t}","F");
    leg.AddEntry(zvv__, "Z invisible","F");
    leg.AddEntry(ewk__, "EWK","F");
    
    if(shouldIplotSignals && drawSignals)
	{

        leg.AddEntry(zp600__, "[m_{Z'},m_{A0},m_{#chi}] = [600,300,100] GeV @ 1pb","L");
        leg.AddEntry(zp1200__, "[m_{Z'},m_{A0},m_{#chi}] = [1200,300,100] GeV @ 1pb","L");
        //legDM.AddEntry(zp1700__, "[m_{Z'},m_{A0},m_{#chi}] = [1700,300,100] GeV @ 1pb","L");
        //legDM.AddEntry(zp2000__, "[m_{Z'},m_{A0},m_{#chi}] = [2000,300,100] GeV @ 1pb","L");
        //legDM.AddEntry(zp2500__, "[m_{Z'},m_{A0},m_{#chi}] = [2500,300,100] GeV @ 1pb","L");
	}

    leg.DrawClone();

    TLatex L1;
    L1.DrawLatexNDC(0.12,0.83,"CMS (#it{Preliminary})");

    TLatex L2;
    //L2.SetFont(42);
    L2.SetTextSize(0.06);
    L2.DrawLatexNDC(0.68,0.925,("#font[42]{#bf{"+lumiAndRtS+"}}").c_str());

    TLatex L3;
    L3.SetTextColor(4);
    L3.SetTextSize(1.2 * L1.GetTextSize());

    L3.DrawLatexNDC(0.12,0.76,chan_label.c_str());

    TLatex L4;
    //L2.SetFont(42);
    L4.SetTextSize(0.06);
    L4.DrawLatexNDC(0.15,0.925,("#font[42]{#bf{"+global_title+"}}").c_str());

    C->cd(2);

    TH1F * ratio__ = (TH1F*)num__->Clone("ratio__");
    ratio__ -> Reset();
    
    //Set denom for data ratio with sqrtN as error
    TH1F * mc = (TH1F*)den__->Clone("mc");
    for(int i=1;i<mc->GetNbinsX();i++)
    {
        float y = TMath::Sqrt(mc->GetBinContent(i));
        mc->SetBinError(i,y);
    }
    //ratio__ -> Divide(num__, mc, 1, 1, "B");
    ratio__ -> Divide(num__, mc);

    ratio__ -> SetTitle("");
    ratio__ -> GetYaxis()->SetTitle("Observed/Model");
    ratio__ -> SetMaximum(1.5);
    ratio__ -> SetMinimum(0.5);

    ////////////////
    
    TH1F * shaded__ = (TH1F*)den__->Clone("shaded__");
	shaded__->Reset();

    TH1F * shadedD__ = (TH1F*)den__->Clone("shadedD__");
	shadedD__->Reset();

    
    for(int k=0; k<shaded__->GetNbinsX()+1; ++k) 
    {

        if(ratio__->GetBinContent(k)!=0) 
        {

            ratio__->SetBinContent(k, ratio__->GetBinContent(k));

        }
        float z = 0;
        if(den__->GetBinContent(k)>0)
        {

            double err_ratio = sqrt(num__->GetBinContent(k))/den__->GetBinContent(k);
            ratio__->SetBinError(k, err_ratio);

        }


        double shaded_err = 0;
        if(den__->GetBinContent(k)!=0)
        {
            shaded_err = (den__->GetBinError(k))/(den__->GetBinContent(k));
            shaded__->SetBinContent(k,2*fabs(shaded_err));
            shadedD__->SetBinContent(k,-1*fabs(shaded_err) + 1.);

        }

    }
    
	ratio__ -> DrawCopy("PE");

	//shaded__->SetFillColor(kGray+2);
    shaded__->SetFillColorAlpha(0, 0.0);
    shadedD__->SetFillColorAlpha(0, 0.0);
	shaded__->SetLineColor(kGray+2);
	shadedD__->SetLineColor(kGray+2);
    shaded__->SetLineWidth(2);
	shadedD__->SetLineWidth(2);
	//gStyle->SetHatchesLineWidth(1);
	shaded__->SetFillStyle(3005);


    THStack * bandStack = new THStack();
    bandStack->Add(shadedD__);
    bandStack->Add(shaded__);
    bandStack->Draw("noclear histsames");
    
    //shadedD__->DrawCopy("hist sames");
	//shadedD__ -> DrawCopy("hist sames");
	
	ratio__ -> DrawCopy("PEsames");

        // /* draw the signals */
     //    TH1F * ZP800__ = (TH1F*) ZP800->Clone("ZP800__");
        

        // for(int k=0; k<ZP800__->GetNbinsX()+1; ++k) 
     //    {

     //    	if(den__->GetBinContent(k)>0)
     //    	{
     //    		double cont = ((ZP800__->GetBinContent(k)+den__->GetBinContent(k))/den__->GetBinContent(k)) - 1;
     //    		double errS = ((ZP800__->GetBinError(k)+den__->GetBinContent(k))/den__->GetBinContent(k)) - 1;

     //    		ZP800__->SetBinContent(k, cont);	
     //    		ZP800__->SetBinError(k, errS);	

     //    	}

     //    }	    

     //    ZP800__->SetMarkerStyle(22);
     //    ZP800__->SetMarkerColor(2);
        // ZP800__->SetLineColor(2);


     //    ZP800__->DrawCopy("PE sames");
     //    ZP800__->DrawCopy("hist sames");

	TLine line1;
	line1.SetLineColor(1);
	line1.SetLineWidth(1);
	//line1.DrawLine(ratio__->GetBinLowEdge(1),0,
	//ratio__->GetBinLowEdge(ratio__->GetNbinsX())+ratio__->GetBinWidth(ratio__->GetNbinsX()),0);

    C->Update();

    C->cd(2);
    
	TLegend leg2(.13,.75,.32,.85);
    leg2.SetTextSize(.04);
    leg2.SetTextFont(22);
    leg2.SetBorderSize(0);
    leg2.SetFillColor(0);	

    leg2.AddEntry(ratio__, "Obs./Bkg.","PE");
    leg2.AddEntry(shaded__, "Bkg. #pm 1 #sigma (syst.)","F");

    leg2.DrawClone();

	gPad->RedrawAxis("XY");
	gPad->SetGridx(1);
	gPad->SetGridy(1);

	gPad->Update();
////////////////


    TLatex L5;
    //L2.SetFont(42);
    L5.SetTextSize(0.06);
    L5.DrawLatexNDC(0.1,0.925,"#font[42]{#bf{Ratio}}");


    C->Update();

    /* delete pointers */
	delete data__;
    delete vv__;
    delete zvv__;
    delete ewk__;
	delete tt__;
	delete ztt__;
    delete zl__;
    delete zj__;
	delete w__;
    delete zhtautau__;
    delete glugluhtautau__;
    delete vbfhtautau__;
    delete smh__;
	delete qcd__;
    delete zp600__;
	delete zp1200__;
	delete zp1700__;
	delete zp2000__;
	delete zp2500__;
	delete ratio__;



}





// ////////////////
//     TH1F * shaded__ = (TH1F*)den__->Clone("shaded__");
// 	shaded__->Reset();

//     TH1F * shadedD__ = (TH1F*)den__->Clone("shadedD__");
// 	shadedD__->Reset();

//     for(int k=0; k<shaded__->GetNbinsX(); ++k) 
//     	{

//     		if(ratio__->GetBinContent(k)!=0) ratio__->SetBinContent(k, ratio__->GetBinContent(k)-1);
//     		float z = 0;
//     		if(den__->GetBinContent(k)>0)
//     		{
//     			shaded__->SetBinContent(k,   0+den__->GetBinError(k)/den__->GetBinContent(k)  );
//     			shadedD__->SetBinContent(k,   0-den__->GetBinError(k)/den__->GetBinContent(k)  );

//     		}
//     	}






// 	ratio__ -> DrawCopy("PE");


// 	shaded__->SetFillColor(3);
// 	shadedD__->SetFillColor(3);

// 	shaded__ -> DrawCopy("hist sames");
// 	shadedD__->DrawCopy("hist sames");
// 	ratio__ -> DrawCopy("PEsames");

// 	TLine line1;
// 	line1.SetLineColor(2);
// 	line1.SetLineWidth(2);
// 	line1.DrawLine(ratio__->GetBinLowEdge(1),0,
// 	ratio__->GetBinLowEdge(ratio__->GetNbinsX())+ratio__->GetBinWidth(ratio__->GetNbinsX()),0);


