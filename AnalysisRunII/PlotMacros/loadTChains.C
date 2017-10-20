#include "TChain.h"
#include <iostream>
#include "TString.h"
#include <stdio.h>      
#include <stdlib.h>     
#include <sstream>
#include <vector>
#include <algorithm>
#include <unistd.h>

bool useHTbinnedW = 1;

void fileLoad();

void setup_MC_files(std::string);

void setup_files_muTau();
void setup_files_eleTau();
void setup_files_eleMu();
void setup_files_TauTau();

void testComplete();

void reset_files();

TChain * DATA = new TChain("TauCheck");

TChain * DY = new TChain("TauCheck");
TChain * W = new TChain("TauCheck");
TChain * VV = new TChain("TauCheck");
TChain * ZVV = new TChain("TauCheck");
TChain * EWK = new TChain("TauCheck");

TChain * TT = new TChain("TauCheck");

TChain * MZP600_MA0300  = new TChain("TauCheck");
TChain * MZP600_MA0400  = new TChain("TauCheck");

TChain * MZP800_MA0300  = new TChain("TauCheck");
TChain * MZP800_MA0400  = new TChain("TauCheck");
TChain * MZP800_MA0500  = new TChain("TauCheck");
TChain * MZP800_MA0600  = new TChain("TauCheck");

TChain * MZP1000_MA0300  = new TChain("TauCheck");
TChain * MZP1000_MA0400  = new TChain("TauCheck");
TChain * MZP1000_MA0500  = new TChain("TauCheck");
TChain * MZP1000_MA0600  = new TChain("TauCheck");
TChain * MZP1000_MA0700  = new TChain("TauCheck");
TChain * MZP1000_MA0800  = new TChain("TauCheck");

TChain * MZP1200_MA0300  = new TChain("TauCheck");
TChain * MZP1200_MA0400  = new TChain("TauCheck");
TChain * MZP1200_MA0500  = new TChain("TauCheck");
TChain * MZP1200_MA0600  = new TChain("TauCheck");
TChain * MZP1200_MA0700  = new TChain("TauCheck");
TChain * MZP1200_MA0800  = new TChain("TauCheck");

TChain * MZP1400_MA0300  = new TChain("TauCheck");
TChain * MZP1400_MA0400  = new TChain("TauCheck");
TChain * MZP1400_MA0500  = new TChain("TauCheck");
TChain * MZP1400_MA0600  = new TChain("TauCheck");
TChain * MZP1400_MA0700  = new TChain("TauCheck");
TChain * MZP1400_MA0800  = new TChain("TauCheck");

TChain * MZP1700_MA0300  = new TChain("TauCheck");
TChain * MZP1700_MA0400  = new TChain("TauCheck");
TChain * MZP1700_MA0500  = new TChain("TauCheck");
TChain * MZP1700_MA0600  = new TChain("TauCheck");
TChain * MZP1700_MA0700  = new TChain("TauCheck");
TChain * MZP1700_MA0800  = new TChain("TauCheck");

TChain * MZP2000_MA0300  = new TChain("TauCheck");
TChain * MZP2000_MA0400  = new TChain("TauCheck");
TChain * MZP2000_MA0500  = new TChain("TauCheck");
TChain * MZP2000_MA0600  = new TChain("TauCheck");
TChain * MZP2000_MA0700  = new TChain("TauCheck");
TChain * MZP2000_MA0800  = new TChain("TauCheck");

TChain * MZP2500_MA0300  = new TChain("TauCheck");
TChain * MZP2500_MA0400  = new TChain("TauCheck");
TChain * MZP2500_MA0500  = new TChain("TauCheck");
TChain * MZP2500_MA0600  = new TChain("TauCheck");
TChain * MZP2500_MA0700  = new TChain("TauCheck");
TChain * MZP2500_MA0800  = new TChain("TauCheck");

TChain * MZP10_MX1 = new TChain("TauCheck");
TChain * MZP20_MX1 = new TChain("TauCheck");
TChain * MZP50_MX1 = new TChain("TauCheck");
TChain * MZP100_MX1 = new TChain("TauCheck");
TChain * MZP300_MX1 = new TChain("TauCheck");
TChain * MZP500_MX1 = new TChain("TauCheck");
TChain * MZP1000_MX1 = new TChain("TauCheck");
TChain * MZP2000_MX1 = new TChain("TauCheck");
TChain * MZP10000_MX1 = new TChain("TauCheck");

TChain * MZP10_MX50 = new TChain("TauCheck");
TChain * MZP50_MX50 = new TChain("TauCheck");
TChain * MZP95_MX50 = new TChain("TauCheck");
TChain * MZP200_MX50 = new TChain("TauCheck");
TChain * MZP300_MX50 = new TChain("TauCheck");
TChain * MZP10000_MX50 = new TChain("TauCheck");

TChain * MZP10_MX150 = new TChain("TauCheck");
TChain * MZP200_MX150 = new TChain("TauCheck");
TChain * MZP295_MX150 = new TChain("TauCheck");
TChain * MZP500_MX150 = new TChain("TauCheck");
TChain * MZP1000_MX150 = new TChain("TauCheck");
TChain * MZP10000_MX150 = new TChain("TauCheck");

TChain * MZP10_MX500 = new TChain("TauCheck");
TChain * MZP500_MX500 = new TChain("TauCheck");
TChain * MZP995_MX500 = new TChain("TauCheck");
TChain * MZP2000_MX500 = new TChain("TauCheck");
TChain * MZP10000_MX500 = new TChain("TauCheck");

TChain * MZP10_MX1000 = new TChain("TauCheck");
TChain * MZP1000_MX1000 = new TChain("TauCheck");
TChain * MZP1995_MX1000 = new TChain("TauCheck");
TChain * MZP10000_MX1000 = new TChain("TauCheck");

//BKG subgroups
//W
TChain * WJetsToLNu = new TChain("TauCheck");

TChain * W1JetsToLNu = new TChain("TauCheck");
TChain * W2JetsToLNu = new TChain("TauCheck");
TChain * W3JetsToLNu = new TChain("TauCheck");
TChain * W4JetsToLNu = new TChain("TauCheck");

TChain * WJetsToLNuHT100To200 = new TChain("TauCheck");
TChain * WJetsToLNuHT200To400 = new TChain("TauCheck");
TChain * WJetsToLNuHT400To600 = new TChain("TauCheck");
TChain * WJetsToLNuHT600To800 = new TChain("TauCheck");
TChain * WJetsToLNuHT800To1200 = new TChain("TauCheck");
TChain * WJetsToLNuHT1200To2500 = new TChain("TauCheck");
TChain * WJetsToLNuHT2500ToInf = new TChain("TauCheck");

//EWK
TChain * EWKWMinus2Jets = new TChain("TauCheck");
TChain * EWKWPlus2Jets = new TChain("TauCheck");
TChain * EWKZ2Jets_ZToLL = new TChain("TauCheck");
TChain * EWKZ2Jets_ZToNuNu = new TChain("TauCheck");

//DY
TChain * DYinc = new TChain("TauCheck");
TChain * DY1Jets = new TChain("TauCheck");
TChain * DY2Jets = new TChain("TauCheck");
TChain * DY3Jets = new TChain("TauCheck");
TChain * DY4Jets = new TChain("TauCheck");

//ZNuNu
TChain * ZJetsToNuNuHT100To200 = new TChain("TauCheck");
TChain * ZJetsToNuNuHT200To400 = new TChain("TauCheck");
TChain * ZJetsToNuNuHT400To600 = new TChain("TauCheck");
TChain * ZJetsToNuNuHT600To800 = new TChain("TauCheck");
TChain * ZJetsToNuNuHT800To1200 = new TChain("TauCheck");
TChain * ZJetsToNuNuHT1200To2500 = new TChain("TauCheck");
TChain * ZJetsToNuNuHT2500ToInf = new TChain("TauCheck");

//DY channels divided using bools
TChain * ZTT = new TChain("TauCheck");
TChain * ZJ = new TChain("TauCheck");
TChain * ZL = new TChain("TauCheck");

//TT subgroups
TChain * TTT = new TChain("TauCheck");
TChain * TTJ = new TChain("TauCheck");

TChain * TTZToLLNuNu = new TChain("TauCheck");
TChain * TTWJetsToLNu = new TChain("TauCheck");

//VV subgroups
TChain * VVT = new TChain("TauCheck");
TChain * VVJ = new TChain("TauCheck");

//SMH subgroups
//H
TChain * GluGluHTauTau = new TChain("TauCheck");
TChain * VBFHTauTau = new TChain("TauCheck");
TChain * ZHTauTau = new TChain("TauCheck");
TChain * WminusHTauTau = new TChain("TauCheck");
TChain * WplusHTauTau = new TChain("TauCheck");
TChain * VBFHToWW = new TChain("TauCheck");
TChain * GluGluHToWW = new TChain("TauCheck");
TChain * ttHJetToNonbb = new TChain("TauCheck");

//VV

TChain * ST_tW_antitop_5f_inclusiveDecays = new TChain("TauCheck");
TChain * ST_tW_top_5f_inclusiveDecays = new TChain("TauCheck");
TChain * ST_t_channel_antitop_4f_inclusiveDecays = new TChain("TauCheck");
TChain * ST_t_channel_top_4f_inclusiveDecays = new TChain("TauCheck");

TChain * WWTo2L2Nu = new TChain("TauCheck");
TChain * ZZTo2L2Nu = new TChain("TauCheck");
TChain * WWTo1L1Nu2Q = new TChain("TauCheck");
TChain * WZTo1L1Nu2Q = new TChain("TauCheck");
TChain * WZTo1L3Nu = new TChain("TauCheck");
TChain * WZTo2L2Q = new TChain("TauCheck");
TChain * ZZTo2L2Q = new TChain("TauCheck");
TChain * ZZTo2Q2Nu = new TChain("TauCheck");
TChain * ZZTo4L = new TChain("TauCheck");
TChain * WWW = new TChain("TauCheck");
TChain * WWZ = new TChain("TauCheck");
TChain * WZZ = new TChain("TauCheck");
TChain * ZZZ = new TChain("TauCheck");

TChain * WWTo2L2Nu_DoubleScattering = new TChain("TauCheck");
TChain * WWToLNuQQ = new TChain("TauCheck");
TChain * WWJJToLNuLNu_EWK = new TChain("TauCheck");
TChain * WpWpJJ_QCD = new TChain("TauCheck");
TChain * WpWpJJ_EWK = new TChain("TauCheck");
TChain * GluGluWWTo2L2Nu = new TChain("TauCheck");
TChain * WZTo3LNu = new TChain("TauCheck");

void fileLoad (TChain * CHAIN, TString PATH_TO_FILES, std::string FILE_NAME)
{

	FILE * CondorJobs;
	CondorJobs = popen("xrdfs root://cmseos.fnal.gov ls "+PATH_TO_FILES,"r");
	std::vector <std::string> condorVector;

	if ( CondorJobs != NULL )
    {
      char aline [ 10000 ]; 
      while ( fgets ( aline, sizeof aline, CondorJobs ) != NULL ) 
      {
         std::string astring(aline);
         astring.erase( std::remove( astring.begin(), astring.end(), '\n' ), astring.end() ) ;
         
         //Check for missing files
         
         std::string checkString = "/eos/uscms" + astring + "/" + FILE_NAME;
         //std::cout << checkString << std::endl;
         const char * checkChar = checkString.c_str();
         if (access(checkChar,F_OK) == -1) { std::cout << "NO FILE " << checkString << std::endl;}
         
         
         astring = "root://cmsxrootd.fnal.gov/" + astring + "/" + FILE_NAME;
         condorVector.push_back(astring);
      }
      fclose ( CondorJobs );
    }
    
    for(std::size_t i = 0; i<condorVector.size(); ++i)
	{
		CHAIN->Add(condorVector[i].c_str());
    }

}

void setup_MC_files(std::string FILE_NAME)
{

    fileLoad(GluGluHTauTau, "/store/user/gfunk/AugustProductionV1_8_0_26p1_V3/SIGNAL_GluGluHTauTau/GluGluHToTauTau_M125_13TeV_powheg_pythia8/", FILE_NAME);
    fileLoad(VBFHTauTau, "/store/user/gfunk/AugustProductionV1_8_0_26p1_V3/SIGNAL_VBFHTauTau/VBFHToTauTau_M125_13TeV_powheg_pythia8/", FILE_NAME);
    fileLoad(ZHTauTau, "/store/user/gfunk/AugustProductionV1_8_0_26p1_V3/SIGNAL_ZHTauTau/ZHToTauTau_M125_13TeV_powheg_pythia8/", FILE_NAME);
    fileLoad(WminusHTauTau, "/store/user/gfunk/AugustProductionV1_8_0_26p1_V3/SIGNAL_WHTauTau/WplusHToTauTau_M125_13TeV_powheg_pythia8/", FILE_NAME);
    fileLoad(WplusHTauTau, "/store/user/gfunk/AugustProductionV1_8_0_26p1_V3/SIGNAL_WHTauTau/WminusHToTauTau_M125_13TeV_powheg_pythia8/", FILE_NAME);
    fileLoad(VBFHToWW, "/store/user/gfunk/AugustProductionV1_8_0_26p1_V3/SIGNAL_HWW2L2Nu/VBFHToWWTo2L2Nu_M125_13TeV_powheg_JHUgen_pythia8/", FILE_NAME);
    fileLoad(GluGluHToWW, "/store/user/gfunk/AugustProductionV1_8_0_26p1_V3/SIGNAL_HWW2L2Nu/GluGluHToWWTo2L2Nu_M125_13TeV_powheg_JHUgen_pythia8/", FILE_NAME);
    fileLoad(ttHJetToNonbb, "/store/user/gfunk/AugustProductionV1_8_0_26p1_V3/SIGNAL_ttH/ttHJetToNonbb_M125_13TeV_amcatnloFXFX_madspin_pythia8_mWCutfix/", FILE_NAME);
    
    fileLoad(TT, "/store/user/gfunk/AugustProductionV1_8_0_26p1_V3/TT/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/", FILE_NAME);
    //fileLoad(TT, "/store/user/gfunk/AugustProductionV1_8_0_26p1_V3/VV/TTZToLLNuNu_M-10_TuneCUETP8M1_13TeV-amcatnlo-pythia8/", FILE_NAME);
    //fileLoad(TT, "/store/user/gfunk/AugustProductionV1_8_0_26p1_V3/VV/TTWJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8/", FILE_NAME);
    
    fileLoad(ST_t_channel_antitop_4f_inclusiveDecays, "/store/user/gfunk/AugustProductionV1_8_0_26p1_V3/VV/ST_t-channel_antitop_4f_inclusiveDecays_13TeV-powhegV2-madspin-pythia8_TuneCUETP8M1/", FILE_NAME);
    fileLoad(ST_t_channel_top_4f_inclusiveDecays, "/store/user/gfunk/AugustProductionV1_8_0_26p1_V3/VV/ST_t-channel_top_4f_inclusiveDecays_13TeV-powhegV2-madspin-pythia8_TuneCUETP8M1/", FILE_NAME);
    fileLoad(ST_tW_antitop_5f_inclusiveDecays, "/store/user/gfunk/AugustProductionV1_8_0_26p1_V3/VV/ST_tW_antitop_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1/", FILE_NAME);
    fileLoad(ST_tW_top_5f_inclusiveDecays, "/store/user/gfunk/AugustProductionV1_8_0_26p1_V3/VV/ST_tW_top_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1/", FILE_NAME);

    fileLoad(WWTo2L2Nu, "/store/user/gfunk/AugustProductionV1_8_0_26p1_V3/VV/WWTo2L2Nu_13TeV-powheg/", FILE_NAME);
    fileLoad(ZZTo2L2Nu, "/store/user/gfunk/AugustProductionV1_8_0_26p1_V3/VV/ZZTo2L2Nu_13TeV_powheg_pythia8/", FILE_NAME);
    
    fileLoad(WWTo1L1Nu2Q, "/store/user/gfunk/AugustProductionV1_8_0_26p1_V3/VV/WWTo1L1Nu2Q_13TeV_amcatnloFXFX_madspin_pythia8/", FILE_NAME);
    fileLoad(WZTo1L1Nu2Q, "/store/user/gfunk/AugustProductionV1_8_0_26p1_V3/VV/WZTo1L1Nu2Q_13TeV_amcatnloFXFX_madspin_pythia8/", FILE_NAME);
    fileLoad(WZTo1L3Nu, "/store/user/gfunk/AugustProductionV1_8_0_26p1_V3/VV/WZTo1L3Nu_13TeV_amcatnloFXFX_madspin_pythia8/", FILE_NAME);
    fileLoad(WZTo2L2Q, "/store/user/gfunk/AugustProductionV1_8_0_26p1_V3/VV/WZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8/", FILE_NAME);
    fileLoad(ZZTo2L2Q, "/store/user/gfunk/AugustProductionV1_8_0_26p1_V3/VV/ZZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8/", FILE_NAME);
    fileLoad(ZZTo2Q2Nu, "/store/user/gfunk/AugustProductionV1_8_0_26p1_V3/VV/ZZTo2Q2Nu_13TeV_amcatnloFXFX_madspin_pythia8/", FILE_NAME);
    fileLoad(ZZTo4L, "/store/user/gfunk/AugustProductionV1_8_0_26p1_V3/VV/ZZTo4L_13TeV-amcatnloFXFX-pythia8/",FILE_NAME);
    fileLoad(ZZZ, "/store/user/gfunk/AugustProductionV1_8_0_26p1_V3/VV/ZZZ_TuneCUETP8M1_13TeV-amcatnlo-pythia8/", FILE_NAME);
    fileLoad(WZZ, "/store/user/gfunk/AugustProductionV1_8_0_26p1_V3/VV/WZZ_TuneCUETP8M1_13TeV-amcatnlo-pythia8/", FILE_NAME);
    fileLoad(WWZ, "/store/user/gfunk/AugustProductionV1_8_0_26p1_V3/VV/WWZ_TuneCUETP8M1_13TeV-amcatnlo-pythia8/", FILE_NAME);
    fileLoad(WWW, "/store/user/gfunk/AugustProductionV1_8_0_26p1_V3/VV/WWW_4F_TuneCUETP8M1_13TeV-amcatnlo-pythia8/", FILE_NAME);
    
    fileLoad(WWTo2L2Nu_DoubleScattering, "/store/user/gfunk/AugustProductionV1_8_0_26p1_V3/VV/WWTo2L2Nu_DoubleScattering_13TeV-pythia8/", FILE_NAME);
    fileLoad(WWToLNuQQ, "/store/user/gfunk/AugustProductionV1_8_0_26p1_V3/VV/WWToLNuQQ_13TeV-powheg/", FILE_NAME);
    fileLoad(WWJJToLNuLNu_EWK, "/store/user/gfunk/AugustProductionV1_8_0_26p1_V3/VV/WWJJToLNuLNu_EWK_noTop_13TeV-madgraph-pythia8/", FILE_NAME);
    fileLoad(WpWpJJ_QCD, "/store/user/gfunk/AugustProductionV1_8_0_26p1_V3/VV/WpWpJJ_QCD_TuneCUETP8M1_13TeV-madgraph-pythia8/", FILE_NAME);
    fileLoad(WpWpJJ_EWK, "/store/user/gfunk/AugustProductionV1_8_0_26p1_V3/VV/WpWpJJ_EWK_TuneCUETP8M1_13TeV-madgraph-pythia8/", FILE_NAME);
    fileLoad(GluGluWWTo2L2Nu, "/store/user/gfunk/AugustProductionV1_8_0_26p1_V3/VV/GluGluWWTo2L2Nu_MCFM_13TeV/", FILE_NAME);
    fileLoad(WZTo3LNu, "/store/user/gfunk/AugustProductionV1_8_0_26p1_V3/VV/WZTo3LNu_TuneCUETP8M1_13TeV-powheg-pythia8/", FILE_NAME);
    
    if(useHTbinnedW==1)
    {
        fileLoad(WJetsToLNu, "/store/user/gfunk/AugustProductionV1_8_0_26p1_V3/W/WJetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/", FILE_NAME);
        fileLoad(WJetsToLNuHT100To200, "/store/user/gfunk/AugustProductionV1_8_0_26p1_V3/W/WJetsToLNu_HT-100To200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/", FILE_NAME);
        fileLoad(WJetsToLNuHT200To400, "/store/user/gfunk/AugustProductionV1_8_0_26p1_V3/W/WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/", FILE_NAME);
        fileLoad(WJetsToLNuHT400To600, "/store/user/gfunk/AugustProductionV1_8_0_26p1_V3/W/WJetsToLNu_HT-400To600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/", FILE_NAME);
        fileLoad(WJetsToLNuHT600To800, "/store/user/gfunk/AugustProductionV1_8_0_26p1_V3/W/WJetsToLNu_HT-600To800_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/", FILE_NAME);
        fileLoad(WJetsToLNuHT800To1200, "/store/user/gfunk/AugustProductionV1_8_0_26p1_V3/W/WJetsToLNu_HT-800To1200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/", FILE_NAME);
        fileLoad(WJetsToLNuHT2500ToInf, "/store/user/gfunk/AugustProductionV1_8_0_26p1_V3/W/WJetsToLNu_HT-1200To2500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/", FILE_NAME);
        fileLoad(WJetsToLNuHT1200To2500, "/store/user/gfunk/AugustProductionV1_8_0_26p1_V3/W/WJetsToLNu_HT-2500ToInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/", FILE_NAME);
    }
    else
    {
        fileLoad(W1JetsToLNu, "/store/user/gfunk/AugustProductionV1_8_0_26p1_V3/W/W1JetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/", FILE_NAME);
        fileLoad(W2JetsToLNu, "/store/user/gfunk/AugustProductionV1_8_0_26p1_V3/W/W2JetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/", FILE_NAME);
        fileLoad(W3JetsToLNu, "/store/user/gfunk/AugustProductionV1_8_0_26p1_V3/W/W3JetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/", FILE_NAME);
        fileLoad(W4JetsToLNu, "/store/user/gfunk/AugustProductionV1_8_0_26p1_V3/W/W4JetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/", FILE_NAME);
    }
    
    fileLoad(DY1Jets, "/store/user/gfunk/AugustProductionV1_8_0_26p1_V3/DY/DY1JetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/", FILE_NAME);
    fileLoad(DY2Jets, "/store/user/gfunk/AugustProductionV1_8_0_26p1_V3/DY/DY2JetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/", FILE_NAME);
    fileLoad(DY3Jets, "/store/user/gfunk/AugustProductionV1_8_0_26p1_V3/DY/DY3JetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/", FILE_NAME);
    fileLoad(DY4Jets, "/store/user/gfunk/AugustProductionV1_8_0_26p1_V3/DY/DY4JetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/", FILE_NAME);
    fileLoad(DYinc, "/store/user/gfunk/AugustProductionV1_8_0_26p1_V3/DY/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/", FILE_NAME);
    
    fileLoad(ZJetsToNuNuHT100To200, "/store/user/gfunk/AugustProductionV1_8_0_26p1_V3/DYinv/ZJetsToNuNu_HT-100To200_13TeV-madgraph/", FILE_NAME);
    fileLoad(ZJetsToNuNuHT200To400, "/store/user/gfunk/AugustProductionV1_8_0_26p1_V3/DYinv/ZJetsToNuNu_HT-200To400_13TeV-madgraph/", FILE_NAME);
    fileLoad(ZJetsToNuNuHT400To600, "/store/user/gfunk/AugustProductionV1_8_0_26p1_V3/DYinv/ZJetsToNuNu_HT-400To600_13TeV-madgraph/", FILE_NAME);
    fileLoad(ZJetsToNuNuHT600To800, "/store/user/gfunk/AugustProductionV1_8_0_26p1_V3/DYinv/ZJetsToNuNu_HT-600To800_13TeV-madgraph/", FILE_NAME);
    fileLoad(ZJetsToNuNuHT800To1200, "/store/user/gfunk/AugustProductionV1_8_0_26p1_V3/DYinv/ZJetsToNuNu_HT-800To1200_13TeV-madgraph/", FILE_NAME);
    fileLoad(ZJetsToNuNuHT1200To2500, "/store/user/gfunk/AugustProductionV1_8_0_26p1_V3/DYinv/ZJetsToNuNu_HT-1200To2500_13TeV-madgraph/", FILE_NAME);
    fileLoad(ZJetsToNuNuHT2500ToInf, "/store/user/gfunk/AugustProductionV1_8_0_26p1_V3/DYinv/ZJetsToNuNu_HT-2500ToInf_13TeV-madgraph/", FILE_NAME);
    
    fileLoad(EWKWMinus2Jets, "/store/user/gfunk/AugustProductionV1_8_0_26p1_V3/EWK/EWKWMinus2Jets_WToLNu_M-50_13TeV-madgraph-pythia8/", FILE_NAME);
    fileLoad(EWKWPlus2Jets, "/store/user/gfunk/AugustProductionV1_8_0_26p1_V3/EWK/EWKWPlus2Jets_WToLNu_M-50_13TeV-madgraph-pythia8/", FILE_NAME);
    fileLoad(EWKZ2Jets_ZToLL, "/store/user/gfunk/AugustProductionV1_8_0_26p1_V3/EWK/EWKZ2Jets_ZToLL_M-50_13TeV-madgraph-pythia8/", FILE_NAME);
    fileLoad(EWKZ2Jets_ZToNuNu, "/store/user/gfunk/AugustProductionV1_8_0_26p1_V3/EWK/EWKZ2Jets_ZToNuNu_13TeV-madgraph-pythia8/", FILE_NAME);
    
    fileLoad(MZP600_MA0300, "/store/user/gfunk/AugustProductionV1_8_0_26p1_V3/SIGNAL_MONO_HIGGS/ZprimeToA0hToA0chichihtautau_2HDM_MZp-600_MA0-300_13TeV-madgraph-pythia8/", FILE_NAME);
    fileLoad(MZP600_MA0400, "/store/user/gfunk/AugustProductionV1_8_0_26p1_V3/SIGNAL_MONO_HIGGS/ZprimeToA0hToA0chichihtautau_2HDM_MZp-600_MA0-400_13TeV-madgraph/", FILE_NAME);
    
    fileLoad(MZP800_MA0300, "/store/user/gfunk/AugustProductionV1_8_0_26p1_V3/SIGNAL_MONO_HIGGS/ZprimeToA0hToA0chichihtautau_2HDM_MZp-800_MA0-300_13TeV-madgraph-pythia8/", FILE_NAME);
    fileLoad(MZP800_MA0400, "/store/user/gfunk/AugustProductionV1_8_0_26p1_V3/SIGNAL_MONO_HIGGS/ZprimeToA0hToA0chichihtautau_2HDM_MZp-800_MA0-400_13TeV-madgraph/", FILE_NAME);
    fileLoad(MZP800_MA0500, "/store/user/gfunk/AugustProductionV1_8_0_26p1_V3/SIGNAL_MONO_HIGGS/ZprimeToA0hToA0chichihtautau_2HDM_MZp-800_MA0-500_13TeV-madgraph/", FILE_NAME);
    fileLoad(MZP800_MA0600, "/store/user/gfunk/AugustProductionV1_8_0_26p1_V3/SIGNAL_MONO_HIGGS/ZprimeToA0hToA0chichihtautau_2HDM_MZp-800_MA0-600_13TeV-madgraph/", FILE_NAME);
    
    fileLoad(MZP1000_MA0300, "/store/user/gfunk/AugustProductionV1_8_0_26p1_V3/SIGNAL_MONO_HIGGS/ZprimeToA0hToA0chichihtautau_2HDM_MZp-1000_MA0-300_13TeV-madgraph-pythia8/", FILE_NAME);
    fileLoad(MZP1000_MA0400, "/store/user/gfunk/AugustProductionV1_8_0_26p1_V3/SIGNAL_MONO_HIGGS/ZprimeToA0hToA0chichihtautau_2HDM_MZp-1000_MA0-400_13TeV-madgraph/", FILE_NAME);
    fileLoad(MZP1000_MA0500, "/store/user/gfunk/AugustProductionV1_8_0_26p1_V3/SIGNAL_MONO_HIGGS/ZprimeToA0hToA0chichihtautau_2HDM_MZp-1000_MA0-500_13TeV-madgraph/", FILE_NAME);
    fileLoad(MZP1000_MA0600, "/store/user/gfunk/AugustProductionV1_8_0_26p1_V3/SIGNAL_MONO_HIGGS/ZprimeToA0hToA0chichihtautau_2HDM_MZp-1000_MA0-600_13TeV-madgraph/", FILE_NAME);
    fileLoad(MZP1000_MA0700, "/store/user/gfunk/AugustProductionV1_8_0_26p1_V3/SIGNAL_MONO_HIGGS/ZprimeToA0hToA0chichihtautau_2HDM_MZp-1000_MA0-700_13TeV-madgraph/", FILE_NAME);
    fileLoad(MZP1000_MA0800, "/store/user/gfunk/AugustProductionV1_8_0_26p1_V3/SIGNAL_MONO_HIGGS/ZprimeToA0hToA0chichihtautau_2HDM_MZp-1000_MA0-800_13TeV-madgraph/", FILE_NAME);
    
    fileLoad(MZP1200_MA0300, "/store/user/gfunk/AugustProductionV1_8_0_26p1_V3/SIGNAL_MONO_HIGGS/ZprimeToA0hToA0chichihtautau_2HDM_MZp-1200_MA0-300_13TeV-madgraph-pythia8/", FILE_NAME);
    fileLoad(MZP1200_MA0400, "/store/user/gfunk/AugustProductionV1_8_0_26p1_V3/SIGNAL_MONO_HIGGS/ZprimeToA0hToA0chichihtautau_2HDM_MZp-1200_MA0-400_13TeV-madgraph/", FILE_NAME);
    fileLoad(MZP1200_MA0500, "/store/user/gfunk/AugustProductionV1_8_0_26p1_V3/SIGNAL_MONO_HIGGS/ZprimeToA0hToA0chichihtautau_2HDM_MZp-1200_MA0-500_13TeV-madgraph/", FILE_NAME);
    fileLoad(MZP1200_MA0600, "/store/user/gfunk/AugustProductionV1_8_0_26p1_V3/SIGNAL_MONO_HIGGS/ZprimeToA0hToA0chichihtautau_2HDM_MZp-1200_MA0-600_13TeV-madgraph/", FILE_NAME);
    fileLoad(MZP1200_MA0700, "/store/user/gfunk/AugustProductionV1_8_0_26p1_V3/SIGNAL_MONO_HIGGS/ZprimeToA0hToA0chichihtautau_2HDM_MZp-1200_MA0-700_13TeV-madgraph/", FILE_NAME);
    fileLoad(MZP1200_MA0800, "/store/user/gfunk/AugustProductionV1_8_0_26p1_V3/SIGNAL_MONO_HIGGS/ZprimeToA0hToA0chichihtautau_2HDM_MZp-1200_MA0-800_13TeV-madgraph/", FILE_NAME);
    
    fileLoad(MZP1400_MA0300, "/store/user/gfunk/AugustProductionV1_8_0_26p1_V3/SIGNAL_MONO_HIGGS/ZprimeToA0hToA0chichihtautau_2HDM_MZp-1400_MA0-300_13TeV-madgraph-pythia8/", FILE_NAME);
    fileLoad(MZP1400_MA0400, "/store/user/gfunk/AugustProductionV1_8_0_26p1_V3/SIGNAL_MONO_HIGGS/ZprimeToA0hToA0chichihtautau_2HDM_MZp-1400_MA0-400_13TeV-madgraph/", FILE_NAME);
    fileLoad(MZP1400_MA0500, "/store/user/gfunk/AugustProductionV1_8_0_26p1_V3/SIGNAL_MONO_HIGGS/ZprimeToA0hToA0chichihtautau_2HDM_MZp-1400_MA0-500_13TeV-madgraph/", FILE_NAME);
    fileLoad(MZP1400_MA0600, "/store/user/gfunk/AugustProductionV1_8_0_26p1_V3/SIGNAL_MONO_HIGGS/ZprimeToA0hToA0chichihtautau_2HDM_MZp-1400_MA0-600_13TeV-madgraph/", FILE_NAME);
    fileLoad(MZP1400_MA0700, "/store/user/gfunk/AugustProductionV1_8_0_26p1_V3/SIGNAL_MONO_HIGGS/ZprimeToA0hToA0chichihtautau_2HDM_MZp-1400_MA0-700_13TeV-madgraph/", FILE_NAME);
    fileLoad(MZP1400_MA0800, "/store/user/gfunk/AugustProductionV1_8_0_26p1_V3/SIGNAL_MONO_HIGGS/ZprimeToA0hToA0chichihtautau_2HDM_MZp-1400_MA0-800_13TeV-madgraph/", FILE_NAME);
    
    fileLoad(MZP1700_MA0300, "/store/user/gfunk/AugustProductionV1_8_0_26p1_V3/SIGNAL_MONO_HIGGS/ZprimeToA0hToA0chichihtautau_2HDM_MZp-1700_MA0-300_13TeV-madgraph-pythia8/", FILE_NAME);
    fileLoad(MZP1700_MA0400, "/store/user/gfunk/AugustProductionV1_8_0_26p1_V3/SIGNAL_MONO_HIGGS/ZprimeToA0hToA0chichihtautau_2HDM_MZp-1700_MA0-400_13TeV-madgraph/", FILE_NAME);
    fileLoad(MZP1700_MA0500, "/store/user/gfunk/AugustProductionV1_8_0_26p1_V3/SIGNAL_MONO_HIGGS/ZprimeToA0hToA0chichihtautau_2HDM_MZp-1700_MA0-500_13TeV-madgraph/", FILE_NAME);
    fileLoad(MZP1700_MA0600, "/store/user/gfunk/AugustProductionV1_8_0_26p1_V3/SIGNAL_MONO_HIGGS/ZprimeToA0hToA0chichihtautau_2HDM_MZp-1700_MA0-600_13TeV-madgraph/", FILE_NAME);
    fileLoad(MZP1700_MA0700, "/store/user/gfunk/AugustProductionV1_8_0_26p1_V3/SIGNAL_MONO_HIGGS/ZprimeToA0hToA0chichihtautau_2HDM_MZp-1700_MA0-700_13TeV-madgraph/", FILE_NAME);
    fileLoad(MZP1700_MA0800, "/store/user/gfunk/AugustProductionV1_8_0_26p1_V3/SIGNAL_MONO_HIGGS/ZprimeToA0hToA0chichihtautau_2HDM_MZp-1700_MA0-800_13TeV-madgraph/", FILE_NAME);
    
    fileLoad(MZP2000_MA0300, "/store/user/gfunk/AugustProductionV1_8_0_26p1_V3/SIGNAL_MONO_HIGGS/ZprimeToA0hToA0chichihtautau_2HDM_MZp-2000_MA0-300_13TeV-madgraph-pythia8/", FILE_NAME);
    fileLoad(MZP2000_MA0400, "/store/user/gfunk/AugustProductionV1_8_0_26p1_V3/SIGNAL_MONO_HIGGS/ZprimeToA0hToA0chichihtautau_2HDM_MZp-2000_MA0-400_13TeV-madgraph/", FILE_NAME);
    fileLoad(MZP2000_MA0500, "/store/user/gfunk/AugustProductionV1_8_0_26p1_V3/SIGNAL_MONO_HIGGS/ZprimeToA0hToA0chichihtautau_2HDM_MZp-2000_MA0-500_13TeV-madgraph/", FILE_NAME);
    fileLoad(MZP2000_MA0600, "/store/user/gfunk/AugustProductionV1_8_0_26p1_V3/SIGNAL_MONO_HIGGS/ZprimeToA0hToA0chichihtautau_2HDM_MZp-2000_MA0-600_13TeV-madgraph/", FILE_NAME);
    fileLoad(MZP2000_MA0700, "/store/user/gfunk/AugustProductionV1_8_0_26p1_V3/SIGNAL_MONO_HIGGS/ZprimeToA0hToA0chichihtautau_2HDM_MZp-2000_MA0-700_13TeV-madgraph/", FILE_NAME);
    fileLoad(MZP2000_MA0800, "/store/user/gfunk/AugustProductionV1_8_0_26p1_V3/SIGNAL_MONO_HIGGS/ZprimeToA0hToA0chichihtautau_2HDM_MZp-2000_MA0-800_13TeV-madgraph/", FILE_NAME);
    
    fileLoad(MZP2500_MA0300, "/store/user/gfunk/AugustProductionV1_8_0_26p1_V3/SIGNAL_MONO_HIGGS/ZprimeToA0hToA0chichihtautau_2HDM_MZp-2500_MA0-300_13TeV-madgraph-pythia8/", FILE_NAME);
    fileLoad(MZP2500_MA0400, "/store/user/gfunk/AugustProductionV1_8_0_26p1_V3/SIGNAL_MONO_HIGGS/ZprimeToA0hToA0chichihtautau_2HDM_MZp-2500_MA0-400_13TeV-madgraph/", FILE_NAME);
    fileLoad(MZP2500_MA0500, "/store/user/gfunk/AugustProductionV1_8_0_26p1_V3/SIGNAL_MONO_HIGGS/ZprimeToA0hToA0chichihtautau_2HDM_MZp-2500_MA0-500_13TeV-madgraph/", FILE_NAME);
    fileLoad(MZP2500_MA0600, "/store/user/gfunk/AugustProductionV1_8_0_26p1_V3/SIGNAL_MONO_HIGGS/ZprimeToA0hToA0chichihtautau_2HDM_MZp-2500_MA0-600_13TeV-madgraph/", FILE_NAME);
    fileLoad(MZP2500_MA0700, "/store/user/gfunk/AugustProductionV1_8_0_26p1_V3/SIGNAL_MONO_HIGGS/ZprimeToA0hToA0chichihtautau_2HDM_MZp-2500_MA0-700_13TeV-madgraph/", FILE_NAME);
    fileLoad(MZP2500_MA0800, "/store/user/gfunk/AugustProductionV1_8_0_26p1_V3/SIGNAL_MONO_HIGGS/ZprimeToA0hToA0chichihtautau_2HDM_MZp-2500_MA0-800_13TeV-madgraph/", FILE_NAME);
    
    fileLoad(MZP10_MX1, "/store/user/gfunk/AugustProductionV1_8_0_26p1_V3/SIGNAL_MONO_HIGGS/MonoHtautau_ZpBaryonic_MZp-10_MChi-1_13TeV-madgraph/", FILE_NAME);
    fileLoad(MZP20_MX1, "/store/user/gfunk/AugustProductionV1_8_0_26p1_V3/SIGNAL_MONO_HIGGS/MonoHtautau_ZpBaryonic_MZp-20_MChi-1_13TeV-madgraph/", FILE_NAME);
    fileLoad(MZP50_MX1, "/store/user/gfunk/AugustProductionV1_8_0_26p1_V3/SIGNAL_MONO_HIGGS/MonoHtautau_ZpBaryonic_MZp-50_MChi-1_13TeV-madgraph/", FILE_NAME);
    fileLoad(MZP100_MX1, "/store/user/gfunk/AugustProductionV1_8_0_26p1_V3/SIGNAL_MONO_HIGGS/MonoHtautau_ZpBaryonic_MZp-100_MChi-1_13TeV-madgraph/", FILE_NAME);
    fileLoad(MZP300_MX1, "/store/user/gfunk/AugustProductionV1_8_0_26p1_V3/SIGNAL_MONO_HIGGS/MonoHtautau_ZpBaryonic_MZp-300_MChi-1_13TeV-madgraph/", FILE_NAME);
    fileLoad(MZP500_MX1, "/store/user/gfunk/AugustProductionV1_8_0_26p1_V3/SIGNAL_MONO_HIGGS/MonoHtautau_ZpBaryonic_MZp-500_MChi-1_13TeV-madgraph/", FILE_NAME);
    fileLoad(MZP1000_MX1, "/store/user/gfunk/AugustProductionV1_8_0_26p1_V3/SIGNAL_MONO_HIGGS/MonoHtautau_ZpBaryonic_MZp-1000_MChi-1_13TeV-madgraph/", FILE_NAME);
    fileLoad(MZP2000_MX1, "/store/user/gfunk/AugustProductionV1_8_0_26p1_V3/SIGNAL_MONO_HIGGS/MonoHtautau_ZpBaryonic_MZp-2000_MChi-1_13TeV-madgraph/", FILE_NAME);
    fileLoad(MZP10000_MX1, "/store/user/gfunk/AugustProductionV1_8_0_26p1_V3/SIGNAL_MONO_HIGGS/MonoHtautau_ZpBaryonic_MZp-10000_MChi-1_13TeV-madgraph/", FILE_NAME);
    
    fileLoad(MZP10_MX50, "/store/user/gfunk/AugustProductionV1_8_0_26p1_V3/SIGNAL_MONO_HIGGS/MonoHtautau_ZpBaryonic_MZp-10_MChi-50_13TeV-madgraph/", FILE_NAME);
    fileLoad(MZP50_MX50, "/store/user/gfunk/AugustProductionV1_8_0_26p1_V3/SIGNAL_MONO_HIGGS/MonoHtautau_ZpBaryonic_MZp-50_MChi-50_13TeV-madgraph/", FILE_NAME);
    fileLoad(MZP95_MX50, "/store/user/gfunk/AugustProductionV1_8_0_26p1_V3/SIGNAL_MONO_HIGGS/MonoHtautau_ZpBaryonic_MZp-95_MChi-50_13TeV-madgraph/", FILE_NAME);
    fileLoad(MZP200_MX50, "/store/user/gfunk/AugustProductionV1_8_0_26p1_V3/SIGNAL_MONO_HIGGS/MonoHtautau_ZpBaryonic_MZp-200_MChi-50_13TeV-madgraph/", FILE_NAME);
    fileLoad(MZP300_MX50, "/store/user/gfunk/AugustProductionV1_8_0_26p1_V3/SIGNAL_MONO_HIGGS/MonoHtautau_ZpBaryonic_MZp-300_MChi-50_13TeV-madgraph/", FILE_NAME);
    fileLoad(MZP10000_MX50, "/store/user/gfunk/AugustProductionV1_8_0_26p1_V3/SIGNAL_MONO_HIGGS/MonoHtautau_ZpBaryonic_MZp-10000_MChi-50_13TeV-madgraph/", FILE_NAME);

    fileLoad(MZP10_MX150, "/store/user/gfunk/AugustProductionV1_8_0_26p1_V3/SIGNAL_MONO_HIGGS/MonoHtautau_ZpBaryonic_MZp-10_MChi-150_13TeV-madgraph/", FILE_NAME);
    fileLoad(MZP200_MX150, "/store/user/gfunk/AugustProductionV1_8_0_26p1_V3/SIGNAL_MONO_HIGGS/MonoHtautau_ZpBaryonic_MZp-200_MChi-150_13TeV-madgraph/", FILE_NAME);
    fileLoad(MZP295_MX150, "/store/user/gfunk/AugustProductionV1_8_0_26p1_V3/SIGNAL_MONO_HIGGS/MonoHtautau_ZpBaryonic_MZp-295_MChi-150_13TeV-madgraph/", FILE_NAME);
    fileLoad(MZP500_MX150, "/store/user/gfunk/AugustProductionV1_8_0_26p1_V3/SIGNAL_MONO_HIGGS/MonoHtautau_ZpBaryonic_MZp-500_MChi-150_13TeV-madgraph/", FILE_NAME);
    fileLoad(MZP1000_MX150, "/store/user/gfunk/AugustProductionV1_8_0_26p1_V3/SIGNAL_MONO_HIGGS/MonoHtautau_ZpBaryonic_MZp-1000_MChi-150_13TeV-madgraph/", FILE_NAME);
    fileLoad(MZP10000_MX150, "/store/user/gfunk/AugustProductionV1_8_0_26p1_V3/SIGNAL_MONO_HIGGS/MonoHtautau_ZpBaryonic_MZp-10000_MChi-150_13TeV-madgraph/", FILE_NAME);
    
    fileLoad(MZP10_MX500, "/store/user/gfunk/AugustProductionV1_8_0_26p1_V3/SIGNAL_MONO_HIGGS/MonoHtautau_ZpBaryonic_MZp-10_MChi-500_13TeV-madgraph/", FILE_NAME);
    fileLoad(MZP500_MX500, "/store/user/gfunk/AugustProductionV1_8_0_26p1_V3/SIGNAL_MONO_HIGGS/MonoHtautau_ZpBaryonic_MZp-500_MChi-500_13TeV-madgraph/", FILE_NAME);
    fileLoad(MZP995_MX500, "/store/user/gfunk/AugustProductionV1_8_0_26p1_V3/SIGNAL_MONO_HIGGS/MonoHtautau_ZpBaryonic_MZp-995_MChi-500_13TeV-madgraph/", FILE_NAME);
    fileLoad(MZP2000_MX500, "/store/user/gfunk/AugustProductionV1_8_0_26p1_V3/SIGNAL_MONO_HIGGS/MonoHtautau_ZpBaryonic_MZp-2000_MChi-500_13TeV-madgraph/", FILE_NAME);
    fileLoad(MZP10000_MX500, "/store/user/gfunk/AugustProductionV1_8_0_26p1_V3/SIGNAL_MONO_HIGGS/MonoHtautau_ZpBaryonic_MZp-10000_MChi-500_13TeV-madgraph/", FILE_NAME);
    
    fileLoad(MZP10_MX1000, "/store/user/gfunk/AugustProductionV1_8_0_26p1_V3/SIGNAL_MONO_HIGGS/MonoHtautau_ZpBaryonic_MZp-10_MChi-1000_13TeV-madgraph/", FILE_NAME);
    fileLoad(MZP1000_MX1000, "/store/user/gfunk/AugustProductionV1_8_0_26p1_V3/SIGNAL_MONO_HIGGS/MonoHtautau_ZpBaryonic_MZp-1000_MChi-1000_13TeV-madgraph/", FILE_NAME);
    fileLoad(MZP1995_MX1000, "/store/user/gfunk/AugustProductionV1_8_0_26p1_V3/SIGNAL_MONO_HIGGS/MonoHtautau_ZpBaryonic_MZp-1995_MChi-1000_13TeV-madgraph/", FILE_NAME);
    fileLoad(MZP10000_MX1000, "/store/user/gfunk/AugustProductionV1_8_0_26p1_V3/SIGNAL_MONO_HIGGS/MonoHtautau_ZpBaryonic_MZp-10000_MChi-1000_13TeV-madgraph/", FILE_NAME);
    
    //Consolidate TChains to main BKG categories
    
    DY->Add(DYinc);
    DY->Add(DY1Jets);
    DY->Add(DY2Jets);
    DY->Add(DY3Jets);
    DY->Add(DY4Jets);
    DY->Add(EWKZ2Jets_ZToLL);
    DY->Add(EWKZ2Jets_ZToNuNu);

    ZTT = (TChain*) DY->CopyTree("IsZTT==1");
    ZJ = (TChain*) DY->CopyTree("IsZJ==1");
    ZL = (TChain*) DY->CopyTree("IsZL==1");
    
    TTT = (TChain*) TT->CopyTree("IsTTT==1");
    TTJ = (TChain*) TT->CopyTree("IsTTT==0");
    
    ZVV->Add(ZJetsToNuNuHT100To200);
    ZVV->Add(ZJetsToNuNuHT200To400);
    ZVV->Add(ZJetsToNuNuHT400To600);
    ZVV->Add(ZJetsToNuNuHT600To800);
    ZVV->Add(ZJetsToNuNuHT800To1200);
    ZVV->Add(ZJetsToNuNuHT1200To2500);
    ZVV->Add(ZJetsToNuNuHT2500ToInf);
    
    //Main W samples, choose 1 binning type
    if(useHTbinnedW==1)
    {
        W->Add(WJetsToLNu);
        W->Add(WJetsToLNuHT100To200);
        W->Add(WJetsToLNuHT200To400);
        W->Add(WJetsToLNuHT400To600);
        W->Add(WJetsToLNuHT600To800);
        W->Add(WJetsToLNuHT800To1200);
        W->Add(WJetsToLNuHT1200To2500);
        W->Add(WJetsToLNuHT2500ToInf);
    }
    else
    {
        W->Add(WJetsToLNu);
        W->Add(W1JetsToLNu);
        W->Add(W2JetsToLNu);
        W->Add(W3JetsToLNu);
        W->Add(W4JetsToLNu);
    }
    
    W->Add(EWKWMinus2Jets);
    W->Add(EWKWPlus2Jets);
    
    VV->Add(ST_t_channel_antitop_4f_inclusiveDecays);
    VV->Add(ST_t_channel_top_4f_inclusiveDecays);
    VV->Add(ST_tW_antitop_5f_inclusiveDecays);
    VV->Add(ST_tW_top_5f_inclusiveDecays);
    
    VV->Add(WWTo2L2Nu);
    VV->Add(ZZTo2L2Nu);
    
    VV->Add(WWTo1L1Nu2Q);
    VV->Add(WZTo1L1Nu2Q);
    VV->Add(WZTo1L3Nu);
    VV->Add(WZTo2L2Q);
    VV->Add(ZZTo2L2Q);
    VV->Add(ZZTo2Q2Nu);
    VV->Add(ZZTo4L);
    VV->Add(WWW);
    VV->Add(WWZ);
    VV->Add(WZZ);
    VV->Add(ZZZ);
    
    VV->Add(WWTo2L2Nu_DoubleScattering);
    VV->Add(WWToLNuQQ);
    VV->Add(WWJJToLNuLNu_EWK);
    VV->Add(WpWpJJ_QCD);
    VV->Add(WpWpJJ_EWK);
    VV->Add(GluGluWWTo2L2Nu);
    VV->Add(WZTo3LNu);
    
    
    GluGluHTauTau->Add(GluGluHToWW);
    GluGluHTauTau->Add(VBFHToWW);
    GluGluHTauTau->Add(WminusHTauTau);
    GluGluHTauTau->Add(WplusHTauTau);
    GluGluHTauTau->Add(ttHJetToNonbb);
    
    VVT = (TChain*) VV->CopyTree("IsZTT==1");
    VVJ = (TChain*) VV->CopyTree("IsZTT==0");
    
    std::cout<<" DATA size "<<DATA->GetEntries()<<"\n";
    std::cout<<" TTT size "<<TTT->GetEntries()<<"\n";
    std::cout<<" TTJ size "<<TTJ->GetEntries()<<"\n";
    std::cout<<" W size "<<W->GetEntries()<<"\n";
    std::cout<<" ZTT size "<<ZTT->GetEntries()<<"\n";
    std::cout<<" ZL size "<<ZL->GetEntries()<<"\n";
    std::cout<<" ZJ size "<<ZJ->GetEntries()<<"\n";
    std::cout<<" EWK size" <<EWK->GetEntries()<<"\n";
    std::cout<<" ZVV size "<<ZVV->GetEntries()<<"\n";
    std::cout<<" VVT size "<<VVT->GetEntries()<<"\n";
    std::cout<<" VVJ size "<<VVJ->GetEntries()<<"\n";
    std::cout<<" ZHTauTau size "<<ZHTauTau->GetEntries()<<"\n";
    std::cout<<" GluGluHTauTau size "<<GluGluHTauTau->GetEntries()<<"\n";
    std::cout<<" VBFHTauTau size "<<VBFHTauTau->GetEntries()<<"\n";

}


void setup_files_eleTau()
{
    std::cout<<" setup_files_eleTau TChains set up .... \n";
    fileLoad(DATA, "/store/user/gfunk/AugustProductionV1_8_0_26p1_V3/DATA/SingleElectron/", "davis_syncTree_BASELINE_EleTau.root");
    setup_MC_files("davis_syncTree_BASELINE_EleTau.root");
}

void setup_files_muTau()
{
    std::cout<<" setup_files_muTau TChains set up .... \n";
    fileLoad(DATA, "/store/user/gfunk/AugustProductionV1_8_0_26p1_V3/DATA/SingleMuon/", "davis_syncTree_BASELINE_MuTau.root");
    setup_MC_files("davis_syncTree_BASELINE_MuTau.root");
}

void setup_files_TauTau()
{
    std::cout<<" setup_files_TauTau TChains set up .... \n";
    fileLoad(DATA, "/store/user/gfunk/AugustProductionV1_8_0_26p1_V3/DATA/Tau/", "davis_syncTree_BASELINE_TauTau.root");
    setup_MC_files("davis_syncTree_BASELINE_TauTau.root");
}

void testComplete()
{

    fileLoad(DATA, "/store/user/gfunk/AugustProductionV1_8_0_26p1_V3/DATA/Tau/", "davis_syncTree_BASELINE_TauTau.root");
    fileLoad(DATA, "/store/user/gfunk/AugustProductionV1_8_0_26p1_V3/DATA/SingleElectron/", "davis_syncTree_BASELINE_EleTau.root");
    fileLoad(DATA, "/store/user/gfunk/AugustProductionV1_8_0_26p1_V3/DATA/SingleMuon/", "davis_syncTree_BASELINE_MuTau.root");

    setup_MC_files("davis_syncTree_BASELINE_TauTau.root");
    setup_MC_files("davis_syncTree_BASELINE_EleTau.root");
    setup_MC_files("davis_syncTree_BASELINE_MuTau.root");

    reset_files();

}


void reset_files()
{
    DATA->Reset();
    DY->Reset();
    W->Reset();
    VV->Reset();
    ZVV->Reset();
    EWK->Reset();
    TT->Reset();

    MZP600_MA0300->Reset();
    MZP600_MA0400->Reset();

    MZP800_MA0300->Reset();
    MZP800_MA0400->Reset();
    MZP800_MA0500->Reset();
    MZP800_MA0600->Reset();

    MZP1000_MA0300->Reset();
    MZP1000_MA0400->Reset();
    MZP1000_MA0500->Reset();
    MZP1000_MA0600->Reset();
    MZP1000_MA0700->Reset();
    MZP1000_MA0800->Reset();

    MZP1200_MA0300->Reset();
    MZP1200_MA0400->Reset();
    MZP1200_MA0500->Reset();
    MZP1200_MA0600->Reset();
    MZP1200_MA0700->Reset();
    MZP1200_MA0800->Reset();

    MZP1400_MA0300->Reset();
    MZP1400_MA0400->Reset();
    MZP1400_MA0500->Reset();
    MZP1400_MA0600->Reset();
    MZP1400_MA0700->Reset();
    MZP1400_MA0800->Reset();

    MZP1700_MA0300->Reset();
    MZP1700_MA0400->Reset();
    MZP1700_MA0500->Reset();
    MZP1700_MA0600->Reset();
    MZP1700_MA0700->Reset();
    MZP1700_MA0800->Reset();

    MZP2000_MA0300->Reset();
    MZP2000_MA0400->Reset();
    MZP2000_MA0500->Reset();
    MZP2000_MA0600->Reset();
    MZP2000_MA0700->Reset();
    MZP2000_MA0800->Reset();

    MZP2500_MA0300->Reset();
    MZP2500_MA0400->Reset();
    MZP2500_MA0500->Reset();
    MZP2500_MA0600->Reset();
    MZP2500_MA0700->Reset();
    MZP2500_MA0800->Reset();
    
    MZP10_MX1->Reset();
    MZP20_MX1->Reset();
    MZP50_MX1->Reset();
    MZP100_MX1->Reset();
    MZP300_MX1->Reset();
    MZP500_MX1->Reset();
    MZP1000_MX1->Reset();
    MZP2000_MX1->Reset();
    MZP10000_MX1->Reset();
    
    MZP10_MX50->Reset();
    MZP50_MX50->Reset();
    MZP95_MX50->Reset();
    MZP200_MX50->Reset();
    MZP300_MX50->Reset();
    MZP10000_MX50->Reset();

    MZP10_MX150->Reset();
    MZP200_MX150->Reset();
    MZP295_MX150->Reset();
    MZP500_MX150->Reset();
    MZP1000_MX150->Reset();
    MZP10000_MX150->Reset();
    
    MZP10_MX500->Reset();
    MZP500_MX500->Reset();
    MZP995_MX500->Reset();
    MZP2000_MX500->Reset();
    MZP10000_MX500->Reset();
    
    MZP10_MX1000->Reset();
    MZP1000_MX1000->Reset();
    MZP1995_MX1000->Reset();
    MZP10000_MX1000->Reset();

    //BKG subgroups
    //W
    WJetsToLNu->Reset();
    W1JetsToLNu->Reset();
    W2JetsToLNu->Reset();
    W3JetsToLNu->Reset();
    W4JetsToLNu->Reset();

    WJetsToLNuHT100To200->Reset();
    WJetsToLNuHT200To400->Reset();
    WJetsToLNuHT400To600->Reset();
    WJetsToLNuHT600To800->Reset();
    WJetsToLNuHT800To1200->Reset();
    WJetsToLNuHT1200To2500->Reset();
    WJetsToLNuHT2500ToInf->Reset();

    //EWK
    EWKWMinus2Jets->Reset();
    EWKWPlus2Jets->Reset();
    EWKZ2Jets_ZToLL->Reset();
    EWKZ2Jets_ZToNuNu->Reset();

    //DY
    DYinc->Reset();
    DY1Jets->Reset();
    DY2Jets->Reset();
    DY3Jets->Reset();
    DY4Jets->Reset();

    //ZNuNu
    ZJetsToNuNuHT100To200->Reset();
    ZJetsToNuNuHT200To400->Reset();
    ZJetsToNuNuHT400To600->Reset();
    ZJetsToNuNuHT600To800->Reset();
    ZJetsToNuNuHT800To1200->Reset();
    ZJetsToNuNuHT1200To2500->Reset();
    ZJetsToNuNuHT2500ToInf->Reset();

    //DY channels divided using bools
    ZTT->Reset();
    ZJ->Reset();
    ZL->Reset();

    //TT subgroups
    TTT->Reset();
    TTJ->Reset();

    //SMH subgroups
    //H
    GluGluHTauTau->Reset();
    VBFHTauTau->Reset();
    WminusHTauTau->Reset();
    WplusHTauTau->Reset();
    ZHTauTau->Reset();
    GluGluHToWW->Reset();
    VBFHToWW->Reset();
    ttHJetToNonbb->Reset();

    //VV
    VVT->Reset();
    VVJ->Reset();
    
    //VV subgroups
    ST_tW_antitop_5f_inclusiveDecays->Reset();
    ST_tW_top_5f_inclusiveDecays->Reset();
    ST_t_channel_antitop_4f_inclusiveDecays->Reset();
    ST_t_channel_top_4f_inclusiveDecays->Reset();
    
    WWTo2L2Nu->Reset();
    ZZTo2L2Nu->Reset();
    WWTo1L1Nu2Q->Reset();
    WZTo1L1Nu2Q->Reset();
    WZTo1L3Nu->Reset();
    WZTo2L2Q->Reset();
    ZZTo2L2Q->Reset();
    ZZTo2Q2Nu->Reset();
    ZZTo4L->Reset();
    WWW->Reset();
    WWZ->Reset();
    WZZ->Reset();
    ZZZ->Reset();
    
    TTZToLLNuNu->Reset();
    TTWJetsToLNu->Reset();
    WWTo2L2Nu_DoubleScattering->Reset();
    WWToLNuQQ->Reset();
    WWJJToLNuLNu_EWK->Reset();
    WpWpJJ_QCD->Reset();
    WpWpJJ_EWK->Reset();
    GluGluWWTo2L2Nu->Reset();
    WZTo3LNu->Reset();
    
}
