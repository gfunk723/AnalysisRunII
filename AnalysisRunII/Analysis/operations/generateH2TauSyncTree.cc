#include "generateH2TauSyncTree.h"
#include <iostream>
#include <fstream>
#include <string>

generateH2TauSyncTree::generateH2TauSyncTree(FlatTreeReader R_, bool run_, std::string fileOutName_)
{
	m_run = run_;
	R = R_;

	// init counters
	num_total = 0;
	num_et = 0;
	num_em = 0;
	num_tt = 0;
	num_mt = 0;

	if(m_run)
	{
	std::string MuTauName = "davis_syncTree_"+fileOutName_+"MuTau.root";
	std::string EleTauName = "davis_syncTree_"+fileOutName_+"EleTau.root";
	std::string TauTauName = "davis_syncTree_"+fileOutName_+"TauTau.root";
	std::string EleMuName = "davis_syncTree_"+fileOutName_+"EleMu.root";

	outFile_MuTau = new TFile(MuTauName.c_str(),"RECREATE");
	outFile_EleTau  = new TFile(EleTauName.c_str(),"RECREATE");
	outFile_TauTau  = new TFile(TauTauName.c_str(),"RECREATE");
	outFile_EleMu  = new TFile(EleMuName.c_str(),"RECREATE");

	tree_MuTau = new TTree("tree_MuTau","tree_MuTau");
	tree_EleTau = new TTree("tree_EleTau","tree_EleTau");
	tree_TauTau = new TTree("tree_TauTau","tree_TauTau");
	tree_EleMu = new TTree("tree_EleMu","tree_EleMu");

	setupBranches(tree_MuTau);
	setupBranches(tree_EleTau);
	setupBranches(tree_TauTau);
	setupBranches(tree_EleMu);
    
	/* init the susy signal pt reweight tool */
	NLO_ReadFile();	

	/* init qcd weight tool */
	// with DZeta cut ->
	//qcdWeights = new QCDModelForEMu("QCD_weight_emu.root");

	// w/o DZeta cut ->
	//qcdWeightsNoDZeta = new QCDModelForEMu("QCD_weight_emu_nodzeta.root");

	/* init HTT lepton sf tool */

	sfTool_Muon_IdIso0p15_eff = new ScaleFactor();
	sfTool_Muon_IdIso0p20_eff = new ScaleFactor();
    
	sfTool_Electron_IdIso0p10_eff = new ScaleFactor();
	sfTool_Electron_SingleEle_eff = new ScaleFactor();
    
	sfTool_Electron_IdIso0p15_eff = new ScaleFactor();
    sfTool_Muon_SingleMu_eff = new ScaleFactor();
    
	//sfTool_Muon_Mu8_eff = new ScaleFactor();
	//sfTool_Muon_Mu17_eff = new ScaleFactor();
	//sfTool_Electron_Ele17_eff = new ScaleFactor();
	//sfTool_Electron_Ele12_eff = new ScaleFactor();

	sfTool_Muon_IdIso0p15_eff->init_ScaleFactor("Muon_IdIso_IsoLt0p15_2016BtoH_eff.root"); // 2016
	sfTool_Muon_IdIso0p20_eff->init_ScaleFactor("Muon_IdIso_IsoLt0p2_2016BtoH_eff.root"); // 2016
	
	sfTool_Electron_IdIso0p10_eff->init_ScaleFactor("Electron_IdIso_IsoLt0p1_eff.root"); // 2016
	sfTool_Electron_IdIso0p15_eff->init_ScaleFactor("Electron_IdIso_IsoLt0p15_eff.root"); // 2016

	sfTool_Muon_SingleMu_eff->init_ScaleFactor("Muon_IsoMu24_OR_TkIsoMu24_2016BtoH_eff.root"); // 2016
    
	sfTool_Electron_SingleEle_eff->init_ScaleFactor("Electron_Ele25WPTight_eff.root"); // 2016
    
	//sfTool_Muon_Mu8_eff->init_ScaleFactor("Muon_Mu8leg_2016BtoH_eff.root"); //2016
	//sfTool_Muon_Mu17_eff->init_ScaleFactor("Muon_Mu17_eff.root");
	//sfTool_Electron_Ele17_eff->init_ScaleFactor("Electron_Ele17_eff.root");
	//sfTool_Electron_Ele12_eff->init_ScaleFactor("Electron_Ele12_eff.root"); // 2016
    
    //Corr Workspace Init
    TFile t("htt_scalefactors_v16_4.root");
    tw = (RooWorkspace*)t.Get("w");
    t.Close();
    
    TFile f("htt_scalefactors_sm_moriond_v2.root");
    w = (RooWorkspace*)f.Get("w");
    f.Close();
    
    /* k factor initialization */

    EWK_Wcorr->Divide(LO_Wcorr);
    EWK_Wcorr_dNLO->Divide(NLO_Wcorr);
    
	}

	else
	{
		outFile_MuTau = nullptr;
		outFile_EleTau  = nullptr; 
		outFile_TauTau  = nullptr;
		outFile_EleMu  = nullptr;

		tree_MuTau = nullptr;
		tree_EleTau = nullptr;
		tree_TauTau = nullptr;
		tree_EleMu = nullptr;
	}


	initScaleFactorParametersRunII();
    
    //mt_MZP600A0400_reader->AddVariable( "read_pt_1", &rpt_1 );
    mt_MZP600A0400_reader->AddVariable( "read_pt_2", &rpt_2 );
    mt_MZP600A0400_reader->AddVariable( "read_mt_tot", &rmt_tot );
    mt_MZP600A0400_reader->AddVariable( "read_pfmt_1", &rpfmt_1 );
    mt_MZP600A0400_reader->AddVariable( "read_pt_tt", &rpt_tt );
    //mt_MZP600A0400_reader->AddVariable( "read_m_vis", &rm_vis );
    mt_MZP600A0400_reader->AddVariable( "read_met", &rmet );
    mt_MZP600A0400_reader->AddVariable( "read_P_chi_pf", &rP_chi_pf );
    mt_MZP600A0400_reader->AddVariable( "read_LPT", &rLPT );
    mt_MZP600A0400_reader->AddVariable( "read_DeltaR_leg1_leg2", &rDeltaR_leg1_leg2 );
    mt_MZP600A0400_reader->AddVariable( "read_cos_DeltaPhi_leg1_leg2", &rcos_DeltaPhi_leg1_leg2 );
    mt_MZP600A0400_reader->AddVariable( "read_cos_DeltaPhi_PFMET_Higgs", &rcos_DeltaPhi_PFMET_Higgs );
    
    mt_MZP600A0400_reader->AddSpectator( "read_ratio_weight", &rratio_weight );
    mt_MZP600A0400_reader->AddSpectator( "read_final_weight", &rfinal_weight );
    mt_MZP600A0400_reader->AddSpectator( "read_npu", &rnpu );
    mt_MZP600A0400_reader->AddSpectator( "read_event", &revent );
    mt_MZP600A0400_reader->AddSpectator( "read_randNum", &rrandNum );
    mt_MZP600A0400_reader->AddSpectator( "read_DataCardInt", &rDataCardInt );
    mt_MZP600A0400_reader->AddSpectator( "read_IsZTT", &rIsZTT );
    mt_MZP600A0400_reader->AddSpectator( "read_IsZJ", &rIsZJ );
    mt_MZP600A0400_reader->AddSpectator( "read_IsZL", &rIsZL );
    mt_MZP600A0400_reader->AddSpectator( "read_IsTTT", &rIsTTT );
    
    mt_MZP600A0400_reader->BookMVA("MLPBNN","TMVAClassification_MLPBNN_MuTau_MZP600_MA0400_MDM100_nodes30_epoch1000.weights.xml");
    
    mt_MZP800A0400_reader->AddVariable( "read_pt_1", &rpt_1 );
    mt_MZP800A0400_reader->AddVariable( "read_pt_2", &rpt_2 );
    mt_MZP800A0400_reader->AddVariable( "read_mt_tot", &rmt_tot );
    mt_MZP800A0400_reader->AddVariable( "read_pfmt_1", &rpfmt_1 );
    mt_MZP800A0400_reader->AddVariable( "read_pt_tt", &rpt_tt );
    //mt_MZP800A0400_reader->AddVariable( "read_m_vis", &rm_vis );
    mt_MZP800A0400_reader->AddVariable( "read_met", &rmet );
    mt_MZP800A0400_reader->AddVariable( "read_P_chi_pf", &rP_chi_pf );
    mt_MZP800A0400_reader->AddVariable( "read_LPT", &rLPT );
    mt_MZP800A0400_reader->AddVariable( "read_DeltaR_leg1_leg2", &rDeltaR_leg1_leg2 );
    mt_MZP800A0400_reader->AddVariable( "read_cos_DeltaPhi_leg1_leg2", &rcos_DeltaPhi_leg1_leg2 );
    mt_MZP800A0400_reader->AddVariable( "read_cos_DeltaPhi_PFMET_Higgs", &rcos_DeltaPhi_PFMET_Higgs );
    
    mt_MZP800A0400_reader->AddSpectator( "read_ratio_weight", &rratio_weight );
    mt_MZP800A0400_reader->AddSpectator( "read_final_weight", &rfinal_weight );
    mt_MZP800A0400_reader->AddSpectator( "read_npu", &rnpu );
    mt_MZP800A0400_reader->AddSpectator( "read_event", &revent );
    mt_MZP800A0400_reader->AddSpectator( "read_randNum", &rrandNum );
    mt_MZP800A0400_reader->AddSpectator( "read_DataCardInt", &rDataCardInt );
    mt_MZP800A0400_reader->AddSpectator( "read_IsZTT", &rIsZTT );
    mt_MZP800A0400_reader->AddSpectator( "read_IsZJ", &rIsZJ );
    mt_MZP800A0400_reader->AddSpectator( "read_IsZL", &rIsZL );
    mt_MZP800A0400_reader->AddSpectator( "read_IsTTT", &rIsTTT );
    
    mt_MZP800A0400_reader->BookMVA("MLPBNN","TMVAClassification_MLPBNN_MuTau_MZP800_MA0400_MDM100_nodes30_epoch1000.weights.xml");

    mt_MZP1000A0400_reader->AddVariable( "read_pt_1", &rpt_1 );
    mt_MZP1000A0400_reader->AddVariable( "read_pt_2", &rpt_2 );
    mt_MZP1000A0400_reader->AddVariable( "read_mt_tot", &rmt_tot );
    mt_MZP1000A0400_reader->AddVariable( "read_pfmt_1", &rpfmt_1 );
    mt_MZP1000A0400_reader->AddVariable( "read_pt_tt", &rpt_tt );
    mt_MZP1000A0400_reader->AddVariable( "read_m_vis", &rm_vis );
    mt_MZP1000A0400_reader->AddVariable( "read_met", &rmet );
    mt_MZP1000A0400_reader->AddVariable( "read_P_chi_pf", &rP_chi_pf );
    mt_MZP1000A0400_reader->AddVariable( "read_LPT", &rLPT );
    mt_MZP1000A0400_reader->AddVariable( "read_DeltaR_leg1_leg2", &rDeltaR_leg1_leg2 );
    mt_MZP1000A0400_reader->AddVariable( "read_cos_DeltaPhi_leg1_leg2", &rcos_DeltaPhi_leg1_leg2 );
    mt_MZP1000A0400_reader->AddVariable( "read_cos_DeltaPhi_PFMET_Higgs", &rcos_DeltaPhi_PFMET_Higgs );
    
    mt_MZP1000A0400_reader->AddSpectator( "read_ratio_weight", &rratio_weight );
    mt_MZP1000A0400_reader->AddSpectator( "read_final_weight", &rfinal_weight );
    mt_MZP1000A0400_reader->AddSpectator( "read_npu", &rnpu );
    mt_MZP1000A0400_reader->AddSpectator( "read_event", &revent );
    mt_MZP1000A0400_reader->AddSpectator( "read_randNum", &rrandNum );
    mt_MZP1000A0400_reader->AddSpectator( "read_DataCardInt", &rDataCardInt );
    mt_MZP1000A0400_reader->AddSpectator( "read_IsZTT", &rIsZTT );
    mt_MZP1000A0400_reader->AddSpectator( "read_IsZJ", &rIsZJ );
    mt_MZP1000A0400_reader->AddSpectator( "read_IsZL", &rIsZL );
    mt_MZP1000A0400_reader->AddSpectator( "read_IsTTT", &rIsTTT );
    
    mt_MZP1000A0400_reader->BookMVA("MLPBNN","TMVAClassification_MLPBNN_MuTau_MZP1000_MA0400_MDM100_nodes30_epoch1000.weights.xml");
    
    mt_MZP1200A0400_reader->AddVariable( "read_pt_1", &rpt_1 );
    mt_MZP1200A0400_reader->AddVariable( "read_pt_2", &rpt_2 );
    mt_MZP1200A0400_reader->AddVariable( "read_mt_tot", &rmt_tot );
    mt_MZP1200A0400_reader->AddVariable( "read_pfmt_1", &rpfmt_1 );
    mt_MZP1200A0400_reader->AddVariable( "read_pt_tt", &rpt_tt );
    mt_MZP1200A0400_reader->AddVariable( "read_m_vis", &rm_vis );
    mt_MZP1200A0400_reader->AddVariable( "read_met", &rmet );
    mt_MZP1200A0400_reader->AddVariable( "read_P_chi_pf", &rP_chi_pf );
    mt_MZP1200A0400_reader->AddVariable( "read_LPT", &rLPT );
    mt_MZP1200A0400_reader->AddVariable( "read_DeltaR_leg1_leg2", &rDeltaR_leg1_leg2 );
    mt_MZP1200A0400_reader->AddVariable( "read_cos_DeltaPhi_leg1_leg2", &rcos_DeltaPhi_leg1_leg2 );
    mt_MZP1200A0400_reader->AddVariable( "read_cos_DeltaPhi_PFMET_Higgs", &rcos_DeltaPhi_PFMET_Higgs );
    
    mt_MZP1200A0400_reader->AddSpectator( "read_ratio_weight", &rratio_weight );
    mt_MZP1200A0400_reader->AddSpectator( "read_final_weight", &rfinal_weight );
    mt_MZP1200A0400_reader->AddSpectator( "read_npu", &rnpu );
    mt_MZP1200A0400_reader->AddSpectator( "read_event", &revent );
    mt_MZP1200A0400_reader->AddSpectator( "read_randNum", &rrandNum );
    mt_MZP1200A0400_reader->AddSpectator( "read_DataCardInt", &rDataCardInt );
    mt_MZP1200A0400_reader->AddSpectator( "read_IsZTT", &rIsZTT );
    mt_MZP1200A0400_reader->AddSpectator( "read_IsZJ", &rIsZJ );
    mt_MZP1200A0400_reader->AddSpectator( "read_IsZL", &rIsZL );
    mt_MZP1200A0400_reader->AddSpectator( "read_IsTTT", &rIsTTT );
    
    mt_MZP1200A0400_reader->BookMVA("MLPBNN","TMVAClassification_MLPBNN_MuTau_MZP1200_MA0400_MDM100_nodes30_epoch1000.weights.xml");
    
    //et_MZP600A0400_reader->AddVariable( "read_pt_1", &rpt_1 );
    et_MZP600A0400_reader->AddVariable( "read_pt_2", &rpt_2 );
    et_MZP600A0400_reader->AddVariable( "read_mt_tot", &rmt_tot );
    et_MZP600A0400_reader->AddVariable( "read_pfmt_1", &rpfmt_1 );
    et_MZP600A0400_reader->AddVariable( "read_pt_tt", &rpt_tt );
    //et_MZP600A0400_reader->AddVariable( "read_m_vis", &rm_vis );
    et_MZP600A0400_reader->AddVariable( "read_met", &rmet );
    et_MZP600A0400_reader->AddVariable( "read_P_chi_pf", &rP_chi_pf );
    et_MZP600A0400_reader->AddVariable( "read_LPT", &rLPT );
    et_MZP600A0400_reader->AddVariable( "read_DeltaR_leg1_leg2", &rDeltaR_leg1_leg2 );
    et_MZP600A0400_reader->AddVariable( "read_cos_DeltaPhi_leg1_leg2", &rcos_DeltaPhi_leg1_leg2 );
    et_MZP600A0400_reader->AddVariable( "read_cos_DeltaPhi_PFMET_Higgs", &rcos_DeltaPhi_PFMET_Higgs );
    
    et_MZP600A0400_reader->AddSpectator( "read_ratio_weight", &rratio_weight );
    et_MZP600A0400_reader->AddSpectator( "read_final_weight", &rfinal_weight );
    et_MZP600A0400_reader->AddSpectator( "read_npu", &rnpu );
    et_MZP600A0400_reader->AddSpectator( "read_event", &revent );
    et_MZP600A0400_reader->AddSpectator( "read_randNum", &rrandNum );
    et_MZP600A0400_reader->AddSpectator( "read_DataCardInt", &rDataCardInt );
    et_MZP600A0400_reader->AddSpectator( "read_IsZTT", &rIsZTT );
    et_MZP600A0400_reader->AddSpectator( "read_IsZJ", &rIsZJ );
    et_MZP600A0400_reader->AddSpectator( "read_IsZL", &rIsZL );
    et_MZP600A0400_reader->AddSpectator( "read_IsTTT", &rIsTTT );

    et_MZP600A0400_reader->BookMVA("MLPBNN","TMVAClassification_MLPBNN_EleTau_MZP600_MA0400_MDM100_nodes30_epoch1000.weights.xml");
    
    et_MZP800A0400_reader->AddVariable( "read_pt_1", &rpt_1 );
    et_MZP800A0400_reader->AddVariable( "read_pt_2", &rpt_2 );
    et_MZP800A0400_reader->AddVariable( "read_mt_tot", &rmt_tot );
    et_MZP800A0400_reader->AddVariable( "read_pfmt_1", &rpfmt_1 );
    et_MZP800A0400_reader->AddVariable( "read_pt_tt", &rpt_tt );
    //et_MZP800A0400_reader->AddVariable( "read_m_vis", &rm_vis );
    et_MZP800A0400_reader->AddVariable( "read_met", &rmet );
    et_MZP800A0400_reader->AddVariable( "read_P_chi_pf", &rP_chi_pf );
    et_MZP800A0400_reader->AddVariable( "read_LPT", &rLPT );
    et_MZP800A0400_reader->AddVariable( "read_DeltaR_leg1_leg2", &rDeltaR_leg1_leg2 );
    et_MZP800A0400_reader->AddVariable( "read_cos_DeltaPhi_leg1_leg2", &rcos_DeltaPhi_leg1_leg2 );
    et_MZP800A0400_reader->AddVariable( "read_cos_DeltaPhi_PFMET_Higgs", &rcos_DeltaPhi_PFMET_Higgs );
    
    et_MZP800A0400_reader->AddSpectator( "read_ratio_weight", &rratio_weight );
    et_MZP800A0400_reader->AddSpectator( "read_final_weight", &rfinal_weight );
    et_MZP800A0400_reader->AddSpectator( "read_npu", &rnpu );
    et_MZP800A0400_reader->AddSpectator( "read_event", &revent );
    et_MZP800A0400_reader->AddSpectator( "read_randNum", &rrandNum );
    et_MZP800A0400_reader->AddSpectator( "read_DataCardInt", &rDataCardInt );
    et_MZP800A0400_reader->AddSpectator( "read_IsZTT", &rIsZTT );
    et_MZP800A0400_reader->AddSpectator( "read_IsZJ", &rIsZJ );
    et_MZP800A0400_reader->AddSpectator( "read_IsZL", &rIsZL );
    et_MZP800A0400_reader->AddSpectator( "read_IsTTT", &rIsTTT );

    et_MZP800A0400_reader->BookMVA("MLPBNN","TMVAClassification_MLPBNN_EleTau_MZP800_MA0400_MDM100_nodes30_epoch1000.weights.xml");
    
    
    et_MZP1000A0400_reader->AddVariable( "read_pt_1", &rpt_1 );
    et_MZP1000A0400_reader->AddVariable( "read_pt_2", &rpt_2 );
    et_MZP1000A0400_reader->AddVariable( "read_mt_tot", &rmt_tot );
    et_MZP1000A0400_reader->AddVariable( "read_pfmt_1", &rpfmt_1 );
    et_MZP1000A0400_reader->AddVariable( "read_pt_tt", &rpt_tt );
    //et_MZP1000A0400_reader->AddVariable( "read_m_vis", &rm_vis );
    et_MZP1000A0400_reader->AddVariable( "read_met", &rmet );
    et_MZP1000A0400_reader->AddVariable( "read_P_chi_pf", &rP_chi_pf );
    et_MZP1000A0400_reader->AddVariable( "read_LPT", &rLPT );
    et_MZP1000A0400_reader->AddVariable( "read_DeltaR_leg1_leg2", &rDeltaR_leg1_leg2 );
    et_MZP1000A0400_reader->AddVariable( "read_cos_DeltaPhi_leg1_leg2", &rcos_DeltaPhi_leg1_leg2 );
    et_MZP1000A0400_reader->AddVariable( "read_cos_DeltaPhi_PFMET_Higgs", &rcos_DeltaPhi_PFMET_Higgs );
    
    et_MZP1000A0400_reader->AddSpectator( "read_ratio_weight", &rratio_weight );
    et_MZP1000A0400_reader->AddSpectator( "read_final_weight", &rfinal_weight );
    et_MZP1000A0400_reader->AddSpectator( "read_npu", &rnpu );
    et_MZP1000A0400_reader->AddSpectator( "read_event", &revent );
    et_MZP1000A0400_reader->AddSpectator( "read_randNum", &rrandNum );
    et_MZP1000A0400_reader->AddSpectator( "read_DataCardInt", &rDataCardInt );
    et_MZP1000A0400_reader->AddSpectator( "read_IsZTT", &rIsZTT );
    et_MZP1000A0400_reader->AddSpectator( "read_IsZJ", &rIsZJ );
    et_MZP1000A0400_reader->AddSpectator( "read_IsZL", &rIsZL );
    et_MZP1000A0400_reader->AddSpectator( "read_IsTTT", &rIsTTT );

    et_MZP1000A0400_reader->BookMVA("MLPBNN","TMVAClassification_MLPBNN_EleTau_MZP1000_MA0400_MDM100_nodes30_epoch1000.weights.xml");
    
    et_MZP1200A0400_reader->AddVariable( "read_pt_1", &rpt_1 );
    et_MZP1200A0400_reader->AddVariable( "read_pt_2", &rpt_2 );
    et_MZP1200A0400_reader->AddVariable( "read_mt_tot", &rmt_tot );
    et_MZP1200A0400_reader->AddVariable( "read_pfmt_1", &rpfmt_1 );
    et_MZP1200A0400_reader->AddVariable( "read_pt_tt", &rpt_tt );
    et_MZP1200A0400_reader->AddVariable( "read_m_vis", &rm_vis );
    et_MZP1200A0400_reader->AddVariable( "read_met", &rmet );
    et_MZP1200A0400_reader->AddVariable( "read_P_chi_pf", &rP_chi_pf );
    et_MZP1200A0400_reader->AddVariable( "read_LPT", &rLPT );
    et_MZP1200A0400_reader->AddVariable( "read_DeltaR_leg1_leg2", &rDeltaR_leg1_leg2 );
    et_MZP1200A0400_reader->AddVariable( "read_cos_DeltaPhi_leg1_leg2", &rcos_DeltaPhi_leg1_leg2 );
    et_MZP1200A0400_reader->AddVariable( "read_cos_DeltaPhi_PFMET_Higgs", &rcos_DeltaPhi_PFMET_Higgs );
    
    et_MZP1200A0400_reader->AddSpectator( "read_ratio_weight", &rratio_weight );
    et_MZP1200A0400_reader->AddSpectator( "read_final_weight", &rfinal_weight );
    et_MZP1200A0400_reader->AddSpectator( "read_npu", &rnpu );
    et_MZP1200A0400_reader->AddSpectator( "read_event", &revent );
    et_MZP1200A0400_reader->AddSpectator( "read_randNum", &rrandNum );
    et_MZP1200A0400_reader->AddSpectator( "read_DataCardInt", &rDataCardInt );
    et_MZP1200A0400_reader->AddSpectator( "read_IsZTT", &rIsZTT );
    et_MZP1200A0400_reader->AddSpectator( "read_IsZJ", &rIsZJ );
    et_MZP1200A0400_reader->AddSpectator( "read_IsZL", &rIsZL );
    et_MZP1200A0400_reader->AddSpectator( "read_IsTTT", &rIsTTT );

    et_MZP1200A0400_reader->BookMVA("MLPBNN","TMVAClassification_MLPBNN_EleTau_MZP1200_MA0400_MDM100_nodes30_epoch1000.weights.xml");

    //tt_MZP600A0400_reader->AddVariable( "read_pt_1", &rpt_1 );
    //tt_MZP600A0400_reader->AddVariable( "read_pt_2", &rpt_2 );
    tt_MZP600A0400_reader->AddVariable( "read_mt_tot", &rmt_tot );
    tt_MZP600A0400_reader->AddVariable( "read_pfmt_1", &rpfmt_1 );
    tt_MZP600A0400_reader->AddVariable( "read_pt_tt", &rpt_tt );
    //tt_MZP600A0400_reader->AddVariable( "read_m_vis", &rm_vis );
    tt_MZP600A0400_reader->AddVariable( "read_met", &rmet );
    tt_MZP600A0400_reader->AddVariable( "read_LPT", &rLPT );
    tt_MZP600A0400_reader->AddVariable( "read_DeltaR_leg1_leg2", &rDeltaR_leg1_leg2 );
    tt_MZP600A0400_reader->AddVariable( "read_cos_DeltaPhi_leg1_leg2", &rcos_DeltaPhi_leg1_leg2 );
    tt_MZP600A0400_reader->AddVariable( "read_cos_DeltaPhi_PFMET_Higgs", &rcos_DeltaPhi_PFMET_Higgs );
    
    tt_MZP600A0400_reader->AddSpectator( "read_ratio_weight", &rratio_weight );
    tt_MZP600A0400_reader->AddSpectator( "read_final_weight", &rfinal_weight );
    tt_MZP600A0400_reader->AddSpectator( "read_npu", &rnpu );
    tt_MZP600A0400_reader->AddSpectator( "read_event", &revent );
    tt_MZP600A0400_reader->AddSpectator( "read_randNum", &rrandNum );
    tt_MZP600A0400_reader->AddSpectator( "read_DataCardInt", &rDataCardInt );
    tt_MZP600A0400_reader->AddSpectator( "read_IsZTT", &rIsZTT );
    tt_MZP600A0400_reader->AddSpectator( "read_IsZJ", &rIsZJ );
    tt_MZP600A0400_reader->AddSpectator( "read_IsZL", &rIsZL );
    tt_MZP600A0400_reader->AddSpectator( "read_IsTTT", &rIsTTT );

    tt_MZP600A0400_reader->BookMVA("MLPBNN","TMVAClassification_MLPBNN_TauTau_MZP600_MA0400_MDM100_nodes30_epoch1000.weights.xml");
    
    tt_MZP800A0400_reader->AddVariable( "read_pt_1", &rpt_1 );
    tt_MZP800A0400_reader->AddVariable( "read_pt_2", &rpt_2 );
    tt_MZP800A0400_reader->AddVariable( "read_mt_tot", &rmt_tot );
    tt_MZP800A0400_reader->AddVariable( "read_pfmt_1", &rpfmt_1 );
    tt_MZP800A0400_reader->AddVariable( "read_pt_tt", &rpt_tt );
    tt_MZP800A0400_reader->AddVariable( "read_m_vis", &rm_vis );
    tt_MZP800A0400_reader->AddVariable( "read_met", &rmet );
    tt_MZP800A0400_reader->AddVariable( "read_LPT", &rLPT );
    tt_MZP800A0400_reader->AddVariable( "read_DeltaR_leg1_leg2", &rDeltaR_leg1_leg2 );
    tt_MZP800A0400_reader->AddVariable( "read_cos_DeltaPhi_leg1_leg2", &rcos_DeltaPhi_leg1_leg2 );
    tt_MZP800A0400_reader->AddVariable( "read_cos_DeltaPhi_PFMET_Higgs", &rcos_DeltaPhi_PFMET_Higgs );
    
    tt_MZP800A0400_reader->AddSpectator( "read_ratio_weight", &rratio_weight );
    tt_MZP800A0400_reader->AddSpectator( "read_final_weight", &rfinal_weight );
    tt_MZP800A0400_reader->AddSpectator( "read_npu", &rnpu );
    tt_MZP800A0400_reader->AddSpectator( "read_event", &revent );
    tt_MZP800A0400_reader->AddSpectator( "read_randNum", &rrandNum );
    tt_MZP800A0400_reader->AddSpectator( "read_DataCardInt", &rDataCardInt );
    tt_MZP800A0400_reader->AddSpectator( "read_IsZTT", &rIsZTT );
    tt_MZP800A0400_reader->AddSpectator( "read_IsZJ", &rIsZJ );
    tt_MZP800A0400_reader->AddSpectator( "read_IsZL", &rIsZL );
    tt_MZP800A0400_reader->AddSpectator( "read_IsTTT", &rIsTTT );

    tt_MZP800A0400_reader->BookMVA("MLPBNN","TMVAClassification_MLPBNN_TauTau_MZP800_MA0400_MDM100_nodes30_epoch1000.weights.xml");

    tt_MZP1000A0400_reader->AddVariable( "read_pt_1", &rpt_1 );
    tt_MZP1000A0400_reader->AddVariable( "read_pt_2", &rpt_2 );
    tt_MZP1000A0400_reader->AddVariable( "read_mt_tot", &rmt_tot );
    tt_MZP1000A0400_reader->AddVariable( "read_pfmt_1", &rpfmt_1 );
    tt_MZP1000A0400_reader->AddVariable( "read_pt_tt", &rpt_tt );
    tt_MZP1000A0400_reader->AddVariable( "read_m_vis", &rm_vis );
    tt_MZP1000A0400_reader->AddVariable( "read_met", &rmet );
    tt_MZP1000A0400_reader->AddVariable( "read_LPT", &rLPT );
    tt_MZP1000A0400_reader->AddVariable( "read_DeltaR_leg1_leg2", &rDeltaR_leg1_leg2 );
    tt_MZP1000A0400_reader->AddVariable( "read_cos_DeltaPhi_leg1_leg2", &rcos_DeltaPhi_leg1_leg2 );
    tt_MZP1000A0400_reader->AddVariable( "read_cos_DeltaPhi_PFMET_Higgs", &rcos_DeltaPhi_PFMET_Higgs );
    
    tt_MZP1000A0400_reader->AddSpectator( "read_ratio_weight", &rratio_weight );
    tt_MZP1000A0400_reader->AddSpectator( "read_final_weight", &rfinal_weight );
    tt_MZP1000A0400_reader->AddSpectator( "read_npu", &rnpu );
    tt_MZP1000A0400_reader->AddSpectator( "read_event", &revent );
    tt_MZP1000A0400_reader->AddSpectator( "read_randNum", &rrandNum );
    tt_MZP1000A0400_reader->AddSpectator( "read_DataCardInt", &rDataCardInt );
    tt_MZP1000A0400_reader->AddSpectator( "read_IsZTT", &rIsZTT );
    tt_MZP1000A0400_reader->AddSpectator( "read_IsZJ", &rIsZJ );
    tt_MZP1000A0400_reader->AddSpectator( "read_IsZL", &rIsZL );
    tt_MZP1000A0400_reader->AddSpectator( "read_IsTTT", &rIsTTT );

    tt_MZP1000A0400_reader->BookMVA("MLPBNN","TMVAClassification_MLPBNN_TauTau_MZP1000_MA0400_MDM100_nodes30_epoch1000.weights.xml");
    
    tt_MZP1200A0400_reader->AddVariable( "read_pt_1", &rpt_1 );
    tt_MZP1200A0400_reader->AddVariable( "read_pt_2", &rpt_2 );
    tt_MZP1200A0400_reader->AddVariable( "read_mt_tot", &rmt_tot );
    tt_MZP1200A0400_reader->AddVariable( "read_pfmt_1", &rpfmt_1 );
    tt_MZP1200A0400_reader->AddVariable( "read_pt_tt", &rpt_tt );
    tt_MZP1200A0400_reader->AddVariable( "read_m_vis", &rm_vis );
    tt_MZP1200A0400_reader->AddVariable( "read_met", &rmet );
    tt_MZP1200A0400_reader->AddVariable( "read_LPT", &rLPT );
    tt_MZP1200A0400_reader->AddVariable( "read_DeltaR_leg1_leg2", &rDeltaR_leg1_leg2 );
    tt_MZP1200A0400_reader->AddVariable( "read_cos_DeltaPhi_leg1_leg2", &rcos_DeltaPhi_leg1_leg2 );
    tt_MZP1200A0400_reader->AddVariable( "read_cos_DeltaPhi_PFMET_Higgs", &rcos_DeltaPhi_PFMET_Higgs );
    
    tt_MZP1200A0400_reader->AddSpectator( "read_ratio_weight", &rratio_weight );
    tt_MZP1200A0400_reader->AddSpectator( "read_final_weight", &rfinal_weight );
    tt_MZP1200A0400_reader->AddSpectator( "read_npu", &rnpu );
    tt_MZP1200A0400_reader->AddSpectator( "read_event", &revent );
    tt_MZP1200A0400_reader->AddSpectator( "read_randNum", &rrandNum );
    tt_MZP1200A0400_reader->AddSpectator( "read_DataCardInt", &rDataCardInt );
    tt_MZP1200A0400_reader->AddSpectator( "read_IsZTT", &rIsZTT );
    tt_MZP1200A0400_reader->AddSpectator( "read_IsZJ", &rIsZJ );
    tt_MZP1200A0400_reader->AddSpectator( "read_IsZL", &rIsZL );
    tt_MZP1200A0400_reader->AddSpectator( "read_IsTTT", &rIsTTT );

    tt_MZP1200A0400_reader->BookMVA("MLPBNN","TMVAClassification_MLPBNN_TauTau_MZP1200_MA0400_MDM100_nodes30_epoch1000.weights.xml");
}



void generateH2TauSyncTree::finish()
{
	if(m_run)
	{
		outFile_MuTau->cd();
		tree_MuTau->CloneTree()->Write("TauCheck");	/* need to clone because of shared addresses */		
		outFile_MuTau->Close();

		outFile_EleTau->cd();
		tree_EleTau->CloneTree()->Write("TauCheck");		
		outFile_EleTau->Close();

		outFile_TauTau->cd();
		tree_TauTau->CloneTree()->Write("TauCheck");		
		outFile_TauTau->Close();		

		outFile_EleMu->cd();
		tree_EleMu->CloneTree()->Write("TauCheck");		
		outFile_EleMu->Close();

	}

	 std::cout<<"final count etau = "<<num_et<<"\n";
	 std::cout<<"final count mtau = "<<num_mt<<"\n";
	 std::cout<<"final count tt = "<<num_tt<<"\n";
	 std::cout<<"final count em = "<<num_em<<"\n";
}

generateH2TauSyncTree::~generateH2TauSyncTree()
{
	if(!m_run) /* Close on the outFile deletes the TH1F pointers */
	{
		delete tree_MuTau;
		delete tree_EleTau;
		delete tree_TauTau;
		delete tree_EleMu;
	}	
	delete outFile_MuTau;
	delete outFile_EleTau;
	delete outFile_TauTau;
	delete outFile_EleMu;
	

};	

void generateH2TauSyncTree::handleEvent()
{

	/* count the total events seen */

	num_total ++; 

	/* reset values before handling each event */

	reset();

    /* set tMVA flag type */
    float FRACTION_  = 0.1; /* this is how much of a mc sample is testing, and training - so 0.1 means 10% test, 10% train, 80% analysis */
    TRandom3 rand_;
    rand_.SetSeed(0); /* unique random pattern, once before the event loop */
    float randVal = rand_.Uniform(1.);
    
    if (randVal < FRACTION_)
    {
        flag_MVAEventType = 0;
    }
    else if ((randVal > FRACTION_) && (randVal < 2.*FRACTION_))
    {
        flag_MVAEventType = 1;
    }
    else {flag_MVAEventType = -1;}

    randNum = randVal;

	/* 
      note: two factors can control what variables are available in the FlatTuple
	  -- (1) both tau ES and electron ES are nominal (it will be nan for channels without e or tau)
	  -- (2) the FlatTuple was not produced under small tree conditions 
	
	  the following booleans will allow us to handle the differences in FlatTuple content without 
	  reader asserts
	*/

	bool eventHasNominalLeptonEnergyScales = (fabs(R.getF("TauEsNumberSigmasShifted")) != 1. && fabs(R.getF("ElectronEsNumberSigmasShifted")) != 1.);
	bool eventIsNotSmallTree = !(R.getB("isSmallTree"));


	/* declare & init 4-vectors for computation of event quantities */
    TLorentzVector l1(0.,0.,0.,0.); /* leg1 */
	TLorentzVector l2(0.,0.,0.,0.); /* leg2 */
    
	TLorentzVector l1NC(0.,0.,0.,0.); /* leg1 */
	TLorentzVector l2NC(0.,0.,0.,0.); /* leg2 */
    
    TLorentzVector l1TESUp(0.,0.,0.,0.);
    TLorentzVector l1TESDown(0.,0.,0.,0.);
    TLorentzVector l2TESUp(0.,0.,0.,0.);
    TLorentzVector l2TESDown(0.,0.,0.,0.);
    
    TLorentzVector l1TESUpNC(0.,0.,0.,0.);
    TLorentzVector l1TESDownNC(0.,0.,0.,0.);
    TLorentzVector l2TESUpNC(0.,0.,0.,0.);
    TLorentzVector l2TESDownNC(0.,0.,0.,0.);
    
    TLorentzVector TSDiff(0.,0.,0.,0.);
    TLorentzVector TSDiffTESUp(0.,0.,0.,0.);
    TLorentzVector TSDiffTESDown(0.,0.,0.,0.);
    TLorentzVector TSDiffL1TESUp(0.,0.,0.,0.);
    TLorentzVector TSDiffL1TESDown(0.,0.,0.,0.);
    TLorentzVector TSDiffL2TESUp(0.,0.,0.,0.);
    TLorentzVector TSDiffL2TESDown(0.,0.,0.,0.);

	TLorentzVector mvaMetVec(0.,0.,0.,0.);	
	TLorentzVector mvaMetVec_uncorr(0.,0.,0.,0.);	
	TLorentzVector mvaMetVec_responseUP(0.,0.,0.,0.);	
	TLorentzVector mvaMetVec_responseDOWN(0.,0.,0.,0.);	
	TLorentzVector mvaMetVec_resolutionUP(0.,0.,0.,0.);	
	TLorentzVector mvaMetVec_resolutionDOWN(0.,0.,0.,0.);	
	TLorentzVector pfMetVec(0.,0.,0.,0.);
    TLorentzVector pfMetVecUESUp(0.,0.,0.,0.);
    TLorentzVector pfMetVecUESDown(0.,0.,0.,0.);
    TLorentzVector pfMetVecJEnUp(0.,0.,0.,0.);
    TLorentzVector pfMetVecJEnDown(0.,0.,0.,0.);
    TLorentzVector pfMetVecTESUp(0.,0.,0.,0.);
    TLorentzVector pfMetVecTESDown(0.,0.,0.,0.);
    TLorentzVector pfMetVecL1TESUp(0.,0.,0.,0.);
    TLorentzVector pfMetVecL1TESDown(0.,0.,0.,0.);
    TLorentzVector pfMetVecL2TESUp(0.,0.,0.,0.);
    TLorentzVector pfMetVecL2TESDown(0.,0.,0.,0.);
	TLorentzVector puppiMetVec(0.,0.,0.,0.);
    
    /* set non corrected leg1 4-vector */
	l1NC.SetPtEtaPhiM(R.getD("leg1_pt"),R.getD("leg1_eta"),R.getD("leg1_phi"),R.getD("leg1_M"));

	/* set non corrected leg2 4-vector */
	l2NC.SetPtEtaPhiM(R.getD("leg2_pt"),R.getD("leg2_eta"),R.getD("leg2_phi"),R.getD("leg2_M"));
    
    std::vector <double> dummyShift;
    
    dummyShift = getTauShift(R.getI("leg1_decayMode"), R.getI("leg1_leptonType"), 0);
    double pxS1 = l1NC.Px()*dummyShift[0];
    double pyS1 = l1NC.Py()*dummyShift[0];
    double pzS1 = l1NC.Pz()*dummyShift[0];
    double massS1 = l1NC.M()*dummyShift[1];
    double enS1 = TMath::Sqrt(pxS1*pxS1 + pyS1*pyS1 + pzS1*pzS1 + massS1*massS1);
    
    dummyShift = getTauShift(R.getI("leg2_decayMode"), R.getI("leg2_leptonType"), 0);
    double pxS2 = l2NC.Px()*dummyShift[0];
    double pyS2 = l2NC.Py()*dummyShift[0];
    double pzS2 = l2NC.Pz()*dummyShift[0];
    double massS2 = l2NC.M()*dummyShift[1];
    double enS2 = TMath::Sqrt(pxS2*pxS2 + pyS2*pyS2 + pzS2*pzS2 + massS2*massS2);

	/* set leg1 4-vector */
	l1.SetPxPyPzE(pxS1,pyS1,pzS1,enS1);

	/* set the leg2 4-vector */
	l2.SetPxPyPzE(pxS2,pyS2,pzS2,enS2);
    
    /* Set original - corrected vectors for MET calculations */
    TSDiff = l1NC + l2NC - l1 - l2;
    
    dummyShift = getTauShift(R.getI("leg1_decayMode"), R.getI("leg1_leptonType"), 1);
    double pxS1Up = l1NC.Px()*dummyShift[0];
    double pyS1Up = l1NC.Py()*dummyShift[0];
    double pzS1Up = l1NC.Pz()*dummyShift[0];
    double massS1Up = l1NC.M()*dummyShift[1];
    double enS1Up = TMath::Sqrt(pxS1Up*pxS1Up + pyS1Up*pyS1Up + pzS1Up*pzS1Up + massS1Up*massS1Up);
    
    dummyShift = getTauShift(R.getI("leg2_decayMode"), R.getI("leg2_leptonType"), 1);
    double pxS2Up = l2NC.Px()*dummyShift[0];
    double pyS2Up = l2NC.Py()*dummyShift[0];
    double pzS2Up = l2NC.Pz()*dummyShift[0];
    double massS2Up = l2NC.M()*dummyShift[1];
    double enS2Up = TMath::Sqrt(pxS2Up*pxS2Up + pyS2Up*pyS2Up + pzS2Up*pzS2Up + massS2Up*massS2Up);
    
    dummyShift = getTauShift(R.getI("leg1_decayMode"), R.getI("leg1_leptonType"), -1);
    double pxS1Down = l1NC.Px()*dummyShift[0];
    double pyS1Down = l1NC.Py()*dummyShift[0];
    double pzS1Down = l1NC.Pz()*dummyShift[0];
    double massS1Down = l1NC.M()*dummyShift[1];
    double enS1Down = TMath::Sqrt(pxS1Down*pxS1Down + pyS1Down*pyS1Down + pzS1Down*pzS1Down + massS1Down*massS1Down);
    
    dummyShift = getTauShift(R.getI("leg2_decayMode"), R.getI("leg2_leptonType"), -1);
    double pxS2Down = l2NC.Px()*dummyShift[0];
    double pyS2Down = l2NC.Py()*dummyShift[0];
    double pzS2Down = l2NC.Pz()*dummyShift[0];
    double massS2Down = l2NC.M()*dummyShift[1];
    double enS2Down = TMath::Sqrt(pxS2Down*pxS2Down + pyS2Down*pyS2Down + pzS2Down*pzS2Down + massS2Down*massS2Down);
    
    /* set leg1 4-vector Variants*/
	l1TESUp.SetPxPyPzE(pxS1Up,pyS1Up,pzS1Up,enS1Up);
    l1TESDown.SetPxPyPzE(pxS1Down,pyS1Down,pzS1Down,enS1Down);
    
	/* set leg1 4-vector Variants*/
	l2TESUp.SetPxPyPzE(pxS2Up,pyS2Up,pzS2Up,enS2Up);
    l2TESDown.SetPxPyPzE(pxS2Down,pyS2Down,pzS2Down,enS2Down);
    
    dummyShift = getTauShift(R.getI("leg1_decayMode"), R.getI("leg1_leptonType"), 1);
    double pxS1UpNC = l1NC.Px()*dummyShift[2];
    double pyS1UpNC = l1NC.Py()*dummyShift[2];
    double pzS1UpNC = l1NC.Pz()*dummyShift[2];
    double massS1UpNC = l1NC.M()*dummyShift[3];
    double enS1UpNC = TMath::Sqrt(pxS1UpNC*pxS1UpNC + pyS1UpNC*pyS1UpNC + pzS1UpNC*pzS1UpNC + massS1UpNC*massS1UpNC);
    
    dummyShift = getTauShift(R.getI("leg2_decayMode"), R.getI("leg2_leptonType"), 1);
    double pxS2UpNC = l2NC.Px()*dummyShift[2];
    double pyS2UpNC = l2NC.Py()*dummyShift[2];
    double pzS2UpNC = l2NC.Pz()*dummyShift[2];
    double massS2UpNC = l2NC.M()*dummyShift[3];
    double enS2UpNC = TMath::Sqrt(pxS2UpNC*pxS2UpNC + pyS2UpNC*pyS2UpNC + pzS2UpNC*pzS2UpNC + massS2UpNC*massS2UpNC);
    
    dummyShift = getTauShift(R.getI("leg1_decayMode"), R.getI("leg1_leptonType"), -1);
    double pxS1DownNC = l1NC.Px()*dummyShift[2];
    double pyS1DownNC = l1NC.Py()*dummyShift[2];
    double pzS1DownNC = l1NC.Pz()*dummyShift[2];
    double massS1DownNC = l1NC.M()*dummyShift[3];
    double enS1DownNC = TMath::Sqrt(pxS1DownNC*pxS1DownNC + pyS1DownNC*pyS1DownNC + pzS1DownNC*pzS1DownNC + massS1DownNC*massS1DownNC);
    
    dummyShift = getTauShift(R.getI("leg2_decayMode"), R.getI("leg2_leptonType"), -1);
    double pxS2DownNC = l2NC.Px()*dummyShift[2];
    double pyS2DownNC = l2NC.Py()*dummyShift[2];
    double pzS2DownNC = l2NC.Pz()*dummyShift[2];
    double massS2DownNC = l2NC.M()*dummyShift[3];
    double enS2DownNC = TMath::Sqrt(pxS2DownNC*pxS2DownNC + pyS2DownNC*pyS2DownNC + pzS2DownNC*pzS2DownNC + massS2DownNC*massS2DownNC);
    
    /* set leg1 4-vector Variants not corrected */
	l1TESUpNC.SetPxPyPzE(pxS1UpNC,pyS1UpNC,pzS1UpNC,enS1UpNC);
    l1TESDownNC.SetPxPyPzE(pxS1DownNC,pyS1DownNC,pzS1DownNC,enS1DownNC);
    
	/* set leg1 4-vector Variants not corrected */
	l2TESUpNC.SetPxPyPzE(pxS2UpNC,pyS2UpNC,pzS2UpNC,enS2UpNC);
    l2TESDownNC.SetPxPyPzE(pxS2DownNC,pyS2DownNC,pzS2DownNC,enS2DownNC);
    
    TSDiffTESUp = l1TESUpNC + l2TESUpNC - l1TESUp - l2TESUp;
    TSDiffTESDown = l1TESDownNC + l2TESDownNC - l1TESDown - l2TESDown;
    
    TSDiffL1TESUp = l1TESUpNC + l2NC - l1TESUp - l2;
    TSDiffL1TESDown = l1TESDownNC + l2NC - l1TESDown - l2;
    
    TSDiffL2TESUp = l1NC + l2TESUpNC - l1 - l2TESUp;
    TSDiffL2TESDown = l1NC + l2TESDownNC - l1 - l2TESDown;
    
	/* set the MET 4-vectors */

    //MVA MET and Variants
	mvaMetVec.SetPtEtaPhiM(R.getD("corr_mvaMET"),0.0,R.getD("corr_mvaMETphi"),0.0); mvaMetVec += TSDiff;
	mvaMetVec_uncorr.SetPtEtaPhiM(R.getD("uncorr_mvaMET"),0.0,R.getD("uncorr_mvaMETphi"),0.0); mvaMetVec_uncorr += TSDiff;

	/* some variants of met are only stored for nominal Energy scale */

   	if(eventHasNominalLeptonEnergyScales)
   	{
		mvaMetVec_responseUP.SetPtEtaPhiM(R.getD("responseUP_mvaMET"),0.0,R.getD("responseUP_mvaMETphi"),0.0); mvaMetVec_responseUP += TSDiff;
		mvaMetVec_responseDOWN.SetPtEtaPhiM(R.getD("responseDOWN_mvaMET"),0.0,R.getD("responseDOWN_mvaMETphi"),0.0); mvaMetVec_responseDOWN += TSDiff;
		mvaMetVec_resolutionUP.SetPtEtaPhiM(R.getD("resolutionUP_mvaMET"),0.0,R.getD("resolutionUP_mvaMETphi"),0.0); mvaMetVec_resolutionUP += TSDiff;
		mvaMetVec_resolutionDOWN.SetPtEtaPhiM(R.getD("resolutionDOWN_mvaMET"),0.0,R.getD("resolutionDOWN_mvaMETphi"),0.0); mvaMetVec_resolutionDOWN += TSDiff;
   	}

    //PFMET TYPE1 and Variants
	pfMetVec.SetPtEtaPhiM(R.getD("pfmet_type1_Pt"),0.0,R.getD("pfmet_type1_Phi"),0.0); pfMetVec += TSDiff;
    pfMetVecUESUp.SetPtEtaPhiM(R.getD("pfmet_type1_UnclusteredEnUp_Pt"),0.0,R.getD("pfmet_type1_UnclusteredEnUp_Phi"),0.0); pfMetVecUESUp += TSDiff;
    pfMetVecUESDown.SetPtEtaPhiM(R.getD("pfmet_type1_UnclusteredEnDown_Pt"),0.0,R.getD("pfmet_type1_UnclusteredEnDown_Phi"),0.0); pfMetVecUESDown += TSDiff;
    pfMetVecJEnUp.SetPtEtaPhiM(R.getD("pfmet_type1_JetEnUp_Pt"),0.0,R.getD("pfmet_type1_JetEnUp_Phi"),0.0); pfMetVecJEnUp += TSDiff;
    pfMetVecJEnDown.SetPtEtaPhiM(R.getD("pfmet_type1_JetEnDown_Pt"),0.0,R.getD("pfmet_type1_JetEnDown_Phi"),0.0); pfMetVecJEnDown += TSDiff;
    
    //correct for TES shifts, only for PF TYPE1 now
    pfMetVecTESUp.SetPtEtaPhiM(R.getD("pfmet_type1_UnclusteredEnUp_Pt"),0.0,R.getD("pfmet_type1_UnclusteredEnUp_Phi"),0.0); pfMetVecTESUp += TSDiffTESUp;
    pfMetVecTESDown.SetPtEtaPhiM(R.getD("pfmet_type1_UnclusteredEnDown_Pt"),0.0,R.getD("pfmet_type1_UnclusteredEnDown_Phi"),0.0); pfMetVecTESDown += TSDiffTESDown;
    
    pfMetVecL1TESUp.SetPtEtaPhiM(R.getD("pfmet_type1_UnclusteredEnUp_Pt"),0.0,R.getD("pfmet_type1_UnclusteredEnUp_Phi"),0.0); pfMetVecL1TESUp += TSDiffL1TESUp;
    pfMetVecL1TESDown.SetPtEtaPhiM(R.getD("pfmet_type1_UnclusteredEnDown_Pt"),0.0,R.getD("pfmet_type1_UnclusteredEnDown_Phi"),0.0); pfMetVecL1TESDown += TSDiffL1TESDown;
    
    pfMetVecL2TESUp.SetPtEtaPhiM(R.getD("pfmet_type1_UnclusteredEnUp_Pt"),0.0,R.getD("pfmet_type1_UnclusteredEnUp_Phi"),0.0); pfMetVecL2TESUp += TSDiffL2TESUp;
    pfMetVecL2TESDown.SetPtEtaPhiM(R.getD("pfmet_type1_UnclusteredEnDown_Pt"),0.0,R.getD("pfmet_type1_UnclusteredEnDown_Phi"),0.0); pfMetVecL2TESDown += TSDiffL2TESDown;
    
    //PUPPI
	puppiMetVec.SetPtEtaPhiM(R.getD("puppiMET"),0.0,R.getD("puppiMETphi"),0.0); puppiMetVec += TSDiff;

    /* The LHE weights & Scale Factor Vector, mapping is available in FlatTuple production log file */
    
    originalXWGTUP = R.getF("originalXWGTUP");
    
    if( R.getS("DataCard") == "SIGNAL_MONO_HIGGS" && eventHasNominalLeptonEnergyScales)
    {
        theory_scale_factors = R.getVF("theory_scale_factors");
    }

    /* basic info */

	pairRank = R.getUI("pairRank");
 	isOsPair = R.getI("isOsPair");
    isBoostedChannelPair = R.getB("isBoostedChannelPair");
    
    if (R.getS("DataCard") == "DY") DataCardInt = 1;
    if (R.getS("DataCard") == "TT") DataCardInt = 2;
    if (R.getS("DataCard") == "VV") DataCardInt = 3;
    if (R.getS("DataCard") == "VVV") DataCardInt = 4;
    if (R.getS("DataCard") == "EWK") DataCardInt = 5;
    if (R.getS("DataCard") == "W") DataCardInt = 6;
    if (R.getS("DataCard") == "SMHIGGS") DataCardInt = 7;
    if (R.getS("DataCard") == "DYINV") DataCardInt = 8;
    if (R.getS("DataCard") == "MonoHiggs_2HDM") DataCardInt = 9;
    if (R.getS("DataCard") == "MonoHiggs_ZpBaryonic") DataCardInt = 10;

 	/* gen boson/top info */

	genBosonTotal_pt = R.getD("genBosonTotal_pt"); 
    genBosonTotal_eta = R.getD("genBosonTotal_eta");
    genBosonTotal_phi = R.getD("genBosonTotal_phi");
    genBosonTotal_M = R.getD("genBosonTotal_M");

	genBosonVisible_pt = R.getD("genBosonVisible_pt");
	genBosonVisible_eta = R.getD("genBosonVisible_eta");
	genBosonVisible_phi = R.getD("genBosonVisible_phi");
	genBosonVisible_M = R.getD("genBosonVisible_M");
    
    leg1_GENMOTHERpdgId = R.getI("leg1_GENMOTHERpdgId");
    leg2_GENMOTHERpdgId = R.getI("leg2_GENMOTHERpdgId");
    
    if(R.getS("KeyName") == "WJetsToLNu" ||\
	   R.getS("KeyName") == "WJetsToLNuext2-v1" ||\
	   R.getS("KeyName") == "WJetsToLNu_HT-70To100" ||\
	   R.getS("KeyName") == "WJetsToLNu_HT-100To200orig" ||\
	   R.getS("KeyName") == "WJetsToLNu_HT-100To200" ||\
       R.getS("KeyName") == "WJetsToLNu_HT-100To200ext2-v1" ||\
       R.getS("KeyName") == "WJetsToLNu_HT-200To400orig" ||\
       R.getS("KeyName") == "WJetsToLNu_HT-200To400ext1-v1" ||\
       R.getS("KeyName") == "WJetsToLNu_HT-200To400" ||\
       R.getS("KeyName") == "WJetsToLNu_HT-400To600orig" ||\
       R.getS("KeyName") == "WJetsToLNu_HT-400To600" ||\
       R.getS("KeyName") == "WJetsToLNu_HT-600To800orig" ||\
       R.getS("KeyName") == "WJetsToLNu_HT-600To800" ||\
       R.getS("KeyName") == "WJetsToLNu_HT-800To1200orig" ||\
       R.getS("KeyName") == "WJetsToLNu_HT-800To1200" ||\
       R.getS("KeyName") == "WJetsToLNu_HT-1200To2500orig" ||\
       R.getS("KeyName") == "WJetsToLNu_HT-1200To2500" ||\
       R.getS("KeyName") == "WJetsToLNu_HT-2500ToInforig" ||\
       R.getS("KeyName") == "WJetsToLNu_HT-2500ToInf" ||\
       R.getS("KeyName") == "ZZTo2L2Nu" ||\
       R.getS("KeyName") == "WWTo2L2Nu")
    {
        genBosonTotal_Wpt = R.getD("MaxPtGenBoson_WisconinStyle_pt");
    }
    
	genTopPt1 = R.getD("genTopPt1");
	genTopPt2 = R.getD("genTopPt2");

	/* event ID */

	run = R.getUI("run");
	event = R.getUI("event");
	evt = R.getUI("event");
	lumi = R.getUI("luminosityBlock");

	/* vertices & pileup */

	npv = R.getI("NumberOfGoodVertices");
 	npu = R.getF("NumTruePileUpInt");
 	rho = R.getD("rho");
 	puweight = R.getD("puWeight");

	/* leg 1 quantities */

 	pt_1 				= l1.Pt();
	phi_1 				= l1.Phi();
	eta_1 				= l1.Eta();
	m_1 				= l1.M();
	iso_1 				= R.getF("leg1_RelIso");
	dZ_1 				= R.getF("leg1_dz");
	d0_1 				= R.getF("leg1_dxy");
	q_1 				= R.getI("leg1_charge");
	id_e_mva_nt_loose_1 = R.getF("leg1_raw_electronMVA");
	tau_decay_mode_1 	= R.getI("leg1_decayMode");
	ZimpactTau_1 		= R.getF("leg1_ZimpactTau");
	dzTauVertex_1 		= R.getF("leg1_dzTauVertex");
	mt_1  				= GetTransverseMass(pfMetVec,l1);
	pfmt_1 				= GetTransverseMass(pfMetVec,l1);
    pfmt_1_UESUp		= GetTransverseMass(pfMetVecUESUp,l1);
    pfmt_1_UESDown 		= GetTransverseMass(pfMetVecUESDown,l1);
    pfmt_1_JEnUp		= GetTransverseMass(pfMetVecJEnUp,l1);
    pfmt_1_JEnDown 		= GetTransverseMass(pfMetVecJEnDown,l1);
	puppimt_1 			= GetTransverseMass(puppiMetVec,l1);
	mt_uncorr_1   		= GetTransverseMass(mvaMetVec_uncorr,l1);
    pt_1_flat 			= l1NC.Pt();
	phi_1_flat 			= l1NC.Phi();
	eta_1_flat 			= l1NC.Eta();
	m_1_flat 			= l1NC.M();
    
    //TES
    pt_1_TESUp          = l1TESUp.Pt();
    mt_1_TESUp  		= GetTransverseMass(pfMetVecTESUp,l1TESUp);
    pfmt_1_TESUp        = GetTransverseMass(pfMetVecTESUp,l1TESUp);
    mt_1_L1TESUp  		= GetTransverseMass(pfMetVecL1TESUp,l1TESUp);
    pfmt_1_L1TESUp      = GetTransverseMass(pfMetVecL1TESUp,l1TESUp);
    mt_1_L2TESUp  		= GetTransverseMass(pfMetVecL2TESUp,l1);
    pfmt_1_L2TESUp      = GetTransverseMass(pfMetVecL2TESUp,l1);
    
    pt_1_TESDown        = l1TESDown.Pt();
    mt_1_TESDown		= GetTransverseMass(pfMetVecTESDown,l1TESDown);
    pfmt_1_TESDown      = GetTransverseMass(pfMetVecTESDown,l1TESDown);
    mt_1_L1TESDown      = GetTransverseMass(pfMetVecL1TESDown,l1TESDown);
    pfmt_1_L1TESDown    = GetTransverseMass(pfMetVecL1TESDown,l1TESDown);
    mt_1_L2TESDown  	= GetTransverseMass(pfMetVecL2TESDown,l1);
    pfmt_1_L2TESDown    = GetTransverseMass(pfMetVecL2TESDown,l1);
    
	if(eventHasNominalLeptonEnergyScales)
	{
		responseUP_MTmvaMET_1		=  GetTransverseMass(mvaMetVec_responseUP,l1);
		responseDOWN_MTmvaMET_1		=  GetTransverseMass(mvaMetVec_responseDOWN,l1);
		resolutionUP_MTmvaMET_1		=  GetTransverseMass(mvaMetVec_resolutionUP,l1);
		resolutionDOWN_MTmvaMET_1	=  GetTransverseMass(mvaMetVec_resolutionDOWN,l1);
	}

	gen_match_1 		= R.getI("leg1_MCMatchType");
	genMCmatch_pt_1		= R.getD("leg1_genMCmatch_pt");
	genMCmatch_eta_1	= R.getD("leg1_genMCmatch_eta");
	genMCmatch_phi_1	= R.getD("leg1_genMCmatch_phi");
	genMCmatch_M_1		= R.getD("leg1_genMCmatch_M");
	MCMatchPdgId_1		= R.getI("leg1_MCMatchPdgId");

	byIsolationMVArun2v1DBdR03oldDMwLTraw_1 = R.getF("leg1_byIsolationMVArun2v1DBdR03oldDMwLTraw");
	byTightIsolationMVArun2v1DBdR03oldDMwLT_1 = R.getF("leg1_byTightIsolationMVArun2v1DBdR03oldDMwLT");
	byVTightIsolationMVArun2v1DBdR03oldDMwLT_1 = R.getF("leg1_byVTightIsolationMVArun2v1DBdR03oldDMwLT");
	byLooseIsolationMVArun2v1DBdR03oldDMwLT_1 = R.getF("leg1_byLooseIsolationMVArun2v1DBdR03oldDMwLT");
	byMediumIsolationMVArun2v1DBdR03oldDMwLT_1 = R.getF("leg1_byMediumIsolationMVArun2v1DBdR03oldDMwLT");
	byVLooseIsolationMVArun2v1DBdR03oldDMwLT_1 = R.getF("leg1_byVLooseIsolationMVArun2v1DBdR03oldDMwLT");
	byVVTightIsolationMVArun2v1DBdR03oldDMwLT_1 = R.getF("leg1_byVVTightIsolationMVArun2v1DBdR03oldDMwLT");
    byLooseIsolationMVArun2v1DBoldDMwLT_1 = R.getF("leg1_byLooseIsolationMVArun2v1DBoldDMwLT");
    byTightIsolationMVArun2v1DBoldDMwLT_1 = R.getF("leg1_byTightIsolationMVArun2v1DBoldDMwLT");
	againstElectronVLooseMVA6_1 = R.getF("leg1_againstElectronVLooseMVA6");
	againstMuonTight3_1 = R.getF("leg1_againstMuonTight3");
	againstElectronTightMVA6_1 = R.getF("leg1_againstElectronTightMVA6");
	againstMuonLoose3_1 = R.getF("leg1_againstMuonLoose3");
	decayModeFinding_1 = R.getF("leg1_decayModeFinding");
	//byIsolationMVA3oldDMwLTraw_1 = R.getF("leg1_byIsolationMVA3oldDMwLTraw");
	byCombinedIsolationDeltaBetaCorrRaw3Hits_1 = R.getF("leg1_byCombinedIsolationDeltaBetaCorrRaw3Hits");
	byIsolationMVArun2v1DBnewDMwLTraw_1 = R.getF("leg1_byIsolationMVArun2v1DBnewDMwLTraw");
	decayModeFindingNewDMs_1 = R.getF("leg1_decayModeFindingNewDMs");

	/* leg 2 quantities */

 	pt_2 				= l2.Pt();
	phi_2 				= l2.Phi();
	eta_2 				= l2.Eta();
	m_2 				= l2.M();
	iso_2 				= R.getF("leg2_RelIso");
	dZ_2 				= R.getF("leg2_dz");
	d0_2 				= R.getF("leg2_dxy");
	q_2 				= R.getI("leg2_charge");
	id_e_mva_nt_loose_2 = R.getF("leg2_raw_electronMVA");
	tau_decay_mode_2 	= R.getI("leg2_decayMode");
	ZimpactTau_2 		= R.getF("leg2_ZimpactTau");
	dzTauVertex_2 		= R.getF("leg2_dzTauVertex");
	mt_2  				= GetTransverseMass(pfMetVec,l2);
	pfmt_2 				= GetTransverseMass(pfMetVec,l2);
	puppimt_2 			= GetTransverseMass(puppiMetVec,l2);
	mt_uncorr_2   		= GetTransverseMass(mvaMetVec_uncorr,l2);
    pt_2_flat 			= l2NC.Pt();
	phi_2_flat 			= l2NC.Phi();
	eta_2_flat 			= l2NC.Eta();
	m_2_flat 			= l2NC.M();
    
    //TES
    pt_2_TESUp          = l2TESUp.Pt();
    mt_2_TESUp  		= GetTransverseMass(pfMetVecTESUp,l2TESUp);
    pfmt_2_TESUp        = GetTransverseMass(pfMetVecTESUp,l2TESUp);
    mt_2_L1TESUp  		= GetTransverseMass(pfMetVecL1TESUp,l2);
    pfmt_2_L1TESUp      = GetTransverseMass(pfMetVecL1TESUp,l2);
    mt_2_L2TESUp  		= GetTransverseMass(pfMetVecL2TESUp,l2TESUp);
    pfmt_2_L2TESUp      = GetTransverseMass(pfMetVecL2TESUp,l2TESUp);
    
    pt_2_TESDown        = l2TESDown.Pt();
    mt_2_TESDown  		= GetTransverseMass(pfMetVecTESDown,l2TESDown);
    pfmt_2_TESDown      = GetTransverseMass(pfMetVecTESDown,l2TESDown);
    mt_2_L1TESDown  	= GetTransverseMass(pfMetVecL1TESDown,l2);
    pfmt_2_L1TESDown    = GetTransverseMass(pfMetVecL1TESDown,l2);
    mt_2_L2TESDown  	= GetTransverseMass(pfMetVecL2TESDown,l2TESDown);
    pfmt_2_L2TESDown    = GetTransverseMass(pfMetVecL2TESDown,l2TESDown);

	if(eventHasNominalLeptonEnergyScales)
	{
		responseUP_MTmvaMET_2		=  GetTransverseMass(mvaMetVec_responseUP,l2);
		responseDOWN_MTmvaMET_2		=  GetTransverseMass(mvaMetVec_responseDOWN,l2);
		resolutionUP_MTmvaMET_2		=  GetTransverseMass(mvaMetVec_resolutionUP,l2);
		resolutionDOWN_MTmvaMET_2	=  GetTransverseMass(mvaMetVec_resolutionDOWN,l2);
	}

	gen_match_2 		= R.getI("leg2_MCMatchType");
	genMCmatch_pt_2		= R.getD("leg2_genMCmatch_pt");
	genMCmatch_eta_2	= R.getD("leg2_genMCmatch_eta");
	genMCmatch_phi_2	= R.getD("leg2_genMCmatch_phi");
	genMCmatch_M_2		= R.getD("leg2_genMCmatch_M");
	MCMatchPdgId_2		= R.getI("leg2_MCMatchPdgId");

	byIsolationMVArun2v1DBdR03oldDMwLTraw_2 = R.getF("leg2_byIsolationMVArun2v1DBdR03oldDMwLTraw");
	byTightIsolationMVArun2v1DBdR03oldDMwLT_2 = R.getF("leg2_byTightIsolationMVArun2v1DBdR03oldDMwLT");
	byVTightIsolationMVArun2v1DBdR03oldDMwLT_2 = R.getF("leg2_byVTightIsolationMVArun2v1DBdR03oldDMwLT");
	byLooseIsolationMVArun2v1DBdR03oldDMwLT_2 = R.getF("leg2_byLooseIsolationMVArun2v1DBdR03oldDMwLT");
	byMediumIsolationMVArun2v1DBdR03oldDMwLT_2 = R.getF("leg2_byMediumIsolationMVArun2v1DBdR03oldDMwLT");
	byVLooseIsolationMVArun2v1DBdR03oldDMwLT_2 = R.getF("leg2_byVLooseIsolationMVArun2v1DBdR03oldDMwLT");
	byVVTightIsolationMVArun2v1DBdR03oldDMwLT_2 = R.getF("leg2_byVVTightIsolationMVArun2v1DBdR03oldDMwLT");
    byLooseIsolationMVArun2v1DBoldDMwLT_2 = R.getF("leg2_byLooseIsolationMVArun2v1DBoldDMwLT");
    byTightIsolationMVArun2v1DBoldDMwLT_2 = R.getF("leg2_byTightIsolationMVArun2v1DBoldDMwLT");
	againstElectronVLooseMVA6_2 = R.getF("leg2_againstElectronVLooseMVA6");
	againstMuonTight3_2 = R.getF("leg2_againstMuonTight3");
	againstElectronTightMVA6_2 = R.getF("leg2_againstElectronTightMVA6");
	againstMuonLoose3_2 = R.getF("leg2_againstMuonLoose3");
	decayModeFinding_2 = R.getF("leg2_decayModeFinding");
	//byIsolationMVA3oldDMwLTraw_2 = R.getF("leg2_byIsolationMVA3oldDMwLTraw");
	byCombinedIsolationDeltaBetaCorrRaw3Hits_2 = R.getF("leg2_byCombinedIsolationDeltaBetaCorrRaw3Hits");
	byIsolationMVArun2v1DBnewDMwLTraw_2 = R.getF("leg2_byIsolationMVArun2v1DBnewDMwLTraw");
	decayModeFindingNewDMs_2 = R.getF("leg2_decayModeFindingNewDMs");

	/* di-tau system */

	pt_tt = (l1+l2).Pt();
    pt_tt_TESUp = (l1TESUp+l2TESUp).Pt();
    pt_tt_TESDown = (l1TESDown+l2TESDown).Pt();
    pt_tt_L1TESUp = (l1TESUp+l2).Pt();
    pt_tt_L1TESDown = (l1TESDown+l2).Pt();
    pt_tt_L2TESUp = (l1+l2TESUp).Pt();
    pt_tt_L2TESDown = (l1+l2TESDown).Pt();
	mt_tot = mtTotCalc(l1, l2, pfMetVec);
    mt_tot_UESUp = mtTotCalc(l1, l2, pfMetVecUESUp);
    mt_tot_UESDown = mtTotCalc(l1, l2, pfMetVecUESDown);
    mt_tot_TESUp = mtTotCalc(l1TESUp, l2TESUp, pfMetVecTESUp);
    mt_tot_TESDown = mtTotCalc(l1TESDown, l2TESDown, pfMetVecTESDown);
    mt_tot_L1TESUp = mtTotCalc(l1TESUp, l2, pfMetVecL1TESUp);
    mt_tot_L1TESDown = mtTotCalc(l1TESDown, l2, pfMetVecL1TESDown);
    mt_tot_L2TESUp = mtTotCalc(l1, l2TESUp, pfMetVecL2TESUp);
    mt_tot_L2TESDown = mtTotCalc(l1, l2TESDown, pfMetVecL2TESDown);
	m_vis = (l1+l2).M();
    m_vis_flat = R.getD("VISMass");
    m_vis_TESUp = (l1TESUp + l2TESUp).M();
    m_vis_TESDown = (l1TESDown + l2TESDown).M();
    m_vis_L1TESUp = (l1TESUp + l2).M();
    m_vis_L1TESDown = (l1TESDown + l2).M();
    m_vis_L2TESUp = (l1 + l2TESUp).M();
    m_vis_L2TESDown = (l1 + l2TESDown).M();
	DeltaR_leg1_leg2 = R.getD("DeltaR_leg1_leg2");
    cos_DeltaPhi_leg1_leg2 = cos(l1.DeltaPhi(l2));
    cos_DeltaPhi_PFMET_Higgs = cos((pfMetVec).DeltaPhi(l1+l2));
    cos_DeltaPhi_PFMET_Higgs_UESUp = cos((pfMetVecUESUp).DeltaPhi(l1+l2));
    cos_DeltaPhi_PFMET_Higgs_UESDown = cos((pfMetVecUESDown).DeltaPhi(l1+l2));
    cos_DeltaPhi_PFMET_Higgs_JEnUp = cos((pfMetVecJEnUp).DeltaPhi(l1+l2));
    cos_DeltaPhi_PFMET_Higgs_JEnDown = cos((pfMetVecJEnDown).DeltaPhi(l1+l2));
    cos_DeltaPhi_PFMET_Higgs_TESUp = cos((pfMetVecTESUp).DeltaPhi(l1TESUp+l2TESUp));
    cos_DeltaPhi_PFMET_Higgs_TESDown = cos((pfMetVecTESDown).DeltaPhi(l1TESDown+l2TESDown));
    cos_DeltaPhi_PFMET_Higgs_L1TESUp = cos((pfMetVecL1TESUp).DeltaPhi(l1TESUp+l2));
    cos_DeltaPhi_PFMET_Higgs_L1TESDown = cos((pfMetVecL1TESDown).DeltaPhi(l1TESDown+l2));
    cos_DeltaPhi_PFMET_Higgs_L2TESUp = cos((pfMetVecL2TESUp).DeltaPhi(l1+l2TESUp));
    cos_DeltaPhi_PFMET_Higgs_L2TESDown = cos((pfMetVecL2TESDown).DeltaPhi(l1+l2TESDown));
    
	/* sv fit -- only keeping variants which use mvaMET, study by H2Tau group showed met variants have minimal impact on SVMass shape
	and so those are omitted */
    /*
	m_sv 						= R.getD("SVMass");
	mt_sv 						= R.getD("SVTransverseMass");
	SVFit_mvaMET_diTau_pt		= R.getD("SVFit_mvaMET_diTau_pt");
	SVFit_mvaMET_diTau_eta		= R.getD("SVFit_mvaMET_diTau_eta");
	SVFit_mvaMET_diTau_phi		= R.getD("SVFit_mvaMET_diTau_phi");
	SVFit_mvaMET_FittedMET		= R.getD("SVFit_mvaMET_FittedMET");
	SVFit_mvaMET_FittedMETphi	= R.getD("SVFit_mvaMET_FittedMETphi");
    */

	/* met-related branches */
	mvamet				= mvaMetVec.Pt();
	mvametphi			= mvaMetVec.Phi();
	met					= pfMetVec.Pt();
	metphi				= pfMetVec.Phi();
	puppimet			= puppiMetVec.Pt();
	puppimetphi			= puppiMetVec.Phi();
	uncorr_mvamet		= mvaMetVec_uncorr.Pt();
	uncorr_mvametphi	= mvaMetVec_uncorr.Phi();
    
    pfmet_type1_Pt = pfMetVec.Pt();
	pfmet_type1_Phi = pfMetVec.Phi();
	pfmet_type1_MT1 = GetTransverseMass(pfMetVec,l1);
	pfmet_type1_MT2 = GetTransverseMass(pfMetVec,l2);
    
    pfmet_type1_UnclusteredEnUp_Pt = pfMetVecUESUp.Pt();
	pfmet_type1_UnclusteredEnUp_Phi = pfMetVecUESUp.Phi();
    pfmet_type1_UnclusteredEnUp_MT1 = GetTransverseMass(pfMetVecUESUp,l1);
	pfmet_type1_UnclusteredEnUp_MT2 = GetTransverseMass(pfMetVecUESUp,l2);
	pfmet_type1_UnclusteredEnDown_Pt = pfMetVecUESDown.Pt();
	pfmet_type1_UnclusteredEnDown_Phi = pfMetVecUESDown.Phi();
	pfmet_type1_UnclusteredEnDown_MT1 = GetTransverseMass(pfMetVecUESDown,l1);
	pfmet_type1_UnclusteredEnDown_MT2 = GetTransverseMass(pfMetVecUESDown,l2);
    
    pfmet_type1_TESUp_Pt = pfMetVecTESUp.Pt();
	pfmet_type1_TESUp_Phi = pfMetVecTESUp.Phi();
	pfmet_type1_TESUp_MT1 = GetTransverseMass(pfMetVecTESUp,l1TESUp);
	pfmet_type1_TESUp_MT2 = GetTransverseMass(pfMetVecTESUp,l2TESUp);
	pfmet_type1_TESDown_Pt = pfMetVecTESDown.Pt();
	pfmet_type1_TESDown_Phi = pfMetVecTESDown.Phi();
	pfmet_type1_TESDown_MT1 = GetTransverseMass(pfMetVecTESDown,l1TESDown);
	pfmet_type1_TESDown_MT2 = GetTransverseMass(pfMetVecTESDown,l2TESDown);
    
    pfmet_type1_L1TESUp_Pt = pfMetVecL1TESUp.Pt();
	pfmet_type1_L1TESDown_Pt = pfMetVecL1TESDown.Pt();
    
    pfmet_type1_L2TESUp_Pt = pfMetVecL2TESUp.Pt();
	pfmet_type1_L2TESDown_Pt = pfMetVecL2TESDown.Pt();
    
    pfmet_type1_JetEnUp_Pt = pfMetVecJEnUp.Pt();
	pfmet_type1_JetEnUp_Phi = pfMetVecJEnUp.Phi();
	pfmet_type1_JetEnUp_MT1 = GetTransverseMass(pfMetVecJEnUp,l1);
	pfmet_type1_JetEnUp_MT2 = GetTransverseMass(pfMetVecJEnUp,l2);
	pfmet_type1_JetEnDown_Pt = pfMetVecJEnDown.Pt();
	pfmet_type1_JetEnDown_Phi = pfMetVecJEnDown.Phi();
	pfmet_type1_JetEnDown_MT1 = GetTransverseMass(pfMetVecJEnDown,l1);
	pfmet_type1_JetEnDown_MT2 = GetTransverseMass(pfMetVecJEnDown,l2);


   	if(eventHasNominalLeptonEnergyScales)
   	{
   		responseUP_mvaMET			= mvaMetVec_responseUP.Pt();
   		responseUP_mvaMETphi		= mvaMetVec_responseUP.Phi();
   		responseDOWN_mvaMET			= mvaMetVec_responseDOWN.Pt();
   		responseDOWN_mvaMETphi		= mvaMetVec_responseDOWN.Phi();
   		resolutionUP_mvaMET			= mvaMetVec_resolutionUP.Pt();
   		resolutionUP_mvaMETphi		= mvaMetVec_resolutionUP.Phi();
   		resolutionDOWN_mvaMET		= mvaMetVec_resolutionDOWN.Pt();
   		resolutionDOWN_mvaMETphi	= mvaMetVec_resolutionDOWN.Phi();
   	}

	/* sig matrix using mva met */

	mvacov00			= R.getD("mvaMET_cov00");     
	mvacov01			= R.getD("mvaMET_cov01");
	mvacov10			= R.getD("mvaMET_cov10");
	mvacov11			= R.getD("mvaMET_cov11");

	/* sig matrix using  met */

	if(eventIsNotSmallTree)
	{
		metcov00			= R.getD("pfMET_cov00");     
		metcov01			= R.getD("pfMET_cov01");
		metcov10			= R.getD("pfMET_cov10");
		metcov11			= R.getD("pfMET_cov11");
	}

    //These are not corrected by shift!
	genMET = R.getD("genMET");
	genMETphi = R.getD("genMETphi");
    genMETeta = R.getD("genMETeta");
	genMETmass = R.getD("genMETmass");
	pfmet_raw_Pt = R.getD("pfmet_raw_Pt");
	pfmet_raw_Phi = R.getD("pfmet_raw_Phi");
	pfmet_raw_MT1 = R.getD("pfmet_raw_MT1");
	pfmet_raw_MT2 = R.getD("pfmet_raw_MT2");
	pfmet_type1_JetResUp_Pt = R.getD("pfmet_type1_JetResUp_Pt");
	pfmet_type1_JetResUp_Phi = R.getD("pfmet_type1_JetResUp_Phi");
	pfmet_type1_JetResUp_MT1 = R.getD("pfmet_type1_JetResUp_MT1");
	pfmet_type1_JetResUp_MT2 = R.getD("pfmet_type1_JetResUp_MT2");
	pfmet_type1_JetResDown_Pt = R.getD("pfmet_type1_JetResDown_Pt");
	pfmet_type1_JetResDown_Phi = R.getD("pfmet_type1_JetResDown_Phi");
	pfmet_type1_JetResDown_MT1 = R.getD("pfmet_type1_JetResDown_MT1");
	pfmet_type1_JetResDown_MT2 = R.getD("pfmet_type1_JetResDown_MT2");
	pfmet_type1_MuonEnUp_Pt = R.getD("pfmet_type1_MuonEnUp_Pt");
	pfmet_type1_MuonEnUp_Phi = R.getD("pfmet_type1_MuonEnUp_Phi");
	pfmet_type1_MuonEnUp_MT1 = R.getD("pfmet_type1_MuonEnUp_MT1");
	pfmet_type1_MuonEnUp_MT2 = R.getD("pfmet_type1_MuonEnUp_MT2");
	pfmet_type1_MuonEnDown_Pt = R.getD("pfmet_type1_MuonEnDown_Pt");
	pfmet_type1_MuonEnDown_Phi = R.getD("pfmet_type1_MuonEnDown_Phi");
	pfmet_type1_MuonEnDown_MT1 = R.getD("pfmet_type1_MuonEnDown_MT1");
	pfmet_type1_MuonEnDown_MT2 = R.getD("pfmet_type1_MuonEnDown_MT2");
	pfmet_type1_ElectronEnUp_Pt = R.getD("pfmet_type1_ElectronEnUp_Pt");
	pfmet_type1_ElectronEnUp_Phi = R.getD("pfmet_type1_ElectronEnUp_Phi");
	pfmet_type1_ElectronEnUp_MT1 = R.getD("pfmet_type1_ElectronEnUp_MT1");
	pfmet_type1_ElectronEnUp_MT2 = R.getD("pfmet_type1_ElectronEnUp_MT2");
	pfmet_type1_ElectronEnDown_Pt = R.getD("pfmet_type1_ElectronEnDown_Pt");
	pfmet_type1_ElectronEnDown_Phi = R.getD("pfmet_type1_ElectronEnDown_Phi");
	pfmet_type1_ElectronEnDown_MT1 = R.getD("pfmet_type1_ElectronEnDown_MT1");
	pfmet_type1_ElectronEnDown_MT2 = R.getD("pfmet_type1_ElectronEnDown_MT2");
	pfmet_type1_PhotonEnUp_Pt = R.getD("pfmet_type1_PhotonEnUp_Pt");
	pfmet_type1_PhotonEnUp_Phi = R.getD("pfmet_type1_PhotonEnUp_Phi");
	pfmet_type1_PhotonEnUp_MT1 = R.getD("pfmet_type1_PhotonEnUp_MT1");
	pfmet_type1_PhotonEnUp_MT2 = R.getD("pfmet_type1_PhotonEnUp_MT2");
	pfmet_type1_PhotonEnDown_Pt = R.getD("pfmet_type1_PhotonEnDown_Pt");
	pfmet_type1_PhotonEnDown_Phi = R.getD("pfmet_type1_PhotonEnDown_Phi");
	pfmet_type1_PhotonEnDown_MT1 = R.getD("pfmet_type1_PhotonEnDown_MT1");
	pfmet_type1_PhotonEnDown_MT2 = R.getD("pfmet_type1_PhotonEnDown_MT2");

	/* p_zeta variables */
    
    mt_tot_JEnUp = mtTotCalc(l1, l2, pfMetVecJEnUp);
    mt_tot_JEnDown = mtTotCalc(l1, l2, pfMetVecJEnDown);
    
    weight_ttPtUp = ttTrigPtShape(1, l1, l2);
    weight_ttPtDown =  ttTrigPtShape(0, l1, l2);

	pzetavis 			= pzetaVisCalc(l1,l2);

	pzetamiss			= pzetaMissCalc(l1,l2,mvaMetVec);
	pfpzetamiss 		= pzetaMissCalc(l1,l2,pfMetVec);
	puppipzetamiss		= pzetaMissCalc(l1,l2,puppiMetVec);
	pzetamiss_uncorr	= pzetaMissCalc(l1,l2,mvaMetVec_uncorr);

   	if(eventHasNominalLeptonEnergyScales)
   	{   
		pzetamiss_responseUP 		= pzetaMissCalc(l1,l2,mvaMetVec_responseUP);
		pzetamiss_responseDOWN 		= pzetaMissCalc(l1,l2,mvaMetVec_responseDOWN);
		pzetamiss_resolutionUP		= pzetaMissCalc(l1,l2,mvaMetVec_resolutionUP);
		pzetamiss_resolutionDOWN	= pzetaMissCalc(l1,l2,mvaMetVec_resolutionDOWN);
   	}

   	/* fill the jet and b-tag info */

   	/* for nominal jets */
   	std::string argString = "";

   	jetINFOstruct.reset();
	fillJetBranches(eventHasNominalLeptonEnergyScales, eventIsNotSmallTree, argString, jetINFOstruct);

	njets = jetINFOstruct.m_njets;
	njetspt20 = jetINFOstruct.m_njetspt20;
	mjj = jetINFOstruct.m_mjj;
	jdeta = jetINFOstruct.m_jdeta;
	njetingap = jetINFOstruct.m_njetingap;
	njetingap20 = jetINFOstruct.m_njetingap20;
	jdphi = jetINFOstruct.m_jdphi;
	jpt_1 = jetINFOstruct.m_jpt_1;
	jeta_1 = jetINFOstruct.m_jeta_1;
	jphi_1 = jetINFOstruct.m_jphi_1;
	jm_1 = jetINFOstruct.m_jm_1;
	jmva_1 = jetINFOstruct.m_jmva_1;
	jpt_2 = jetINFOstruct.m_jpt_2;
	jeta_2 = jetINFOstruct.m_jeta_2;
	jphi_2 = jetINFOstruct.m_jphi_2;
	jm_2 = jetINFOstruct.m_jm_2;
	jmva_2 = jetINFOstruct.m_jmva_2;
	nbtag = jetINFOstruct.m_nbtag;
	nbtag_oneSigmaUp = jetINFOstruct.m_nbtag_oneSigmaUp;
	nbtag_oneSigmaDown = jetINFOstruct.m_nbtag_oneSigmaDown;
	bpt_1 = jetINFOstruct.m_bpt_1;
	beta_1 = jetINFOstruct.m_beta_1;
	bphi_1 = jetINFOstruct.m_bphi_1;
	bm_1 = jetINFOstruct.m_bm_1;
	bmva_1 = jetINFOstruct.m_bmva_1;
	bcsv_1 = jetINFOstruct.m_bcsv_1;
	bpt_2 = jetINFOstruct.m_bpt_2;
	beta_2 = jetINFOstruct.m_beta_2;
	bphi_2 = jetINFOstruct.m_bphi_2;
	bm_2 = jetINFOstruct.m_bm_2;
	bmva_2 = jetINFOstruct.m_bmva_2;
	bcsv_2 = jetINFOstruct.m_bcsv_2;
	nbtag_LooseWp = jetINFOstruct.m_nbtag_LooseWp;
	nbtag_LooseWp_oneSigmaUp = jetINFOstruct.m_nbtag_LooseWp_oneSigmaUp;
	nbtag_LooseWp_oneSigmaDown = jetINFOstruct.m_nbtag_LooseWp_oneSigmaDown;
	bpt_1_LooseWp = jetINFOstruct.m_bpt_1_LooseWp;
	beta_1_LooseWp = jetINFOstruct.m_beta_1_LooseWp;
	bphi_1_LooseWp = jetINFOstruct.m_bphi_1_LooseWp;
	bm_1_LooseWp = jetINFOstruct.m_bm_1_LooseWp;
	bmva_1_LooseWp = jetINFOstruct.m_bmva_1_LooseWp;
	bcsv_1_LooseWp = jetINFOstruct.m_bcsv_1_LooseWp;
	bpt_2_LooseWp = jetINFOstruct.m_bpt_2_LooseWp;
	beta_2_LooseWp = jetINFOstruct.m_beta_2_LooseWp;
	bphi_2_LooseWp = jetINFOstruct.m_bphi_2_LooseWp;
	bm_2_LooseWp = jetINFOstruct.m_bm_2_LooseWp;
	bmva_2_LooseWp = jetINFOstruct.m_bmva_2_LooseWp;
	bcsv_2_LooseWp = jetINFOstruct.m_bcsv_2_LooseWp;
	nbtag_TightWp = jetINFOstruct.m_nbtag_TightWp;
	nbtag_TightWp_oneSigmaUp = jetINFOstruct.m_nbtag_TightWp_oneSigmaUp;
	nbtag_TightWp_oneSigmaDown = jetINFOstruct.m_nbtag_TightWp_oneSigmaDown;
	bpt_1_TightWp = jetINFOstruct.m_bpt_1_TightWp;
	beta_1_TightWp = jetINFOstruct.m_beta_1_TightWp;
	bphi_1_TightWp = jetINFOstruct.m_bphi_1_TightWp;
	bm_1_TightWp = jetINFOstruct.m_bm_1_TightWp;
	bmva_1_TightWp = jetINFOstruct.m_bmva_1_TightWp;
	bcsv_1_TightWp = jetINFOstruct.m_bcsv_1_TightWp;
	bpt_2_TightWp = jetINFOstruct.m_bpt_2_TightWp;
	beta_2_TightWp = jetINFOstruct.m_beta_2_TightWp;
	bphi_2_TightWp = jetINFOstruct.m_bphi_2_TightWp;
	bm_2_TightWp = jetINFOstruct.m_bm_2_TightWp;
	bmva_2_TightWp = jetINFOstruct.m_bmva_2_TightWp;
	bcsv_2_TightWp = jetINFOstruct.m_bcsv_2_TightWp;

   	/* for _JECshiftedUp jets */
   	argString = "_JECshiftedUp";

   	jetINFOstruct.reset();
	fillJetBranches(eventHasNominalLeptonEnergyScales, eventIsNotSmallTree, argString, jetINFOstruct);

    /* the event-based mono-H style btag scale factors, necessary values fixed, different branch name format */
   
   BtagEventSFproduct_looseWpDown = R.getD("BtagEventSFproduct_looseWpDown");
   BtagEventSFproduct_looseWpCentral = R.getD("BtagEventSFproduct_looseWpCentral");
   BtagEventSFproduct_looseWpUp = R.getD("BtagEventSFproduct_looseWpUp");
   BtagEventSFproduct_mediumWpDown = R.getD("jets_zero_btag_event_weight_down");
   BtagEventSFproduct_mediumWpCentral = R.getD("jets_zero_btag_event_weight");
   BtagEventSFproduct_mediumWpUp = R.getD("jets_zero_btag_event_weight_up");
   BtagEventSFproduct_tightWpDown = R.getD("BtagEventSFproduct_tightWpDown");
   BtagEventSFproduct_tightWpCentral = R.getD("BtagEventSFproduct_tightWpCentral");
   BtagEventSFproduct_tightWpUp = R.getD("BtagEventSFproduct_tightWpUp");
   BtagEventSFproduct_looseWpCentral_JECshiftedUp = R.getD("BtagEventSFproduct_looseWpCentral_JECshiftedUp");
   BtagEventSFproduct_mediumWpCentral_JECshiftedUp = R.getD("jets_JECshiftedUp_zero_btag_event_weight");
   BtagEventSFproduct_tightWpCentral_JECshiftedUp = R.getD("BtagEventSFproduct_tightWpCentral_JECshiftedUp");
   BtagEventSFproduct_looseWpCentral_JECshiftedDown = R.getD("BtagEventSFproduct_looseWpCentral_JECshiftedDown");
   BtagEventSFproduct_mediumWpCentral_JECshiftedDown = R.getD("jets_JECshiftedDown_zero_btag_event_weight");
   BtagEventSFproduct_tightWpCentral_JECshiftedDown = R.getD("BtagEventSFproduct_tightWpCentral_JECshiftedDown");
   BtagEventSFproduct_looseWpCentral_JERup = R.getD("BtagEventSFproduct_looseWpCentral_JERup");
   BtagEventSFproduct_mediumWpCentral_JERup = R.getD("jets_JERup_zero_btag_event_weight");
   BtagEventSFproduct_tightWpCentral_JERup = R.getD("BtagEventSFproduct_tightWpCentral_JERup");
   BtagEventSFproduct_looseWpCentral_JERdown = R.getD("BtagEventSFproduct_looseWpCentral_JERdown");
   BtagEventSFproduct_mediumWpCentral_JERdown = R.getD("jets_JERdown_zero_btag_event_weight");
   BtagEventSFproduct_tightWpCentral_JERdown = R.getD("BtagEventSFproduct_tightWpCentral_JERdown");
   
   /*
   BtagEventSFproduct_looseWpDown = R.getD("BtagEventSFproduct_looseWpDown");
   BtagEventSFproduct_looseWpCentral = R.getD("BtagEventSFproduct_looseWpCentral");
   BtagEventSFproduct_looseWpUp = R.getD("BtagEventSFproduct_looseWpUp");
   BtagEventSFproduct_mediumWpDown = R.getD("BtagEventSFproduct_mediumWpDown");
   BtagEventSFproduct_mediumWpCentral = R.getD("BtagEventSFproduct_mediumWpCentral");
   BtagEventSFproduct_mediumWpUp = R.getD("BtagEventSFproduct_mediumWpUp");
   BtagEventSFproduct_tightWpDown = R.getD("BtagEventSFproduct_tightWpDown");
   BtagEventSFproduct_tightWpCentral = R.getD("BtagEventSFproduct_tightWpCentral");
   BtagEventSFproduct_tightWpUp = R.getD("BtagEventSFproduct_tightWpUp");
   BtagEventSFproduct_looseWpCentral_JECshiftedUp = R.getD("BtagEventSFproduct_looseWpCentral_JECshiftedUp");
   BtagEventSFproduct_mediumWpCentral_JECshiftedUp = R.getD("BtagEventSFproduct_mediumWpCentral_JECshiftedUp");
   BtagEventSFproduct_tightWpCentral_JECshiftedUp = R.getD("BtagEventSFproduct_tightWpCentral_JECshiftedUp");
   BtagEventSFproduct_looseWpCentral_JECshiftedDown = R.getD("BtagEventSFproduct_looseWpCentral_JECshiftedDown");
   BtagEventSFproduct_mediumWpCentral_JECshiftedDown = R.getD("BtagEventSFproduct_mediumWpCentral_JECshiftedDown");
   BtagEventSFproduct_tightWpCentral_JECshiftedDown = R.getD("BtagEventSFproduct_tightWpCentral_JECshiftedDown");
   BtagEventSFproduct_looseWpCentral_JERup = R.getD("BtagEventSFproduct_looseWpCentral_JERup");
   BtagEventSFproduct_mediumWpCentral_JERup = R.getD("BtagEventSFproduct_mediumWpCentral_JERup");
   BtagEventSFproduct_tightWpCentral_JERup = R.getD("BtagEventSFproduct_tightWpCentral_JERup");
   BtagEventSFproduct_looseWpCentral_JERdown = R.getD("BtagEventSFproduct_looseWpCentral_JERdown");
   BtagEventSFproduct_mediumWpCentral_JERdown = R.getD("BtagEventSFproduct_mediumWpCentral_JERdown");
   BtagEventSFproduct_tightWpCentral_JERdown = R.getD("BtagEventSFproduct_tightWpCentral_JERdown");
   */
    
	njets_JECshiftedUp = jetINFOstruct.m_njets;
	njetspt20_JECshiftedUp = jetINFOstruct.m_njetspt20;
	mjj_JECshiftedUp = jetINFOstruct.m_mjj;
	jdeta_JECshiftedUp = jetINFOstruct.m_jdeta;
	njetingap_JECshiftedUp = jetINFOstruct.m_njetingap;
	njetingap20_JECshiftedUp = jetINFOstruct.m_njetingap20;
	jdphi_JECshiftedUp = jetINFOstruct.m_jdphi;
	jpt_1_JECshiftedUp = jetINFOstruct.m_jpt_1;
	jeta_1_JECshiftedUp = jetINFOstruct.m_jeta_1;
	jphi_1_JECshiftedUp = jetINFOstruct.m_jphi_1;
	jm_1_JECshiftedUp = jetINFOstruct.m_jm_1;
	jmva_1_JECshiftedUp = jetINFOstruct.m_jmva_1;
	jpt_2_JECshiftedUp = jetINFOstruct.m_jpt_2;
	jeta_2_JECshiftedUp = jetINFOstruct.m_jeta_2;
	jphi_2_JECshiftedUp = jetINFOstruct.m_jphi_2;
	jm_2_JECshiftedUp = jetINFOstruct.m_jm_2;
	jmva_2_JECshiftedUp = jetINFOstruct.m_jmva_2;
	nbtag_JECshiftedUp = jetINFOstruct.m_nbtag;
	bpt_1_JECshiftedUp = jetINFOstruct.m_bpt_1;
	beta_1_JECshiftedUp = jetINFOstruct.m_beta_1;
	bphi_1_JECshiftedUp = jetINFOstruct.m_bphi_1;
	bm_1_JECshiftedUp = jetINFOstruct.m_bm_1;
	bmva_1_JECshiftedUp = jetINFOstruct.m_bmva_1;
	bcsv_1_JECshiftedUp = jetINFOstruct.m_bcsv_1;
	bpt_2_JECshiftedUp = jetINFOstruct.m_bpt_2;
	beta_2_JECshiftedUp = jetINFOstruct.m_beta_2;
	bphi_2_JECshiftedUp = jetINFOstruct.m_bphi_2;
	bm_2_JECshiftedUp = jetINFOstruct.m_bm_2;
	bmva_2_JECshiftedUp = jetINFOstruct.m_bmva_2;
	bcsv_2_JECshiftedUp = jetINFOstruct.m_bcsv_2;
	nbtag_LooseWp_JECshiftedUp = jetINFOstruct.m_nbtag_LooseWp;
	bpt_1_LooseWp_JECshiftedUp = jetINFOstruct.m_bpt_1_LooseWp;
	beta_1_LooseWp_JECshiftedUp = jetINFOstruct.m_beta_1_LooseWp;
	bphi_1_LooseWp_JECshiftedUp = jetINFOstruct.m_bphi_1_LooseWp;
	bm_1_LooseWp_JECshiftedUp = jetINFOstruct.m_bm_1_LooseWp;
	bmva_1_LooseWp_JECshiftedUp = jetINFOstruct.m_bmva_1_LooseWp;
	bcsv_1_LooseWp_JECshiftedUp = jetINFOstruct.m_bcsv_1_LooseWp;
	bpt_2_LooseWp_JECshiftedUp = jetINFOstruct.m_bpt_2_LooseWp;
	beta_2_LooseWp_JECshiftedUp = jetINFOstruct.m_beta_2_LooseWp;
	bphi_2_LooseWp_JECshiftedUp = jetINFOstruct.m_bphi_2_LooseWp;
	bm_2_LooseWp_JECshiftedUp = jetINFOstruct.m_bm_2_LooseWp;
	bmva_2_LooseWp_JECshiftedUp = jetINFOstruct.m_bmva_2_LooseWp;
	bcsv_2_LooseWp_JECshiftedUp = jetINFOstruct.m_bcsv_2_LooseWp;
	nbtag_TightWp_JECshiftedUp = jetINFOstruct.m_nbtag_TightWp;
	bpt_1_TightWp_JECshiftedUp = jetINFOstruct.m_bpt_1_TightWp;
	beta_1_TightWp_JECshiftedUp = jetINFOstruct.m_beta_1_TightWp;
	bphi_1_TightWp_JECshiftedUp = jetINFOstruct.m_bphi_1_TightWp;
	bm_1_TightWp_JECshiftedUp = jetINFOstruct.m_bm_1_TightWp;
	bmva_1_TightWp_JECshiftedUp = jetINFOstruct.m_bmva_1_TightWp;
	bcsv_1_TightWp_JECshiftedUp = jetINFOstruct.m_bcsv_1_TightWp;
	bpt_2_TightWp_JECshiftedUp = jetINFOstruct.m_bpt_2_TightWp;
	beta_2_TightWp_JECshiftedUp = jetINFOstruct.m_beta_2_TightWp;
	bphi_2_TightWp_JECshiftedUp = jetINFOstruct.m_bphi_2_TightWp;
	bm_2_TightWp_JECshiftedUp = jetINFOstruct.m_bm_2_TightWp;
	bmva_2_TightWp_JECshiftedUp = jetINFOstruct.m_bmva_2_TightWp;
	bcsv_2_TightWp_JECshiftedUp = jetINFOstruct.m_bcsv_2_TightWp;
    
    if (R.getB("isRealData")==0)
   {
        BtagEventSFproduct_or_DataTag_Central = BtagEventSFproduct_mediumWpCentral;
        BtagEventSFproduct_or_DataTag_Up = BtagEventSFproduct_mediumWpUp;
        BtagEventSFproduct_or_DataTag_Down = BtagEventSFproduct_mediumWpDown;
        BtagEventSFproduct_or_DataTag_Central_JECshiftedUp = BtagEventSFproduct_mediumWpCentral_JECshiftedUp;
        BtagEventSFproduct_or_DataTag_Central_JECshiftedDown = BtagEventSFproduct_mediumWpCentral_JECshiftedDown;
   }
   else if (R.getB("isRealData")==1)
   {
        if (nbtag==0) BtagEventSFproduct_or_DataTag_Central = 1.0;
        if (nbtag_oneSigmaUp==0) BtagEventSFproduct_or_DataTag_Up = 1.0;
        if (nbtag_oneSigmaDown==0) BtagEventSFproduct_or_DataTag_Down = 1.0;
        if (nbtag_JECshiftedUp==0) BtagEventSFproduct_or_DataTag_Central_JECshiftedUp = 1.0;
        if (nbtag_JECshiftedDown==0) BtagEventSFproduct_or_DataTag_Central_JECshiftedDown = 1.0;
   }
   
   	/* for _JECshiftedDown jets */
   	argString = "_JECshiftedDown";

   	jetINFOstruct.reset();
	fillJetBranches(eventHasNominalLeptonEnergyScales, eventIsNotSmallTree, argString, jetINFOstruct);

	njets_JECshiftedDown = jetINFOstruct.m_njets;
	njetspt20_JECshiftedDown = jetINFOstruct.m_njetspt20;
	mjj_JECshiftedDown = jetINFOstruct.m_mjj;
	jdeta_JECshiftedDown = jetINFOstruct.m_jdeta;
	njetingap_JECshiftedDown = jetINFOstruct.m_njetingap;
	njetingap20_JECshiftedDown = jetINFOstruct.m_njetingap20;
	jdphi_JECshiftedDown = jetINFOstruct.m_jdphi;
	jpt_1_JECshiftedDown = jetINFOstruct.m_jpt_1;
	jeta_1_JECshiftedDown = jetINFOstruct.m_jeta_1;
	jphi_1_JECshiftedDown = jetINFOstruct.m_jphi_1;
	jm_1_JECshiftedDown = jetINFOstruct.m_jm_1;
	jmva_1_JECshiftedDown = jetINFOstruct.m_jmva_1;
	jpt_2_JECshiftedDown = jetINFOstruct.m_jpt_2;
	jeta_2_JECshiftedDown = jetINFOstruct.m_jeta_2;
	jphi_2_JECshiftedDown = jetINFOstruct.m_jphi_2;
	jm_2_JECshiftedDown = jetINFOstruct.m_jm_2;
	jmva_2_JECshiftedDown = jetINFOstruct.m_jmva_2;
	nbtag_JECshiftedDown = jetINFOstruct.m_nbtag;
	bpt_1_JECshiftedDown = jetINFOstruct.m_bpt_1;
	beta_1_JECshiftedDown = jetINFOstruct.m_beta_1;
	bphi_1_JECshiftedDown = jetINFOstruct.m_bphi_1;
	bm_1_JECshiftedDown = jetINFOstruct.m_bm_1;
	bmva_1_JECshiftedDown = jetINFOstruct.m_bmva_1;
	bcsv_1_JECshiftedDown = jetINFOstruct.m_bcsv_1;
	bpt_2_JECshiftedDown = jetINFOstruct.m_bpt_2;
	beta_2_JECshiftedDown = jetINFOstruct.m_beta_2;
	bphi_2_JECshiftedDown = jetINFOstruct.m_bphi_2;
	bm_2_JECshiftedDown = jetINFOstruct.m_bm_2;
	bmva_2_JECshiftedDown = jetINFOstruct.m_bmva_2;
	bcsv_2_JECshiftedDown = jetINFOstruct.m_bcsv_2;
	nbtag_LooseWp_JECshiftedDown = jetINFOstruct.m_nbtag_LooseWp;
	bpt_1_LooseWp_JECshiftedDown = jetINFOstruct.m_bpt_1_LooseWp;
	beta_1_LooseWp_JECshiftedDown = jetINFOstruct.m_beta_1_LooseWp;
	bphi_1_LooseWp_JECshiftedDown = jetINFOstruct.m_bphi_1_LooseWp;
	bm_1_LooseWp_JECshiftedDown = jetINFOstruct.m_bm_1_LooseWp;
	bmva_1_LooseWp_JECshiftedDown = jetINFOstruct.m_bmva_1_LooseWp;
	bcsv_1_LooseWp_JECshiftedDown = jetINFOstruct.m_bcsv_1_LooseWp;
	bpt_2_LooseWp_JECshiftedDown = jetINFOstruct.m_bpt_2_LooseWp;
	beta_2_LooseWp_JECshiftedDown = jetINFOstruct.m_beta_2_LooseWp;
	bphi_2_LooseWp_JECshiftedDown = jetINFOstruct.m_bphi_2_LooseWp;
	bm_2_LooseWp_JECshiftedDown = jetINFOstruct.m_bm_2_LooseWp;
	bmva_2_LooseWp_JECshiftedDown = jetINFOstruct.m_bmva_2_LooseWp;
	bcsv_2_LooseWp_JECshiftedDown = jetINFOstruct.m_bcsv_2_LooseWp;
	nbtag_TightWp_JECshiftedDown = jetINFOstruct.m_nbtag_TightWp;
	bpt_1_TightWp_JECshiftedDown = jetINFOstruct.m_bpt_1_TightWp;
	beta_1_TightWp_JECshiftedDown = jetINFOstruct.m_beta_1_TightWp;
	bphi_1_TightWp_JECshiftedDown = jetINFOstruct.m_bphi_1_TightWp;
	bm_1_TightWp_JECshiftedDown = jetINFOstruct.m_bm_1_TightWp;
	bmva_1_TightWp_JECshiftedDown = jetINFOstruct.m_bmva_1_TightWp;
	bcsv_1_TightWp_JECshiftedDown = jetINFOstruct.m_bcsv_1_TightWp;
	bpt_2_TightWp_JECshiftedDown = jetINFOstruct.m_bpt_2_TightWp;
	beta_2_TightWp_JECshiftedDown = jetINFOstruct.m_beta_2_TightWp;
	bphi_2_TightWp_JECshiftedDown = jetINFOstruct.m_bphi_2_TightWp;
	bm_2_TightWp_JECshiftedDown = jetINFOstruct.m_bm_2_TightWp;
	bmva_2_TightWp_JECshiftedDown = jetINFOstruct.m_bmva_2_TightWp;
	bcsv_2_TightWp_JECshiftedDown = jetINFOstruct.m_bcsv_2_TightWp;
	
   	/* for _JERup jets */
   	argString = "_JERup";

   	jetINFOstruct.reset();
	fillJetBranches(eventHasNominalLeptonEnergyScales, eventIsNotSmallTree, argString, jetINFOstruct);

	njets_JERup = jetINFOstruct.m_njets;
	njetspt20_JERup = jetINFOstruct.m_njetspt20;
	mjj_JERup = jetINFOstruct.m_mjj;
	jdeta_JERup = jetINFOstruct.m_jdeta;
	njetingap_JERup = jetINFOstruct.m_njetingap;
	njetingap20_JERup = jetINFOstruct.m_njetingap20;
	jdphi_JERup = jetINFOstruct.m_jdphi;
	jpt_1_JERup = jetINFOstruct.m_jpt_1;
	jeta_1_JERup = jetINFOstruct.m_jeta_1;
	jphi_1_JERup = jetINFOstruct.m_jphi_1;
	jm_1_JERup = jetINFOstruct.m_jm_1;
	jmva_1_JERup = jetINFOstruct.m_jmva_1;
	jpt_2_JERup = jetINFOstruct.m_jpt_2;
	jeta_2_JERup = jetINFOstruct.m_jeta_2;
	jphi_2_JERup = jetINFOstruct.m_jphi_2;
	jm_2_JERup = jetINFOstruct.m_jm_2;
	jmva_2_JERup = jetINFOstruct.m_jmva_2;
	nbtag_JERup = jetINFOstruct.m_nbtag;
	bpt_1_JERup = jetINFOstruct.m_bpt_1;
	beta_1_JERup = jetINFOstruct.m_beta_1;
	bphi_1_JERup = jetINFOstruct.m_bphi_1;
	bm_1_JERup = jetINFOstruct.m_bm_1;
	bmva_1_JERup = jetINFOstruct.m_bmva_1;
	bcsv_1_JERup = jetINFOstruct.m_bcsv_1;
	bpt_2_JERup = jetINFOstruct.m_bpt_2;
	beta_2_JERup = jetINFOstruct.m_beta_2;
	bphi_2_JERup = jetINFOstruct.m_bphi_2;
	bm_2_JERup = jetINFOstruct.m_bm_2;
	bmva_2_JERup = jetINFOstruct.m_bmva_2;
	bcsv_2_JERup = jetINFOstruct.m_bcsv_2;
	nbtag_LooseWp_JERup = jetINFOstruct.m_nbtag_LooseWp;
	bpt_1_LooseWp_JERup = jetINFOstruct.m_bpt_1_LooseWp;
	beta_1_LooseWp_JERup = jetINFOstruct.m_beta_1_LooseWp;
	bphi_1_LooseWp_JERup = jetINFOstruct.m_bphi_1_LooseWp;
	bm_1_LooseWp_JERup = jetINFOstruct.m_bm_1_LooseWp;
	bmva_1_LooseWp_JERup = jetINFOstruct.m_bmva_1_LooseWp;
	bcsv_1_LooseWp_JERup = jetINFOstruct.m_bcsv_1_LooseWp;
	bpt_2_LooseWp_JERup = jetINFOstruct.m_bpt_2_LooseWp;
	beta_2_LooseWp_JERup = jetINFOstruct.m_beta_2_LooseWp;
	bphi_2_LooseWp_JERup = jetINFOstruct.m_bphi_2_LooseWp;
	bm_2_LooseWp_JERup = jetINFOstruct.m_bm_2_LooseWp;
	bmva_2_LooseWp_JERup = jetINFOstruct.m_bmva_2_LooseWp;
	bcsv_2_LooseWp_JERup = jetINFOstruct.m_bcsv_2_LooseWp;
	nbtag_TightWp_JERup = jetINFOstruct.m_nbtag_TightWp;
	bpt_1_TightWp_JERup = jetINFOstruct.m_bpt_1_TightWp;
	beta_1_TightWp_JERup = jetINFOstruct.m_beta_1_TightWp;
	bphi_1_TightWp_JERup = jetINFOstruct.m_bphi_1_TightWp;
	bm_1_TightWp_JERup = jetINFOstruct.m_bm_1_TightWp;
	bmva_1_TightWp_JERup = jetINFOstruct.m_bmva_1_TightWp;
	bcsv_1_TightWp_JERup = jetINFOstruct.m_bcsv_1_TightWp;
	bpt_2_TightWp_JERup = jetINFOstruct.m_bpt_2_TightWp;
	beta_2_TightWp_JERup = jetINFOstruct.m_beta_2_TightWp;
	bphi_2_TightWp_JERup = jetINFOstruct.m_bphi_2_TightWp;
	bm_2_TightWp_JERup = jetINFOstruct.m_bm_2_TightWp;
	bmva_2_TightWp_JERup = jetINFOstruct.m_bmva_2_TightWp;
	bcsv_2_TightWp_JERup = jetINFOstruct.m_bcsv_2_TightWp;

   	/* for _JERdown jets */
   	argString = "_JERdown";

   	jetINFOstruct.reset();
	fillJetBranches(eventHasNominalLeptonEnergyScales, eventIsNotSmallTree, argString, jetINFOstruct);

	njets_JERdown = jetINFOstruct.m_njets;
	njetspt20_JERdown = jetINFOstruct.m_njetspt20;
	mjj_JERdown = jetINFOstruct.m_mjj;
	jdeta_JERdown = jetINFOstruct.m_jdeta;
	njetingap_JERdown = jetINFOstruct.m_njetingap;
	njetingap20_JERdown = jetINFOstruct.m_njetingap20;
	jdphi_JERdown = jetINFOstruct.m_jdphi;
	jpt_1_JERdown = jetINFOstruct.m_jpt_1;
	jeta_1_JERdown = jetINFOstruct.m_jeta_1;
	jphi_1_JERdown = jetINFOstruct.m_jphi_1;
	jm_1_JERdown = jetINFOstruct.m_jm_1;
	jmva_1_JERdown = jetINFOstruct.m_jmva_1;
	jpt_2_JERdown = jetINFOstruct.m_jpt_2;
	jeta_2_JERdown = jetINFOstruct.m_jeta_2;
	jphi_2_JERdown = jetINFOstruct.m_jphi_2;
	jm_2_JERdown = jetINFOstruct.m_jm_2;
	jmva_2_JERdown = jetINFOstruct.m_jmva_2;
	nbtag_JERdown = jetINFOstruct.m_nbtag;
	bpt_1_JERdown = jetINFOstruct.m_bpt_1;
	beta_1_JERdown = jetINFOstruct.m_beta_1;
	bphi_1_JERdown = jetINFOstruct.m_bphi_1;
	bm_1_JERdown = jetINFOstruct.m_bm_1;
	bmva_1_JERdown = jetINFOstruct.m_bmva_1;
	bcsv_1_JERdown = jetINFOstruct.m_bcsv_1;
	bpt_2_JERdown = jetINFOstruct.m_bpt_2;
	beta_2_JERdown = jetINFOstruct.m_beta_2;
	bphi_2_JERdown = jetINFOstruct.m_bphi_2;
	bm_2_JERdown = jetINFOstruct.m_bm_2;
	bmva_2_JERdown = jetINFOstruct.m_bmva_2;
	bcsv_2_JERdown = jetINFOstruct.m_bcsv_2;
	nbtag_LooseWp_JERdown = jetINFOstruct.m_nbtag_LooseWp;
	bpt_1_LooseWp_JERdown = jetINFOstruct.m_bpt_1_LooseWp;
	beta_1_LooseWp_JERdown = jetINFOstruct.m_beta_1_LooseWp;
	bphi_1_LooseWp_JERdown = jetINFOstruct.m_bphi_1_LooseWp;
	bm_1_LooseWp_JERdown = jetINFOstruct.m_bm_1_LooseWp;
	bmva_1_LooseWp_JERdown = jetINFOstruct.m_bmva_1_LooseWp;
	bcsv_1_LooseWp_JERdown = jetINFOstruct.m_bcsv_1_LooseWp;
	bpt_2_LooseWp_JERdown = jetINFOstruct.m_bpt_2_LooseWp;
	beta_2_LooseWp_JERdown = jetINFOstruct.m_beta_2_LooseWp;
	bphi_2_LooseWp_JERdown = jetINFOstruct.m_bphi_2_LooseWp;
	bm_2_LooseWp_JERdown = jetINFOstruct.m_bm_2_LooseWp;
	bmva_2_LooseWp_JERdown = jetINFOstruct.m_bmva_2_LooseWp;
	bcsv_2_LooseWp_JERdown = jetINFOstruct.m_bcsv_2_LooseWp;
	nbtag_TightWp_JERdown = jetINFOstruct.m_nbtag_TightWp;
	bpt_1_TightWp_JERdown = jetINFOstruct.m_bpt_1_TightWp;
	beta_1_TightWp_JERdown = jetINFOstruct.m_beta_1_TightWp;
	bphi_1_TightWp_JERdown = jetINFOstruct.m_bphi_1_TightWp;
	bm_1_TightWp_JERdown = jetINFOstruct.m_bm_1_TightWp;
	bmva_1_TightWp_JERdown = jetINFOstruct.m_bmva_1_TightWp;
	bcsv_1_TightWp_JERdown = jetINFOstruct.m_bcsv_1_TightWp;
	bpt_2_TightWp_JERdown = jetINFOstruct.m_bpt_2_TightWp;
	beta_2_TightWp_JERdown = jetINFOstruct.m_beta_2_TightWp;
	bphi_2_TightWp_JERdown = jetINFOstruct.m_bphi_2_TightWp;
	bm_2_TightWp_JERdown = jetINFOstruct.m_bm_2_TightWp;
	bmva_2_TightWp_JERdown = jetINFOstruct.m_bmva_2_TightWp;
	bcsv_2_TightWp_JERdown = jetINFOstruct.m_bcsv_2_TightWp;
	
	// extra lepton + dilepton vetoes

	dilepton_veto = 0.0;
	extraelec_veto = 0.0;
	extramuon_veto = 0.0;

	if(R.getF("DiMuon_Flag") > 0.5 || R.getF("DiElectron_Flag") > 0.5) dilepton_veto = 1.0;
	if(R.getF("ThirdElectron_Flag") > 0.5) extraelec_veto = 1.0;
	if(R.getF("ThirdMuon_Flag") > 0.5) extramuon_veto = 1.0;

	/* veto leptons needed as is from FlatTuple for VH channels of SM Higgs analysis */

	veto_leptonType = R.getVI("veto_leptonType");
	veto_pt = R.getVD("veto_pt");
	veto_eta = R.getVD("veto_eta");
	veto_phi = R.getVD("veto_phi");
	veto_M = R.getVD("veto_M");
	veto_charge = R.getVI("veto_charge");
	veto_dxy = R.getVF("veto_dxy");
	veto_dz = R.getVF("veto_dz");
	veto_RelIso = R.getVF("veto_RelIso");
    veto_passesLooseMuonId = R.getVF("veto_passesLooseMuonId");
    veto_passesMediumMuonId_ICHEP16 = R.getVF("veto_passesMediumMuonId_ICHEP16");
    veto_passesMediumMuonId_Moriond17 = R.getVF("veto_passesMediumMuonId_Moriond17");
    veto_passesTightMuonId = R.getVF("veto_passesTightMuonId");

	veto_passesMediumMuonId = R.getVF("veto_passesMediumMuonId");
	veto_passElectronMVA80 = R.getVF("veto_passElectronMVA80");
	veto_passElectronMVA90 = R.getVF("veto_passElectronMVA90");
	veto_passVetoElectronCutBased = R.getVF("veto_passVetoElectronCutBased");
    veto_passTightElectronCutBased = R.getVF("veto_passTightElectronCutBased");
	veto_isTrackerGlobalPFMuon = R.getVF("veto_isTrackerGlobalPFMuon");
	veto_numberOfMissingInnerHits = R.getVF("veto_numberOfMissingInnerHits");
	veto_numberOfMissingOuterHits = R.getVF("veto_numberOfMissingOuterHits");
	veto_passConversionVeto = R.getVF("veto_passConversionVeto");
	veto_LeptonPassesThirdElectronVetoCuts = R.getVI("veto_LeptonPassesThirdElectronVetoCuts");
	veto_LeptonPassesThirdMuonVetoCuts = R.getVI("veto_LeptonPassesThirdMuonVetoCuts");
	veto_LeptonPassesDiElectronVetoCuts = R.getVI("veto_LeptonPassesDiElectronVetoCuts");
	veto_LeptonPassesDiMuonVetoCuts = R.getVI("veto_LeptonPassesDiMuonVetoCuts");

	/* MET FILTERS */
	HBHENoiseFilter 					= R.getB("HBHENoiseFilter");
	HBHENoiseIsoFilter 				  	= R.getB("HBHENoiseIsoFilter");
	CSCTightHalo2015Filter 			   	= R.getB("CSCTightHalo2015Filter");
	EcalDeadCellTriggerPrimitiveFilter 	= R.getB("EcalDeadCellTriggerPrimitiveFilter");
	goodVerticesFilter 					= R.getB("goodVerticesFilter");
	eeBadScFilter 						= R.getB("eeBadScFilter");
	chargedHadronTrackResolutionFilter 	= R.getB("chargedHadronTrackResolutionFilter");
	muonBadTrackFilter 					= R.getB("muonBadTrackFilter");
	globalTightHalo2016Filter    		= R.getB("globalTightHalo2016Filter");
	BadChargedCandidateFilter	    	= R.getB("BadChargedCandidateFilter");
	BadPFMuonFilter	    				= R.getB("BadPFMuonFilter");

    /* bad/duplicate muon taggers */

    BadMuonTaggedMoriond17              = R.getB("BadMuonTaggedMoriond17");
    DuplicateMuonTaggedMoriond17        = R.getB("DuplicateMuonTaggedMoriond17");

	/* lhe and gen info */

	NUP = R.getI("hepNUP");
	weight = 1.0; //weight  to be computed
	generatorEventWeight = R.getD("generatorEventWeight");
	lheHT = R.getD("lheHT");
	lheOutGoingPartons = R.getI("lheOutGoingPartons");
	lheZmass = R.getD("lheZmass");

	/* dy classification */
	IsZTT = R.getI("IsZTT");
	IsZL = R.getI("IsZL");
	IsZJ = R.getI("IsZJ");
	IsZLL = R.getI("IsZLL");
    
    bool ltauTTTcheck = ((R.getI("leg2_MCMatchType") == 5) && (((R.getI("CandidateEventType") == 5)) || (R.getI("CandidateEventType") == 3)));
    bool tautauTTTcheck = ((R.getI("leg1_MCMatchType") == 5) && (R.getI("leg2_MCMatchType") == 5) && (R.getI("CandidateEventType") == 6));
    
    /* tt classification */
    if(ltauTTTcheck || tautauTTTcheck)
    {
        IsTTT = 1;
	}
    
	/* summary info */
	DataSet = R.getS("DataSet");
	EventTotal = R.getI("EventTotal");
	EventType = R.getS("EventType");
	KeyName = R.getS("KeyName");
	DataCard = R.getS("DataCard");
	CrossSection = R.getD("CrossSection");
	FilterEff = R.getD("FilterEff");
	isSmallTree = R.getB("isSmallTree");
	TauEsNumberSigmasShifted = R.getF("TauEsNumberSigmasShifted");
	ElectronEsNumberSigmasShifted = R.getF("ElectronEsNumberSigmasShifted");

	// 2016 triggers -- will cause crash if applied different samples
	
	leg1_HLT_DoubleMediumIsoPFTau35_Trk1_eta2p1_Reg = R.getF("leg1_HLT_DoubleMediumIsoPFTau35_Trk1_eta2p1_Reg");
    leg1_HLT_DoubleMediumCombinedIsoPFTau35_Trk1_eta2p1_Reg = R.getF("leg1_HLT_DoubleMediumCombinedIsoPFTau35_Trk1_eta2p1_Reg");
	leg1_HLT_Ele25_eta2p1_WPTight_Gsf = R.getF("leg1_HLT_Ele25_eta2p1_WPTight_Gsf");
    leg1_HLT_Ele45_WPLoose_Gsf_L1JetTauSeeded = R.getF("leg1_HLT_Ele45_WPLoose_Gsf_L1JetTauSeeded");
	leg1_HLT_IsoMu22 = R.getF("leg1_HLT_IsoMu22");
	leg1_HLT_IsoTkMu22 = R.getF("leg1_HLT_IsoTkMu22");
    leg1_HLT_IsoMu22_eta2p1 = R.getF("leg1_HLT_IsoMu22_eta2p1");
    leg1_HLT_IsoTkMu22_eta2p1 = R.getF("leg1_HLT_IsoTkMu22_eta2p1");
    leg1_HLT_IsoMu24 = R.getF("leg1_HLT_IsoMu24");
    leg1_HLT_IsoTkMu24 = R.getF("leg1_HLT_IsoTkMu24");
    
	leg2_HLT_DoubleMediumIsoPFTau35_Trk1_eta2p1_Reg = R.getF("leg2_HLT_DoubleMediumIsoPFTau35_Trk1_eta2p1_Reg");
    leg2_HLT_DoubleMediumCombinedIsoPFTau35_Trk1_eta2p1_Reg = R.getF("leg2_HLT_DoubleMediumCombinedIsoPFTau35_Trk1_eta2p1_Reg");
	leg2_HLT_Ele25_eta2p1_WPTight_Gsf = R.getF("leg2_HLT_Ele25_eta2p1_WPTight_Gsf");
    leg2_HLT_Ele45_WPLoose_Gsf_L1JetTauSeeded = R.getF("leg2_HLT_Ele45_WPLoose_Gsf_L1JetTauSeeded");
	leg2_HLT_IsoMu22 = R.getF("leg2_HLT_IsoMu22");
	leg2_HLT_IsoTkMu22 = R.getF("leg2_HLT_IsoTkMu22");
    leg2_HLT_IsoMu22_eta2p1 = R.getF("leg2_HLT_IsoMu22_eta2p1");
    leg2_HLT_IsoTkMu22_eta2p1 = R.getF("leg2_HLT_IsoTkMu22_eta2p1");
    leg2_HLT_IsoMu24 = R.getF("leg2_HLT_IsoMu24");
    leg2_HLT_IsoTkMu24 = R.getF("leg2_HLT_IsoTkMu24");

    //use if NO triggers in MC
    /*
    if (R.getB("isRealData") == 1)
    {
        if(R.getI("CandidateEventType")==3 && (leg1_HLT_Ele25_eta2p1_WPTight_Gsf>0.5))
        {
            pairGoodForTrigger = 1;
        }
        else if(R.getI("CandidateEventType")==5 && (leg1_HLT_IsoMu24>0.5 || leg1_HLT_IsoTkMu24>0.5))
        {
            pairGoodForTrigger = 1;
        }
        else if(R.getI("CandidateEventType")==6 && ((leg1_HLT_DoubleMediumIsoPFTau35_Trk1_eta2p1_Reg>0.5 && leg2_HLT_DoubleMediumIsoPFTau35_Trk1_eta2p1_Reg>0.5) || (leg1_HLT_DoubleMediumCombinedIsoPFTau35_Trk1_eta2p1_Reg>0.5 && leg2_HLT_DoubleMediumCombinedIsoPFTau35_Trk1_eta2p1_Reg>0.5)))
        {
            pairGoodForTrigger = 1;
        }
        else pairGoodForTrigger = 0;
    
    }
    else pairGoodForTrigger = 1;
    */
    
    if(R.getI("CandidateEventType")==3 && (leg1_HLT_Ele25_eta2p1_WPTight_Gsf>0.5))
    {
        pairGoodForTrigger = 1;
    }
    else if(R.getI("CandidateEventType")==5 && (leg1_HLT_IsoMu24>0.5 || leg1_HLT_IsoTkMu24>0.5))
    {
        pairGoodForTrigger = 1;
    }
    else if(R.getI("CandidateEventType")==6 && ((leg1_HLT_DoubleMediumIsoPFTau35_Trk1_eta2p1_Reg>0.5 && leg2_HLT_DoubleMediumIsoPFTau35_Trk1_eta2p1_Reg>0.5) || (leg1_HLT_DoubleMediumCombinedIsoPFTau35_Trk1_eta2p1_Reg>0.5 && leg2_HLT_DoubleMediumCombinedIsoPFTau35_Trk1_eta2p1_Reg>0.5)))
    {
        pairGoodForTrigger = 1;
    }
    else {pairGoodForTrigger = 0;}
    
    //Use instead of above elseif IF not using MC triggers in Di-Tau
    
    /*
    if(R.getI("CandidateEventType")==6)
    {
        if (R.getB("isRealData") == 0)
        {
            pairGoodForTrigger = 1;
        }
        else if(R.getB("isRealData") == 1 && ((leg1_HLT_DoubleMediumIsoPFTau35_Trk1_eta2p1_Reg>0.5 && leg2_HLT_DoubleMediumIsoPFTau35_Trk1_eta2p1_Reg>0.5) || (leg1_HLT_DoubleMediumCombinedIsoPFTau35_Trk1_eta2p1_Reg>0.5 && leg2_HLT_DoubleMediumCombinedIsoPFTau35_Trk1_eta2p1_Reg>0.5)))
        {
            pairGoodForTrigger = 1;
        }
    }
    */
    
    //get LPT value
    LPT = computeLPT(0,l1,l2);
    LPT_TESUp = computeLPT(0,l1TESUp,l2TESUp);
    LPT_TESDown = computeLPT(0,l1TESDown,l2TESDown);
    LPT_L1TESUp = computeLPT(0,l1TESUp,l2);
    LPT_L1TESDown = computeLPT(0,l1TESDown,l2);
    LPT_L2TESUp = computeLPT(0,l1,l2TESUp);
    LPT_L2TESDown = computeLPT(0,l1,l2TESDown);
    
	// Pchi and Mmin Calculations

	std::vector <double> dummyPchiMmin;

	dummyPchiMmin = computePchi_and_Mmin(0, met, metphi, l1);
    P_chi = dummyPchiMmin[0];
    M_min = dummyPchiMmin[1];

	dummyPchiMmin = computePchi_and_Mmin(0, uncorr_mvamet, uncorr_mvametphi, l1);
	P_chi_uncorr = dummyPchiMmin[0];
	M_min_uncorr= dummyPchiMmin[1];
	
	dummyPchiMmin = computePchi_and_Mmin(0, met, metphi, l1);
	P_chi_pf = dummyPchiMmin[0];
	M_min_pf = dummyPchiMmin[1];
	
	dummyPchiMmin = computePchi_and_Mmin(0, puppimet, puppimetphi, l1);
	P_chi_puppi = dummyPchiMmin[0];
	M_min_puppi = dummyPchiMmin[1];
	

	if(eventHasNominalLeptonEnergyScales)
	{
    
        //MVAMET
        
		dummyPchiMmin = computePchi_and_Mmin(0, responseUP_mvaMET, responseUP_mvaMETphi, l1);
		P_chi_responseUP = dummyPchiMmin[0];
		M_min_responseUP = dummyPchiMmin[1];
		
		dummyPchiMmin = computePchi_and_Mmin(0, responseDOWN_mvaMET, responseDOWN_mvaMETphi, l1);
		P_chi_responseDOWN = dummyPchiMmin[0];
		M_min_responseDOWN = dummyPchiMmin[1];
		
		dummyPchiMmin = computePchi_and_Mmin(0, resolutionUP_mvaMET, resolutionUP_mvaMETphi, l1);
		P_chi_resolutionUP = dummyPchiMmin[0];
		M_min_resolutionUP = dummyPchiMmin[1];
		
		dummyPchiMmin = computePchi_and_Mmin(0, resolutionDOWN_mvaMET, resolutionDOWN_mvaMETphi, l1);
		P_chi_resolutionDOWN = dummyPchiMmin[0];
		M_min_resolutionDOWN = dummyPchiMmin[1];
        
        //PFMET
        
        dummyPchiMmin = computePchi_and_Mmin(0, pfmet_type1_UnclusteredEnUp_Pt, pfmet_type1_UnclusteredEnUp_Phi, l1);
		P_chi_pf_UESUp = dummyPchiMmin[0];
        
        dummyPchiMmin = computePchi_and_Mmin(0, pfmet_type1_UnclusteredEnDown_Pt, pfmet_type1_UnclusteredEnDown_Phi, l1);
		P_chi_pf_UESDown = dummyPchiMmin[0];
        
        dummyPchiMmin = computePchi_and_Mmin(0, pfmet_type1_TESUp_Pt, pfmet_type1_TESUp_Phi, l1TESUp);
		P_chi_pf_TESUp = dummyPchiMmin[0];
        
        dummyPchiMmin = computePchi_and_Mmin(0, pfmet_type1_TESDown_Pt, pfmet_type1_TESDown_Phi, l1TESDown);
		P_chi_pf_TESDown = dummyPchiMmin[0];
        
        dummyPchiMmin = computePchi_and_Mmin(0, pfmet_type1_JetEnUp_Pt, pfmet_type1_JetEnUp_Phi, l1);
		P_chi_pf_JEnUp = dummyPchiMmin[0];
        
        dummyPchiMmin = computePchi_and_Mmin(0, pfmet_type1_JetEnDown_Pt, pfmet_type1_JetEnDown_Phi, l1);
		P_chi_pf_JEnDown = dummyPchiMmin[0];
        
	}
    
	//////////////////////////////////////////////////////////////////////////////

	// Nominal TMVA variable calculations
    
    rpt_1 = pt_1;
    rpt_2 = pt_2;
    rpfmt_1 = pfmt_1;
    rmt_tot = mt_tot;
    rpt_tt = pt_tt;
    rm_vis = m_vis;
    rmet = met;
    rP_chi_pf = P_chi_pf;
    rLPT = LPT;
    rDeltaR_leg1_leg2 = DeltaR_leg1_leg2;
    rcos_DeltaPhi_leg1_leg2 = cos_DeltaPhi_leg1_leg2;
    rcos_DeltaPhi_PFMET_Higgs = cos_DeltaPhi_PFMET_Higgs;
    
    rratio_weight = -999.0;
    rfinal_weight = final_weight;
    rnpu = npu;
    revent = evt;
    rrandNum = randNum;
    rDataCardInt = DataCardInt;
    rIsZTT = IsZTT;
    rIsZJ = IsZJ;
    rIsZL = IsZL;
    rIsTTT = IsTTT;
    
    mvaVar_mt_MZP600A0400 = mt_MZP600A0400_reader->EvaluateMVA( "MLPBNN" );
    mvaVar_mt_MZP800A0400 = mt_MZP800A0400_reader->EvaluateMVA( "MLPBNN" );
    mvaVar_mt_MZP1000A0400 = mt_MZP1000A0400_reader->EvaluateMVA( "MLPBNN" );
    mvaVar_mt_MZP1200A0400 = mt_MZP1200A0400_reader->EvaluateMVA( "MLPBNN" );
    
    mvaVar_et_MZP600A0400 = et_MZP600A0400_reader->EvaluateMVA( "MLPBNN" );
    mvaVar_et_MZP800A0400 = et_MZP800A0400_reader->EvaluateMVA( "MLPBNN" );
    mvaVar_et_MZP1000A0400 = et_MZP1000A0400_reader->EvaluateMVA("MLPBNN");
    mvaVar_et_MZP1200A0400 = et_MZP1200A0400_reader->EvaluateMVA( "MLPBNN" );
    
    mvaVar_tt_MZP600A0400 = tt_MZP600A0400_reader->EvaluateMVA( "MLPBNN" );
    mvaVar_tt_MZP800A0400 = tt_MZP800A0400_reader->EvaluateMVA( "MLPBNN" );
    mvaVar_tt_MZP1000A0400 = tt_MZP1000A0400_reader->EvaluateMVA( "MLPBNN" );
    mvaVar_tt_MZP1200A0400 = tt_MZP1200A0400_reader->EvaluateMVA( "MLPBNN" );
    
    rpfmt_1 = pfmt_1_UESUp;
    rmt_tot = mt_tot_UESUp;
    rmet = pfmet_type1_UnclusteredEnUp_Pt;
    rP_chi_pf = P_chi_pf_UESUp;
    rcos_DeltaPhi_PFMET_Higgs = cos_DeltaPhi_PFMET_Higgs_UESUp;
    
    mvaVar_mt_MZP600A0400_UESUp = mt_MZP600A0400_reader->EvaluateMVA( "MLPBNN" );
    mvaVar_mt_MZP800A0400_UESUp = mt_MZP800A0400_reader->EvaluateMVA( "MLPBNN" );
    mvaVar_mt_MZP1000A0400_UESUp = mt_MZP1000A0400_reader->EvaluateMVA( "MLPBNN" );
    mvaVar_mt_MZP1200A0400_UESUp = mt_MZP1200A0400_reader->EvaluateMVA( "MLPBNN" );
    
    mvaVar_et_MZP600A0400_UESUp = et_MZP600A0400_reader->EvaluateMVA( "MLPBNN" );
    mvaVar_et_MZP800A0400_UESUp = et_MZP800A0400_reader->EvaluateMVA( "MLPBNN" );
    mvaVar_et_MZP1000A0400_UESUp = et_MZP1000A0400_reader->EvaluateMVA("MLPBNN");
    mvaVar_et_MZP1200A0400_UESUp = et_MZP1200A0400_reader->EvaluateMVA( "MLPBNN" );
    
    mvaVar_tt_MZP600A0400_UESUp = tt_MZP600A0400_reader->EvaluateMVA( "MLPBNN" );
    mvaVar_tt_MZP800A0400_UESUp = tt_MZP800A0400_reader->EvaluateMVA( "MLPBNN" );
    mvaVar_tt_MZP1000A0400_UESUp = tt_MZP1000A0400_reader->EvaluateMVA( "MLPBNN" );
    mvaVar_tt_MZP1200A0400_UESUp = tt_MZP1200A0400_reader->EvaluateMVA( "MLPBNN" );
    
    rpfmt_1 = pfmt_1_UESDown;
    rmt_tot = mt_tot_UESDown;
    rmet = pfmet_type1_UnclusteredEnDown_Pt;
    rP_chi_pf = P_chi_pf_UESDown;
    rcos_DeltaPhi_PFMET_Higgs = cos_DeltaPhi_PFMET_Higgs_UESDown;
    
    mvaVar_mt_MZP600A0400_UESDown = mt_MZP600A0400_reader->EvaluateMVA( "MLPBNN" );
    mvaVar_mt_MZP800A0400_UESDown = mt_MZP800A0400_reader->EvaluateMVA( "MLPBNN" );
    mvaVar_mt_MZP1000A0400_UESDown = mt_MZP1000A0400_reader->EvaluateMVA( "MLPBNN" );
    mvaVar_mt_MZP1200A0400_UESDown = mt_MZP1200A0400_reader->EvaluateMVA( "MLPBNN" );
    
    mvaVar_et_MZP600A0400_UESDown = et_MZP600A0400_reader->EvaluateMVA( "MLPBNN" );
    mvaVar_et_MZP800A0400_UESDown = et_MZP800A0400_reader->EvaluateMVA( "MLPBNN" );
    mvaVar_et_MZP1000A0400_UESDown = et_MZP1000A0400_reader->EvaluateMVA("MLPBNN");
    mvaVar_et_MZP1200A0400_UESDown = et_MZP1200A0400_reader->EvaluateMVA( "MLPBNN" );
    
    mvaVar_tt_MZP600A0400_UESDown = tt_MZP600A0400_reader->EvaluateMVA( "MLPBNN" );
    mvaVar_tt_MZP800A0400_UESDown = tt_MZP800A0400_reader->EvaluateMVA( "MLPBNN" );
    mvaVar_tt_MZP1000A0400_UESDown = tt_MZP1000A0400_reader->EvaluateMVA( "MLPBNN" );
    mvaVar_tt_MZP1200A0400_UESDown = tt_MZP1200A0400_reader->EvaluateMVA( "MLPBNN" );
    
    rpfmt_1 = pfmt_1_JEnUp;
    rmt_tot = mt_tot_JEnUp;
    rmet = pfmet_type1_JetEnUp_Pt;
    rP_chi_pf = P_chi_pf_JEnUp;
    rcos_DeltaPhi_PFMET_Higgs = cos_DeltaPhi_PFMET_Higgs_JEnUp;
    
    mvaVar_mt_MZP600A0400_JEnUp = mt_MZP600A0400_reader->EvaluateMVA( "MLPBNN" );
    mvaVar_mt_MZP800A0400_JEnUp = mt_MZP800A0400_reader->EvaluateMVA( "MLPBNN" );
    mvaVar_mt_MZP1000A0400_JEnUp = mt_MZP1000A0400_reader->EvaluateMVA( "MLPBNN" );
    mvaVar_mt_MZP1200A0400_JEnUp = mt_MZP1200A0400_reader->EvaluateMVA( "MLPBNN" );
    
    mvaVar_et_MZP600A0400_JEnUp = et_MZP600A0400_reader->EvaluateMVA( "MLPBNN" );
    mvaVar_et_MZP800A0400_JEnUp = et_MZP800A0400_reader->EvaluateMVA( "MLPBNN" );
    mvaVar_et_MZP1000A0400_JEnUp = et_MZP1000A0400_reader->EvaluateMVA("MLPBNN");
    mvaVar_et_MZP1200A0400_JEnUp = et_MZP1200A0400_reader->EvaluateMVA( "MLPBNN" );
    
    mvaVar_tt_MZP600A0400_JEnUp = tt_MZP600A0400_reader->EvaluateMVA( "MLPBNN" );
    mvaVar_tt_MZP800A0400_JEnUp = tt_MZP800A0400_reader->EvaluateMVA( "MLPBNN" );
    mvaVar_tt_MZP1000A0400_JEnUp = tt_MZP1000A0400_reader->EvaluateMVA( "MLPBNN" );
    mvaVar_tt_MZP1200A0400_JEnUp = tt_MZP1200A0400_reader->EvaluateMVA( "MLPBNN" );
    
    rpfmt_1 = pfmt_1_JEnDown;
    rmt_tot = mt_tot_JEnDown;
    rmet = pfmet_type1_JetEnDown_Pt;
    rP_chi_pf = P_chi_pf_JEnDown;
    rcos_DeltaPhi_PFMET_Higgs = cos_DeltaPhi_PFMET_Higgs_JEnDown;
    
    mvaVar_mt_MZP600A0400_JEnDown = mt_MZP600A0400_reader->EvaluateMVA( "MLPBNN" );
    mvaVar_mt_MZP800A0400_JEnDown = mt_MZP800A0400_reader->EvaluateMVA( "MLPBNN" );
    mvaVar_mt_MZP1000A0400_JEnDown = mt_MZP1000A0400_reader->EvaluateMVA( "MLPBNN" );
    mvaVar_mt_MZP1200A0400_JEnDown = mt_MZP1200A0400_reader->EvaluateMVA( "MLPBNN" );
    
    mvaVar_et_MZP600A0400_JEnDown = et_MZP600A0400_reader->EvaluateMVA( "MLPBNN" );
    mvaVar_et_MZP800A0400_JEnDown = et_MZP800A0400_reader->EvaluateMVA( "MLPBNN" );
    mvaVar_et_MZP1000A0400_JEnDown = et_MZP1000A0400_reader->EvaluateMVA("MLPBNN");
    mvaVar_et_MZP1200A0400_JEnDown = et_MZP1200A0400_reader->EvaluateMVA( "MLPBNN" );
    
    mvaVar_tt_MZP600A0400_JEnDown = tt_MZP600A0400_reader->EvaluateMVA( "MLPBNN" );
    mvaVar_tt_MZP800A0400_JEnDown = tt_MZP800A0400_reader->EvaluateMVA( "MLPBNN" );
    mvaVar_tt_MZP1000A0400_JEnDown = tt_MZP1000A0400_reader->EvaluateMVA( "MLPBNN" );
    mvaVar_tt_MZP1200A0400_JEnDown = tt_MZP1200A0400_reader->EvaluateMVA( "MLPBNN" );
    
    rpt_1 = pt_1_TESUp;
    rpt_2 = pt_2_TESUp;
    rpfmt_1 = pfmt_1_TESUp;
    rmt_tot = mt_tot_TESUp;
    rpt_tt = pt_tt_TESUp;
    rm_vis = m_vis_TESUp;
    rmet = pfmet_type1_TESUp_Pt;
    rP_chi_pf = P_chi_pf_TESUp;
    rLPT = LPT_TESUp;
    rDeltaR_leg1_leg2 = DeltaR_leg1_leg2;
    rcos_DeltaPhi_leg1_leg2 = cos_DeltaPhi_leg1_leg2;
    rcos_DeltaPhi_PFMET_Higgs = cos_DeltaPhi_PFMET_Higgs_TESUp;
    
    mvaVar_mt_MZP600A0400_TESUp = mt_MZP600A0400_reader->EvaluateMVA( "MLPBNN" );
    mvaVar_mt_MZP800A0400_TESUp = mt_MZP800A0400_reader->EvaluateMVA( "MLPBNN" );
    mvaVar_mt_MZP1000A0400_TESUp = mt_MZP1000A0400_reader->EvaluateMVA( "MLPBNN" );
    mvaVar_mt_MZP1200A0400_TESUp = mt_MZP1200A0400_reader->EvaluateMVA( "MLPBNN" );
    
    mvaVar_et_MZP600A0400_TESUp = et_MZP600A0400_reader->EvaluateMVA( "MLPBNN" );
    mvaVar_et_MZP800A0400_TESUp = et_MZP800A0400_reader->EvaluateMVA( "MLPBNN" );
    mvaVar_et_MZP1000A0400_TESUp = et_MZP1000A0400_reader->EvaluateMVA("MLPBNN");
    mvaVar_et_MZP1200A0400_TESUp = et_MZP1200A0400_reader->EvaluateMVA( "MLPBNN" );
    
    mvaVar_tt_MZP600A0400_TESUp = tt_MZP600A0400_reader->EvaluateMVA( "MLPBNN" );
    mvaVar_tt_MZP800A0400_TESUp = tt_MZP800A0400_reader->EvaluateMVA( "MLPBNN" );
    mvaVar_tt_MZP1000A0400_TESUp = tt_MZP1000A0400_reader->EvaluateMVA( "MLPBNN" );
    mvaVar_tt_MZP1200A0400_TESUp = tt_MZP1200A0400_reader->EvaluateMVA( "MLPBNN" );
    
    rpt_1 = pt_1_TESDown;
    rpt_2 = pt_2_TESDown;
    rpfmt_1 = pfmt_1_TESDown;
    rmt_tot = mt_tot_TESDown;
    rpt_tt = pt_tt_TESDown;
    rm_vis = m_vis_TESDown;
    rmet = pfmet_type1_TESDown_Pt;
    rP_chi_pf = P_chi_pf_TESDown;
    rLPT = LPT_TESDown;
    rDeltaR_leg1_leg2 = DeltaR_leg1_leg2;
    rcos_DeltaPhi_leg1_leg2 = cos_DeltaPhi_leg1_leg2;
    rcos_DeltaPhi_PFMET_Higgs = cos_DeltaPhi_PFMET_Higgs_TESDown;
    
    mvaVar_mt_MZP600A0400_TESDown = mt_MZP600A0400_reader->EvaluateMVA( "MLPBNN" );
    mvaVar_mt_MZP800A0400_TESDown = mt_MZP800A0400_reader->EvaluateMVA( "MLPBNN" );
    mvaVar_mt_MZP1000A0400_TESDown = mt_MZP1000A0400_reader->EvaluateMVA( "MLPBNN" );
    mvaVar_mt_MZP1200A0400_TESDown = mt_MZP1200A0400_reader->EvaluateMVA( "MLPBNN" );
    
    mvaVar_et_MZP600A0400_TESDown = et_MZP600A0400_reader->EvaluateMVA( "MLPBNN" );
    mvaVar_et_MZP800A0400_TESDown = et_MZP800A0400_reader->EvaluateMVA( "MLPBNN" );
    mvaVar_et_MZP1000A0400_TESDown = et_MZP1000A0400_reader->EvaluateMVA("MLPBNN");
    mvaVar_et_MZP1200A0400_TESDown = et_MZP1200A0400_reader->EvaluateMVA( "MLPBNN" );
    
    mvaVar_tt_MZP600A0400_TESDown = tt_MZP600A0400_reader->EvaluateMVA( "MLPBNN" );
    mvaVar_tt_MZP800A0400_TESDown = tt_MZP800A0400_reader->EvaluateMVA( "MLPBNN" );
    mvaVar_tt_MZP1000A0400_TESDown = tt_MZP1000A0400_reader->EvaluateMVA( "MLPBNN" );
    mvaVar_tt_MZP1200A0400_TESDown = tt_MZP1200A0400_reader->EvaluateMVA( "MLPBNN" );
    
    /* handle separate Leg DM shifts */
    
    pt_1_dm0TESUp = pt_1;
    mt_1_dm0TESUp = mt_1;
    pfmt_1_dm0TESUp = pfmt_1;

    pt_2_dm0TESUp = pt_2;
    mt_2_dm0TESUp = mt_2;
    pfmt_2_dm0TESUp = pfmt_2;

    mt_tot_dm0TESUp = mt_tot;
    pt_tt_dm0TESUp = pt_tt;
    m_vis_dm0TESUp = m_vis;
    pfmet_type1_dm0TESUp_Pt = pfmet_type1_Pt;
    P_chi_pf_dm0TESUp = P_chi_pf;
    LPT_dm0TESUp = LPT;
    cos_DeltaPhi_PFMET_Higgs_dm0TESUp = cos_DeltaPhi_PFMET_Higgs;

    pt_1_dm0TESDown = pt_1;
    mt_1_dm0TESDown = mt_1;
    pfmt_1_dm0TESDown = pfmt_1;

    pt_2_dm0TESDown = pt_2;
    mt_2_dm0TESDown = mt_2;
    pfmt_2_dm0TESDown = pfmt_2;

    mt_tot_dm0TESDown = mt_tot;
    pt_tt_dm0TESDown = pt_tt;
    m_vis_dm0TESDown = m_vis;
    pfmet_type1_dm0TESDown_Pt = pfmet_type1_Pt;
    P_chi_pf_dm0TESDown = P_chi_pf;
    LPT_dm0TESDown = LPT;
    cos_DeltaPhi_PFMET_Higgs_dm0TESDown = cos_DeltaPhi_PFMET_Higgs;

    mvaVar_mt_MZP600A0400_dm0TESUp = mvaVar_mt_MZP600A0400;
    mvaVar_mt_MZP800A0400_dm0TESUp = mvaVar_mt_MZP800A0400;
    mvaVar_mt_MZP1000A0400_dm0TESUp = mvaVar_mt_MZP1000A0400;
    mvaVar_mt_MZP1200A0400_dm0TESUp = mvaVar_mt_MZP1200A0400;

    mvaVar_et_MZP600A0400_dm0TESUp = mvaVar_et_MZP600A0400;
    mvaVar_et_MZP800A0400_dm0TESUp = mvaVar_et_MZP800A0400;
    mvaVar_et_MZP1000A0400_dm0TESUp = mvaVar_et_MZP1000A0400;
    mvaVar_et_MZP1200A0400_dm0TESUp = mvaVar_et_MZP1200A0400;
        
    mvaVar_tt_MZP600A0400_dm0TESUp = mvaVar_tt_MZP600A0400;
    mvaVar_tt_MZP800A0400_dm0TESUp = mvaVar_tt_MZP800A0400;
    mvaVar_tt_MZP1000A0400_dm0TESUp = mvaVar_tt_MZP1000A0400;
    mvaVar_tt_MZP1200A0400_dm0TESUp = mvaVar_tt_MZP1200A0400;

    mvaVar_mt_MZP600A0400_dm0TESDown = mvaVar_mt_MZP600A0400;
    mvaVar_mt_MZP800A0400_dm0TESDown = mvaVar_mt_MZP800A0400;
    mvaVar_mt_MZP1000A0400_dm0TESDown = mvaVar_mt_MZP1000A0400;
    mvaVar_mt_MZP1200A0400_dm0TESDown = mvaVar_mt_MZP1200A0400;

    mvaVar_et_MZP600A0400_dm0TESDown = mvaVar_et_MZP600A0400;
    mvaVar_et_MZP800A0400_dm0TESDown = mvaVar_et_MZP800A0400;
    mvaVar_et_MZP1000A0400_dm0TESDown = mvaVar_et_MZP1000A0400;
    mvaVar_et_MZP1200A0400_dm0TESDown = mvaVar_et_MZP1200A0400;
        
    mvaVar_tt_MZP600A0400_dm0TESDown = mvaVar_tt_MZP600A0400;
    mvaVar_tt_MZP800A0400_dm0TESDown = mvaVar_tt_MZP800A0400;
    mvaVar_tt_MZP1000A0400_dm0TESDown = mvaVar_tt_MZP1000A0400;
    mvaVar_tt_MZP1200A0400_dm0TESDown = mvaVar_tt_MZP1200A0400;
    
    pt_1_dm1TESUp = pt_1;
    mt_1_dm1TESUp = mt_1;
    pfmt_1_dm1TESUp = pfmt_1;

    pt_2_dm1TESUp = pt_2;
    mt_2_dm1TESUp = mt_2;
    pfmt_2_dm1TESUp = pfmt_2;

    mt_tot_dm1TESUp = mt_tot;
    pt_tt_dm1TESUp = pt_tt;
    m_vis_dm1TESUp = m_vis;
    pfmet_type1_dm1TESUp_Pt = pfmet_type1_Pt;
    P_chi_pf_dm1TESUp = P_chi_pf;
    LPT_dm1TESUp = LPT;
    cos_DeltaPhi_PFMET_Higgs_dm1TESUp = cos_DeltaPhi_PFMET_Higgs;

    pt_1_dm1TESDown = pt_1;
    mt_1_dm1TESDown = mt_1;
    pfmt_1_dm1TESDown = pfmt_1;

    pt_2_dm1TESDown = pt_2;
    mt_2_dm1TESDown = mt_2;
    pfmt_2_dm1TESDown = pfmt_2;

    mt_tot_dm1TESDown = mt_tot;
    pt_tt_dm1TESDown = pt_tt;
    m_vis_dm1TESDown = m_vis;
    pfmet_type1_dm1TESDown_Pt = pfmet_type1_Pt;
    P_chi_pf_dm1TESDown = P_chi_pf;
    LPT_dm1TESDown = LPT;
    cos_DeltaPhi_PFMET_Higgs_dm1TESDown = cos_DeltaPhi_PFMET_Higgs;

    mvaVar_mt_MZP600A0400_dm1TESUp = mvaVar_mt_MZP600A0400;
    mvaVar_mt_MZP800A0400_dm1TESUp = mvaVar_mt_MZP800A0400;
    mvaVar_mt_MZP1000A0400_dm1TESUp = mvaVar_mt_MZP1000A0400;
    mvaVar_mt_MZP1200A0400_dm1TESUp = mvaVar_mt_MZP1200A0400;

    mvaVar_et_MZP600A0400_dm1TESUp = mvaVar_et_MZP600A0400;
    mvaVar_et_MZP800A0400_dm1TESUp = mvaVar_et_MZP800A0400;
    mvaVar_et_MZP1000A0400_dm1TESUp = mvaVar_et_MZP1000A0400;
    mvaVar_et_MZP1200A0400_dm1TESUp = mvaVar_et_MZP1200A0400;
        
    mvaVar_tt_MZP600A0400_dm1TESUp = mvaVar_tt_MZP600A0400;
    mvaVar_tt_MZP800A0400_dm1TESUp = mvaVar_tt_MZP800A0400;
    mvaVar_tt_MZP1000A0400_dm1TESUp = mvaVar_tt_MZP1000A0400;
    mvaVar_tt_MZP1200A0400_dm1TESUp = mvaVar_tt_MZP1200A0400;

    mvaVar_mt_MZP600A0400_dm1TESDown = mvaVar_mt_MZP600A0400;
    mvaVar_mt_MZP800A0400_dm1TESDown = mvaVar_mt_MZP800A0400;
    mvaVar_mt_MZP1000A0400_dm1TESDown = mvaVar_mt_MZP1000A0400;
    mvaVar_mt_MZP1200A0400_dm1TESDown = mvaVar_mt_MZP1200A0400;

    mvaVar_et_MZP600A0400_dm1TESDown = mvaVar_et_MZP600A0400;
    mvaVar_et_MZP800A0400_dm1TESDown = mvaVar_et_MZP800A0400;
    mvaVar_et_MZP1000A0400_dm1TESDown = mvaVar_et_MZP1000A0400;
    mvaVar_et_MZP1200A0400_dm1TESDown = mvaVar_et_MZP1200A0400;
        
    mvaVar_tt_MZP600A0400_dm1TESDown = mvaVar_tt_MZP600A0400;
    mvaVar_tt_MZP800A0400_dm1TESDown = mvaVar_tt_MZP800A0400;
    mvaVar_tt_MZP1000A0400_dm1TESDown = mvaVar_tt_MZP1000A0400;
    mvaVar_tt_MZP1200A0400_dm1TESDown = mvaVar_tt_MZP1200A0400;
    
    pt_1_dm10TESUp = pt_1;
    mt_1_dm10TESUp = mt_1;
    pfmt_1_dm10TESUp = pfmt_1;

    pt_2_dm10TESUp = pt_2;
    mt_2_dm10TESUp = mt_2;
    pfmt_2_dm10TESUp = pfmt_2;

    mt_tot_dm10TESUp = mt_tot;
    pt_tt_dm10TESUp = pt_tt;
    m_vis_dm10TESUp = m_vis;
    pfmet_type1_dm10TESUp_Pt = pfmet_type1_Pt;
    P_chi_pf_dm10TESUp = P_chi_pf;
    LPT_dm10TESUp = LPT;
    cos_DeltaPhi_PFMET_Higgs_dm10TESUp = cos_DeltaPhi_PFMET_Higgs;

    pt_1_dm10TESDown = pt_1;
    mt_1_dm10TESDown = mt_1;
    pfmt_1_dm10TESDown = pfmt_1;

    pt_2_dm10TESDown = pt_2;
    mt_2_dm10TESDown = mt_2;
    pfmt_2_dm10TESDown = pfmt_2;

    mt_tot_dm10TESDown = mt_tot;
    pt_tt_dm10TESDown = pt_tt;
    m_vis_dm10TESDown = m_vis;
    pfmet_type1_dm10TESDown_Pt = pfmet_type1_Pt;
    P_chi_pf_dm10TESDown = P_chi_pf;
    LPT_dm10TESDown = LPT;
    cos_DeltaPhi_PFMET_Higgs_dm10TESDown = cos_DeltaPhi_PFMET_Higgs;

    mvaVar_mt_MZP600A0400_dm10TESUp = mvaVar_mt_MZP600A0400;
    mvaVar_mt_MZP800A0400_dm10TESUp = mvaVar_mt_MZP800A0400;
    mvaVar_mt_MZP1000A0400_dm10TESUp = mvaVar_mt_MZP1000A0400;
    mvaVar_mt_MZP1200A0400_dm10TESUp = mvaVar_mt_MZP1200A0400;

    mvaVar_et_MZP600A0400_dm10TESUp = mvaVar_et_MZP600A0400;
    mvaVar_et_MZP800A0400_dm10TESUp = mvaVar_et_MZP800A0400;
    mvaVar_et_MZP1000A0400_dm10TESUp = mvaVar_et_MZP1000A0400;
    mvaVar_et_MZP1200A0400_dm10TESUp = mvaVar_et_MZP1200A0400;
        
    mvaVar_tt_MZP600A0400_dm10TESUp = mvaVar_tt_MZP600A0400;
    mvaVar_tt_MZP800A0400_dm10TESUp = mvaVar_tt_MZP800A0400;
    mvaVar_tt_MZP1000A0400_dm10TESUp = mvaVar_tt_MZP1000A0400;
    mvaVar_tt_MZP1200A0400_dm10TESUp = mvaVar_tt_MZP1200A0400;

    mvaVar_mt_MZP600A0400_dm10TESDown = mvaVar_mt_MZP600A0400;
    mvaVar_mt_MZP800A0400_dm10TESDown = mvaVar_mt_MZP800A0400;
    mvaVar_mt_MZP1000A0400_dm10TESDown = mvaVar_mt_MZP1000A0400;
    mvaVar_mt_MZP1200A0400_dm10TESDown = mvaVar_mt_MZP1200A0400;

    mvaVar_et_MZP600A0400_dm10TESDown = mvaVar_et_MZP600A0400;
    mvaVar_et_MZP800A0400_dm10TESDown = mvaVar_et_MZP800A0400;
    mvaVar_et_MZP1000A0400_dm10TESDown = mvaVar_et_MZP1000A0400;
    mvaVar_et_MZP1200A0400_dm10TESDown = mvaVar_et_MZP1200A0400;
        
    mvaVar_tt_MZP600A0400_dm10TESDown = mvaVar_tt_MZP600A0400;
    mvaVar_tt_MZP800A0400_dm10TESDown = mvaVar_tt_MZP800A0400;
    mvaVar_tt_MZP1000A0400_dm10TESDown = mvaVar_tt_MZP1000A0400;
    mvaVar_tt_MZP1200A0400_dm10TESDown = mvaVar_tt_MZP1200A0400;
    
    if (R.getI("CandidateEventType")==3 || R.getI("CandidateEventType")==5)
    {
    
        if (R.getI("leg2_decayMode")==0)
        {
        
            pt_1_dm0TESUp = pt_1_TESUp;
            mt_1_dm0TESUp = mt_1_TESUp;
            pfmt_1_dm0TESUp = pfmt_1_TESUp;

            pt_2_dm0TESUp = pt_2_TESUp;
            mt_2_dm0TESUp = mt_2_TESUp;
            pfmt_2_dm0TESUp = pfmt_2_TESUp;

            mt_tot_dm0TESUp = mt_tot_TESUp;
            pt_tt_dm0TESUp = pt_tt_TESUp;
            m_vis_dm0TESUp = m_vis_TESUp;
            pfmet_type1_dm0TESUp_Pt = pfmet_type1_TESUp_Pt;
            P_chi_pf_dm0TESUp = P_chi_pf_TESUp;
            LPT_dm0TESUp = LPT_TESUp;
            cos_DeltaPhi_PFMET_Higgs_dm0TESUp = cos_DeltaPhi_PFMET_Higgs_TESUp;

            pt_1_dm0TESDown = pt_1_TESDown;
            mt_1_dm0TESDown = mt_1_TESDown;
            pfmt_1_dm0TESDown = pfmt_1_TESDown;

            pt_2_dm0TESDown = pt_2_TESDown;
            mt_2_dm0TESDown = mt_2_TESDown;
            pfmt_2_dm0TESDown = pfmt_2_TESDown;

            mt_tot_dm0TESDown = mt_tot_TESDown;
            pt_tt_dm0TESDown = pt_tt_TESDown;
            m_vis_dm0TESDown = m_vis_TESDown;
            pfmet_type1_dm0TESDown_Pt = pfmet_type1_TESDown_Pt;
            P_chi_pf_dm0TESDown = P_chi_pf_TESDown;
            LPT_dm0TESDown = LPT_TESDown;
            cos_DeltaPhi_PFMET_Higgs_dm0TESDown = cos_DeltaPhi_PFMET_Higgs_TESDown;
            
            rpt_1 = pt_1_dm0TESUp;
            rpt_2 = pt_2_dm0TESUp;
            rpfmt_1 = pfmt_1_dm0TESUp;
            rmt_tot = mt_tot_dm0TESUp;
            rpt_tt = pt_tt_dm0TESUp;
            rm_vis = m_vis_dm0TESUp;
            rmet = pfmet_type1_dm0TESUp_Pt;
            rP_chi_pf = P_chi_pf_dm0TESUp;
            rLPT = LPT_dm0TESUp;
            rDeltaR_leg1_leg2 = DeltaR_leg1_leg2;
            rcos_DeltaPhi_leg1_leg2 = cos_DeltaPhi_leg1_leg2;
            rcos_DeltaPhi_PFMET_Higgs = cos_DeltaPhi_PFMET_Higgs_dm0TESUp;
            
            mvaVar_mt_MZP600A0400_dm0TESUp = mt_MZP600A0400_reader->EvaluateMVA( "MLPBNN" );
            mvaVar_mt_MZP800A0400_dm0TESUp = mt_MZP800A0400_reader->EvaluateMVA( "MLPBNN" );
            mvaVar_mt_MZP1000A0400_dm0TESUp = mt_MZP1000A0400_reader->EvaluateMVA( "MLPBNN" );
            mvaVar_mt_MZP1200A0400_dm0TESUp = mt_MZP1200A0400_reader->EvaluateMVA( "MLPBNN" );
            
            mvaVar_et_MZP600A0400_dm0TESUp = et_MZP600A0400_reader->EvaluateMVA( "MLPBNN" );
            mvaVar_et_MZP800A0400_dm0TESUp = et_MZP800A0400_reader->EvaluateMVA( "MLPBNN" );
            mvaVar_et_MZP1000A0400_dm0TESUp = et_MZP1000A0400_reader->EvaluateMVA("MLPBNN");
            mvaVar_et_MZP1200A0400_dm0TESUp = et_MZP1200A0400_reader->EvaluateMVA( "MLPBNN" );
            
            mvaVar_tt_MZP600A0400_dm0TESUp = tt_MZP600A0400_reader->EvaluateMVA( "MLPBNN" );
            mvaVar_tt_MZP800A0400_dm0TESUp = tt_MZP800A0400_reader->EvaluateMVA( "MLPBNN" );
            mvaVar_tt_MZP1000A0400_dm0TESUp = tt_MZP1000A0400_reader->EvaluateMVA( "MLPBNN" );
            mvaVar_tt_MZP1200A0400_dm0TESUp = tt_MZP1200A0400_reader->EvaluateMVA( "MLPBNN" );
            
            rpt_1 = pt_1_dm0TESDown;
            rpt_2 = pt_2_dm0TESDown;
            rpfmt_1 = pfmt_1_dm0TESDown;
            rmt_tot = mt_tot_dm0TESDown;
            rpt_tt = pt_tt_dm0TESDown;
            rm_vis = m_vis_dm0TESDown;
            rmet = pfmet_type1_dm0TESDown_Pt;
            rP_chi_pf = P_chi_pf_dm0TESDown;
            rLPT = LPT_dm0TESDown;
            rDeltaR_leg1_leg2 = DeltaR_leg1_leg2;
            rcos_DeltaPhi_leg1_leg2 = cos_DeltaPhi_leg1_leg2;
            rcos_DeltaPhi_PFMET_Higgs = cos_DeltaPhi_PFMET_Higgs_dm0TESDown;
            
            mvaVar_mt_MZP600A0400_dm0TESDown = mt_MZP600A0400_reader->EvaluateMVA( "MLPBNN" );
            mvaVar_mt_MZP800A0400_dm0TESDown = mt_MZP800A0400_reader->EvaluateMVA( "MLPBNN" );
            mvaVar_mt_MZP1000A0400_dm0TESDown = mt_MZP1000A0400_reader->EvaluateMVA( "MLPBNN" );
            mvaVar_mt_MZP1200A0400_dm0TESDown = mt_MZP1200A0400_reader->EvaluateMVA( "MLPBNN" );
            
            mvaVar_et_MZP600A0400_dm0TESDown = et_MZP600A0400_reader->EvaluateMVA( "MLPBNN" );
            mvaVar_et_MZP800A0400_dm0TESDown = et_MZP800A0400_reader->EvaluateMVA( "MLPBNN" );
            mvaVar_et_MZP1000A0400_dm0TESDown = et_MZP1000A0400_reader->EvaluateMVA("MLPBNN");
            mvaVar_et_MZP1200A0400_dm0TESDown = et_MZP1200A0400_reader->EvaluateMVA( "MLPBNN" );
            
            mvaVar_tt_MZP600A0400_dm0TESDown = tt_MZP600A0400_reader->EvaluateMVA( "MLPBNN" );
            mvaVar_tt_MZP800A0400_dm0TESDown = tt_MZP800A0400_reader->EvaluateMVA( "MLPBNN" );
            mvaVar_tt_MZP1000A0400_dm0TESDown = tt_MZP1000A0400_reader->EvaluateMVA( "MLPBNN" );
            mvaVar_tt_MZP1200A0400_dm0TESDown = tt_MZP1200A0400_reader->EvaluateMVA( "MLPBNN" );

        }
        
        if (R.getI("leg2_decayMode")==1)
        {
        
            pt_1_dm1TESUp = pt_1_TESUp;
            mt_1_dm1TESUp = mt_1_TESUp;
            pfmt_1_dm1TESUp = pfmt_1_TESUp;

            pt_2_dm1TESUp = pt_2_TESUp;
            mt_2_dm1TESUp = mt_2_TESUp;
            pfmt_2_dm1TESUp = pfmt_2_TESUp;

            mt_tot_dm1TESUp = mt_tot_TESUp;
            pt_tt_dm1TESUp = pt_tt_TESUp;
            m_vis_dm1TESUp = m_vis_TESUp;
            pfmet_type1_dm1TESUp_Pt = pfmet_type1_TESUp_Pt;
            P_chi_pf_dm1TESUp = P_chi_pf_TESUp;
            LPT_dm1TESUp = LPT_TESUp;
            cos_DeltaPhi_PFMET_Higgs_dm1TESUp = cos_DeltaPhi_PFMET_Higgs_TESUp;

            pt_1_dm1TESDown = pt_1_TESDown;
            mt_1_dm1TESDown = mt_1_TESDown;
            pfmt_1_dm1TESDown = pfmt_1_TESDown;

            pt_2_dm1TESDown = pt_2_TESDown;
            mt_2_dm1TESDown = mt_2_TESDown;
            pfmt_2_dm1TESDown = pfmt_2_TESDown;

            mt_tot_dm1TESDown = mt_tot_TESDown;
            pt_tt_dm1TESDown = pt_tt_TESDown;
            m_vis_dm1TESDown = m_vis_TESDown;
            pfmet_type1_dm1TESDown_Pt = pfmet_type1_TESDown_Pt;
            P_chi_pf_dm1TESDown = P_chi_pf_TESDown;
            LPT_dm1TESDown = LPT_TESDown;
            cos_DeltaPhi_PFMET_Higgs_dm1TESDown = cos_DeltaPhi_PFMET_Higgs_TESDown;

            rpt_1 = pt_1_dm1TESUp;
            rpt_2 = pt_2_dm1TESUp;
            rpfmt_1 = pfmt_1_dm1TESUp;
            rmt_tot = mt_tot_dm1TESUp;
            rpt_tt = pt_tt_dm1TESUp;
            rm_vis = m_vis_dm1TESUp;
            rmet = pfmet_type1_dm1TESUp_Pt;
            rP_chi_pf = P_chi_pf_dm1TESUp;
            rLPT = LPT_dm1TESUp;
            rDeltaR_leg1_leg2 = DeltaR_leg1_leg2;
            rcos_DeltaPhi_leg1_leg2 = cos_DeltaPhi_leg1_leg2;
            rcos_DeltaPhi_PFMET_Higgs = cos_DeltaPhi_PFMET_Higgs_dm1TESUp;
            
            mvaVar_mt_MZP600A0400_dm1TESUp = mt_MZP600A0400_reader->EvaluateMVA( "MLPBNN" );
            mvaVar_mt_MZP800A0400_dm1TESUp = mt_MZP800A0400_reader->EvaluateMVA( "MLPBNN" );
            mvaVar_mt_MZP1000A0400_dm1TESUp = mt_MZP1000A0400_reader->EvaluateMVA( "MLPBNN" );
            mvaVar_mt_MZP1200A0400_dm1TESUp = mt_MZP1200A0400_reader->EvaluateMVA( "MLPBNN" );
            
            mvaVar_et_MZP600A0400_dm1TESUp = et_MZP600A0400_reader->EvaluateMVA( "MLPBNN" );
            mvaVar_et_MZP800A0400_dm1TESUp = et_MZP800A0400_reader->EvaluateMVA( "MLPBNN" );
            mvaVar_et_MZP1000A0400_dm1TESUp = et_MZP1000A0400_reader->EvaluateMVA("MLPBNN");
            mvaVar_et_MZP1200A0400_dm1TESUp = et_MZP1200A0400_reader->EvaluateMVA( "MLPBNN" );
            
            mvaVar_tt_MZP600A0400_dm1TESUp = tt_MZP600A0400_reader->EvaluateMVA( "MLPBNN" );
            mvaVar_tt_MZP800A0400_dm1TESUp = tt_MZP800A0400_reader->EvaluateMVA( "MLPBNN" );
            mvaVar_tt_MZP1000A0400_dm1TESUp = tt_MZP1000A0400_reader->EvaluateMVA( "MLPBNN" );
            mvaVar_tt_MZP1200A0400_dm1TESUp = tt_MZP1200A0400_reader->EvaluateMVA( "MLPBNN" );
            
            rpt_1 = pt_1_dm1TESDown;
            rpt_2 = pt_2_dm1TESDown;
            rpfmt_1 = pfmt_1_dm1TESDown;
            rmt_tot = mt_tot_dm1TESDown;
            rpt_tt = pt_tt_dm1TESDown;
            rm_vis = m_vis_dm1TESDown;
            rmet = pfmet_type1_dm1TESDown_Pt;
            rP_chi_pf = P_chi_pf_dm1TESDown;
            rLPT = LPT_dm1TESDown;
            rDeltaR_leg1_leg2 = DeltaR_leg1_leg2;
            rcos_DeltaPhi_leg1_leg2 = cos_DeltaPhi_leg1_leg2;
            rcos_DeltaPhi_PFMET_Higgs = cos_DeltaPhi_PFMET_Higgs_dm1TESDown;
            
            mvaVar_mt_MZP600A0400_dm1TESDown = mt_MZP600A0400_reader->EvaluateMVA( "MLPBNN" );
            mvaVar_mt_MZP800A0400_dm1TESDown = mt_MZP800A0400_reader->EvaluateMVA( "MLPBNN" );
            mvaVar_mt_MZP1000A0400_dm1TESDown = mt_MZP1000A0400_reader->EvaluateMVA( "MLPBNN" );
            mvaVar_mt_MZP1200A0400_dm1TESDown = mt_MZP1200A0400_reader->EvaluateMVA( "MLPBNN" );
            
            mvaVar_et_MZP600A0400_dm1TESDown = et_MZP600A0400_reader->EvaluateMVA( "MLPBNN" );
            mvaVar_et_MZP800A0400_dm1TESDown = et_MZP800A0400_reader->EvaluateMVA( "MLPBNN" );
            mvaVar_et_MZP1000A0400_dm1TESDown = et_MZP1000A0400_reader->EvaluateMVA("MLPBNN");
            mvaVar_et_MZP1200A0400_dm1TESDown = et_MZP1200A0400_reader->EvaluateMVA( "MLPBNN" );
            
            mvaVar_tt_MZP600A0400_dm1TESDown = tt_MZP600A0400_reader->EvaluateMVA( "MLPBNN" );
            mvaVar_tt_MZP800A0400_dm1TESDown = tt_MZP800A0400_reader->EvaluateMVA( "MLPBNN" );
            mvaVar_tt_MZP1000A0400_dm1TESDown = tt_MZP1000A0400_reader->EvaluateMVA( "MLPBNN" );
            mvaVar_tt_MZP1200A0400_dm1TESDown = tt_MZP1200A0400_reader->EvaluateMVA( "MLPBNN" );
        }
        
        if (R.getI("leg2_decayMode")==10)
        {
        
            pt_1_dm10TESUp = pt_1_TESUp;
            mt_1_dm10TESUp = mt_1_TESUp;
            pfmt_1_dm10TESUp = pfmt_1_TESUp;

            pt_2_dm10TESUp = pt_2_TESUp;
            mt_2_dm10TESUp = mt_2_TESUp;
            pfmt_2_dm10TESUp = pfmt_2_TESUp;

            mt_tot_dm10TESUp = mt_tot_TESUp;
            pt_tt_dm10TESUp = pt_tt_TESUp;
            m_vis_dm10TESUp = m_vis_TESUp;
            pfmet_type1_dm10TESUp_Pt = pfmet_type1_TESUp_Pt;
            P_chi_pf_dm10TESUp = P_chi_pf_TESUp;
            LPT_dm10TESUp = LPT_TESUp;
            cos_DeltaPhi_PFMET_Higgs_dm10TESUp = cos_DeltaPhi_PFMET_Higgs_TESUp;

            pt_1_dm10TESDown = pt_1_TESDown;
            mt_1_dm10TESDown = mt_1_TESDown;
            pfmt_1_dm10TESDown = pfmt_1_TESDown;

            pt_2_dm10TESDown = pt_2_TESDown;
            mt_2_dm10TESDown = mt_2_TESDown;
            pfmt_2_dm10TESDown = pfmt_2_TESDown;

            mt_tot_dm10TESDown = mt_tot_TESDown;
            pt_tt_dm10TESDown = pt_tt_TESDown;
            m_vis_dm10TESDown = m_vis_TESDown;
            pfmet_type1_dm10TESDown_Pt = pfmet_type1_TESDown_Pt;
            P_chi_pf_dm10TESDown = P_chi_pf_TESDown;
            LPT_dm10TESDown = LPT_TESDown;
            cos_DeltaPhi_PFMET_Higgs_dm10TESDown = cos_DeltaPhi_PFMET_Higgs_TESDown;

            mvaVar_mt_MZP600A0400_dm10TESUp = mt_MZP600A0400_reader->EvaluateMVA( "MLPBNN" );
            mvaVar_mt_MZP800A0400_dm10TESUp = mt_MZP800A0400_reader->EvaluateMVA( "MLPBNN" );
            mvaVar_mt_MZP1000A0400_dm10TESUp = mt_MZP1000A0400_reader->EvaluateMVA( "MLPBNN" );
            mvaVar_mt_MZP1200A0400_dm10TESUp = mt_MZP1200A0400_reader->EvaluateMVA( "MLPBNN" );
            
            mvaVar_et_MZP600A0400_dm10TESUp = et_MZP600A0400_reader->EvaluateMVA( "MLPBNN" );
            mvaVar_et_MZP800A0400_dm10TESUp = et_MZP800A0400_reader->EvaluateMVA( "MLPBNN" );
            mvaVar_et_MZP1000A0400_dm10TESUp = et_MZP1000A0400_reader->EvaluateMVA("MLPBNN");
            mvaVar_et_MZP1200A0400_dm10TESUp = et_MZP1200A0400_reader->EvaluateMVA( "MLPBNN" );
            
            mvaVar_tt_MZP600A0400_dm10TESUp = tt_MZP600A0400_reader->EvaluateMVA( "MLPBNN" );
            mvaVar_tt_MZP800A0400_dm10TESUp = tt_MZP800A0400_reader->EvaluateMVA( "MLPBNN" );
            mvaVar_tt_MZP1000A0400_dm10TESUp = tt_MZP1000A0400_reader->EvaluateMVA( "MLPBNN" );
            mvaVar_tt_MZP1200A0400_dm10TESUp = tt_MZP1200A0400_reader->EvaluateMVA( "MLPBNN" );
            
            rpt_1 = pt_1_dm10TESDown;
            rpt_2 = pt_2_dm10TESDown;
            rpfmt_1 = pfmt_1_dm10TESDown;
            rmt_tot = mt_tot_dm10TESDown;
            rpt_tt = pt_tt_dm10TESDown;
            rm_vis = m_vis_dm10TESDown;
            rmet = pfmet_type1_dm10TESDown_Pt;
            rP_chi_pf = P_chi_pf_dm10TESDown;
            rLPT = LPT_dm10TESDown;
            rDeltaR_leg1_leg2 = DeltaR_leg1_leg2;
            rcos_DeltaPhi_leg1_leg2 = cos_DeltaPhi_leg1_leg2;
            rcos_DeltaPhi_PFMET_Higgs = cos_DeltaPhi_PFMET_Higgs_dm10TESDown;
            
            mvaVar_mt_MZP600A0400_dm10TESDown = mt_MZP600A0400_reader->EvaluateMVA( "MLPBNN" );
            mvaVar_mt_MZP800A0400_dm10TESDown = mt_MZP800A0400_reader->EvaluateMVA( "MLPBNN" );
            mvaVar_mt_MZP1000A0400_dm10TESDown = mt_MZP1000A0400_reader->EvaluateMVA( "MLPBNN" );
            mvaVar_mt_MZP1200A0400_dm10TESDown = mt_MZP1200A0400_reader->EvaluateMVA( "MLPBNN" );
            
            mvaVar_et_MZP600A0400_dm10TESDown = et_MZP600A0400_reader->EvaluateMVA( "MLPBNN" );
            mvaVar_et_MZP800A0400_dm10TESDown = et_MZP800A0400_reader->EvaluateMVA( "MLPBNN" );
            mvaVar_et_MZP1000A0400_dm10TESDown = et_MZP1000A0400_reader->EvaluateMVA("MLPBNN");
            mvaVar_et_MZP1200A0400_dm10TESDown = et_MZP1200A0400_reader->EvaluateMVA( "MLPBNN" );
            
            mvaVar_tt_MZP600A0400_dm10TESDown = tt_MZP600A0400_reader->EvaluateMVA( "MLPBNN" );
            mvaVar_tt_MZP800A0400_dm10TESDown = tt_MZP800A0400_reader->EvaluateMVA( "MLPBNN" );
            mvaVar_tt_MZP1000A0400_dm10TESDown = tt_MZP1000A0400_reader->EvaluateMVA( "MLPBNN" );
            mvaVar_tt_MZP1200A0400_dm10TESDown = tt_MZP1200A0400_reader->EvaluateMVA( "MLPBNN" );
        }
    }
    
    else if (R.getI("CandidateEventType")==6)
    {
    
        if (R.getI("leg1_decayMode")==0 && R.getI("leg2_decayMode")==0)
        {
        
            pt_1_dm0TESUp = pt_1_TESUp;
            mt_1_dm0TESUp = mt_1_TESUp;
            pfmt_1_dm0TESUp = pfmt_1_TESUp;

            pt_2_dm0TESUp = pt_2_TESUp;
            mt_2_dm0TESUp = mt_2_TESUp;
            pfmt_2_dm0TESUp = pfmt_2_TESUp;

            mt_tot_dm0TESUp = mt_tot_TESUp;
            pt_tt_dm0TESUp = pt_tt_TESUp;
            m_vis_dm0TESUp = m_vis_TESUp;
            pfmet_type1_dm0TESUp_Pt = pfmet_type1_TESUp_Pt;
            P_chi_pf_dm0TESUp = P_chi_pf_TESUp;
            LPT_dm0TESUp = LPT_TESUp;
            cos_DeltaPhi_PFMET_Higgs_dm0TESUp = cos_DeltaPhi_PFMET_Higgs_TESUp;

            pt_1_dm0TESDown = pt_1_TESDown;
            mt_1_dm0TESDown = mt_1_TESDown;
            pfmt_1_dm0TESDown = pfmt_1_TESDown;

            pt_2_dm0TESDown = pt_2_TESDown;
            mt_2_dm0TESDown = mt_2_TESDown;
            pfmt_2_dm0TESDown = pfmt_2_TESDown;

            mt_tot_dm0TESDown = mt_tot_TESDown;
            pt_tt_dm0TESDown = pt_tt_TESDown;
            m_vis_dm0TESDown = m_vis_TESDown;
            pfmet_type1_dm0TESDown_Pt = pfmet_type1_TESDown_Pt;
            P_chi_pf_dm0TESDown = P_chi_pf_TESDown;
            LPT_dm0TESDown = LPT_TESDown;
            cos_DeltaPhi_PFMET_Higgs_dm0TESDown = cos_DeltaPhi_PFMET_Higgs_TESDown;

            rpt_1 = pt_1_dm0TESUp;
            rpt_2 = pt_2_dm0TESUp;
            rpfmt_1 = pfmt_1_dm0TESUp;
            rmt_tot = mt_tot_dm0TESUp;
            rpt_tt = pt_tt_dm0TESUp;
            rm_vis = m_vis_dm0TESUp;
            rmet = pfmet_type1_dm0TESUp_Pt;
            rP_chi_pf = P_chi_pf_dm0TESUp;
            rLPT = LPT_dm0TESUp;
            rDeltaR_leg1_leg2 = DeltaR_leg1_leg2;
            rcos_DeltaPhi_leg1_leg2 = cos_DeltaPhi_leg1_leg2;
            rcos_DeltaPhi_PFMET_Higgs = cos_DeltaPhi_PFMET_Higgs_dm0TESUp;
            
            mvaVar_mt_MZP600A0400_dm0TESUp = mt_MZP600A0400_reader->EvaluateMVA( "MLPBNN" );
            mvaVar_mt_MZP800A0400_dm0TESUp = mt_MZP800A0400_reader->EvaluateMVA( "MLPBNN" );
            mvaVar_mt_MZP1000A0400_dm0TESUp = mt_MZP1000A0400_reader->EvaluateMVA( "MLPBNN" );
            mvaVar_mt_MZP1200A0400_dm0TESUp = mt_MZP1200A0400_reader->EvaluateMVA( "MLPBNN" );
            
            mvaVar_et_MZP600A0400_dm0TESUp = et_MZP600A0400_reader->EvaluateMVA( "MLPBNN" );
            mvaVar_et_MZP800A0400_dm0TESUp = et_MZP800A0400_reader->EvaluateMVA( "MLPBNN" );
            mvaVar_et_MZP1000A0400_dm0TESUp = et_MZP1000A0400_reader->EvaluateMVA("MLPBNN");
            mvaVar_et_MZP1200A0400_dm0TESUp = et_MZP1200A0400_reader->EvaluateMVA( "MLPBNN" );
            
            mvaVar_tt_MZP600A0400_dm0TESUp = tt_MZP600A0400_reader->EvaluateMVA( "MLPBNN" );
            mvaVar_tt_MZP800A0400_dm0TESUp = tt_MZP800A0400_reader->EvaluateMVA( "MLPBNN" );
            mvaVar_tt_MZP1000A0400_dm0TESUp = tt_MZP1000A0400_reader->EvaluateMVA( "MLPBNN" );
            mvaVar_tt_MZP1200A0400_dm0TESUp = tt_MZP1200A0400_reader->EvaluateMVA( "MLPBNN" );
            
            rpt_1 = pt_1_dm0TESDown;
            rpt_2 = pt_2_dm0TESDown;
            rpfmt_1 = pfmt_1_dm0TESDown;
            rmt_tot = mt_tot_dm0TESDown;
            rpt_tt = pt_tt_dm0TESDown;
            rm_vis = m_vis_dm0TESDown;
            rmet = pfmet_type1_dm0TESDown_Pt;
            rP_chi_pf = P_chi_pf_dm0TESDown;
            rLPT = LPT_dm0TESDown;
            rDeltaR_leg1_leg2 = DeltaR_leg1_leg2;
            rcos_DeltaPhi_leg1_leg2 = cos_DeltaPhi_leg1_leg2;
            rcos_DeltaPhi_PFMET_Higgs = cos_DeltaPhi_PFMET_Higgs_dm0TESDown;
            
            mvaVar_mt_MZP600A0400_dm0TESDown = mt_MZP600A0400_reader->EvaluateMVA( "MLPBNN" );
            mvaVar_mt_MZP800A0400_dm0TESDown = mt_MZP800A0400_reader->EvaluateMVA( "MLPBNN" );
            mvaVar_mt_MZP1000A0400_dm0TESDown = mt_MZP1000A0400_reader->EvaluateMVA( "MLPBNN" );
            mvaVar_mt_MZP1200A0400_dm0TESDown = mt_MZP1200A0400_reader->EvaluateMVA( "MLPBNN" );
            
            mvaVar_et_MZP600A0400_dm0TESDown = et_MZP600A0400_reader->EvaluateMVA( "MLPBNN" );
            mvaVar_et_MZP800A0400_dm0TESDown = et_MZP800A0400_reader->EvaluateMVA( "MLPBNN" );
            mvaVar_et_MZP1000A0400_dm0TESDown = et_MZP1000A0400_reader->EvaluateMVA("MLPBNN");
            mvaVar_et_MZP1200A0400_dm0TESDown = et_MZP1200A0400_reader->EvaluateMVA( "MLPBNN" );
            
            mvaVar_tt_MZP600A0400_dm0TESDown = tt_MZP600A0400_reader->EvaluateMVA( "MLPBNN" );
            mvaVar_tt_MZP800A0400_dm0TESDown = tt_MZP800A0400_reader->EvaluateMVA( "MLPBNN" );
            mvaVar_tt_MZP1000A0400_dm0TESDown = tt_MZP1000A0400_reader->EvaluateMVA( "MLPBNN" );
            mvaVar_tt_MZP1200A0400_dm0TESDown = tt_MZP1200A0400_reader->EvaluateMVA( "MLPBNN" );
            
        }
        if (R.getI("leg1_decayMode")==0 && R.getI("leg2_decayMode")!=0)
        {
        
            pt_1_dm0TESUp = pt_1_L1TESUp;
            mt_1_dm0TESUp = mt_1_L1TESUp;
            pfmt_1_dm0TESUp = pfmt_1_L1TESUp;

            pt_2_dm0TESUp = pt_2_L1TESUp;
            mt_2_dm0TESUp = mt_2_L1TESUp;
            pfmt_2_dm0TESUp = pfmt_2_L1TESUp;

            mt_tot_dm0TESUp = mt_tot_L1TESUp;
            pt_tt_dm0TESUp = pt_tt_L1TESUp;
            m_vis_dm0TESUp = m_vis_L1TESUp;
            pfmet_type1_dm0TESUp_Pt = pfmet_type1_L1TESUp_Pt;
            P_chi_pf_dm0TESUp = P_chi_pf_L1TESUp;
            LPT_dm0TESUp = LPT_L1TESUp;
            cos_DeltaPhi_PFMET_Higgs_dm0TESUp = cos_DeltaPhi_PFMET_Higgs_L1TESUp;

            pt_1_dm0TESDown = pt_1_L1TESDown;
            mt_1_dm0TESDown = mt_1_L1TESDown;
            pfmt_1_dm0TESDown = pfmt_1_L1TESDown;

            pt_2_dm0TESDown = pt_2_L1TESDown;
            mt_2_dm0TESDown = mt_2_L1TESDown;
            pfmt_2_dm0TESDown = pfmt_2_L1TESDown;

            mt_tot_dm0TESDown = mt_tot_L1TESDown;
            pt_tt_dm0TESDown = pt_tt_L1TESDown;
            m_vis_dm0TESDown = m_vis_L1TESDown;
            pfmet_type1_dm0TESDown_Pt = pfmet_type1_L1TESDown_Pt;
            P_chi_pf_dm0TESDown = P_chi_pf_L1TESDown;
            LPT_dm0TESDown = LPT_L1TESDown;
            cos_DeltaPhi_PFMET_Higgs_dm0TESDown = cos_DeltaPhi_PFMET_Higgs_L1TESDown;

            rpt_1 = pt_1_L1TESUp;
            rpt_2 = pt_2_L1TESUp;
            rpfmt_1 = pfmt_1_L1TESUp;
            rmt_tot = mt_tot_L1TESUp;
            rpt_tt = pt_tt_L1TESUp;
            rm_vis = m_vis_L1TESUp;
            rmet = pfmet_type1_L1TESUp_Pt;
            rP_chi_pf = P_chi_pf_L1TESUp;
            rLPT = LPT_L1TESUp;
            rDeltaR_leg1_leg2 = DeltaR_leg1_leg2;
            rcos_DeltaPhi_leg1_leg2 = cos_DeltaPhi_leg1_leg2;
            rcos_DeltaPhi_PFMET_Higgs = cos_DeltaPhi_PFMET_Higgs_L1TESUp;
            
            mvaVar_mt_MZP600A0400_dm0TESUp = mt_MZP600A0400_reader->EvaluateMVA( "MLPBNN" );
            mvaVar_mt_MZP800A0400_dm0TESUp = mt_MZP800A0400_reader->EvaluateMVA( "MLPBNN" );
            mvaVar_mt_MZP1000A0400_dm0TESUp = mt_MZP1000A0400_reader->EvaluateMVA( "MLPBNN" );
            mvaVar_mt_MZP1200A0400_dm0TESUp = mt_MZP1200A0400_reader->EvaluateMVA( "MLPBNN" );
            
            mvaVar_et_MZP600A0400_dm0TESUp = et_MZP600A0400_reader->EvaluateMVA( "MLPBNN" );
            mvaVar_et_MZP800A0400_dm0TESUp = et_MZP800A0400_reader->EvaluateMVA( "MLPBNN" );
            mvaVar_et_MZP1000A0400_dm0TESUp = et_MZP1000A0400_reader->EvaluateMVA("MLPBNN");
            mvaVar_et_MZP1200A0400_dm0TESUp = et_MZP1200A0400_reader->EvaluateMVA( "MLPBNN" );
            
            mvaVar_tt_MZP600A0400_dm0TESUp = tt_MZP600A0400_reader->EvaluateMVA( "MLPBNN" );
            mvaVar_tt_MZP800A0400_dm0TESUp = tt_MZP800A0400_reader->EvaluateMVA( "MLPBNN" );
            mvaVar_tt_MZP1000A0400_dm0TESUp = tt_MZP1000A0400_reader->EvaluateMVA( "MLPBNN" );
            mvaVar_tt_MZP1200A0400_dm0TESUp = tt_MZP1200A0400_reader->EvaluateMVA( "MLPBNN" );
            
            rpt_1 = pt_1_L1TESDown;
            rpt_2 = pt_2_L1TESDown;
            rpfmt_1 = pfmt_1_L1TESDown;
            rmt_tot = mt_tot_L1TESDown;
            rpt_tt = pt_tt_L1TESDown;
            rm_vis = m_vis_L1TESDown;
            rmet = pfmet_type1_L1TESDown_Pt;
            rP_chi_pf = P_chi_pf_L1TESDown;
            rLPT = LPT_L1TESDown;
            rDeltaR_leg1_leg2 = DeltaR_leg1_leg2;
            rcos_DeltaPhi_leg1_leg2 = cos_DeltaPhi_leg1_leg2;
            rcos_DeltaPhi_PFMET_Higgs = cos_DeltaPhi_PFMET_Higgs_L1TESDown;
            
            mvaVar_mt_MZP600A0400_dm0TESDown = mt_MZP600A0400_reader->EvaluateMVA( "MLPBNN" );
            mvaVar_mt_MZP800A0400_dm0TESDown = mt_MZP800A0400_reader->EvaluateMVA( "MLPBNN" );
            mvaVar_mt_MZP1000A0400_dm0TESDown = mt_MZP1000A0400_reader->EvaluateMVA( "MLPBNN" );
            mvaVar_mt_MZP1200A0400_dm0TESDown = mt_MZP1200A0400_reader->EvaluateMVA( "MLPBNN" );
            
            mvaVar_et_MZP600A0400_dm0TESDown = et_MZP600A0400_reader->EvaluateMVA( "MLPBNN" );
            mvaVar_et_MZP800A0400_dm0TESDown = et_MZP800A0400_reader->EvaluateMVA( "MLPBNN" );
            mvaVar_et_MZP1000A0400_dm0TESDown = et_MZP1000A0400_reader->EvaluateMVA("MLPBNN");
            mvaVar_et_MZP1200A0400_dm0TESDown = et_MZP1200A0400_reader->EvaluateMVA( "MLPBNN" );
            
            mvaVar_tt_MZP600A0400_dm0TESDown = tt_MZP600A0400_reader->EvaluateMVA( "MLPBNN" );
            mvaVar_tt_MZP800A0400_dm0TESDown = tt_MZP800A0400_reader->EvaluateMVA( "MLPBNN" );
            mvaVar_tt_MZP1000A0400_dm0TESDown = tt_MZP1000A0400_reader->EvaluateMVA( "MLPBNN" );
            mvaVar_tt_MZP1200A0400_dm0TESDown = tt_MZP1200A0400_reader->EvaluateMVA( "MLPBNN" );
        }
        if (R.getI("leg1_decayMode")!=0 && R.getI("leg2_decayMode")==0)
        {
        
            pt_1_dm0TESUp = pt_1_L2TESUp;
            mt_1_dm0TESUp = mt_1_L2TESUp;
            pfmt_1_dm0TESUp = pfmt_1_L2TESUp;

            pt_2_dm0TESUp = pt_2_L2TESUp;
            mt_2_dm0TESUp = mt_2_L2TESUp;
            pfmt_2_dm0TESUp = pfmt_2_L2TESUp;

            mt_tot_dm0TESUp = mt_tot_L2TESUp;
            pt_tt_dm0TESUp = pt_tt_L2TESUp;
            m_vis_dm0TESUp = m_vis_L2TESUp;
            pfmet_type1_dm0TESUp_Pt = pfmet_type1_L2TESUp_Pt;
            P_chi_pf_dm0TESUp = P_chi_pf_L2TESUp;
            LPT_dm0TESUp = LPT_L2TESUp;
            cos_DeltaPhi_PFMET_Higgs_dm0TESUp = cos_DeltaPhi_PFMET_Higgs_L2TESUp;

            pt_1_dm0TESDown = pt_1_L2TESDown;
            mt_1_dm0TESDown = mt_1_L2TESDown;
            pfmt_1_dm0TESDown = pfmt_1_L2TESDown;

            pt_2_dm0TESDown = pt_2_L2TESDown;
            mt_2_dm0TESDown = mt_2_L2TESDown;
            pfmt_2_dm0TESDown = pfmt_2_L2TESDown;

            mt_tot_dm0TESDown = mt_tot_L2TESDown;
            pt_tt_dm0TESDown = pt_tt_L2TESDown;
            m_vis_dm0TESDown = m_vis_L2TESDown;
            pfmet_type1_dm0TESDown_Pt = pfmet_type1_L2TESDown_Pt;
            P_chi_pf_dm0TESDown = P_chi_pf_L2TESDown;
            LPT_dm0TESDown = LPT_L2TESDown;
            cos_DeltaPhi_PFMET_Higgs_dm0TESDown = cos_DeltaPhi_PFMET_Higgs_L2TESDown;

            rpt_1 = pt_1_L2TESUp;
            rpt_2 = pt_2_L2TESUp;
            rpfmt_1 = pfmt_1_L2TESUp;
            rmt_tot = mt_tot_L2TESUp;
            rpt_tt = pt_tt_L2TESUp;
            rm_vis = m_vis_L2TESUp;
            rmet = pfmet_type1_L2TESUp_Pt;
            rP_chi_pf = P_chi_pf_L2TESUp;
            rLPT = LPT_L2TESUp;
            rDeltaR_leg1_leg2 = DeltaR_leg1_leg2;
            rcos_DeltaPhi_leg1_leg2 = cos_DeltaPhi_leg1_leg2;
            rcos_DeltaPhi_PFMET_Higgs = cos_DeltaPhi_PFMET_Higgs_L2TESUp;
            
            mvaVar_mt_MZP600A0400_dm0TESUp = mt_MZP600A0400_reader->EvaluateMVA( "MLPBNN" );
            mvaVar_mt_MZP800A0400_dm0TESUp = mt_MZP800A0400_reader->EvaluateMVA( "MLPBNN" );
            mvaVar_mt_MZP1000A0400_dm0TESUp = mt_MZP1000A0400_reader->EvaluateMVA( "MLPBNN" );
            mvaVar_mt_MZP1200A0400_dm0TESUp = mt_MZP1200A0400_reader->EvaluateMVA( "MLPBNN" );
            
            mvaVar_et_MZP600A0400_dm0TESUp = et_MZP600A0400_reader->EvaluateMVA( "MLPBNN" );
            mvaVar_et_MZP800A0400_dm0TESUp = et_MZP800A0400_reader->EvaluateMVA( "MLPBNN" );
            mvaVar_et_MZP1000A0400_dm0TESUp = et_MZP1000A0400_reader->EvaluateMVA("MLPBNN");
            mvaVar_et_MZP1200A0400_dm0TESUp = et_MZP1200A0400_reader->EvaluateMVA( "MLPBNN" );
            
            mvaVar_tt_MZP600A0400_dm0TESUp = tt_MZP600A0400_reader->EvaluateMVA( "MLPBNN" );
            mvaVar_tt_MZP800A0400_dm0TESUp = tt_MZP800A0400_reader->EvaluateMVA( "MLPBNN" );
            mvaVar_tt_MZP1000A0400_dm0TESUp = tt_MZP1000A0400_reader->EvaluateMVA( "MLPBNN" );
            mvaVar_tt_MZP1200A0400_dm0TESUp = tt_MZP1200A0400_reader->EvaluateMVA( "MLPBNN" );
            
            rpt_1 = pt_1_L2TESDown;
            rpt_2 = pt_2_L2TESDown;
            rpfmt_1 = pfmt_1_L2TESDown;
            rmt_tot = mt_tot_L2TESDown;
            rpt_tt = pt_tt_L2TESDown;
            rm_vis = m_vis_L2TESDown;
            rmet = pfmet_type1_L2TESDown_Pt;
            rP_chi_pf = P_chi_pf_L2TESDown;
            rLPT = LPT_L2TESDown;
            rDeltaR_leg1_leg2 = DeltaR_leg1_leg2;
            rcos_DeltaPhi_leg1_leg2 = cos_DeltaPhi_leg1_leg2;
            rcos_DeltaPhi_PFMET_Higgs = cos_DeltaPhi_PFMET_Higgs_L2TESDown;
            
            mvaVar_mt_MZP600A0400_dm0TESDown = mt_MZP600A0400_reader->EvaluateMVA( "MLPBNN" );
            mvaVar_mt_MZP800A0400_dm0TESDown = mt_MZP800A0400_reader->EvaluateMVA( "MLPBNN" );
            mvaVar_mt_MZP1000A0400_dm0TESDown = mt_MZP1000A0400_reader->EvaluateMVA( "MLPBNN" );
            mvaVar_mt_MZP1200A0400_dm0TESDown = mt_MZP1200A0400_reader->EvaluateMVA( "MLPBNN" );
            
            mvaVar_et_MZP600A0400_dm0TESDown = et_MZP600A0400_reader->EvaluateMVA( "MLPBNN" );
            mvaVar_et_MZP800A0400_dm0TESDown = et_MZP800A0400_reader->EvaluateMVA( "MLPBNN" );
            mvaVar_et_MZP1000A0400_dm0TESDown = et_MZP1000A0400_reader->EvaluateMVA("MLPBNN");
            mvaVar_et_MZP1200A0400_dm0TESDown = et_MZP1200A0400_reader->EvaluateMVA( "MLPBNN" );
            
            mvaVar_tt_MZP600A0400_dm0TESDown = tt_MZP600A0400_reader->EvaluateMVA( "MLPBNN" );
            mvaVar_tt_MZP800A0400_dm0TESDown = tt_MZP800A0400_reader->EvaluateMVA( "MLPBNN" );
            mvaVar_tt_MZP1000A0400_dm0TESDown = tt_MZP1000A0400_reader->EvaluateMVA( "MLPBNN" );
            mvaVar_tt_MZP1200A0400_dm0TESDown = tt_MZP1200A0400_reader->EvaluateMVA( "MLPBNN" );
        }
        
        if (R.getI("leg1_decayMode")==1 && R.getI("leg2_decayMode")==1)
        {
        
            pt_1_dm1TESUp = pt_1_TESUp;
            mt_1_dm1TESUp = mt_1_TESUp;
            pfmt_1_dm1TESUp = pfmt_1_TESUp;

            pt_2_dm1TESUp = pt_2_TESUp;
            mt_2_dm1TESUp = mt_2_TESUp;
            pfmt_2_dm1TESUp = pfmt_2_TESUp;

            mt_tot_dm1TESUp = mt_tot_TESUp;
            pt_tt_dm1TESUp = pt_tt_TESUp;
            m_vis_dm1TESUp = m_vis_TESUp;
            pfmet_type1_dm1TESUp_Pt = pfmet_type1_TESUp_Pt;
            P_chi_pf_dm1TESUp = P_chi_pf_TESUp;
            LPT_dm1TESUp = LPT_TESUp;
            cos_DeltaPhi_PFMET_Higgs_dm1TESUp = cos_DeltaPhi_PFMET_Higgs_TESUp;

            pt_1_dm1TESDown = pt_1_TESDown;
            mt_1_dm1TESDown = mt_1_TESDown;
            pfmt_1_dm1TESDown = pfmt_1_TESDown;

            pt_2_dm1TESDown = pt_2_TESDown;
            mt_2_dm1TESDown = mt_2_TESDown;
            pfmt_2_dm1TESDown = pfmt_2_TESDown;

            mt_tot_dm1TESDown = mt_tot_TESDown;
            pt_tt_dm1TESDown = pt_tt_TESDown;
            m_vis_dm1TESDown = m_vis_TESDown;
            pfmet_type1_dm1TESDown_Pt = pfmet_type1_TESDown_Pt;
            P_chi_pf_dm1TESDown = P_chi_pf_TESDown;
            LPT_dm1TESDown = LPT_TESDown;
            cos_DeltaPhi_PFMET_Higgs_dm1TESDown = cos_DeltaPhi_PFMET_Higgs_TESDown;

            rpt_1 = pt_1_dm1TESUp;
            rpt_2 = pt_2_dm1TESUp;
            rpfmt_1 = pfmt_1_dm1TESUp;
            rmt_tot = mt_tot_dm1TESUp;
            rpt_tt = pt_tt_dm1TESUp;
            rm_vis = m_vis_dm1TESUp;
            rmet = pfmet_type1_dm1TESUp_Pt;
            rP_chi_pf = P_chi_pf_dm1TESUp;
            rLPT = LPT_dm1TESUp;
            rDeltaR_leg1_leg2 = DeltaR_leg1_leg2;
            rcos_DeltaPhi_leg1_leg2 = cos_DeltaPhi_leg1_leg2;
            rcos_DeltaPhi_PFMET_Higgs = cos_DeltaPhi_PFMET_Higgs_dm1TESUp;
            
            mvaVar_mt_MZP600A0400_dm1TESUp = mt_MZP600A0400_reader->EvaluateMVA( "MLPBNN" );
            mvaVar_mt_MZP800A0400_dm1TESUp = mt_MZP800A0400_reader->EvaluateMVA( "MLPBNN" );
            mvaVar_mt_MZP1000A0400_dm1TESUp = mt_MZP1000A0400_reader->EvaluateMVA( "MLPBNN" );
            mvaVar_mt_MZP1200A0400_dm1TESUp = mt_MZP1200A0400_reader->EvaluateMVA( "MLPBNN" );
            
            mvaVar_et_MZP600A0400_dm1TESUp = et_MZP600A0400_reader->EvaluateMVA( "MLPBNN" );
            mvaVar_et_MZP800A0400_dm1TESUp = et_MZP800A0400_reader->EvaluateMVA( "MLPBNN" );
            mvaVar_et_MZP1000A0400_dm1TESUp = et_MZP1000A0400_reader->EvaluateMVA("MLPBNN");
            mvaVar_et_MZP1200A0400_dm1TESUp = et_MZP1200A0400_reader->EvaluateMVA( "MLPBNN" );
            
            mvaVar_tt_MZP600A0400_dm1TESUp = tt_MZP600A0400_reader->EvaluateMVA( "MLPBNN" );
            mvaVar_tt_MZP800A0400_dm1TESUp = tt_MZP800A0400_reader->EvaluateMVA( "MLPBNN" );
            mvaVar_tt_MZP1000A0400_dm1TESUp = tt_MZP1000A0400_reader->EvaluateMVA( "MLPBNN" );
            mvaVar_tt_MZP1200A0400_dm1TESUp = tt_MZP1200A0400_reader->EvaluateMVA( "MLPBNN" );
            
            rpt_1 = pt_1_dm1TESDown;
            rpt_2 = pt_2_dm1TESDown;
            rpfmt_1 = pfmt_1_dm1TESDown;
            rmt_tot = mt_tot_dm1TESDown;
            rpt_tt = pt_tt_dm1TESDown;
            rm_vis = m_vis_dm1TESDown;
            rmet = pfmet_type1_dm1TESDown_Pt;
            rP_chi_pf = P_chi_pf_dm1TESDown;
            rLPT = LPT_dm1TESDown;
            rDeltaR_leg1_leg2 = DeltaR_leg1_leg2;
            rcos_DeltaPhi_leg1_leg2 = cos_DeltaPhi_leg1_leg2;
            rcos_DeltaPhi_PFMET_Higgs = cos_DeltaPhi_PFMET_Higgs_dm1TESDown;
            
            mvaVar_mt_MZP600A0400_dm1TESDown = mt_MZP600A0400_reader->EvaluateMVA( "MLPBNN" );
            mvaVar_mt_MZP800A0400_dm1TESDown = mt_MZP800A0400_reader->EvaluateMVA( "MLPBNN" );
            mvaVar_mt_MZP1000A0400_dm1TESDown = mt_MZP1000A0400_reader->EvaluateMVA( "MLPBNN" );
            mvaVar_mt_MZP1200A0400_dm1TESDown = mt_MZP1200A0400_reader->EvaluateMVA( "MLPBNN" );
            
            mvaVar_et_MZP600A0400_dm1TESDown = et_MZP600A0400_reader->EvaluateMVA( "MLPBNN" );
            mvaVar_et_MZP800A0400_dm1TESDown = et_MZP800A0400_reader->EvaluateMVA( "MLPBNN" );
            mvaVar_et_MZP1000A0400_dm1TESDown = et_MZP1000A0400_reader->EvaluateMVA("MLPBNN");
            mvaVar_et_MZP1200A0400_dm1TESDown = et_MZP1200A0400_reader->EvaluateMVA( "MLPBNN" );
            
            mvaVar_tt_MZP600A0400_dm1TESDown = tt_MZP600A0400_reader->EvaluateMVA( "MLPBNN" );
            mvaVar_tt_MZP800A0400_dm1TESDown = tt_MZP800A0400_reader->EvaluateMVA( "MLPBNN" );
            mvaVar_tt_MZP1000A0400_dm1TESDown = tt_MZP1000A0400_reader->EvaluateMVA( "MLPBNN" );
            mvaVar_tt_MZP1200A0400_dm1TESDown = tt_MZP1200A0400_reader->EvaluateMVA( "MLPBNN" );

        }
        if (R.getI("leg1_decayMode")==1 && R.getI("leg2_decayMode")!=1)
        {
        
            pt_1_dm1TESUp = pt_1_L1TESUp;
            mt_1_dm1TESUp = mt_1_L1TESUp;
            pfmt_1_dm1TESUp = pfmt_1_L1TESUp;

            pt_2_dm1TESUp = pt_2_L1TESUp;
            mt_2_dm1TESUp = mt_2_L1TESUp;
            pfmt_2_dm1TESUp = pfmt_2_L1TESUp;

            mt_tot_dm1TESUp = mt_tot_L1TESUp;
            pt_tt_dm1TESUp = pt_tt_L1TESUp;
            m_vis_dm1TESUp = m_vis_L1TESUp;
            pfmet_type1_dm1TESUp_Pt = pfmet_type1_L1TESUp_Pt;
            P_chi_pf_dm1TESUp = P_chi_pf_L1TESUp;
            LPT_dm1TESUp = LPT_L1TESUp;
            cos_DeltaPhi_PFMET_Higgs_dm1TESUp = cos_DeltaPhi_PFMET_Higgs_L1TESUp;

            pt_1_dm1TESDown = pt_1_L1TESDown;
            mt_1_dm1TESDown = mt_1_L1TESDown;
            pfmt_1_dm1TESDown = pfmt_1_L1TESDown;

            pt_2_dm1TESDown = pt_2_L1TESDown;
            mt_2_dm1TESDown = mt_2_L1TESDown;
            pfmt_2_dm1TESDown = pfmt_2_L1TESDown;

            mt_tot_dm1TESDown = mt_tot_L1TESDown;
            pt_tt_dm1TESDown = pt_tt_L1TESDown;
            m_vis_dm1TESDown = m_vis_L1TESDown;
            pfmet_type1_dm1TESDown_Pt = pfmet_type1_L1TESDown_Pt;
            P_chi_pf_dm1TESDown = P_chi_pf_L1TESDown;
            LPT_dm1TESDown = LPT_L1TESDown;
            cos_DeltaPhi_PFMET_Higgs_dm1TESDown = cos_DeltaPhi_PFMET_Higgs_L1TESDown;

            rpt_1 = pt_1_L1TESUp;
            rpt_2 = pt_2_L1TESUp;
            rpfmt_1 = pfmt_1_L1TESUp;
            rmt_tot = mt_tot_L1TESUp;
            rpt_tt = pt_tt_L1TESUp;
            rm_vis = m_vis_L1TESUp;
            rmet = pfmet_type1_L1TESUp_Pt;
            rP_chi_pf = P_chi_pf_L1TESUp;
            rLPT = LPT_L1TESUp;
            rDeltaR_leg1_leg2 = DeltaR_leg1_leg2;
            rcos_DeltaPhi_leg1_leg2 = cos_DeltaPhi_leg1_leg2;
            rcos_DeltaPhi_PFMET_Higgs = cos_DeltaPhi_PFMET_Higgs_L1TESUp;
            
            mvaVar_mt_MZP600A0400_dm1TESUp = mt_MZP600A0400_reader->EvaluateMVA( "MLPBNN" );
            mvaVar_mt_MZP800A0400_dm1TESUp = mt_MZP800A0400_reader->EvaluateMVA( "MLPBNN" );
            mvaVar_mt_MZP1000A0400_dm1TESUp = mt_MZP1000A0400_reader->EvaluateMVA( "MLPBNN" );
            mvaVar_mt_MZP1200A0400_dm1TESUp = mt_MZP1200A0400_reader->EvaluateMVA( "MLPBNN" );
            
            mvaVar_et_MZP600A0400_dm1TESUp = et_MZP600A0400_reader->EvaluateMVA( "MLPBNN" );
            mvaVar_et_MZP800A0400_dm1TESUp = et_MZP800A0400_reader->EvaluateMVA( "MLPBNN" );
            mvaVar_et_MZP1000A0400_dm1TESUp = et_MZP1000A0400_reader->EvaluateMVA("MLPBNN");
            mvaVar_et_MZP1200A0400_dm1TESUp = et_MZP1200A0400_reader->EvaluateMVA( "MLPBNN" );
            
            mvaVar_tt_MZP600A0400_dm1TESUp = tt_MZP600A0400_reader->EvaluateMVA( "MLPBNN" );
            mvaVar_tt_MZP800A0400_dm1TESUp = tt_MZP800A0400_reader->EvaluateMVA( "MLPBNN" );
            mvaVar_tt_MZP1000A0400_dm1TESUp = tt_MZP1000A0400_reader->EvaluateMVA( "MLPBNN" );
            mvaVar_tt_MZP1200A0400_dm1TESUp = tt_MZP1200A0400_reader->EvaluateMVA( "MLPBNN" );
            
            rpt_1 = pt_1_L1TESDown;
            rpt_2 = pt_2_L1TESDown;
            rpfmt_1 = pfmt_1_L1TESDown;
            rmt_tot = mt_tot_L1TESDown;
            rpt_tt = pt_tt_L1TESDown;
            rm_vis = m_vis_L1TESDown;
            rmet = pfmet_type1_L1TESDown_Pt;
            rP_chi_pf = P_chi_pf_L1TESDown;
            rLPT = LPT_L1TESDown;
            rDeltaR_leg1_leg2 = DeltaR_leg1_leg2;
            rcos_DeltaPhi_leg1_leg2 = cos_DeltaPhi_leg1_leg2;
            rcos_DeltaPhi_PFMET_Higgs = cos_DeltaPhi_PFMET_Higgs_L1TESDown;
            
            mvaVar_mt_MZP600A0400_dm1TESDown = mt_MZP600A0400_reader->EvaluateMVA( "MLPBNN" );
            mvaVar_mt_MZP800A0400_dm1TESDown = mt_MZP800A0400_reader->EvaluateMVA( "MLPBNN" );
            mvaVar_mt_MZP1000A0400_dm1TESDown = mt_MZP1000A0400_reader->EvaluateMVA( "MLPBNN" );
            mvaVar_mt_MZP1200A0400_dm1TESDown = mt_MZP1200A0400_reader->EvaluateMVA( "MLPBNN" );
            
            mvaVar_et_MZP600A0400_dm1TESDown = et_MZP600A0400_reader->EvaluateMVA( "MLPBNN" );
            mvaVar_et_MZP800A0400_dm1TESDown = et_MZP800A0400_reader->EvaluateMVA( "MLPBNN" );
            mvaVar_et_MZP1000A0400_dm1TESDown = et_MZP1000A0400_reader->EvaluateMVA("MLPBNN");
            mvaVar_et_MZP1200A0400_dm1TESDown = et_MZP1200A0400_reader->EvaluateMVA( "MLPBNN" );
            
            mvaVar_tt_MZP600A0400_dm1TESDown = tt_MZP600A0400_reader->EvaluateMVA( "MLPBNN" );
            mvaVar_tt_MZP800A0400_dm1TESDown = tt_MZP800A0400_reader->EvaluateMVA( "MLPBNN" );
            mvaVar_tt_MZP1000A0400_dm1TESDown = tt_MZP1000A0400_reader->EvaluateMVA( "MLPBNN" );
            mvaVar_tt_MZP1200A0400_dm1TESDown = tt_MZP1200A0400_reader->EvaluateMVA( "MLPBNN" );
            
        }
        if (R.getI("leg1_decayMode")!=1 && R.getI("leg2_decayMode")==1)
        {
        
            pt_1_dm1TESUp = pt_1_L2TESUp;
            mt_1_dm1TESUp = mt_1_L2TESUp;
            pfmt_1_dm1TESUp = pfmt_1_L2TESUp;

            pt_2_dm1TESUp = pt_2_L2TESUp;
            mt_2_dm1TESUp = mt_2_L2TESUp;
            pfmt_2_dm1TESUp = pfmt_2_L2TESUp;

            mt_tot_dm1TESUp = mt_tot_L2TESUp;
            pt_tt_dm1TESUp = pt_tt_L2TESUp;
            m_vis_dm1TESUp = m_vis_L2TESUp;
            pfmet_type1_dm1TESUp_Pt = pfmet_type1_L2TESUp_Pt;
            P_chi_pf_dm1TESUp = P_chi_pf_L2TESUp;
            LPT_dm1TESUp = LPT_L2TESUp;
            cos_DeltaPhi_PFMET_Higgs_dm1TESUp = cos_DeltaPhi_PFMET_Higgs_L2TESUp;

            pt_1_dm1TESDown = pt_1_L2TESDown;
            mt_1_dm1TESDown = mt_1_L2TESDown;
            pfmt_1_dm1TESDown = pfmt_1_L2TESDown;

            pt_2_dm1TESDown = pt_2_L2TESDown;
            mt_2_dm1TESDown = mt_2_L2TESDown;
            pfmt_2_dm1TESDown = pfmt_2_L2TESDown;

            mt_tot_dm1TESDown = mt_tot_L2TESDown;
            pt_tt_dm1TESDown = pt_tt_L2TESDown;
            m_vis_dm1TESDown = m_vis_L2TESDown;
            pfmet_type1_dm1TESDown_Pt = pfmet_type1_L2TESDown_Pt;
            P_chi_pf_dm1TESDown = P_chi_pf_L2TESDown;
            LPT_dm1TESDown = LPT_L2TESDown;
            cos_DeltaPhi_PFMET_Higgs_dm1TESDown = cos_DeltaPhi_PFMET_Higgs_L2TESDown;

            rpt_1 = pt_1_L2TESUp;
            rpt_2 = pt_2_L2TESUp;
            rpfmt_1 = pfmt_1_L2TESUp;
            rmt_tot = mt_tot_L2TESUp;
            rpt_tt = pt_tt_L2TESUp;
            rm_vis = m_vis_L2TESUp;
            rmet = pfmet_type1_L2TESUp_Pt;
            rP_chi_pf = P_chi_pf_L2TESUp;
            rLPT = LPT_L2TESUp;
            rDeltaR_leg1_leg2 = DeltaR_leg1_leg2;
            rcos_DeltaPhi_leg1_leg2 = cos_DeltaPhi_leg1_leg2;
            rcos_DeltaPhi_PFMET_Higgs = cos_DeltaPhi_PFMET_Higgs_L2TESUp;
            
            mvaVar_mt_MZP600A0400_dm1TESUp = mt_MZP600A0400_reader->EvaluateMVA( "MLPBNN" );
            mvaVar_mt_MZP800A0400_dm1TESUp = mt_MZP800A0400_reader->EvaluateMVA( "MLPBNN" );
            mvaVar_mt_MZP1000A0400_dm1TESUp = mt_MZP1000A0400_reader->EvaluateMVA( "MLPBNN" );
            mvaVar_mt_MZP1200A0400_dm1TESUp = mt_MZP1200A0400_reader->EvaluateMVA( "MLPBNN" );
            
            mvaVar_et_MZP600A0400_dm1TESUp = et_MZP600A0400_reader->EvaluateMVA( "MLPBNN" );
            mvaVar_et_MZP800A0400_dm1TESUp = et_MZP800A0400_reader->EvaluateMVA( "MLPBNN" );
            mvaVar_et_MZP1000A0400_dm1TESUp = et_MZP1000A0400_reader->EvaluateMVA("MLPBNN");
            mvaVar_et_MZP1200A0400_dm1TESUp = et_MZP1200A0400_reader->EvaluateMVA( "MLPBNN" );
            
            mvaVar_tt_MZP600A0400_dm1TESUp = tt_MZP600A0400_reader->EvaluateMVA( "MLPBNN" );
            mvaVar_tt_MZP800A0400_dm1TESUp = tt_MZP800A0400_reader->EvaluateMVA( "MLPBNN" );
            mvaVar_tt_MZP1000A0400_dm1TESUp = tt_MZP1000A0400_reader->EvaluateMVA( "MLPBNN" );
            mvaVar_tt_MZP1200A0400_dm1TESUp = tt_MZP1200A0400_reader->EvaluateMVA( "MLPBNN" );
            
            rpt_1 = pt_1_L2TESDown;
            rpt_2 = pt_2_L2TESDown;
            rpfmt_1 = pfmt_1_L2TESDown;
            rmt_tot = mt_tot_L2TESDown;
            rpt_tt = pt_tt_L2TESDown;
            rm_vis = m_vis_L2TESDown;
            rmet = pfmet_type1_L2TESDown_Pt;
            rP_chi_pf = P_chi_pf_L2TESDown;
            rLPT = LPT_L2TESDown;
            rDeltaR_leg1_leg2 = DeltaR_leg1_leg2;
            rcos_DeltaPhi_leg1_leg2 = cos_DeltaPhi_leg1_leg2;
            rcos_DeltaPhi_PFMET_Higgs = cos_DeltaPhi_PFMET_Higgs_L2TESDown;
            
            mvaVar_mt_MZP600A0400_dm1TESDown = mt_MZP600A0400_reader->EvaluateMVA( "MLPBNN" );
            mvaVar_mt_MZP800A0400_dm1TESDown = mt_MZP800A0400_reader->EvaluateMVA( "MLPBNN" );
            mvaVar_mt_MZP1000A0400_dm1TESDown = mt_MZP1000A0400_reader->EvaluateMVA( "MLPBNN" );
            mvaVar_mt_MZP1200A0400_dm1TESDown = mt_MZP1200A0400_reader->EvaluateMVA( "MLPBNN" );
            
            mvaVar_et_MZP600A0400_dm1TESDown = et_MZP600A0400_reader->EvaluateMVA( "MLPBNN" );
            mvaVar_et_MZP800A0400_dm1TESDown = et_MZP800A0400_reader->EvaluateMVA( "MLPBNN" );
            mvaVar_et_MZP1000A0400_dm1TESDown = et_MZP1000A0400_reader->EvaluateMVA("MLPBNN");
            mvaVar_et_MZP1200A0400_dm1TESDown = et_MZP1200A0400_reader->EvaluateMVA( "MLPBNN" );
            
            mvaVar_tt_MZP600A0400_dm1TESDown = tt_MZP600A0400_reader->EvaluateMVA( "MLPBNN" );
            mvaVar_tt_MZP800A0400_dm1TESDown = tt_MZP800A0400_reader->EvaluateMVA( "MLPBNN" );
            mvaVar_tt_MZP1000A0400_dm1TESDown = tt_MZP1000A0400_reader->EvaluateMVA( "MLPBNN" );
            mvaVar_tt_MZP1200A0400_dm1TESDown = tt_MZP1200A0400_reader->EvaluateMVA( "MLPBNN" );
        }
        
        if (R.getI("leg1_decayMode")==10 && R.getI("leg2_decayMode")==10)
        {
        
            pt_1_dm10TESUp = pt_1_TESUp;
            mt_1_dm10TESUp = mt_1_TESUp;
            pfmt_1_dm10TESUp = pfmt_1_TESUp;

            pt_2_dm10TESUp = pt_2_TESUp;
            mt_2_dm10TESUp = mt_2_TESUp;
            pfmt_2_dm10TESUp = pfmt_2_TESUp;

            mt_tot_dm10TESUp = mt_tot_TESUp;
            pt_tt_dm10TESUp = pt_tt_TESUp;
            m_vis_dm10TESUp = m_vis_TESUp;
            pfmet_type1_dm10TESUp_Pt = pfmet_type1_TESUp_Pt;
            P_chi_pf_dm10TESUp = P_chi_pf_TESUp;
            LPT_dm10TESUp = LPT_TESUp;
            cos_DeltaPhi_PFMET_Higgs_dm10TESUp = cos_DeltaPhi_PFMET_Higgs_TESUp;

            pt_1_dm10TESDown = pt_1_TESDown;
            mt_1_dm10TESDown = mt_1_TESDown;
            pfmt_1_dm10TESDown = pfmt_1_TESDown;

            pt_2_dm10TESDown = pt_2_TESDown;
            mt_2_dm10TESDown = mt_2_TESDown;
            pfmt_2_dm10TESDown = pfmt_2_TESDown;

            mt_tot_dm10TESDown = mt_tot_TESDown;
            pt_tt_dm10TESDown = pt_tt_TESDown;
            m_vis_dm10TESDown = m_vis_TESDown;
            pfmet_type1_dm10TESDown_Pt = pfmet_type1_TESDown_Pt;
            P_chi_pf_dm10TESDown = P_chi_pf_TESDown;
            LPT_dm10TESDown = LPT_TESDown;
            cos_DeltaPhi_PFMET_Higgs_dm10TESDown = cos_DeltaPhi_PFMET_Higgs_TESDown;

            pt_1_dm10TESDown = pt_1_TESDown;
            mt_1_dm10TESDown = mt_1_TESDown;
            pfmt_1_dm10TESDown = pfmt_1_TESDown;

            pt_2_dm10TESDown = pt_2_TESDown;
            mt_2_dm10TESDown = mt_2_TESDown;
            pfmt_2_dm10TESDown = pfmt_2_TESDown;

            mt_tot_dm10TESDown = mt_tot_TESDown;
            pt_tt_dm10TESDown = pt_tt_TESDown;
            m_vis_dm10TESDown = m_vis_TESDown;
            pfmet_type1_dm10TESDown_Pt = pfmet_type1_TESDown_Pt;
            P_chi_pf_dm10TESDown = P_chi_pf_TESDown;
            LPT_dm10TESDown = LPT_TESDown;
            cos_DeltaPhi_PFMET_Higgs_dm10TESDown = cos_DeltaPhi_PFMET_Higgs_TESDown;

            mvaVar_mt_MZP600A0400_dm10TESUp = mt_MZP600A0400_reader->EvaluateMVA( "MLPBNN" );
            mvaVar_mt_MZP800A0400_dm10TESUp = mt_MZP800A0400_reader->EvaluateMVA( "MLPBNN" );
            mvaVar_mt_MZP1000A0400_dm10TESUp = mt_MZP1000A0400_reader->EvaluateMVA( "MLPBNN" );
            mvaVar_mt_MZP1200A0400_dm10TESUp = mt_MZP1200A0400_reader->EvaluateMVA( "MLPBNN" );
            
            mvaVar_et_MZP600A0400_dm10TESUp = et_MZP600A0400_reader->EvaluateMVA( "MLPBNN" );
            mvaVar_et_MZP800A0400_dm10TESUp = et_MZP800A0400_reader->EvaluateMVA( "MLPBNN" );
            mvaVar_et_MZP1000A0400_dm10TESUp = et_MZP1000A0400_reader->EvaluateMVA("MLPBNN");
            mvaVar_et_MZP1200A0400_dm10TESUp = et_MZP1200A0400_reader->EvaluateMVA( "MLPBNN" );
            
            mvaVar_tt_MZP600A0400_dm10TESUp = tt_MZP600A0400_reader->EvaluateMVA( "MLPBNN" );
            mvaVar_tt_MZP800A0400_dm10TESUp = tt_MZP800A0400_reader->EvaluateMVA( "MLPBNN" );
            mvaVar_tt_MZP1000A0400_dm10TESUp = tt_MZP1000A0400_reader->EvaluateMVA( "MLPBNN" );
            mvaVar_tt_MZP1200A0400_dm10TESUp = tt_MZP1200A0400_reader->EvaluateMVA( "MLPBNN" );
            
            rpt_1 = pt_1_dm10TESDown;
            rpt_2 = pt_2_dm10TESDown;
            rpfmt_1 = pfmt_1_dm10TESDown;
            rmt_tot = mt_tot_dm10TESDown;
            rpt_tt = pt_tt_dm10TESDown;
            rm_vis = m_vis_dm10TESDown;
            rmet = pfmet_type1_dm10TESDown_Pt;
            rP_chi_pf = P_chi_pf_dm10TESDown;
            rLPT = LPT_dm10TESDown;
            rDeltaR_leg1_leg2 = DeltaR_leg1_leg2;
            rcos_DeltaPhi_leg1_leg2 = cos_DeltaPhi_leg1_leg2;
            rcos_DeltaPhi_PFMET_Higgs = cos_DeltaPhi_PFMET_Higgs_dm10TESDown;
            
            mvaVar_mt_MZP600A0400_dm10TESDown = mt_MZP600A0400_reader->EvaluateMVA( "MLPBNN" );
            mvaVar_mt_MZP800A0400_dm10TESDown = mt_MZP800A0400_reader->EvaluateMVA( "MLPBNN" );
            mvaVar_mt_MZP1000A0400_dm10TESDown = mt_MZP1000A0400_reader->EvaluateMVA( "MLPBNN" );
            mvaVar_mt_MZP1200A0400_dm10TESDown = mt_MZP1200A0400_reader->EvaluateMVA( "MLPBNN" );
            
            mvaVar_et_MZP600A0400_dm10TESDown = et_MZP600A0400_reader->EvaluateMVA( "MLPBNN" );
            mvaVar_et_MZP800A0400_dm10TESDown = et_MZP800A0400_reader->EvaluateMVA( "MLPBNN" );
            mvaVar_et_MZP1000A0400_dm10TESDown = et_MZP1000A0400_reader->EvaluateMVA("MLPBNN");
            mvaVar_et_MZP1200A0400_dm10TESDown = et_MZP1200A0400_reader->EvaluateMVA( "MLPBNN" );
            
            mvaVar_tt_MZP600A0400_dm10TESDown = tt_MZP600A0400_reader->EvaluateMVA( "MLPBNN" );
            mvaVar_tt_MZP800A0400_dm10TESDown = tt_MZP800A0400_reader->EvaluateMVA( "MLPBNN" );
            mvaVar_tt_MZP1000A0400_dm10TESDown = tt_MZP1000A0400_reader->EvaluateMVA( "MLPBNN" );
            mvaVar_tt_MZP1200A0400_dm10TESDown = tt_MZP1200A0400_reader->EvaluateMVA( "MLPBNN" );
            
        }
        if (R.getI("leg1_decayMode")==10 && R.getI("leg2_decayMode")!=10)
        {
        
            pt_1_dm10TESUp = pt_1_L1TESUp;
            mt_1_dm10TESUp = mt_1_L1TESUp;
            pfmt_1_dm10TESUp = pfmt_1_L1TESUp;

            pt_2_dm10TESUp = pt_2_L1TESUp;
            mt_2_dm10TESUp = mt_2_L1TESUp;
            pfmt_2_dm10TESUp = pfmt_2_L1TESUp;

            mt_tot_dm10TESUp = mt_tot_L1TESUp;
            pt_tt_dm10TESUp = pt_tt_L1TESUp;
            m_vis_dm10TESUp = m_vis_L1TESUp;
            pfmet_type1_dm10TESUp_Pt = pfmet_type1_L1TESUp_Pt;
            P_chi_pf_dm10TESUp = P_chi_pf_L1TESUp;
            LPT_dm10TESUp = LPT_L1TESUp;
            cos_DeltaPhi_PFMET_Higgs_dm10TESUp = cos_DeltaPhi_PFMET_Higgs_L1TESUp;

            pt_1_dm10TESDown = pt_1_L1TESDown;
            mt_1_dm10TESDown = mt_1_L1TESDown;
            pfmt_1_dm10TESDown = pfmt_1_L1TESDown;

            pt_2_dm10TESDown = pt_2_L1TESDown;
            mt_2_dm10TESDown = mt_2_L1TESDown;
            pfmt_2_dm10TESDown = pfmt_2_L1TESDown;

            mt_tot_dm10TESDown = mt_tot_L1TESDown;
            pt_tt_dm10TESDown = pt_tt_L1TESDown;
            m_vis_dm10TESDown = m_vis_L1TESDown;
            pfmet_type1_dm10TESDown_Pt = pfmet_type1_L1TESDown_Pt;
            P_chi_pf_dm10TESDown = P_chi_pf_L1TESDown;
            LPT_dm10TESDown = LPT_L1TESDown;
            cos_DeltaPhi_PFMET_Higgs_dm10TESDown = cos_DeltaPhi_PFMET_Higgs_L1TESDown;

            rpt_1 = pt_1_L1TESUp;
            rpt_2 = pt_2_L1TESUp;
            rpfmt_1 = pfmt_1_L1TESUp;
            rmt_tot = mt_tot_L1TESUp;
            rpt_tt = pt_tt_L1TESUp;
            rm_vis = m_vis_L1TESUp;
            rmet = pfmet_type1_L1TESUp_Pt;
            rP_chi_pf = P_chi_pf_L1TESUp;
            rLPT = LPT_L1TESUp;
            rDeltaR_leg1_leg2 = DeltaR_leg1_leg2;
            rcos_DeltaPhi_leg1_leg2 = cos_DeltaPhi_leg1_leg2;
            rcos_DeltaPhi_PFMET_Higgs = cos_DeltaPhi_PFMET_Higgs_L1TESUp;
            
            mvaVar_mt_MZP600A0400_dm10TESUp = mt_MZP600A0400_reader->EvaluateMVA( "MLPBNN" );
            mvaVar_mt_MZP800A0400_dm10TESUp = mt_MZP800A0400_reader->EvaluateMVA( "MLPBNN" );
            mvaVar_mt_MZP1000A0400_dm10TESUp = mt_MZP1000A0400_reader->EvaluateMVA( "MLPBNN" );
            mvaVar_mt_MZP1200A0400_dm10TESUp = mt_MZP1200A0400_reader->EvaluateMVA( "MLPBNN" );
            
            mvaVar_et_MZP600A0400_dm10TESUp = et_MZP600A0400_reader->EvaluateMVA( "MLPBNN" );
            mvaVar_et_MZP800A0400_dm10TESUp = et_MZP800A0400_reader->EvaluateMVA( "MLPBNN" );
            mvaVar_et_MZP1000A0400_dm10TESUp = et_MZP1000A0400_reader->EvaluateMVA("MLPBNN");
            mvaVar_et_MZP1200A0400_dm10TESUp = et_MZP1200A0400_reader->EvaluateMVA( "MLPBNN" );
            
            mvaVar_tt_MZP600A0400_dm10TESUp = tt_MZP600A0400_reader->EvaluateMVA( "MLPBNN" );
            mvaVar_tt_MZP800A0400_dm10TESUp = tt_MZP800A0400_reader->EvaluateMVA( "MLPBNN" );
            mvaVar_tt_MZP1000A0400_dm10TESUp = tt_MZP1000A0400_reader->EvaluateMVA( "MLPBNN" );
            mvaVar_tt_MZP1200A0400_dm10TESUp = tt_MZP1200A0400_reader->EvaluateMVA( "MLPBNN" );
            
            rpt_1 = pt_1_L1TESDown;
            rpt_2 = pt_2_L1TESDown;
            rpfmt_1 = pfmt_1_L1TESDown;
            rmt_tot = mt_tot_L1TESDown;
            rpt_tt = pt_tt_L1TESDown;
            rm_vis = m_vis_L1TESDown;
            rmet = pfmet_type1_L1TESDown_Pt;
            rP_chi_pf = P_chi_pf_L1TESDown;
            rLPT = LPT_L1TESDown;
            rDeltaR_leg1_leg2 = DeltaR_leg1_leg2;
            rcos_DeltaPhi_leg1_leg2 = cos_DeltaPhi_leg1_leg2;
            rcos_DeltaPhi_PFMET_Higgs = cos_DeltaPhi_PFMET_Higgs_L1TESDown;
            
            mvaVar_mt_MZP600A0400_dm10TESDown = mt_MZP600A0400_reader->EvaluateMVA( "MLPBNN" );
            mvaVar_mt_MZP800A0400_dm10TESDown = mt_MZP800A0400_reader->EvaluateMVA( "MLPBNN" );
            mvaVar_mt_MZP1000A0400_dm10TESDown = mt_MZP1000A0400_reader->EvaluateMVA( "MLPBNN" );
            mvaVar_mt_MZP1200A0400_dm10TESDown = mt_MZP1200A0400_reader->EvaluateMVA( "MLPBNN" );
            
            mvaVar_et_MZP600A0400_dm10TESDown = et_MZP600A0400_reader->EvaluateMVA( "MLPBNN" );
            mvaVar_et_MZP800A0400_dm10TESDown = et_MZP800A0400_reader->EvaluateMVA( "MLPBNN" );
            mvaVar_et_MZP1000A0400_dm10TESDown = et_MZP1000A0400_reader->EvaluateMVA("MLPBNN");
            mvaVar_et_MZP1200A0400_dm10TESDown = et_MZP1200A0400_reader->EvaluateMVA( "MLPBNN" );
            
            mvaVar_tt_MZP600A0400_dm10TESDown = tt_MZP600A0400_reader->EvaluateMVA( "MLPBNN" );
            mvaVar_tt_MZP800A0400_dm10TESDown = tt_MZP800A0400_reader->EvaluateMVA( "MLPBNN" );
            mvaVar_tt_MZP1000A0400_dm10TESDown = tt_MZP1000A0400_reader->EvaluateMVA( "MLPBNN" );
            mvaVar_tt_MZP1200A0400_dm10TESDown = tt_MZP1200A0400_reader->EvaluateMVA( "MLPBNN" );
        }
        if (R.getI("leg1_decayMode")!=10 && R.getI("leg2_decayMode")==10)
        {
        
            pt_1_dm10TESUp = pt_1_L2TESUp;
            mt_1_dm10TESUp = mt_1_L2TESUp;
            pfmt_1_dm10TESUp = pfmt_1_L2TESUp;

            pt_2_dm10TESUp = pt_2_L2TESUp;
            mt_2_dm10TESUp = mt_2_L2TESUp;
            pfmt_2_dm10TESUp = pfmt_2_L2TESUp;

            mt_tot_dm10TESUp = mt_tot_L2TESUp;
            pt_tt_dm10TESUp = pt_tt_L2TESUp;
            m_vis_dm10TESUp = m_vis_L2TESUp;
            pfmet_type1_dm10TESUp_Pt = pfmet_type1_L2TESUp_Pt;
            P_chi_pf_dm10TESUp = P_chi_pf_L2TESUp;
            LPT_dm10TESUp = LPT_L2TESUp;
            cos_DeltaPhi_PFMET_Higgs_dm10TESUp = cos_DeltaPhi_PFMET_Higgs_L2TESUp;

            pt_1_dm10TESDown = pt_1_L2TESDown;
            mt_1_dm10TESDown = mt_1_L2TESDown;
            pfmt_1_dm10TESDown = pfmt_1_L2TESDown;

            pt_2_dm10TESDown = pt_2_L2TESDown;
            mt_2_dm10TESDown = mt_2_L2TESDown;
            pfmt_2_dm10TESDown = pfmt_2_L2TESDown;

            mt_tot_dm10TESDown = mt_tot_L2TESDown;
            pt_tt_dm10TESDown = pt_tt_L2TESDown;
            m_vis_dm10TESDown = m_vis_L2TESDown;
            pfmet_type1_dm10TESDown_Pt = pfmet_type1_L2TESDown_Pt;
            P_chi_pf_dm10TESDown = P_chi_pf_L2TESDown;
            LPT_dm10TESDown = LPT_L2TESDown;
            cos_DeltaPhi_PFMET_Higgs_dm10TESDown = cos_DeltaPhi_PFMET_Higgs_L2TESDown;

            rpt_1 = pt_1_L2TESUp;
            rpt_2 = pt_2_L2TESUp;
            rpfmt_1 = pfmt_1_L2TESUp;
            rmt_tot = mt_tot_L2TESUp;
            rpt_tt = pt_tt_L2TESUp;
            rm_vis = m_vis_L2TESUp;
            rmet = pfmet_type1_L2TESUp_Pt;
            rP_chi_pf = P_chi_pf_L2TESUp;
            rLPT = LPT_L2TESUp;
            rDeltaR_leg1_leg2 = DeltaR_leg1_leg2;
            rcos_DeltaPhi_leg1_leg2 = cos_DeltaPhi_leg1_leg2;
            rcos_DeltaPhi_PFMET_Higgs = cos_DeltaPhi_PFMET_Higgs_L2TESUp;
            
            mvaVar_mt_MZP600A0400_dm10TESUp = mt_MZP600A0400_reader->EvaluateMVA( "MLPBNN" );
            mvaVar_mt_MZP800A0400_dm10TESUp = mt_MZP800A0400_reader->EvaluateMVA( "MLPBNN" );
            mvaVar_mt_MZP1000A0400_dm10TESUp = mt_MZP1000A0400_reader->EvaluateMVA( "MLPBNN" );
            mvaVar_mt_MZP1200A0400_dm10TESUp = mt_MZP1200A0400_reader->EvaluateMVA( "MLPBNN" );
            
            mvaVar_et_MZP600A0400_dm10TESUp = et_MZP600A0400_reader->EvaluateMVA( "MLPBNN" );
            mvaVar_et_MZP800A0400_dm10TESUp = et_MZP800A0400_reader->EvaluateMVA( "MLPBNN" );
            mvaVar_et_MZP1000A0400_dm10TESUp = et_MZP1000A0400_reader->EvaluateMVA("MLPBNN");
            mvaVar_et_MZP1200A0400_dm10TESUp = et_MZP1200A0400_reader->EvaluateMVA( "MLPBNN" );
            
            mvaVar_tt_MZP600A0400_dm10TESUp = tt_MZP600A0400_reader->EvaluateMVA( "MLPBNN" );
            mvaVar_tt_MZP800A0400_dm10TESUp = tt_MZP800A0400_reader->EvaluateMVA( "MLPBNN" );
            mvaVar_tt_MZP1000A0400_dm10TESUp = tt_MZP1000A0400_reader->EvaluateMVA( "MLPBNN" );
            mvaVar_tt_MZP1200A0400_dm10TESUp = tt_MZP1200A0400_reader->EvaluateMVA( "MLPBNN" );
            
            rpt_1 = pt_1_L2TESDown;
            rpt_2 = pt_2_L2TESDown;
            rpfmt_1 = pfmt_1_L2TESDown;
            rmt_tot = mt_tot_L2TESDown;
            rpt_tt = pt_tt_L2TESDown;
            rm_vis = m_vis_L2TESDown;
            rmet = pfmet_type1_L2TESDown_Pt;
            rP_chi_pf = P_chi_pf_L2TESDown;
            rLPT = LPT_L2TESDown;
            rDeltaR_leg1_leg2 = DeltaR_leg1_leg2;
            rcos_DeltaPhi_leg1_leg2 = cos_DeltaPhi_leg1_leg2;
            rcos_DeltaPhi_PFMET_Higgs = cos_DeltaPhi_PFMET_Higgs_L2TESDown;
            
            mvaVar_mt_MZP600A0400_dm10TESDown = mt_MZP600A0400_reader->EvaluateMVA( "MLPBNN" );
            mvaVar_mt_MZP800A0400_dm10TESDown = mt_MZP800A0400_reader->EvaluateMVA( "MLPBNN" );
            mvaVar_mt_MZP1000A0400_dm10TESDown = mt_MZP1000A0400_reader->EvaluateMVA( "MLPBNN" );
            mvaVar_mt_MZP1200A0400_dm10TESDown = mt_MZP1200A0400_reader->EvaluateMVA( "MLPBNN" );
            
            mvaVar_et_MZP600A0400_dm10TESDown = et_MZP600A0400_reader->EvaluateMVA( "MLPBNN" );
            mvaVar_et_MZP800A0400_dm10TESDown = et_MZP800A0400_reader->EvaluateMVA( "MLPBNN" );
            mvaVar_et_MZP1000A0400_dm10TESDown = et_MZP1000A0400_reader->EvaluateMVA("MLPBNN");
            mvaVar_et_MZP1200A0400_dm10TESDown = et_MZP1200A0400_reader->EvaluateMVA( "MLPBNN" );
            
            mvaVar_tt_MZP600A0400_dm10TESDown = tt_MZP600A0400_reader->EvaluateMVA( "MLPBNN" );
            mvaVar_tt_MZP800A0400_dm10TESDown = tt_MZP800A0400_reader->EvaluateMVA( "MLPBNN" );
            mvaVar_tt_MZP1000A0400_dm10TESDown = tt_MZP1000A0400_reader->EvaluateMVA( "MLPBNN" );
            mvaVar_tt_MZP1200A0400_dm10TESDown = tt_MZP1200A0400_reader->EvaluateMVA( "MLPBNN" );
        }
    }

	/* handle event weights */

    final_weight = getFinalWeight(0, l1, l2);
    nominalCrossSection_Weight = getNominalWeight(0);
    puWeight_Weight = R.getD("puWeight");
    TopQuarkPtWeight_Weight = getTopQuarkPtWeight(0);
    ZReWeight_Weight = getZReWeight(0);
    ZReWeight_WeightUp = getZReWeight(0);
    if (ZReWeight_Weight != 0.) ZReWeight_WeightDown = (1./ZReWeight_Weight);
    KReWeight_Weight = getKFactor(0);
    KReWeight_WeightUp = getKFactorSyst(0,0);
    KReWeight_WeightDown = getKFactorSyst(0,1);
    ZZReWeight_Weight = getZZFactor(0);
    ZZReWeight_WeightUp = getZZFactor(1);
    ZZReWeight_WeightDown = getZZFactor(-1);
    WWReWeight_WeightUp = getWWFactor(0);
    JTF_WeightUp = getJetTauFakeFactor(0,1, l1, l2);
    JTF_WeightDown = getJetTauFakeFactor(0,-1, l1, l2);
    NLOReWeight_Weight = getNLOReWeight(0,10);
    sf_ALD = getALDScaleFactors(0);

    sf_IDISO = getFinalScaleFactorsForPair(0,0,1,1,0,0, l1, l2);
    sf_TRIG = getFinalScaleFactorsForPair(0,0,1,0,1,0, l1, l2);
    sf_TRACK = getFinalScaleFactorsForPair(0,0,1,0,0,1, l1, l2);
    
    ScaleFactorsForPair_Weight = getFinalScaleFactorsForPair(0,0,1,0,0,0, l1, l2);
    ScaleFactorsForPair_WeightUp  = getFinalScaleFactorsForPair(0,1,1,0,0,0, l1, l2);
    ScaleFactorsForPair_WeightDown = getFinalScaleFactorsForPair(0,-1,1,0,0,0, l1, l2);

	std::vector <double> qcd_eleMu = getQCDWeightForEleMuChannel(0, l1, l2);

    QCDWeightForEleMuChannel_Weight = qcd_eleMu[0];
    QCDWeightForEleMuChannel_WeightUp = qcd_eleMu[1];
    QCDWeightForEleMuChannel_WeightDown = qcd_eleMu[2];
    QCDWeightForEleMuChannelNoPZetaCut_Weight = qcd_eleMu[3];
	QCDWeightForEleMuChannelNoPZetaCut_WeightUp = qcd_eleMu[4];
    QCDWeightForEleMuChannelNoPZetaCut_WeightDown = qcd_eleMu[5];

	std::vector<double> highPtTauEff = getHighPtTauUncertainty(0);
    highPtTauEff_WeightUp = highPtTauEff[0];
    highPtTauEff_WeightDown = highPtTauEff[1];

	///////// cutoff ---- XXXXX

    //Final cuts plus loose isolation taus, abreviated sync trees make plotting fast
    
    //bool passMuonFilters = (BadMuonTaggedMoriond17==0 && DuplicateMuonTaggedMoriond17==0); Not used with re-miniAOD
    
    //bool passMetFiltersData = (R.getB("isRealData")==1 && HBHENoiseFilter==1 && HBHENoiseIsoFilter==1 && EcalDeadCellTriggerPrimitiveFilter==1 && goodVerticesFilter==1 && eeBadScFilter==1 && chargedHadronTrackResolutionFilter==1 && globalTightHalo2016Filter==1 && BadChargedCandidateFilter==1 && BadPFMuonFilter==1);
    //bool passMetFiltersMC = (R.getB("isRealData")==0 && HBHENoiseFilter==1 && HBHENoiseIsoFilter==1 && EcalDeadCellTriggerPrimitiveFilter==1 && goodVerticesFilter==1 && chargedHadronTrackResolutionFilter==1 && globalTightHalo2016Filter==1 && BadChargedCandidateFilter==1 && BadPFMuonFilter==1);
    
    //With Data Muon Filters
    //bool passMetFiltersData = (R.getB("isRealData")==1 && HBHENoiseFilter==1 && HBHENoiseIsoFilter==1 && EcalDeadCellTriggerPrimitiveFilter==1 && goodVerticesFilter==1 && eeBadScFilter==1 && chargedHadronTrackResolutionFilter==1 && globalTightHalo2016Filter==1 && BadChargedCandidateFilter==1 && BadPFMuonFilter==1 && BadMuonTaggedMoriond17==0 && DuplicateMuonTaggedMoriond17==0);
    
    bool passMetFiltersData = (R.getB("isRealData")==1 && HBHENoiseFilter==1 && HBHENoiseIsoFilter==1 && EcalDeadCellTriggerPrimitiveFilter==1 && goodVerticesFilter==1 && eeBadScFilter==1 && chargedHadronTrackResolutionFilter==1 && globalTightHalo2016Filter==1 && BadChargedCandidateFilter==1 && BadPFMuonFilter==1);
    
    bool passMetFiltersMC = (R.getB("isRealData")==0 && HBHENoiseFilter==1 && HBHENoiseIsoFilter==1 && EcalDeadCellTriggerPrimitiveFilter==1 && goodVerticesFilter==1 && chargedHadronTrackResolutionFilter==1 && globalTightHalo2016Filter==1 && BadChargedCandidateFilter==1 && BadPFMuonFilter==1);
    
    bool passFilters = ((passMetFiltersData || passMetFiltersMC) && isBoostedChannelPair==0);

    //abrev cuts with lower edge TES cuts, should switch
    //bool AbrevCutsTT = (passFilters && pt_1_TESDown > 40. && pt_2_TESDown > 40. && DeltaR_leg1_leg2 > 0.3 && DeltaR_leg1_leg2 < 2.0  && pairGoodForTrigger==1 && extramuon_veto==0  && extraelec_veto==0 && againstElectronVLooseMVA6_1 > 0.5 && againstMuonLoose3_1 > 0.5 && againstElectronVLooseMVA6_2 > 0.5 && againstMuonLoose3_2 > 0.5 && decayModeFinding_1 > 0.5 && decayModeFinding_2 > 0.5 && ((byLooseIsolationMVArun2v1DBdR03oldDMwLT_1 > 0.5 && byLooseIsolationMVArun2v1DBdR03oldDMwLT_2 > 0.5) || (byLooseIsolationMVArun2v1DBoldDMwLT_1 > 0.5 && byLooseIsolationMVArun2v1DBoldDMwLT_2 > 0.5)));
    
    //bool AbrevCutsET = (passFilters && pt_1 > 26. && pt_2_TESDown > 20. && DeltaR_leg1_leg2 > 0.3 && DeltaR_leg1_leg2 < 2.0 && pairGoodForTrigger==1 && extramuon_veto==0  && extraelec_veto==0 && iso_1 < 0.1 && againstElectronTightMVA6_2 > 0.5 && againstMuonLoose3_2 > 0.5 && (byLooseIsolationMVArun2v1DBdR03oldDMwLT_2 > 0.5 || byLooseIsolationMVArun2v1DBoldDMwLT_2 > 0.5) && decayModeFinding_2 > 0.5);
    
    //bool AbrevCutsMT = (passFilters && pt_1 > 26. && pt_2_TESDown > 20. && DeltaR_leg1_leg2 > 0.3 && DeltaR_leg1_leg2 < 2.0 && pairGoodForTrigger==1 && extramuon_veto==0  && extraelec_veto==0 && iso_1 < 0.15 && againstElectronVLooseMVA6_2 > 0.5 && againstMuonTight3_2 > 0.5 && (byLooseIsolationMVArun2v1DBdR03oldDMwLT_2 > 0.5 || byLooseIsolationMVArun2v1DBoldDMwLT_2 > 0.5) && decayModeFinding_2 > 0.5 );
    
    bool AbrevCutsTT = (passFilters && pt_1 > 40. && pt_2 > 40. && DeltaR_leg1_leg2 > 0.3  && pairGoodForTrigger==1 && extramuon_veto==0  && extraelec_veto==0 && againstElectronVLooseMVA6_1 > 0.5 && againstMuonLoose3_1 > 0.5 && againstElectronVLooseMVA6_2 > 0.5 && againstMuonLoose3_2 > 0.5 && decayModeFinding_1 > 0.5 && decayModeFinding_2 > 0.5 && ((byLooseIsolationMVArun2v1DBdR03oldDMwLT_1 > 0.5 && byLooseIsolationMVArun2v1DBdR03oldDMwLT_2 > 0.5) || (byLooseIsolationMVArun2v1DBoldDMwLT_1 > 0.5 && byLooseIsolationMVArun2v1DBoldDMwLT_2 > 0.5)));
    
    bool AbrevCutsET = (passFilters && pt_1 > 26. && pt_2 > 20. && DeltaR_leg1_leg2 > 0.3 && pairGoodForTrigger==1 && extramuon_veto==0  && extraelec_veto==0 && iso_1 < 0.3 && againstElectronTightMVA6_2 > 0.5 && againstMuonLoose3_2 > 0.5 && (byLooseIsolationMVArun2v1DBdR03oldDMwLT_2 > 0.5 || byLooseIsolationMVArun2v1DBoldDMwLT_2 > 0.5) && decayModeFinding_2 > 0.5);
    
    bool AbrevCutsMT = (passFilters && pt_1 > 26. && pt_2 > 20. && DeltaR_leg1_leg2 > 0.3 && pairGoodForTrigger==1 && extramuon_veto==0  && extraelec_veto==0 && iso_1 < 0.3 && againstElectronVLooseMVA6_2 > 0.5 && againstMuonTight3_2 > 0.5 && (byLooseIsolationMVArun2v1DBdR03oldDMwLT_2 > 0.5 || byLooseIsolationMVArun2v1DBoldDMwLT_2 > 0.5) && decayModeFinding_2 > 0.5 );

	if(AbrevCutsET && R.getI("CandidateEventType")==3) {num_et++; tree_EleTau->Fill();}
	else if(R.getI("CandidateEventType")==2) {num_em++; tree_EleMu->Fill();}
	else if(AbrevCutsMT && R.getI("CandidateEventType")==5) {num_mt++; tree_MuTau->Fill();}
	else if(AbrevCutsTT && R.getI("CandidateEventType")==6) {num_tt++; tree_TauTau->Fill();}
	
	if(num_total%1000==0){
	 std::cout<<" etau = "<<num_et<<"\n";
	 std::cout<<" mtau = "<<num_mt<<"\n";
	 std::cout<<" tt = "<<num_tt<<"\n";
	 std::cout<<" em = "<<num_em<<"\n";
	}

}

/* given met and met phi return a size 2 vector with P_chi at [0] and M_min at [1] */
std::vector <double> generateH2TauSyncTree::computePchi_and_Mmin(bool verbose_, double ME_T, double ME_T_phi, TLorentzVector l1)
{

	if(verbose_)
	{
		std::cout<<" calling computePchi_and_Mmin function \n";
	}

	std::vector<double> returnVector;
	returnVector.clear();

    double P_chi_ =  -999.0;
  	double M_min_ =  -999.0;

    double p_lt = 0;
    double eta_l = 0;
    double phi_l = 0;
    double M_l = 0;

    double p_lx = 0;
    double p_ly = 0;
    double p_lz = 0;
    double E_l = 0;

    double p_vx = 0;
    double p_vy = 0;
    double E_v = 0;
    
    //Pchi and Mmin Calculations
    
    if ((R.getI("leg1_leptonType") == 1 && R.getI("leg2_leptonType") == 3) || (R.getI("leg1_leptonType") == 2 && R.getI("leg2_leptonType") == 3))
    {


		if(verbose_)
		{
			std::cout<<" valid lepton types for Pchi and Mmin computation ";
		}
    
        p_lt = l1.Pt();
        eta_l = l1.Eta();
        phi_l = l1.Phi();;
        M_l = l1.M();
    
       //Calculate p_lx and p_ly   
        p_lx = p_lt*cos(phi_l);
        p_ly = p_lt*sin(phi_l);
       
       //Calculate p_lz
        p_lz = p_lt*TMath::SinH(eta_l);

       //Calculate E_l
        E_l = TMath::Sqrt(M_l*M_l + p_lt*p_lt*TMath::CosH(eta_l)*TMath::CosH(eta_l));
       
       //Calculate p_vx and p_vy 
        p_vx = ME_T*TMath::Cos(ME_T_phi);
        p_vy = ME_T*TMath::Sin(ME_T_phi);
       
       //Calculate E_v
        E_v = (ME_T*E_l)/(TMath::Sqrt(E_l*E_l - p_lz*p_lz));
       

       //Calculate P_chi
        P_chi_ = ME_T*(TMath::Cos(phi_l)*TMath::Cos(ME_T_phi) + TMath::Sin(phi_l)*TMath::Sin(ME_T_phi));

       //Calculate M_min
        M_min_ = TMath::Sqrt(((E_l + E_v)*(E_l + E_v)) - ((p_lx + p_vx)*(p_lx + p_vx)) - ((p_ly + p_vy)*(p_ly + p_vy)) - ((1 + (E_v/E_l))*(1 + (E_v/E_l))*p_lz*p_lz));

    

        if(verbose_)
		{
			std::cout<<" returning [ P_chi_ = "<< P_chi_<< " , M_min_ = "<< M_min_ <<" ] \n";
		}

    }

    else
    {
		if(verbose_)
		{
			std::cout<<" invalid lepton types for Pchi and Mmin computation, returning [-999.0, -999.0] \n";
		}

    }

    returnVector.clear();
    returnVector.push_back(P_chi_);
    returnVector.push_back(M_min_);

	return returnVector;
}






double generateH2TauSyncTree::computeLPT( bool verbose_, TLorentzVector l1, TLorentzVector l2)
{

    // LPT Calculation 


	double LPT_ = 0.0;    
    double visFrac1 = 0.;
    double visFrac1Prob = 0.;
    double visFrac2 = 0.;
    double visFrac2Prob = 0.;
    double currentFracProb = 0.;
    double bestFracProb = 0.;
    double bestP = 0.;


    TVector3 tauVisVec1;
    TVector3 tauVisVec2;
    TVector3 totalVisVec;
    tauVisVec1.SetPtEtaPhi(l1.Pt(),l1.Eta(),l1.Phi());
    tauVisVec2.SetPtEtaPhi(l2.Pt(),l2.Eta(),l2.Phi());
    totalVisVec = tauVisVec1 + tauVisVec2;
    //double totalVisPt = totalVisVec.Pt();
    //assign leg PDF histograms
    
    if(verbose_)
    {
    	std::cout << "Start PDF Loops " << R.getI("CandidateEventType") << std::endl;
	    std::cout << "Leg One Decay Mode: " << R.getI("leg1_decayMode") << std::endl;
    	std::cout << "Leg One Num Hads: " << R.getF("leg1_numHadrons") << std::endl;
    	std::cout << "Leg One Num Strips: " << R.getF("leg1_numStrips") << std::endl;
    	std::cout << "Leg Two Decay Mode: " << R.getI("leg2_decayMode") << std::endl;
    	std::cout << "Leg Two Num Hads: " << R.getF("leg2_numHadrons") << std::endl;
    	std::cout << "Leg Two Num Strips: " << R.getF("leg2_numStrips") << std::endl;
    }


    for (int i = 0; i < eHistoFrac->GetNbinsX(); i++)
    {
        if((R.getI("CandidateEventType")==1) || (R.getI("CandidateEventType")==2) || (R.getI("CandidateEventType")==3))
        {
            visFrac1 = eHistoFrac->GetXaxis()->GetBinCenter(i);
            visFrac1Prob = eHistoFrac->GetBinContent(i);
        }
        else if ((R.getI("CandidateEventType")==4) || (R.getI("CandidateEventType")==5))
        {
            visFrac1 = muHistoFrac->GetXaxis()->GetBinCenter(i);
            visFrac1Prob = muHistoFrac->GetBinContent(i);
        }
        else
        {
            if (R.getI("leg1_decayMode")==0)
            {
                visFrac1 = h1p0sHistoFrac->GetXaxis()->GetBinCenter(i);
                visFrac1Prob = h1p0sHistoFrac->GetBinContent(i);
            }
            else if (R.getI("leg1_decayMode")==1)
            {
                visFrac1 = h1p1sHistoFrac->GetXaxis()->GetBinCenter(i);
                visFrac1Prob = h1p1sHistoFrac->GetBinContent(i);
            }
            else if (R.getI("leg1_decayMode")==10)
            {
                visFrac1 = h3p0sHistoFrac->GetXaxis()->GetBinCenter(i);
                visFrac1Prob = h3p0sHistoFrac->GetBinContent(i);
            }
            
        }
        
        //find constrained fraction
        visFrac2 = 2 * (1/visFrac1) * ((tauVisVec1.Mag()*tauVisVec2.Mag())-(tauVisVec1.Dot(tauVisVec2)))/(pow(125.,2));
        
        //find PDF value from constrained fraction
        if(R.getI("CandidateEventType")==1)
        {
            visFrac2Prob = eHistoFrac->GetBinContent(eHistoFrac->FindBin(visFrac2));
        }
        else if ((R.getI("CandidateEventType")==2) || (R.getI("CandidateEventType")==4))
        {
            visFrac2Prob = muHistoFrac->GetBinContent(muHistoFrac->FindBin(visFrac2));
        }
        else
        {
            if (R.getI("leg2_decayMode")==0)
            {
                visFrac2Prob = h1p0sHistoFrac->GetBinContent(h1p0sHistoFrac->FindBin(visFrac2));
            }
            else if (R.getI("leg2_decayMode")==1)
            {
                visFrac2Prob = h1p1sHistoFrac->GetBinContent(h1p1sHistoFrac->FindBin(visFrac2));
            }
            else if (R.getI("leg2_decayMode")==10)
            {
                visFrac2Prob = h3p0sHistoFrac->GetBinContent(h3p0sHistoFrac->FindBin(visFrac2));
            }

        }
        if (visFrac2 > 1.03)
        {
            visFrac2Prob = 0;
        }
        TVector3 tauVec1;
        TVector3 tauVec2;
        TVector3 totalVec;
        tauVec1.SetPtEtaPhi(l1.Pt() * (1/visFrac1),l1.Eta(),l1.Phi());
        tauVec2.SetPtEtaPhi(l2.Pt() * (1/visFrac2),l2.Eta(),l2.Phi());
        totalVec = tauVec1 + tauVec2;
        double totalP = totalVec.Pt();
        currentFracProb = visFrac1Prob * visFrac2Prob;
        if (currentFracProb > bestFracProb)
        {
            bestFracProb = currentFracProb;
            bestP = totalP;
        }
    }
    

    LPT_ = bestP;

    if(verbose_)
    {
    	std::cout << "best P :" << bestP << std::endl;
	    std::cout << "best Prob: " << bestFracProb << std::endl;
    }

    return LPT_;
}

double generateH2TauSyncTree::ttTrigPtShape(bool up, TLorentzVector l1, TLorentzVector l2)
{

    double ttPt1Var = 0.01*(-0.5*(std::min(60.0,l1.Pt())) + 32.);
    double ttPt2Var = 0.01*(-0.5*(std::min(60.0,l2.Pt())) + 32.);
    
    double weightUp = 1.0 + ttPt1Var + ttPt2Var;
    double weightDown =  1.0 - ttPt1Var - ttPt2Var;
    if (up) {return weightUp;}
    else {return weightDown;}

}



void generateH2TauSyncTree::setupBranches(TTree * T)
{
	/* all trees get same branch address, but will only write event
	 to one tree per event based on CandidateEventType */

	T->Branch("final_weight", &final_weight);
	T->Branch("nominalCrossSection_Weight", &nominalCrossSection_Weight);
	T->Branch("puWeight_Weight", &puWeight_Weight);
    T->Branch("sf_IDISO", &sf_IDISO);
    T->Branch("sf_TRIG", &sf_TRIG);
    T->Branch("sf_TRACK", &sf_TRACK);
    T->Branch("sf_ALD", &sf_ALD);
	T->Branch("TopQuarkPtWeight_Weight", &TopQuarkPtWeight_Weight);
	T->Branch("ZReWeight_Weight", &ZReWeight_Weight);
    T->Branch("ZReWeight_WeightUp", &ZReWeight_WeightUp);
    T->Branch("ZReWeight_WeightDown", &ZReWeight_WeightDown);
    T->Branch("KReWeight_Weight", &KReWeight_Weight);
    T->Branch("KReWeight_WeightUp", &KReWeight_WeightUp);
    T->Branch("KReWeight_WeightDown", &KReWeight_WeightDown);
    T->Branch("ZZReWeight_Weight", &ZZReWeight_Weight);
    T->Branch("ZZReWeight_WeightUp", &ZZReWeight_WeightUp);
    T->Branch("ZZReWeight_WeightDown", &ZZReWeight_WeightDown);
    T->Branch("WWReWeight_WeightUp", &WWReWeight_WeightUp);
    T->Branch("JTF_WeightUp", &JTF_WeightUp);
    T->Branch("JTF_WeightDown", &JTF_WeightDown);
	T->Branch("NLOReWeight_Weight", &NLOReWeight_Weight);
	T->Branch("ScaleFactorsForPair_Weight", &ScaleFactorsForPair_Weight);
	T->Branch("ScaleFactorsForPair_WeightUp", &ScaleFactorsForPair_WeightUp);
	T->Branch("ScaleFactorsForPair_WeightDown", &ScaleFactorsForPair_WeightDown);
	T->Branch("QCDWeightForEleMuChannel_Weight", &QCDWeightForEleMuChannel_Weight);
	T->Branch("QCDWeightForEleMuChannel_WeightUp", &QCDWeightForEleMuChannel_WeightUp);
	T->Branch("QCDWeightForEleMuChannel_WeightDown", &QCDWeightForEleMuChannel_WeightDown);
	T->Branch("QCDWeightForEleMuChannelNoPZetaCut_Weight", &QCDWeightForEleMuChannelNoPZetaCut_Weight);
	T->Branch("QCDWeightForEleMuChannelNoPZetaCut_WeightUp", &QCDWeightForEleMuChannelNoPZetaCut_WeightUp);
	T->Branch("QCDWeightForEleMuChannelNoPZetaCut_WeightDown", &QCDWeightForEleMuChannelNoPZetaCut_WeightDown);
	T->Branch("highPtTauEff_WeightUp", &highPtTauEff_WeightUp);
	T->Branch("highPtTauEff_WeightDown", &highPtTauEff_WeightDown);

    T->Branch("flag_MVAEventType", &flag_MVAEventType);
    T->Branch("randNum", &randNum);

	T->Branch("originalXWGTUP", &originalXWGTUP);
	T->Branch("theory_scale_factors", &theory_scale_factors);
	T->Branch("pairRank", &pairRank);
	T->Branch("isOsPair", &isOsPair);
    T->Branch("isBoostedChannelPair", &isBoostedChannelPair);
    T->Branch("DataCardInt", &DataCardInt);
	T->Branch("genBosonTotal_pt", &genBosonTotal_pt);
	T->Branch("genBosonTotal_eta", &genBosonTotal_eta);
	T->Branch("genBosonTotal_phi", &genBosonTotal_phi);
	T->Branch("genBosonTotal_M", &genBosonTotal_M);
	T->Branch("genBosonVisible_pt", &genBosonVisible_pt);
	T->Branch("genBosonVisible_eta", &genBosonVisible_eta);
	T->Branch("genBosonVisible_phi", &genBosonVisible_phi);
	T->Branch("genBosonVisible_M", &genBosonVisible_M);
    T->Branch("genBosonTotal_Wpt", &genBosonTotal_Wpt);
    T->Branch("leg1_GENMOTHERpdgId", &leg1_GENMOTHERpdgId);
    T->Branch("leg2_GENMOTHERpdgId", &leg2_GENMOTHERpdgId);
	T->Branch("genTopPt1", &genTopPt1);
	T->Branch("genTopPt2", &genTopPt2);
	T->Branch("run", &run);
	T->Branch("event", &event);
	T->Branch("evt", &evt);
	T->Branch("lumi", &lumi);
	T->Branch("npv", &npv);
	T->Branch("npu", &npu);
	T->Branch("rho", &rho);
	T->Branch("puweight", &puweight);
	T->Branch("pt_1", &pt_1);
	T->Branch("phi_1", &phi_1);
	T->Branch("eta_1", &eta_1);
	T->Branch("m_1", &m_1);
    T->Branch("pt_1_flat", &pt_1_flat);
	T->Branch("phi_1_flat", &phi_1_flat);
	T->Branch("eta_1_flat", &eta_1_flat);
	T->Branch("m_1_flat", &m_1_flat);
	T->Branch("iso_1", &iso_1);
	T->Branch("dZ_1", &dZ_1);
	T->Branch("dzTauVertex_1",&dzTauVertex_1);
	T->Branch("d0_1", &d0_1);
	T->Branch("q_1", &q_1);
	T->Branch("id_e_mva_nt_loose_1", &id_e_mva_nt_loose_1);
	T->Branch("tau_decay_mode_1", &tau_decay_mode_1);
	T->Branch("ZimpactTau_1", &ZimpactTau_1);
	T->Branch("mt_1", &mt_1);
    T->Branch("pfmt_1", &pfmt_1);
    T->Branch("pfmt_1_UESUp", &pfmt_1_UESUp);
    T->Branch("pfmt_1_UESDown", &pfmt_1_UESDown);
	T->Branch("puppimt_1", &puppimt_1);
	T->Branch("mt_uncorr_1", &mt_uncorr_1);
	T->Branch("responseUP_MTmvaMET_1", &responseUP_MTmvaMET_1);
	T->Branch("responseDOWN_MTmvaMET_1", &responseDOWN_MTmvaMET_1);
	T->Branch("resolutionUP_MTmvaMET_1", &resolutionUP_MTmvaMET_1);
	T->Branch("resolutionDOWN_MTmvaMET_1", &resolutionDOWN_MTmvaMET_1);
	T->Branch("gen_match_1", &gen_match_1);
	T->Branch("genMCmatch_pt_1", &genMCmatch_pt_1);
	T->Branch("genMCmatch_eta_1", &genMCmatch_eta_1);
	T->Branch("genMCmatch_phi_1", &genMCmatch_phi_1);
	T->Branch("genMCmatch_M_1", &genMCmatch_M_1);
	T->Branch("MCMatchPdgId_1", &MCMatchPdgId_1);
	T->Branch("byIsolationMVArun2v1DBdR03oldDMwLTraw_1", &byIsolationMVArun2v1DBdR03oldDMwLTraw_1);
	T->Branch("byTightIsolationMVArun2v1DBdR03oldDMwLT_1", &byTightIsolationMVArun2v1DBdR03oldDMwLT_1);
	T->Branch("byVTightIsolationMVArun2v1DBdR03oldDMwLT_1", &byVTightIsolationMVArun2v1DBdR03oldDMwLT_1);
	T->Branch("byLooseIsolationMVArun2v1DBdR03oldDMwLT_1", &byLooseIsolationMVArun2v1DBdR03oldDMwLT_1);
	T->Branch("byMediumIsolationMVArun2v1DBdR03oldDMwLT_1", &byMediumIsolationMVArun2v1DBdR03oldDMwLT_1);
	T->Branch("byVLooseIsolationMVArun2v1DBdR03oldDMwLT_1", &byVLooseIsolationMVArun2v1DBdR03oldDMwLT_1);
	T->Branch("byVVTightIsolationMVArun2v1DBdR03oldDMwLT_1", &byVVTightIsolationMVArun2v1DBdR03oldDMwLT_1);
    T->Branch("byLooseIsolationMVArun2v1DBoldDMwLT_1", &byLooseIsolationMVArun2v1DBoldDMwLT_1);
    T->Branch("byTightIsolationMVArun2v1DBoldDMwLT_1", &byTightIsolationMVArun2v1DBoldDMwLT_1);
	T->Branch("againstElectronVLooseMVA6_1", &againstElectronVLooseMVA6_1);
	T->Branch("againstMuonTight3_1", &againstMuonTight3_1);
	T->Branch("againstElectronTightMVA6_1", &againstElectronTightMVA6_1);
	T->Branch("againstMuonLoose3_1", &againstMuonLoose3_1);
	T->Branch("decayModeFinding_1", &decayModeFinding_1);
	//T->Branch("byIsolationMVA3oldDMwLTraw_1", &byIsolationMVA3oldDMwLTraw_1);
	T->Branch("byCombinedIsolationDeltaBetaCorrRaw3Hits_1", &byCombinedIsolationDeltaBetaCorrRaw3Hits_1);
	T->Branch("byIsolationMVArun2v1DBnewDMwLTraw_1", &byIsolationMVArun2v1DBnewDMwLTraw_1);
	T->Branch("decayModeFindingNewDMs_1", &decayModeFindingNewDMs_1);
	T->Branch("pt_2", &pt_2);
	T->Branch("phi_2", &phi_2);
	T->Branch("eta_2", &eta_2);
	T->Branch("m_2", &m_2);
    T->Branch("pt_2_flat", &pt_2_flat);
	T->Branch("phi_2_flat", &phi_2_flat);
	T->Branch("eta_2_flat", &eta_2_flat);
	T->Branch("m_2_flat", &m_2_flat);
	T->Branch("iso_2", &iso_2);
	T->Branch("dZ_2", &dZ_2);
	T->Branch("dzTauVertex_2", &dzTauVertex_2);
	T->Branch("d0_2", &d0_2);
	T->Branch("q_2", &q_2);
	T->Branch("id_e_mva_nt_loose_2", &id_e_mva_nt_loose_2);
	T->Branch("tau_decay_mode_2", &tau_decay_mode_2);
	T->Branch("ZimpactTau_2", &ZimpactTau_2);
	T->Branch("mt_2", &mt_2);
	T->Branch("pfmt_2", &pfmt_2);
	T->Branch("puppimt_2", &puppimt_2);
	T->Branch("mt_uncorr_2", &mt_uncorr_2);
	T->Branch("responseUP_MTmvaMET_2", &responseUP_MTmvaMET_2);
	T->Branch("responseDOWN_MTmvaMET_2", &responseDOWN_MTmvaMET_2);
	T->Branch("resolutionUP_MTmvaMET_2", &resolutionUP_MTmvaMET_2);
	T->Branch("resolutionDOWN_MTmvaMET_2", &resolutionDOWN_MTmvaMET_2);
	T->Branch("gen_match_2", &gen_match_2);
	T->Branch("genMCmatch_pt_2", &genMCmatch_pt_2);
	T->Branch("genMCmatch_eta_2", &genMCmatch_eta_2);
	T->Branch("genMCmatch_phi_2", &genMCmatch_phi_2);
	T->Branch("genMCmatch_M_2", &genMCmatch_M_2);
	T->Branch("MCMatchPdgId_2", &MCMatchPdgId_2);
	T->Branch("byIsolationMVArun2v1DBdR03oldDMwLTraw_2", &byIsolationMVArun2v1DBdR03oldDMwLTraw_2);
	T->Branch("byTightIsolationMVArun2v1DBdR03oldDMwLT_2", &byTightIsolationMVArun2v1DBdR03oldDMwLT_2);
	T->Branch("byVTightIsolationMVArun2v1DBdR03oldDMwLT_2", &byVTightIsolationMVArun2v1DBdR03oldDMwLT_2);
	T->Branch("byLooseIsolationMVArun2v1DBdR03oldDMwLT_2", &byLooseIsolationMVArun2v1DBdR03oldDMwLT_2);
	T->Branch("byMediumIsolationMVArun2v1DBdR03oldDMwLT_2", &byMediumIsolationMVArun2v1DBdR03oldDMwLT_2);
	T->Branch("byVLooseIsolationMVArun2v1DBdR03oldDMwLT_2", &byVLooseIsolationMVArun2v1DBdR03oldDMwLT_2);
	T->Branch("byVVTightIsolationMVArun2v1DBdR03oldDMwLT_2", &byVVTightIsolationMVArun2v1DBdR03oldDMwLT_2);
    T->Branch("byLooseIsolationMVArun2v1DBoldDMwLT_2", &byLooseIsolationMVArun2v1DBoldDMwLT_2);
    T->Branch("byTightIsolationMVArun2v1DBoldDMwLT_2", &byTightIsolationMVArun2v1DBoldDMwLT_2);
	T->Branch("againstElectronVLooseMVA6_2", &againstElectronVLooseMVA6_2);
	T->Branch("againstMuonTight3_2", &againstMuonTight3_2);
	T->Branch("againstElectronTightMVA6_2", &againstElectronTightMVA6_2);
	T->Branch("againstMuonLoose3_2", &againstMuonLoose3_2);
	T->Branch("decayModeFinding_2", &decayModeFinding_2);
	//T->Branch("byIsolationMVA3oldDMwLTraw_2", &byIsolationMVA3oldDMwLTraw_2);
	T->Branch("byCombinedIsolationDeltaBetaCorrRaw3Hits_2", &byCombinedIsolationDeltaBetaCorrRaw3Hits_2);
	T->Branch("byIsolationMVArun2v1DBnewDMwLTraw_2", &byIsolationMVArun2v1DBnewDMwLTraw_2);
	T->Branch("decayModeFindingNewDMs_2", &decayModeFindingNewDMs_2);
	T->Branch("pt_tt", &pt_tt);
	T->Branch("DeltaR_leg1_leg2", &DeltaR_leg1_leg2);
    T->Branch("cos_DeltaPhi_leg1_leg2", &cos_DeltaPhi_leg1_leg2);
    T->Branch("cos_DeltaPhi_PFMET_Higgs", &cos_DeltaPhi_PFMET_Higgs);
    T->Branch("cos_DeltaPhi_PFMET_Higgs_UESUp", &cos_DeltaPhi_PFMET_Higgs_UESUp);
    T->Branch("cos_DeltaPhi_PFMET_Higgs_UESDown", &cos_DeltaPhi_PFMET_Higgs_UESDown);
	T->Branch("mt_tot", &mt_tot);
    T->Branch("weight_ttPtUp", &weight_ttPtUp);
    T->Branch("weight_ttPtDown", &weight_ttPtDown);
    T->Branch("mt_tot_UESUp", &mt_tot_UESUp);
    T->Branch("mt_tot_UESDown", &mt_tot_UESDown);
    T->Branch("mt_tot_JEnUp", &mt_tot_JEnUp);
    T->Branch("mt_tot_JEnDown", &mt_tot_JEnDown);
	T->Branch("m_vis", &m_vis);
    T->Branch("m_vis_flat", &m_vis_flat);
	T->Branch("m_sv", &m_sv);
	T->Branch("mt_sv", &mt_sv);
	T->Branch("SVFit_mvaMET_diTau_pt", &SVFit_mvaMET_diTau_pt);
	T->Branch("SVFit_mvaMET_diTau_eta", &SVFit_mvaMET_diTau_eta);
	T->Branch("SVFit_mvaMET_diTau_phi", &SVFit_mvaMET_diTau_phi);
	T->Branch("SVFit_mvaMET_FittedMET", &SVFit_mvaMET_FittedMET);
	T->Branch("SVFit_mvaMET_FittedMETphi", &SVFit_mvaMET_FittedMETphi);
	T->Branch("mvamet", &mvamet);
	T->Branch("mvametphi", &mvametphi);
	T->Branch("met", &met);
	T->Branch("metphi", &metphi);
	T->Branch("puppimet", &puppimet);
	T->Branch("puppimetphi", &puppimetphi);
	T->Branch("uncorr_mvamet", &uncorr_mvamet);
	T->Branch("uncorr_mvametphi", &uncorr_mvametphi);
	T->Branch("responseUP_mvaMET", &responseUP_mvaMET);
	T->Branch("responseUP_mvaMETphi", &responseUP_mvaMETphi);
	T->Branch("responseDOWN_mvaMET", &responseDOWN_mvaMET);
	T->Branch("responseDOWN_mvaMETphi", &responseDOWN_mvaMETphi);
	T->Branch("resolutionUP_mvaMET", &resolutionUP_mvaMET);
	T->Branch("resolutionUP_mvaMETphi", &resolutionUP_mvaMETphi);
	T->Branch("resolutionDOWN_mvaMET", &resolutionDOWN_mvaMET);
	T->Branch("resolutionDOWN_mvaMETphi", &resolutionDOWN_mvaMETphi);
	T->Branch("mvacov00", &mvacov00);
	T->Branch("mvacov01", &mvacov01);
	T->Branch("mvacov10", &mvacov10);
	T->Branch("mvacov11", &mvacov11);
	T->Branch("metcov00", &metcov00);
	T->Branch("metcov01", &metcov01);
	T->Branch("metcov10", &metcov10);
	T->Branch("metcov11", &metcov11);
	T->Branch("genMET", &genMET);
	T->Branch("genMETphi", &genMETphi);
	T->Branch("genMETeta", &genMETeta);
	T->Branch("genMETmass", &genMETmass);

	T->Branch("pfmet_raw_Pt" , &pfmet_raw_Pt);
	T->Branch("pfmet_raw_Phi" , &pfmet_raw_Phi);
	T->Branch("pfmet_raw_MT1" , &pfmet_raw_MT1);
	T->Branch("pfmet_raw_MT2" , &pfmet_raw_MT2);
	T->Branch("pfmet_type1_Pt" , &pfmet_type1_Pt);
	T->Branch("pfmet_type1_Phi" , &pfmet_type1_Phi);
	T->Branch("pfmet_type1_MT1" , &pfmet_type1_MT1);
	T->Branch("pfmet_type1_MT2" , &pfmet_type1_MT2);
	T->Branch("pfmet_type1_JetResUp_Pt" , &pfmet_type1_JetResUp_Pt);
	T->Branch("pfmet_type1_JetResUp_Phi" , &pfmet_type1_JetResUp_Phi);
	T->Branch("pfmet_type1_JetResUp_MT1" , &pfmet_type1_JetResUp_MT1);
	T->Branch("pfmet_type1_JetResUp_MT2" , &pfmet_type1_JetResUp_MT2);
	T->Branch("pfmet_type1_JetResDown_Pt" , &pfmet_type1_JetResDown_Pt);
	T->Branch("pfmet_type1_JetResDown_Phi" , &pfmet_type1_JetResDown_Phi);
	T->Branch("pfmet_type1_JetResDown_MT1" , &pfmet_type1_JetResDown_MT1);
	T->Branch("pfmet_type1_JetResDown_MT2" , &pfmet_type1_JetResDown_MT2);
	T->Branch("pfmet_type1_JetEnUp_Pt" , &pfmet_type1_JetEnUp_Pt);
	T->Branch("pfmet_type1_JetEnUp_Phi" , &pfmet_type1_JetEnUp_Phi);
	T->Branch("pfmet_type1_JetEnUp_MT1" , &pfmet_type1_JetEnUp_MT1);
	T->Branch("pfmet_type1_JetEnUp_MT2" , &pfmet_type1_JetEnUp_MT2);
	T->Branch("pfmet_type1_JetEnDown_Pt" , &pfmet_type1_JetEnDown_Pt);
	T->Branch("pfmet_type1_JetEnDown_Phi" , &pfmet_type1_JetEnDown_Phi);
	T->Branch("pfmet_type1_JetEnDown_MT1" , &pfmet_type1_JetEnDown_MT1);
	T->Branch("pfmet_type1_JetEnDown_MT2" , &pfmet_type1_JetEnDown_MT2);
	T->Branch("pfmet_type1_MuonEnUp_Pt" , &pfmet_type1_MuonEnUp_Pt);
	T->Branch("pfmet_type1_MuonEnUp_Phi" , &pfmet_type1_MuonEnUp_Phi);
	T->Branch("pfmet_type1_MuonEnUp_MT1" , &pfmet_type1_MuonEnUp_MT1);
	T->Branch("pfmet_type1_MuonEnUp_MT2" , &pfmet_type1_MuonEnUp_MT2);
	T->Branch("pfmet_type1_MuonEnDown_Pt" , &pfmet_type1_MuonEnDown_Pt);
	T->Branch("pfmet_type1_MuonEnDown_Phi" , &pfmet_type1_MuonEnDown_Phi);
	T->Branch("pfmet_type1_MuonEnDown_MT1" , &pfmet_type1_MuonEnDown_MT1);
	T->Branch("pfmet_type1_MuonEnDown_MT2" , &pfmet_type1_MuonEnDown_MT2);
	T->Branch("pfmet_type1_ElectronEnUp_Pt" , &pfmet_type1_ElectronEnUp_Pt);
	T->Branch("pfmet_type1_ElectronEnUp_Phi" , &pfmet_type1_ElectronEnUp_Phi);
	T->Branch("pfmet_type1_ElectronEnUp_MT1" , &pfmet_type1_ElectronEnUp_MT1);
	T->Branch("pfmet_type1_ElectronEnUp_MT2" , &pfmet_type1_ElectronEnUp_MT2);
	T->Branch("pfmet_type1_ElectronEnDown_Pt" , &pfmet_type1_ElectronEnDown_Pt);
	T->Branch("pfmet_type1_ElectronEnDown_Phi" , &pfmet_type1_ElectronEnDown_Phi);
	T->Branch("pfmet_type1_ElectronEnDown_MT1" , &pfmet_type1_ElectronEnDown_MT1);
	T->Branch("pfmet_type1_ElectronEnDown_MT2" , &pfmet_type1_ElectronEnDown_MT2);
	T->Branch("pfmet_type1_UnclusteredEnUp_Pt" , &pfmet_type1_UnclusteredEnUp_Pt);
	T->Branch("pfmet_type1_UnclusteredEnUp_Phi" , &pfmet_type1_UnclusteredEnUp_Phi);
	T->Branch("pfmet_type1_UnclusteredEnUp_MT1" , &pfmet_type1_UnclusteredEnUp_MT1);
	T->Branch("pfmet_type1_UnclusteredEnUp_MT2" , &pfmet_type1_UnclusteredEnUp_MT2);
	T->Branch("pfmet_type1_UnclusteredEnDown_Pt" , &pfmet_type1_UnclusteredEnDown_Pt);
	T->Branch("pfmet_type1_UnclusteredEnDown_Phi" , &pfmet_type1_UnclusteredEnDown_Phi);
	T->Branch("pfmet_type1_UnclusteredEnDown_MT1" , &pfmet_type1_UnclusteredEnDown_MT1);
	T->Branch("pfmet_type1_UnclusteredEnDown_MT2" , &pfmet_type1_UnclusteredEnDown_MT2);
	T->Branch("pfmet_type1_PhotonEnUp_Pt" , &pfmet_type1_PhotonEnUp_Pt);
	T->Branch("pfmet_type1_PhotonEnUp_Phi" , &pfmet_type1_PhotonEnUp_Phi);
	T->Branch("pfmet_type1_PhotonEnUp_MT1" , &pfmet_type1_PhotonEnUp_MT1);
	T->Branch("pfmet_type1_PhotonEnUp_MT2" , &pfmet_type1_PhotonEnUp_MT2);
	T->Branch("pfmet_type1_PhotonEnDown_Pt" , &pfmet_type1_PhotonEnDown_Pt);
	T->Branch("pfmet_type1_PhotonEnDown_Phi" , &pfmet_type1_PhotonEnDown_Phi);
	T->Branch("pfmet_type1_PhotonEnDown_MT1" , &pfmet_type1_PhotonEnDown_MT1);
	T->Branch("pfmet_type1_PhotonEnDown_MT2" , &pfmet_type1_PhotonEnDown_MT2);

	
	T->Branch("pzetavis", &pzetavis);
	T->Branch("pzetamiss", &pzetamiss);
	T->Branch("pfpzetamiss", &pfpzetamiss);
	T->Branch("puppipzetamiss", &puppipzetamiss);
	T->Branch("pzetamiss_responseUP", &pzetamiss_responseUP);
	T->Branch("pzetamiss_responseDOWN", &pzetamiss_responseDOWN);
	T->Branch("pzetamiss_resolutionUP", &pzetamiss_resolutionUP);
	T->Branch("pzetamiss_resolutionDOWN", &pzetamiss_resolutionDOWN);
	T->Branch("pzetamiss_uncorr", &pzetamiss_uncorr);
	T->Branch("njets", &njets);
	T->Branch("njetspt20", &njetspt20);
	T->Branch("mjj", &mjj);
	T->Branch("jdeta", &jdeta);
	T->Branch("njetingap", &njetingap);
	T->Branch("njetingap20", &njetingap20);
	T->Branch("jdphi", &jdphi);
	T->Branch("jpt_1", &jpt_1);
	T->Branch("jeta_1", &jeta_1);
	T->Branch("jphi_1", &jphi_1);
	T->Branch("jm_1", &jm_1);
	T->Branch("jmva_1", &jmva_1);
	T->Branch("jpt_2", &jpt_2);
	T->Branch("jeta_2", &jeta_2);
	T->Branch("jphi_2", &jphi_2);
	T->Branch("jm_2", &jm_2);
	T->Branch("jmva_2", &jmva_2);
	T->Branch("nbtag", &nbtag);
	T->Branch("nbtag_oneSigmaUp", &nbtag_oneSigmaUp);
	T->Branch("nbtag_oneSigmaDown", &nbtag_oneSigmaDown);
	T->Branch("bpt_1", &bpt_1);
	T->Branch("beta_1", &beta_1);
	T->Branch("bphi_1", &bphi_1);
	T->Branch("bm_1", &bm_1);
	T->Branch("bmva_1", &bmva_1);
	T->Branch("bcsv_1", &bcsv_1);
	T->Branch("bpt_2", &bpt_2);
	T->Branch("beta_2", &beta_2);
	T->Branch("bphi_2", &bphi_2);
	T->Branch("bm_2", &bm_2);
	T->Branch("bmva_2", &bmva_2);
	T->Branch("bcsv_2", &bcsv_2);
	T->Branch("nbtag_LooseWp", &nbtag_LooseWp);
	T->Branch("nbtag_LooseWp_oneSigmaUp", &nbtag_LooseWp_oneSigmaUp);
	T->Branch("nbtag_LooseWp_oneSigmaDown", &nbtag_LooseWp_oneSigmaDown);
	T->Branch("bpt_1_LooseWp", &bpt_1_LooseWp);
	T->Branch("beta_1_LooseWp", &beta_1_LooseWp);
	T->Branch("bphi_1_LooseWp", &bphi_1_LooseWp);
	T->Branch("bm_1_LooseWp", &bm_1_LooseWp);
	T->Branch("bmva_1_LooseWp", &bmva_1_LooseWp);
	T->Branch("bcsv_1_LooseWp", &bcsv_1_LooseWp);
	T->Branch("bpt_2_LooseWp", &bpt_2_LooseWp);
	T->Branch("beta_2_LooseWp", &beta_2_LooseWp);
	T->Branch("bphi_2_LooseWp", &bphi_2_LooseWp);
	T->Branch("bm_2_LooseWp", &bm_2_LooseWp);
	T->Branch("bmva_2_LooseWp", &bmva_2_LooseWp);
	T->Branch("bcsv_2_LooseWp", &bcsv_2_LooseWp);
	T->Branch("nbtag_TightWp", &nbtag_TightWp);
	T->Branch("nbtag_TightWp_oneSigmaUp", &nbtag_TightWp_oneSigmaUp);
	T->Branch("nbtag_TightWp_oneSigmaDown", &nbtag_TightWp_oneSigmaDown);
	T->Branch("bpt_1_TightWp", &bpt_1_TightWp);
	T->Branch("beta_1_TightWp", &beta_1_TightWp);
	T->Branch("bphi_1_TightWp", &bphi_1_TightWp);
	T->Branch("bm_1_TightWp", &bm_1_TightWp);
	T->Branch("bmva_1_TightWp", &bmva_1_TightWp);
	T->Branch("bcsv_1_TightWp", &bcsv_1_TightWp);
	T->Branch("bpt_2_TightWp", &bpt_2_TightWp);
	T->Branch("beta_2_TightWp", &beta_2_TightWp);
	T->Branch("bphi_2_TightWp", &bphi_2_TightWp);
	T->Branch("bm_2_TightWp", &bm_2_TightWp);
	T->Branch("bmva_2_TightWp", &bmva_2_TightWp);
	T->Branch("bcsv_2_TightWp", &bcsv_2_TightWp);
	T->Branch("njets_JECshiftedUp", &njets_JECshiftedUp);
	T->Branch("njetspt20_JECshiftedUp", &njetspt20_JECshiftedUp);
	T->Branch("mjj_JECshiftedUp", &mjj_JECshiftedUp);
	T->Branch("jdeta_JECshiftedUp", &jdeta_JECshiftedUp);
	T->Branch("njetingap_JECshiftedUp", &njetingap_JECshiftedUp);
	T->Branch("njetingap20_JECshiftedUp", &njetingap20_JECshiftedUp);
	T->Branch("jdphi_JECshiftedUp", &jdphi_JECshiftedUp);
	T->Branch("jpt_1_JECshiftedUp", &jpt_1_JECshiftedUp);
	T->Branch("jeta_1_JECshiftedUp", &jeta_1_JECshiftedUp);
	T->Branch("jphi_1_JECshiftedUp", &jphi_1_JECshiftedUp);
	T->Branch("jm_1_JECshiftedUp", &jm_1_JECshiftedUp);
	T->Branch("jmva_1_JECshiftedUp", &jmva_1_JECshiftedUp);
	T->Branch("jpt_2_JECshiftedUp", &jpt_2_JECshiftedUp);
	T->Branch("jeta_2_JECshiftedUp", &jeta_2_JECshiftedUp);
	T->Branch("jphi_2_JECshiftedUp", &jphi_2_JECshiftedUp);
	T->Branch("jm_2_JECshiftedUp", &jm_2_JECshiftedUp);
	T->Branch("jmva_2_JECshiftedUp", &jmva_2_JECshiftedUp);
	T->Branch("nbtag_JECshiftedUp", &nbtag_JECshiftedUp);
	T->Branch("bpt_1_JECshiftedUp", &bpt_1_JECshiftedUp);
	T->Branch("beta_1_JECshiftedUp", &beta_1_JECshiftedUp);
	T->Branch("bphi_1_JECshiftedUp", &bphi_1_JECshiftedUp);
	T->Branch("bm_1_JECshiftedUp", &bm_1_JECshiftedUp);
	T->Branch("bmva_1_JECshiftedUp", &bmva_1_JECshiftedUp);
	T->Branch("bcsv_1_JECshiftedUp", &bcsv_1_JECshiftedUp);
	T->Branch("bpt_2_JECshiftedUp", &bpt_2_JECshiftedUp);
	T->Branch("beta_2_JECshiftedUp", &beta_2_JECshiftedUp);
	T->Branch("bphi_2_JECshiftedUp", &bphi_2_JECshiftedUp);
	T->Branch("bm_2_JECshiftedUp", &bm_2_JECshiftedUp);
	T->Branch("bmva_2_JECshiftedUp", &bmva_2_JECshiftedUp);
	T->Branch("bcsv_2_JECshiftedUp", &bcsv_2_JECshiftedUp);
	T->Branch("nbtag_LooseWp_JECshiftedUp", &nbtag_LooseWp_JECshiftedUp);
	T->Branch("bpt_1_LooseWp_JECshiftedUp", &bpt_1_LooseWp_JECshiftedUp);
	T->Branch("beta_1_LooseWp_JECshiftedUp", &beta_1_LooseWp_JECshiftedUp);
	T->Branch("bphi_1_LooseWp_JECshiftedUp", &bphi_1_LooseWp_JECshiftedUp);
	T->Branch("bm_1_LooseWp_JECshiftedUp", &bm_1_LooseWp_JECshiftedUp);
	T->Branch("bmva_1_LooseWp_JECshiftedUp", &bmva_1_LooseWp_JECshiftedUp);
	T->Branch("bcsv_1_LooseWp_JECshiftedUp", &bcsv_1_LooseWp_JECshiftedUp);
	T->Branch("bpt_2_LooseWp_JECshiftedUp", &bpt_2_LooseWp_JECshiftedUp);
	T->Branch("beta_2_LooseWp_JECshiftedUp", &beta_2_LooseWp_JECshiftedUp);
	T->Branch("bphi_2_LooseWp_JECshiftedUp", &bphi_2_LooseWp_JECshiftedUp);
	T->Branch("bm_2_LooseWp_JECshiftedUp", &bm_2_LooseWp_JECshiftedUp);
	T->Branch("bmva_2_LooseWp_JECshiftedUp", &bmva_2_LooseWp_JECshiftedUp);
	T->Branch("bcsv_2_LooseWp_JECshiftedUp", &bcsv_2_LooseWp_JECshiftedUp);
	T->Branch("nbtag_TightWp_JECshiftedUp", &nbtag_TightWp_JECshiftedUp);
	T->Branch("bpt_1_TightWp_JECshiftedUp", &bpt_1_TightWp_JECshiftedUp);
	T->Branch("beta_1_TightWp_JECshiftedUp", &beta_1_TightWp_JECshiftedUp);
	T->Branch("bphi_1_TightWp_JECshiftedUp", &bphi_1_TightWp_JECshiftedUp);
	T->Branch("bm_1_TightWp_JECshiftedUp", &bm_1_TightWp_JECshiftedUp);
	T->Branch("bmva_1_TightWp_JECshiftedUp", &bmva_1_TightWp_JECshiftedUp);
	T->Branch("bcsv_1_TightWp_JECshiftedUp", &bcsv_1_TightWp_JECshiftedUp);
	T->Branch("bpt_2_TightWp_JECshiftedUp", &bpt_2_TightWp_JECshiftedUp);
	T->Branch("beta_2_TightWp_JECshiftedUp", &beta_2_TightWp_JECshiftedUp);
	T->Branch("bphi_2_TightWp_JECshiftedUp", &bphi_2_TightWp_JECshiftedUp);
	T->Branch("bm_2_TightWp_JECshiftedUp", &bm_2_TightWp_JECshiftedUp);
	T->Branch("bmva_2_TightWp_JECshiftedUp", &bmva_2_TightWp_JECshiftedUp);
	T->Branch("bcsv_2_TightWp_JECshiftedUp", &bcsv_2_TightWp_JECshiftedUp);
	T->Branch("njets_JECshiftedDown", &njets_JECshiftedDown);
	T->Branch("njetspt20_JECshiftedDown", &njetspt20_JECshiftedDown);
	T->Branch("mjj_JECshiftedDown", &mjj_JECshiftedDown);
	T->Branch("jdeta_JECshiftedDown", &jdeta_JECshiftedDown);
	T->Branch("njetingap_JECshiftedDown", &njetingap_JECshiftedDown);
	T->Branch("njetingap20_JECshiftedDown", &njetingap20_JECshiftedDown);
	T->Branch("jdphi_JECshiftedDown", &jdphi_JECshiftedDown);
	T->Branch("jpt_1_JECshiftedDown", &jpt_1_JECshiftedDown);
	T->Branch("jeta_1_JECshiftedDown", &jeta_1_JECshiftedDown);
	T->Branch("jphi_1_JECshiftedDown", &jphi_1_JECshiftedDown);
	T->Branch("jm_1_JECshiftedDown", &jm_1_JECshiftedDown);
	T->Branch("jmva_1_JECshiftedDown", &jmva_1_JECshiftedDown);
	T->Branch("jpt_2_JECshiftedDown", &jpt_2_JECshiftedDown);
	T->Branch("jeta_2_JECshiftedDown", &jeta_2_JECshiftedDown);
	T->Branch("jphi_2_JECshiftedDown", &jphi_2_JECshiftedDown);
	T->Branch("jm_2_JECshiftedDown", &jm_2_JECshiftedDown);
	T->Branch("jmva_2_JECshiftedDown", &jmva_2_JECshiftedDown);
	T->Branch("nbtag_JECshiftedDown", &nbtag_JECshiftedDown);
	T->Branch("bpt_1_JECshiftedDown", &bpt_1_JECshiftedDown);
	T->Branch("beta_1_JECshiftedDown", &beta_1_JECshiftedDown);
	T->Branch("bphi_1_JECshiftedDown", &bphi_1_JECshiftedDown);
	T->Branch("bm_1_JECshiftedDown", &bm_1_JECshiftedDown);
	T->Branch("bmva_1_JECshiftedDown", &bmva_1_JECshiftedDown);
	T->Branch("bcsv_1_JECshiftedDown", &bcsv_1_JECshiftedDown);
	T->Branch("bpt_2_JECshiftedDown", &bpt_2_JECshiftedDown);
	T->Branch("beta_2_JECshiftedDown", &beta_2_JECshiftedDown);
	T->Branch("bphi_2_JECshiftedDown", &bphi_2_JECshiftedDown);
	T->Branch("bm_2_JECshiftedDown", &bm_2_JECshiftedDown);
	T->Branch("bmva_2_JECshiftedDown", &bmva_2_JECshiftedDown);
	T->Branch("bcsv_2_JECshiftedDown", &bcsv_2_JECshiftedDown);
	T->Branch("nbtag_LooseWp_JECshiftedDown", &nbtag_LooseWp_JECshiftedDown);
	T->Branch("bpt_1_LooseWp_JECshiftedDown", &bpt_1_LooseWp_JECshiftedDown);
	T->Branch("beta_1_LooseWp_JECshiftedDown", &beta_1_LooseWp_JECshiftedDown);
	T->Branch("bphi_1_LooseWp_JECshiftedDown", &bphi_1_LooseWp_JECshiftedDown);
	T->Branch("bm_1_LooseWp_JECshiftedDown", &bm_1_LooseWp_JECshiftedDown);
	T->Branch("bmva_1_LooseWp_JECshiftedDown", &bmva_1_LooseWp_JECshiftedDown);
	T->Branch("bcsv_1_LooseWp_JECshiftedDown", &bcsv_1_LooseWp_JECshiftedDown);
	T->Branch("bpt_2_LooseWp_JECshiftedDown", &bpt_2_LooseWp_JECshiftedDown);
	T->Branch("beta_2_LooseWp_JECshiftedDown", &beta_2_LooseWp_JECshiftedDown);
	T->Branch("bphi_2_LooseWp_JECshiftedDown", &bphi_2_LooseWp_JECshiftedDown);
	T->Branch("bm_2_LooseWp_JECshiftedDown", &bm_2_LooseWp_JECshiftedDown);
	T->Branch("bmva_2_LooseWp_JECshiftedDown", &bmva_2_LooseWp_JECshiftedDown);
	T->Branch("bcsv_2_LooseWp_JECshiftedDown", &bcsv_2_LooseWp_JECshiftedDown);
	T->Branch("nbtag_TightWp_JECshiftedDown", &nbtag_TightWp_JECshiftedDown);
	T->Branch("bpt_1_TightWp_JECshiftedDown", &bpt_1_TightWp_JECshiftedDown);
	T->Branch("beta_1_TightWp_JECshiftedDown", &beta_1_TightWp_JECshiftedDown);
	T->Branch("bphi_1_TightWp_JECshiftedDown", &bphi_1_TightWp_JECshiftedDown);
	T->Branch("bm_1_TightWp_JECshiftedDown", &bm_1_TightWp_JECshiftedDown);
	T->Branch("bmva_1_TightWp_JECshiftedDown", &bmva_1_TightWp_JECshiftedDown);
	T->Branch("bcsv_1_TightWp_JECshiftedDown", &bcsv_1_TightWp_JECshiftedDown);
	T->Branch("bpt_2_TightWp_JECshiftedDown", &bpt_2_TightWp_JECshiftedDown);
	T->Branch("beta_2_TightWp_JECshiftedDown", &beta_2_TightWp_JECshiftedDown);
	T->Branch("bphi_2_TightWp_JECshiftedDown", &bphi_2_TightWp_JECshiftedDown);
	T->Branch("bm_2_TightWp_JECshiftedDown", &bm_2_TightWp_JECshiftedDown);
	T->Branch("bmva_2_TightWp_JECshiftedDown", &bmva_2_TightWp_JECshiftedDown);
	T->Branch("bcsv_2_TightWp_JECshiftedDown", &bcsv_2_TightWp_JECshiftedDown);
	T->Branch("njets_JERup", &njets_JERup);
	T->Branch("njetspt20_JERup", &njetspt20_JERup);
	T->Branch("mjj_JERup", &mjj_JERup);
	T->Branch("jdeta_JERup", &jdeta_JERup);
	T->Branch("njetingap_JERup", &njetingap_JERup);
	T->Branch("njetingap20_JERup", &njetingap20_JERup);
	T->Branch("jdphi_JERup", &jdphi_JERup);
	T->Branch("jpt_1_JERup", &jpt_1_JERup);
	T->Branch("jeta_1_JERup", &jeta_1_JERup);
	T->Branch("jphi_1_JERup", &jphi_1_JERup);
	T->Branch("jm_1_JERup", &jm_1_JERup);
	T->Branch("jmva_1_JERup", &jmva_1_JERup);
	T->Branch("jpt_2_JERup", &jpt_2_JERup);
	T->Branch("jeta_2_JERup", &jeta_2_JERup);
	T->Branch("jphi_2_JERup", &jphi_2_JERup);
	T->Branch("jm_2_JERup", &jm_2_JERup);
	T->Branch("jmva_2_JERup", &jmva_2_JERup);
	T->Branch("nbtag_JERup", &nbtag_JERup);
	T->Branch("bpt_1_JERup", &bpt_1_JERup);
	T->Branch("beta_1_JERup", &beta_1_JERup);
	T->Branch("bphi_1_JERup", &bphi_1_JERup);
	T->Branch("bm_1_JERup", &bm_1_JERup);
	T->Branch("bmva_1_JERup", &bmva_1_JERup);
	T->Branch("bcsv_1_JERup", &bcsv_1_JERup);
	T->Branch("bpt_2_JERup", &bpt_2_JERup);
	T->Branch("beta_2_JERup", &beta_2_JERup);
	T->Branch("bphi_2_JERup", &bphi_2_JERup);
	T->Branch("bm_2_JERup", &bm_2_JERup);
	T->Branch("bmva_2_JERup", &bmva_2_JERup);
	T->Branch("bcsv_2_JERup", &bcsv_2_JERup);
	T->Branch("nbtag_LooseWp_JERup", &nbtag_LooseWp_JERup);
	T->Branch("bpt_1_LooseWp_JERup", &bpt_1_LooseWp_JERup);
	T->Branch("beta_1_LooseWp_JERup", &beta_1_LooseWp_JERup);
	T->Branch("bphi_1_LooseWp_JERup", &bphi_1_LooseWp_JERup);
	T->Branch("bm_1_LooseWp_JERup", &bm_1_LooseWp_JERup);
	T->Branch("bmva_1_LooseWp_JERup", &bmva_1_LooseWp_JERup);
	T->Branch("bcsv_1_LooseWp_JERup", &bcsv_1_LooseWp_JERup);
	T->Branch("bpt_2_LooseWp_JERup", &bpt_2_LooseWp_JERup);
	T->Branch("beta_2_LooseWp_JERup", &beta_2_LooseWp_JERup);
	T->Branch("bphi_2_LooseWp_JERup", &bphi_2_LooseWp_JERup);
	T->Branch("bm_2_LooseWp_JERup", &bm_2_LooseWp_JERup);
	T->Branch("bmva_2_LooseWp_JERup", &bmva_2_LooseWp_JERup);
	T->Branch("bcsv_2_LooseWp_JERup", &bcsv_2_LooseWp_JERup);
	T->Branch("nbtag_TightWp_JERup", &nbtag_TightWp_JERup);
	T->Branch("bpt_1_TightWp_JERup", &bpt_1_TightWp_JERup);
	T->Branch("beta_1_TightWp_JERup", &beta_1_TightWp_JERup);
	T->Branch("bphi_1_TightWp_JERup", &bphi_1_TightWp_JERup);
	T->Branch("bm_1_TightWp_JERup", &bm_1_TightWp_JERup);
	T->Branch("bmva_1_TightWp_JERup", &bmva_1_TightWp_JERup);
	T->Branch("bcsv_1_TightWp_JERup", &bcsv_1_TightWp_JERup);
	T->Branch("bpt_2_TightWp_JERup", &bpt_2_TightWp_JERup);
	T->Branch("beta_2_TightWp_JERup", &beta_2_TightWp_JERup);
	T->Branch("bphi_2_TightWp_JERup", &bphi_2_TightWp_JERup);
	T->Branch("bm_2_TightWp_JERup", &bm_2_TightWp_JERup);
	T->Branch("bmva_2_TightWp_JERup", &bmva_2_TightWp_JERup);
	T->Branch("bcsv_2_TightWp_JERup", &bcsv_2_TightWp_JERup);
	T->Branch("njets_JERdown", &njets_JERdown);
	T->Branch("njetspt20_JERdown", &njetspt20_JERdown);
	T->Branch("mjj_JERdown", &mjj_JERdown);
	T->Branch("jdeta_JERdown", &jdeta_JERdown);
	T->Branch("njetingap_JERdown", &njetingap_JERdown);
	T->Branch("njetingap20_JERdown", &njetingap20_JERdown);
	T->Branch("jdphi_JERdown", &jdphi_JERdown);
	T->Branch("jpt_1_JERdown", &jpt_1_JERdown);
	T->Branch("jeta_1_JERdown", &jeta_1_JERdown);
	T->Branch("jphi_1_JERdown", &jphi_1_JERdown);
	T->Branch("jm_1_JERdown", &jm_1_JERdown);
	T->Branch("jmva_1_JERdown", &jmva_1_JERdown);
	T->Branch("jpt_2_JERdown", &jpt_2_JERdown);
	T->Branch("jeta_2_JERdown", &jeta_2_JERdown);
	T->Branch("jphi_2_JERdown", &jphi_2_JERdown);
	T->Branch("jm_2_JERdown", &jm_2_JERdown);
	T->Branch("jmva_2_JERdown", &jmva_2_JERdown);
	T->Branch("nbtag_JERdown", &nbtag_JERdown);
	T->Branch("bpt_1_JERdown", &bpt_1_JERdown);
	T->Branch("beta_1_JERdown", &beta_1_JERdown);
	T->Branch("bphi_1_JERdown", &bphi_1_JERdown);
	T->Branch("bm_1_JERdown", &bm_1_JERdown);
	T->Branch("bmva_1_JERdown", &bmva_1_JERdown);
	T->Branch("bcsv_1_JERdown", &bcsv_1_JERdown);
	T->Branch("bpt_2_JERdown", &bpt_2_JERdown);
	T->Branch("beta_2_JERdown", &beta_2_JERdown);
	T->Branch("bphi_2_JERdown", &bphi_2_JERdown);
	T->Branch("bm_2_JERdown", &bm_2_JERdown);
	T->Branch("bmva_2_JERdown", &bmva_2_JERdown);
	T->Branch("bcsv_2_JERdown", &bcsv_2_JERdown);
	T->Branch("nbtag_LooseWp_JERdown", &nbtag_LooseWp_JERdown);
	T->Branch("bpt_1_LooseWp_JERdown", &bpt_1_LooseWp_JERdown);
	T->Branch("beta_1_LooseWp_JERdown", &beta_1_LooseWp_JERdown);
	T->Branch("bphi_1_LooseWp_JERdown", &bphi_1_LooseWp_JERdown);
	T->Branch("bm_1_LooseWp_JERdown", &bm_1_LooseWp_JERdown);
	T->Branch("bmva_1_LooseWp_JERdown", &bmva_1_LooseWp_JERdown);
	T->Branch("bcsv_1_LooseWp_JERdown", &bcsv_1_LooseWp_JERdown);
	T->Branch("bpt_2_LooseWp_JERdown", &bpt_2_LooseWp_JERdown);
	T->Branch("beta_2_LooseWp_JERdown", &beta_2_LooseWp_JERdown);
	T->Branch("bphi_2_LooseWp_JERdown", &bphi_2_LooseWp_JERdown);
	T->Branch("bm_2_LooseWp_JERdown", &bm_2_LooseWp_JERdown);
	T->Branch("bmva_2_LooseWp_JERdown", &bmva_2_LooseWp_JERdown);
	T->Branch("bcsv_2_LooseWp_JERdown", &bcsv_2_LooseWp_JERdown);
	T->Branch("nbtag_TightWp_JERdown", &nbtag_TightWp_JERdown);
	T->Branch("bpt_1_TightWp_JERdown", &bpt_1_TightWp_JERdown);
	T->Branch("beta_1_TightWp_JERdown", &beta_1_TightWp_JERdown);
	T->Branch("bphi_1_TightWp_JERdown", &bphi_1_TightWp_JERdown);
	T->Branch("bm_1_TightWp_JERdown", &bm_1_TightWp_JERdown);
	T->Branch("bmva_1_TightWp_JERdown", &bmva_1_TightWp_JERdown);
	T->Branch("bcsv_1_TightWp_JERdown", &bcsv_1_TightWp_JERdown);
	T->Branch("bpt_2_TightWp_JERdown", &bpt_2_TightWp_JERdown);
	T->Branch("beta_2_TightWp_JERdown", &beta_2_TightWp_JERdown);
	T->Branch("bphi_2_TightWp_JERdown", &bphi_2_TightWp_JERdown);
	T->Branch("bm_2_TightWp_JERdown", &bm_2_TightWp_JERdown);
	T->Branch("bmva_2_TightWp_JERdown", &bmva_2_TightWp_JERdown);
	T->Branch("bcsv_2_TightWp_JERdown", &bcsv_2_TightWp_JERdown);
	T->Branch("dilepton_veto", &dilepton_veto);
	T->Branch("extraelec_veto", &extraelec_veto);
	T->Branch("extramuon_veto", &extramuon_veto);
	T->Branch("HBHENoiseFilter", &HBHENoiseFilter);
	T->Branch("HBHENoiseIsoFilter", &HBHENoiseIsoFilter);
	T->Branch("CSCTightHalo2015Filter", &CSCTightHalo2015Filter);
	T->Branch("EcalDeadCellTriggerPrimitiveFilter", &EcalDeadCellTriggerPrimitiveFilter);
	T->Branch("goodVerticesFilter", &goodVerticesFilter);
	T->Branch("eeBadScFilter", &eeBadScFilter);
	T->Branch("chargedHadronTrackResolutionFilter", &chargedHadronTrackResolutionFilter);
	T->Branch("muonBadTrackFilter", &muonBadTrackFilter);
	T->Branch("globalTightHalo2016Filter", &globalTightHalo2016Filter);
	T->Branch("BadChargedCandidateFilter", &BadChargedCandidateFilter);
	T->Branch("BadPFMuonFilter", &BadPFMuonFilter);

    T->Branch("BadMuonTaggedMoriond17", &BadMuonTaggedMoriond17);
    T->Branch("DuplicateMuonTaggedMoriond17", &DuplicateMuonTaggedMoriond17);

    T->Branch("BtagEventSFproduct_looseWpDown",  &BtagEventSFproduct_looseWpDown);
    T->Branch("BtagEventSFproduct_looseWpCentral",  &BtagEventSFproduct_looseWpCentral);
    T->Branch("BtagEventSFproduct_looseWpUp",  &BtagEventSFproduct_looseWpUp);
    T->Branch("BtagEventSFproduct_mediumWpDown",  &BtagEventSFproduct_mediumWpDown);
    T->Branch("BtagEventSFproduct_mediumWpCentral",  &BtagEventSFproduct_mediumWpCentral);
    T->Branch("BtagEventSFproduct_mediumWpUp",  &BtagEventSFproduct_mediumWpUp);
    T->Branch("BtagEventSFproduct_tightWpDown",  &BtagEventSFproduct_tightWpDown);
    T->Branch("BtagEventSFproduct_tightWpCentral",  &BtagEventSFproduct_tightWpCentral);
    T->Branch("BtagEventSFproduct_tightWpUp",  &BtagEventSFproduct_tightWpUp);
    T->Branch("BtagEventSFproduct_looseWpCentral_JECshiftedUp",  &BtagEventSFproduct_looseWpCentral_JECshiftedUp);
    T->Branch("BtagEventSFproduct_mediumWpCentral_JECshiftedUp",  &BtagEventSFproduct_mediumWpCentral_JECshiftedUp);
    T->Branch("BtagEventSFproduct_tightWpCentral_JECshiftedUp",  &BtagEventSFproduct_tightWpCentral_JECshiftedUp);
    T->Branch("BtagEventSFproduct_looseWpCentral_JECshiftedDown",  &BtagEventSFproduct_looseWpCentral_JECshiftedDown);
    T->Branch("BtagEventSFproduct_mediumWpCentral_JECshiftedDown",  &BtagEventSFproduct_mediumWpCentral_JECshiftedDown);
    T->Branch("BtagEventSFproduct_tightWpCentral_JECshiftedDown",  &BtagEventSFproduct_tightWpCentral_JECshiftedDown);
    T->Branch("BtagEventSFproduct_looseWpCentral_JERup",  &BtagEventSFproduct_looseWpCentral_JERup);
    T->Branch("BtagEventSFproduct_mediumWpCentral_JERup",  &BtagEventSFproduct_mediumWpCentral_JERup);
    T->Branch("BtagEventSFproduct_tightWpCentral_JERup",  &BtagEventSFproduct_tightWpCentral_JERup);
    T->Branch("BtagEventSFproduct_looseWpCentral_JERdown",  &BtagEventSFproduct_looseWpCentral_JERdown);
    T->Branch("BtagEventSFproduct_mediumWpCentral_JERdown",  &BtagEventSFproduct_mediumWpCentral_JERdown);
    T->Branch("BtagEventSFproduct_tightWpCentral_JERdown",  &BtagEventSFproduct_tightWpCentral_JERdown);
    
    T->Branch("BtagEventSFproduct_or_DataTag_Central",  &BtagEventSFproduct_or_DataTag_Central);
    T->Branch("BtagEventSFproduct_or_DataTag_Up",  &BtagEventSFproduct_or_DataTag_Up);
    T->Branch("BtagEventSFproduct_or_DataTag_Down",  &BtagEventSFproduct_or_DataTag_Down);
    T->Branch("BtagEventSFproduct_or_DataTag_Central_JECshiftedUp",  &BtagEventSFproduct_or_DataTag_Central_JECshiftedUp);
    T->Branch("BtagEventSFproduct_or_DataTag_Central_JECshiftedDown",  &BtagEventSFproduct_or_DataTag_Central_JECshiftedDown);

	T->Branch("NUP", &NUP);
	T->Branch("weight", &weight);
	T->Branch("generatorEventWeight", &generatorEventWeight);
	T->Branch("lheHT", &lheHT);
	T->Branch("lheOutGoingPartons", &lheOutGoingPartons);
	T->Branch("lheZmass", &lheZmass);
	T->Branch("IsZTT", &IsZTT);
	T->Branch("IsZL", &IsZL);
	T->Branch("IsZJ", &IsZJ);
	T->Branch("IsZLL", &IsZLL);
    T->Branch("IsTTT",&IsTTT);
	T->Branch("DataSet=", &DataSet);
	T->Branch("EventTotal", &EventTotal);
	T->Branch("EventType=", &EventType);
	T->Branch("KeyName=", &KeyName);
	T->Branch("DataCard=", &DataCard);
	T->Branch("CrossSection", &CrossSection);
	T->Branch("FilterEff", &FilterEff);
	T->Branch("isSmallTree", &isSmallTree);
	T->Branch("TauEsNumberSigmasShifted", &TauEsNumberSigmasShifted);
	T->Branch("ElectronEsNumberSigmasShifted", &ElectronEsNumberSigmasShifted);

	T->Branch("leg1_HLT_DoubleMediumIsoPFTau35_Trk1_eta2p1_Reg", &leg1_HLT_DoubleMediumIsoPFTau35_Trk1_eta2p1_Reg);
    T->Branch("leg1_HLT_DoubleMediumCombinedIsoPFTau35_Trk1_eta2p1_Reg", &leg1_HLT_DoubleMediumCombinedIsoPFTau35_Trk1_eta2p1_Reg);
	T->Branch("leg1_HLT_Ele25_eta2p1_WPTight_Gsf", &leg1_HLT_Ele25_eta2p1_WPTight_Gsf);
    T->Branch("leg1_HLT_Ele45_WPLoose_Gsf_L1JetTauSeeded", &leg1_HLT_Ele45_WPLoose_Gsf_L1JetTauSeeded);
	T->Branch("leg1_HLT_IsoMu22", &leg1_HLT_IsoMu22);
	T->Branch("leg1_HLT_IsoTkMu22", &leg1_HLT_IsoTkMu22);
    T->Branch("leg1_HLT_IsoMu22_eta2p1", &leg1_HLT_IsoMu22_eta2p1);
    T->Branch("leg1_HLT_IsoTkMu22_eta2p1", &leg1_HLT_IsoTkMu22_eta2p1);
    T->Branch("leg1_HLT_IsoMu24", &leg1_HLT_IsoMu24);
    T->Branch("leg1_HLT_IsoTkMu24", &leg1_HLT_IsoTkMu24);
    
	T->Branch("leg2_HLT_DoubleMediumIsoPFTau35_Trk1_eta2p1_Reg", &leg2_HLT_DoubleMediumIsoPFTau35_Trk1_eta2p1_Reg);
    T->Branch("leg2_HLT_DoubleMediumCombinedIsoPFTau35_Trk1_eta2p1_Reg", &leg2_HLT_DoubleMediumCombinedIsoPFTau35_Trk1_eta2p1_Reg);
	T->Branch("leg2_HLT_Ele25_eta2p1_WPTight_Gsf", &leg2_HLT_Ele25_eta2p1_WPTight_Gsf);
    T->Branch("leg2_HLT_Ele45_WPLoose_Gsf_L1JetTauSeeded", &leg2_HLT_Ele45_WPLoose_Gsf_L1JetTauSeeded);
	T->Branch("leg2_HLT_IsoMu22", &leg2_HLT_IsoMu22);
	T->Branch("leg2_HLT_IsoTkMu22", &leg2_HLT_IsoTkMu22);
    T->Branch("leg2_HLT_IsoMu22_eta2p1", &leg2_HLT_IsoMu22_eta2p1);
    T->Branch("leg2_HLT_IsoTkMu22_eta2p1", &leg2_HLT_IsoTkMu22_eta2p1);
    T->Branch("leg2_HLT_IsoMu24", &leg2_HLT_IsoMu24);
    T->Branch("leg2_HLT_IsoTkMu24", &leg2_HLT_IsoTkMu24);

    T->Branch("pairGoodForTrigger", &pairGoodForTrigger);

	T->Branch("veto_leptonType", &veto_leptonType);
	T->Branch("veto_pt", &veto_pt);
	T->Branch("veto_eta", &veto_eta);
	T->Branch("veto_phi", &veto_phi);
	T->Branch("veto_M", &veto_M);
	T->Branch("veto_charge", &veto_charge);
	T->Branch("veto_dxy", &veto_dxy);
	T->Branch("veto_dz", &veto_dz);
	T->Branch("veto_RelIso", &veto_RelIso);
    T->Branch("veto_passesLooseMuonId", &veto_passesLooseMuonId);
	T->Branch("veto_passesMediumMuonId", &veto_passesMediumMuonId);
    T->Branch("veto_passesMediumMuonId_ICHEP16", &veto_passesMediumMuonId_ICHEP16);
    T->Branch("veto_passesMediumMuonId_Moriond17", &veto_passesMediumMuonId_Moriond17);
    T->Branch("veto_passesTightMuonId", &veto_passesTightMuonId);
	T->Branch("veto_passElectronMVA80", &veto_passElectronMVA80);
	T->Branch("veto_passElectronMVA90", &veto_passElectronMVA90);
	T->Branch("veto_passVetoElectronCutBased", &veto_passVetoElectronCutBased);
    T->Branch("veto_passTightElectronCutBased", &veto_passTightElectronCutBased);
	T->Branch("veto_isTrackerGlobalPFMuon", &veto_isTrackerGlobalPFMuon);
	T->Branch("veto_numberOfMissingInnerHits", &veto_numberOfMissingInnerHits);
	T->Branch("veto_numberOfMissingOuterHits", &veto_numberOfMissingOuterHits);
	T->Branch("veto_passConversionVeto", &veto_passConversionVeto);
	T->Branch("veto_LeptonPassesThirdElectronVetoCuts", &veto_LeptonPassesThirdElectronVetoCuts);
	T->Branch("veto_LeptonPassesThirdMuonVetoCuts", &veto_LeptonPassesThirdMuonVetoCuts);
	T->Branch("veto_LeptonPassesDiElectronVetoCuts", &veto_LeptonPassesDiElectronVetoCuts);
	T->Branch("veto_LeptonPassesDiMuonVetoCuts", &veto_LeptonPassesDiMuonVetoCuts);
	T->Branch("LPT", &LPT);
	T->Branch("P_chi", &P_chi);
	T->Branch("M_min", &M_min);
	T->Branch("P_chi_pf", &P_chi_pf);
    T->Branch("P_chi_pf_UESUp", &P_chi_pf_UESUp);
    T->Branch("P_chi_pf_UESDown", &P_chi_pf_UESDown);
	T->Branch("M_min_pf", &M_min_pf);
	T->Branch("P_chi_puppi", &P_chi_puppi);
	T->Branch("M_min_puppi", &M_min_puppi);
	T->Branch("P_chi_uncorr", &P_chi_uncorr);
	T->Branch("M_min_uncorr", &M_min_uncorr);
	T->Branch("P_chi_responseUP", &P_chi_responseUP);
	T->Branch("M_min_responseUP", &M_min_responseUP);
	T->Branch("P_chi_responseDOWN", &P_chi_responseDOWN);
	T->Branch("M_min_responseDOWN", &M_min_responseDOWN);
	T->Branch("P_chi_resolutionUP", &P_chi_resolutionUP);
	T->Branch("M_min_resolutionUP", &M_min_resolutionUP);
	T->Branch("P_chi_resolutionDOWN", &P_chi_resolutionDOWN);
	T->Branch("M_min_resolutionDOWN", &M_min_resolutionDOWN);

    T->Branch("mvaVar_mt_MZP600A0400", &mvaVar_mt_MZP600A0400);
    T->Branch("mvaVar_mt_MZP800A0400", &mvaVar_mt_MZP800A0400);
    T->Branch("mvaVar_mt_MZP1000A0400", &mvaVar_mt_MZP1000A0400);
    T->Branch("mvaVar_mt_MZP1200A0400", &mvaVar_mt_MZP1200A0400);
    
    T->Branch("mvaVar_et_MZP600A0400", &mvaVar_et_MZP600A0400);
    T->Branch("mvaVar_et_MZP800A0400", &mvaVar_et_MZP800A0400);
    T->Branch("mvaVar_et_MZP1000A0400", &mvaVar_et_MZP1000A0400);
    T->Branch("mvaVar_et_MZP1200A0400", &mvaVar_et_MZP1200A0400);
    
    T->Branch("mvaVar_tt_MZP600A0400", &mvaVar_tt_MZP600A0400);
    T->Branch("mvaVar_tt_MZP800A0400", &mvaVar_tt_MZP800A0400);
    T->Branch("mvaVar_tt_MZP1000A0400", &mvaVar_tt_MZP1000A0400);
    T->Branch("mvaVar_tt_MZP1200A0400", &mvaVar_tt_MZP1200A0400);
    
    T->Branch("mvaVar_mt_MZP600A0400_UESUp", &mvaVar_mt_MZP600A0400_UESUp);
    T->Branch("mvaVar_mt_MZP800A0400_UESUp", &mvaVar_mt_MZP800A0400_UESUp);
    T->Branch("mvaVar_mt_MZP1000A0400_UESUp", &mvaVar_mt_MZP1000A0400_UESUp);
    T->Branch("mvaVar_mt_MZP1200A0400_UESUp", &mvaVar_mt_MZP1200A0400_UESUp);
    
    T->Branch("mvaVar_et_MZP600A0400_UESUp", &mvaVar_et_MZP600A0400_UESUp);
    T->Branch("mvaVar_et_MZP800A0400_UESUp", &mvaVar_et_MZP800A0400_UESUp);
    T->Branch("mvaVar_et_MZP1000A0400_UESUp", &mvaVar_et_MZP1000A0400_UESUp);
    T->Branch("mvaVar_et_MZP1200A0400_UESUp", &mvaVar_et_MZP1200A0400_UESUp);
    
    T->Branch("mvaVar_tt_MZP600A0400_UESUp", &mvaVar_tt_MZP600A0400_UESUp);
    T->Branch("mvaVar_tt_MZP800A0400_UESUp", &mvaVar_tt_MZP800A0400_UESUp);
    T->Branch("mvaVar_tt_MZP1000A0400_UESUp", &mvaVar_tt_MZP1000A0400_UESUp);
    T->Branch("mvaVar_tt_MZP1200A0400_UESUp", &mvaVar_tt_MZP1200A0400_UESUp);
    
    T->Branch("mvaVar_mt_MZP600A0400_UESDown", &mvaVar_mt_MZP600A0400_UESDown);
    T->Branch("mvaVar_mt_MZP800A0400_UESDown", &mvaVar_mt_MZP800A0400_UESDown);
    T->Branch("mvaVar_mt_MZP1000A0400_UESDown", &mvaVar_mt_MZP1000A0400_UESDown);
    T->Branch("mvaVar_mt_MZP1200A0400_UESDown", &mvaVar_mt_MZP1200A0400_UESDown);
    
    T->Branch("mvaVar_et_MZP600A0400_UESDown", &mvaVar_et_MZP600A0400_UESDown);
    T->Branch("mvaVar_et_MZP800A0400_UESDown", &mvaVar_et_MZP800A0400_UESDown);
    T->Branch("mvaVar_et_MZP1000A0400_UESDown", &mvaVar_et_MZP1000A0400_UESDown);
    T->Branch("mvaVar_et_MZP1200A0400_UESDown", &mvaVar_et_MZP1200A0400_UESDown);
    
    T->Branch("mvaVar_tt_MZP600A0400_UESDown", &mvaVar_tt_MZP600A0400_UESDown);
    T->Branch("mvaVar_tt_MZP800A0400_UESDown", &mvaVar_tt_MZP800A0400_UESDown);
    T->Branch("mvaVar_tt_MZP1000A0400_UESDown", &mvaVar_tt_MZP1000A0400_UESDown);
    T->Branch("mvaVar_tt_MZP1200A0400_UESDown", &mvaVar_tt_MZP1200A0400_UESDown);
    
    T->Branch("pt_1_TESUp", &pt_1_TESUp);
    T->Branch("mt_1_TESUp", &mt_1_TESUp);
    T->Branch("pfmt_1_TESUp", &pfmt_1_TESUp);

    T->Branch("pt_2_TESUp", &pt_2_TESUp);
    T->Branch("mt_2_TESUp", &mt_2_TESUp);
    T->Branch("pfmt_2_TESUp", &pfmt_2_TESUp);

    T->Branch("mt_tot_TESUp", &mt_tot_TESUp);
    T->Branch("pt_tt_TESUp", &pt_tt_TESUp);
    T->Branch("m_vis_TESUp", &m_vis_TESUp);
    T->Branch("pfmet_type1_TESUp_Pt", &pfmet_type1_TESUp_Pt);
    T->Branch("P_chi_pf_TESUp", &P_chi_pf_TESUp);
    T->Branch("LPT_TESUp", &LPT_TESUp);
    T->Branch("DeltaR_leg1_leg2", &DeltaR_leg1_leg2);
    T->Branch("cos_DeltaPhi_leg1_leg2", &cos_DeltaPhi_leg1_leg2);
    T->Branch("cos_DeltaPhi_PFMET_Higgs_TESUp", &cos_DeltaPhi_PFMET_Higgs_TESUp);

    T->Branch("pt_1_TESDown", &pt_1_TESDown);
    T->Branch("mt_1_TESDown", &mt_1_TESDown);
    T->Branch("pfmt_1_TESDown", &pfmt_1_TESDown);

    T->Branch("pt_2_TESDown", &pt_2_TESDown);
    T->Branch("mt_2_TESDown", &mt_2_TESDown);
    T->Branch("pfmt_2_TESDown", &pfmt_2_TESDown);

    T->Branch("mt_tot_TESDown", &mt_tot_TESDown);
    T->Branch("pt_tt_TESDown", &pt_tt_TESDown);
    T->Branch("m_vis_TESDown", &m_vis_TESDown);
    T->Branch("pfmet_type1_TESDown_Pt", &pfmet_type1_TESDown_Pt);
    T->Branch("P_chi_pf_TESDown", &P_chi_pf_TESDown);
    T->Branch("LPT_TESDown", &LPT_TESDown);
    T->Branch("DeltaR_leg1_leg2", &DeltaR_leg1_leg2);
    T->Branch("cos_DeltaPhi_leg1_leg2", &cos_DeltaPhi_leg1_leg2);
    T->Branch("cos_DeltaPhi_PFMET_Higgs_TESDown", &cos_DeltaPhi_PFMET_Higgs_TESDown);

    T->Branch("mvaVar_mt_MZP600A0400_TESUp", &mvaVar_mt_MZP600A0400_TESUp);
    T->Branch("mvaVar_mt_MZP800A0400_TESUp", &mvaVar_mt_MZP800A0400_TESUp);
    T->Branch("mvaVar_mt_MZP1000A0400_TESUp", &mvaVar_mt_MZP1000A0400_TESUp);
    T->Branch("mvaVar_mt_MZP1200A0400_TESUp", &mvaVar_mt_MZP1200A0400_TESUp);

    T->Branch("mvaVar_et_MZP600A0400_TESUp", &mvaVar_et_MZP600A0400_TESUp);
    T->Branch("mvaVar_et_MZP800A0400_TESUp", &mvaVar_et_MZP800A0400_TESUp); 
    T->Branch("mvaVar_et_MZP1000A0400_TESUp", &mvaVar_et_MZP1000A0400_TESUp);
    T->Branch("mvaVar_et_MZP1200A0400_TESUp", &mvaVar_et_MZP1200A0400_TESUp);
        
    T->Branch("mvaVar_tt_MZP600A0400_TESUp", &mvaVar_tt_MZP600A0400_TESUp);
    T->Branch("mvaVar_tt_MZP800A0400_TESUp", &mvaVar_tt_MZP800A0400_TESUp);
    T->Branch("mvaVar_tt_MZP1000A0400_TESUp", &mvaVar_tt_MZP1000A0400_TESUp);
    T->Branch("mvaVar_tt_MZP1200A0400_TESUp", &mvaVar_tt_MZP1200A0400_TESUp);

    T->Branch("mvaVar_mt_MZP600A0400_TESDown", &mvaVar_mt_MZP600A0400_TESDown);
    T->Branch("mvaVar_mt_MZP800A0400_TESDown", &mvaVar_mt_MZP800A0400_TESDown);
    T->Branch("mvaVar_mt_MZP1000A0400_TESDown", &mvaVar_mt_MZP1000A0400_TESDown);
    T->Branch("mvaVar_mt_MZP1200A0400_TESDown", &mvaVar_mt_MZP1200A0400_TESDown);

    T->Branch("mvaVar_et_MZP600A0400_TESDown", &mvaVar_et_MZP600A0400_TESDown);
    T->Branch("mvaVar_et_MZP800A0400_TESDown", &mvaVar_et_MZP800A0400_TESDown); 
    T->Branch("mvaVar_et_MZP1000A0400_TESDown", &mvaVar_et_MZP1000A0400_TESDown);
    T->Branch("mvaVar_et_MZP1200A0400_TESDown", &mvaVar_et_MZP1200A0400_TESDown);
        
    T->Branch("mvaVar_tt_MZP600A0400_TESDown", &mvaVar_tt_MZP600A0400_TESDown);
    T->Branch("mvaVar_tt_MZP800A0400_TESDown", &mvaVar_tt_MZP800A0400_TESDown);
    T->Branch("mvaVar_tt_MZP1000A0400_TESDown", &mvaVar_tt_MZP1000A0400_TESDown);
    T->Branch("mvaVar_tt_MZP1200A0400_TESDown", &mvaVar_tt_MZP1200A0400_TESDown);
    
    T->Branch("pt_1_dm0TESUp", &pt_1_dm0TESUp);
    T->Branch("mt_1_dm0TESUp", &mt_1_dm0TESUp);
    T->Branch("pfmt_1_dm0TESUp", &pfmt_1_dm0TESUp);

    T->Branch("pt_2_dm0TESUp", &pt_2_dm0TESUp);
    T->Branch("mt_2_dm0TESUp", &mt_2_dm0TESUp);
    T->Branch("pfmt_2_dm0TESUp", &pfmt_2_dm0TESUp);

    T->Branch("mt_tot_dm0TESUp", &mt_tot_dm0TESUp);
    T->Branch("pt_tt_dm0TESUp", &pt_tt_dm0TESUp);
    T->Branch("m_vis_dm0TESUp", &m_vis_dm0TESUp);
    T->Branch("pfmet_type1_dm0TESUp_Pt", &pfmet_type1_dm0TESUp_Pt);
    T->Branch("P_chi_pf_dm0TESUp", &P_chi_pf_dm0TESUp);
    T->Branch("LPT_dm0TESUp", &LPT_dm0TESUp);
    T->Branch("DeltaR_leg1_leg2", &DeltaR_leg1_leg2);
    T->Branch("cos_DeltaPhi_leg1_leg2", &cos_DeltaPhi_leg1_leg2);
    T->Branch("cos_DeltaPhi_PFMET_Higgs_dm0TESUp", &cos_DeltaPhi_PFMET_Higgs_dm0TESUp);

    T->Branch("pt_1_dm0TESDown", &pt_1_dm0TESDown);
    T->Branch("mt_1_dm0TESDown", &mt_1_dm0TESDown);
    T->Branch("pfmt_1_dm0TESDown", &pfmt_1_dm0TESDown);

    T->Branch("pt_2_dm0TESDown", &pt_2_dm0TESDown);
    T->Branch("mt_2_dm0TESDown", &mt_2_dm0TESDown);
    T->Branch("pfmt_2_dm0TESDown", &pfmt_2_dm0TESDown);

    T->Branch("mt_tot_dm0TESDown", &mt_tot_dm0TESDown);
    T->Branch("pt_tt_dm0TESDown", &pt_tt_dm0TESDown);
    T->Branch("m_vis_dm0TESDown", &m_vis_dm0TESDown);
    T->Branch("pfmet_type1_dm0TESDown_Pt", &pfmet_type1_dm0TESDown_Pt);
    T->Branch("P_chi_pf_dm0TESDown", &P_chi_pf_dm0TESDown);
    T->Branch("LPT_dm0TESDown", &LPT_dm0TESDown);
    T->Branch("DeltaR_leg1_leg2", &DeltaR_leg1_leg2);
    T->Branch("cos_DeltaPhi_leg1_leg2", &cos_DeltaPhi_leg1_leg2);
    T->Branch("cos_DeltaPhi_PFMET_Higgs_dm0TESDown", &cos_DeltaPhi_PFMET_Higgs_dm0TESDown);

    T->Branch("mvaVar_mt_MZP600A0400_dm0TESUp", &mvaVar_mt_MZP600A0400_dm0TESUp);
    T->Branch("mvaVar_mt_MZP800A0400_dm0TESUp", &mvaVar_mt_MZP800A0400_dm0TESUp);
    T->Branch("mvaVar_mt_MZP1000A0400_dm0TESUp", &mvaVar_mt_MZP1000A0400_dm0TESUp);
    T->Branch("mvaVar_mt_MZP1200A0400_dm0TESUp", &mvaVar_mt_MZP1200A0400_dm0TESUp);

    T->Branch("mvaVar_et_MZP600A0400_dm0TESUp", &mvaVar_et_MZP600A0400_dm0TESUp);
    T->Branch("mvaVar_et_MZP800A0400_dm0TESUp", &mvaVar_et_MZP800A0400_dm0TESUp); 
    T->Branch("mvaVar_et_MZP1000A0400_dm0TESUp", &mvaVar_et_MZP1000A0400_dm0TESUp);
    T->Branch("mvaVar_et_MZP1200A0400_dm0TESUp", &mvaVar_et_MZP1200A0400_dm0TESUp);
        
    T->Branch("mvaVar_tt_MZP600A0400_dm0TESUp", &mvaVar_tt_MZP600A0400_dm0TESUp);
    T->Branch("mvaVar_tt_MZP800A0400_dm0TESUp", &mvaVar_tt_MZP800A0400_dm0TESUp);
    T->Branch("mvaVar_tt_MZP1000A0400_dm0TESUp", &mvaVar_tt_MZP1000A0400_dm0TESUp);
    T->Branch("mvaVar_tt_MZP1200A0400_dm0TESUp", &mvaVar_tt_MZP1200A0400_dm0TESUp);

    T->Branch("mvaVar_mt_MZP600A0400_dm0TESDown", &mvaVar_mt_MZP600A0400_dm0TESDown);
    T->Branch("mvaVar_mt_MZP800A0400_dm0TESDown", &mvaVar_mt_MZP800A0400_dm0TESDown);
    T->Branch("mvaVar_mt_MZP1000A0400_dm0TESDown", &mvaVar_mt_MZP1000A0400_dm0TESDown);
    T->Branch("mvaVar_mt_MZP1200A0400_dm0TESDown", &mvaVar_mt_MZP1200A0400_dm0TESDown);

    T->Branch("mvaVar_et_MZP600A0400_dm0TESDown", &mvaVar_et_MZP600A0400_dm0TESDown);
    T->Branch("mvaVar_et_MZP800A0400_dm0TESDown", &mvaVar_et_MZP800A0400_dm0TESDown); 
    T->Branch("mvaVar_et_MZP1000A0400_dm0TESDown", &mvaVar_et_MZP1000A0400_dm0TESDown);
    T->Branch("mvaVar_et_MZP1200A0400_dm0TESDown", &mvaVar_et_MZP1200A0400_dm0TESDown);
        
    T->Branch("mvaVar_tt_MZP600A0400_dm0TESDown", &mvaVar_tt_MZP600A0400_dm0TESDown);
    T->Branch("mvaVar_tt_MZP800A0400_dm0TESDown", &mvaVar_tt_MZP800A0400_dm0TESDown);
    T->Branch("mvaVar_tt_MZP1000A0400_dm0TESDown", &mvaVar_tt_MZP1000A0400_dm0TESDown);
    T->Branch("mvaVar_tt_MZP1200A0400_dm0TESDown", &mvaVar_tt_MZP1200A0400_dm0TESDown);
    
    T->Branch("pt_1_dm1TESUp", &pt_1_dm1TESUp);
    T->Branch("mt_1_dm1TESUp", &mt_1_dm1TESUp);
    T->Branch("pfmt_1_dm1TESUp", &pfmt_1_dm1TESUp);

    T->Branch("pt_2_dm1TESUp", &pt_2_dm1TESUp);
    T->Branch("mt_2_dm1TESUp", &mt_2_dm1TESUp);
    T->Branch("pfmt_2_dm1TESUp", &pfmt_2_dm1TESUp);

    T->Branch("mt_tot_dm1TESUp", &mt_tot_dm1TESUp);
    T->Branch("pt_tt_dm1TESUp", &pt_tt_dm1TESUp);
    T->Branch("m_vis_dm1TESUp", &m_vis_dm1TESUp);
    T->Branch("pfmet_type1_dm1TESUp_Pt", &pfmet_type1_dm1TESUp_Pt);
    T->Branch("P_chi_pf_dm1TESUp", &P_chi_pf_dm1TESUp);
    T->Branch("LPT_dm1TESUp", &LPT_dm1TESUp);
    T->Branch("DeltaR_leg1_leg2", &DeltaR_leg1_leg2);
    T->Branch("cos_DeltaPhi_leg1_leg2", &cos_DeltaPhi_leg1_leg2);
    T->Branch("cos_DeltaPhi_PFMET_Higgs_dm1TESUp", &cos_DeltaPhi_PFMET_Higgs_dm1TESUp);

    T->Branch("pt_1_dm1TESDown", &pt_1_dm1TESDown);
    T->Branch("mt_1_dm1TESDown", &mt_1_dm1TESDown);
    T->Branch("pfmt_1_dm1TESDown", &pfmt_1_dm1TESDown);

    T->Branch("pt_2_dm1TESDown", &pt_2_dm1TESDown);
    T->Branch("mt_2_dm1TESDown", &mt_2_dm1TESDown);
    T->Branch("pfmt_2_dm1TESDown", &pfmt_2_dm1TESDown);

    T->Branch("mt_tot_dm1TESDown", &mt_tot_dm1TESDown);
    T->Branch("pt_tt_dm1TESDown", &pt_tt_dm1TESDown);
    T->Branch("m_vis_dm1TESDown", &m_vis_dm1TESDown);
    T->Branch("pfmet_type1_dm1TESDown_Pt", &pfmet_type1_dm1TESDown_Pt);
    T->Branch("P_chi_pf_dm1TESDown", &P_chi_pf_dm1TESDown);
    T->Branch("LPT_dm1TESDown", &LPT_dm1TESDown);
    T->Branch("DeltaR_leg1_leg2", &DeltaR_leg1_leg2);
    T->Branch("cos_DeltaPhi_leg1_leg2", &cos_DeltaPhi_leg1_leg2);
    T->Branch("cos_DeltaPhi_PFMET_Higgs_dm1TESDown", &cos_DeltaPhi_PFMET_Higgs_dm1TESDown);

    T->Branch("mvaVar_mt_MZP600A0400_dm1TESUp", &mvaVar_mt_MZP600A0400_dm1TESUp);
    T->Branch("mvaVar_mt_MZP800A0400_dm1TESUp", &mvaVar_mt_MZP800A0400_dm1TESUp);
    T->Branch("mvaVar_mt_MZP1000A0400_dm1TESUp", &mvaVar_mt_MZP1000A0400_dm1TESUp);
    T->Branch("mvaVar_mt_MZP1200A0400_dm1TESUp", &mvaVar_mt_MZP1200A0400_dm1TESUp);

    T->Branch("mvaVar_et_MZP600A0400_dm1TESUp", &mvaVar_et_MZP600A0400_dm1TESUp);
    T->Branch("mvaVar_et_MZP800A0400_dm1TESUp", &mvaVar_et_MZP800A0400_dm1TESUp); 
    T->Branch("mvaVar_et_MZP1000A0400_dm1TESUp", &mvaVar_et_MZP1000A0400_dm1TESUp);
    T->Branch("mvaVar_et_MZP1200A0400_dm1TESUp", &mvaVar_et_MZP1200A0400_dm1TESUp);
        
    T->Branch("mvaVar_tt_MZP600A0400_dm1TESUp", &mvaVar_tt_MZP600A0400_dm1TESUp);
    T->Branch("mvaVar_tt_MZP800A0400_dm1TESUp", &mvaVar_tt_MZP800A0400_dm1TESUp);
    T->Branch("mvaVar_tt_MZP1000A0400_dm1TESUp", &mvaVar_tt_MZP1000A0400_dm1TESUp);
    T->Branch("mvaVar_tt_MZP1200A0400_dm1TESUp", &mvaVar_tt_MZP1200A0400_dm1TESUp);

    T->Branch("mvaVar_mt_MZP600A0400_dm1TESDown", &mvaVar_mt_MZP600A0400_dm1TESDown);
    T->Branch("mvaVar_mt_MZP800A0400_dm1TESDown", &mvaVar_mt_MZP800A0400_dm1TESDown);
    T->Branch("mvaVar_mt_MZP1000A0400_dm1TESDown", &mvaVar_mt_MZP1000A0400_dm1TESDown);
    T->Branch("mvaVar_mt_MZP1200A0400_dm1TESDown", &mvaVar_mt_MZP1200A0400_dm1TESDown);

    T->Branch("mvaVar_et_MZP600A0400_dm1TESDown", &mvaVar_et_MZP600A0400_dm1TESDown);
    T->Branch("mvaVar_et_MZP800A0400_dm1TESDown", &mvaVar_et_MZP800A0400_dm1TESDown); 
    T->Branch("mvaVar_et_MZP1000A0400_dm1TESDown", &mvaVar_et_MZP1000A0400_dm1TESDown);
    T->Branch("mvaVar_et_MZP1200A0400_dm1TESDown", &mvaVar_et_MZP1200A0400_dm1TESDown);
        
    T->Branch("mvaVar_tt_MZP600A0400_dm1TESDown", &mvaVar_tt_MZP600A0400_dm1TESDown);
    T->Branch("mvaVar_tt_MZP800A0400_dm1TESDown", &mvaVar_tt_MZP800A0400_dm1TESDown);
    T->Branch("mvaVar_tt_MZP1000A0400_dm1TESDown", &mvaVar_tt_MZP1000A0400_dm1TESDown);
    T->Branch("mvaVar_tt_MZP1200A0400_dm1TESDown", &mvaVar_tt_MZP1200A0400_dm1TESDown);
    
    T->Branch("pt_1_dm10TESUp", &pt_1_dm10TESUp);
    T->Branch("mt_1_dm10TESUp", &mt_1_dm10TESUp);
    T->Branch("pfmt_1_dm10TESUp", &pfmt_1_dm10TESUp);

    T->Branch("pt_2_dm10TESUp", &pt_2_dm10TESUp);
    T->Branch("mt_2_dm10TESUp", &mt_2_dm10TESUp);
    T->Branch("pfmt_2_dm10TESUp", &pfmt_2_dm10TESUp);
    
    T->Branch("pfmt_1_JEnUp", &pfmt_1_JEnUp);
    T->Branch("pfmt_1_JEnDown", &pfmt_1_JEnDown);

    T->Branch("mt_tot_dm10TESUp", &mt_tot_dm10TESUp);
    T->Branch("pt_tt_dm10TESUp", &pt_tt_dm10TESUp);
    T->Branch("m_vis_dm10TESUp", &m_vis_dm10TESUp);
    T->Branch("pfmet_type1_dm10TESUp_Pt", &pfmet_type1_dm10TESUp_Pt);
    T->Branch("P_chi_pf_dm10TESUp", &P_chi_pf_dm10TESUp);
    T->Branch("LPT_dm10TESUp", &LPT_dm10TESUp);
    T->Branch("DeltaR_leg1_leg2", &DeltaR_leg1_leg2);
    T->Branch("cos_DeltaPhi_leg1_leg2", &cos_DeltaPhi_leg1_leg2);
    T->Branch("cos_DeltaPhi_PFMET_Higgs_dm10TESUp", &cos_DeltaPhi_PFMET_Higgs_dm10TESUp);

    T->Branch("pt_1_dm10TESDown", &pt_1_dm10TESDown);
    T->Branch("mt_1_dm10TESDown", &mt_1_dm10TESDown);
    T->Branch("pfmt_1_dm10TESDown", &pfmt_1_dm10TESDown);

    T->Branch("pt_2_dm10TESDown", &pt_2_dm10TESDown);
    T->Branch("mt_2_dm10TESDown", &mt_2_dm10TESDown);
    T->Branch("pfmt_2_dm10TESDown", &pfmt_2_dm10TESDown);

    T->Branch("mt_tot_dm10TESDown", &mt_tot_dm10TESDown);
    T->Branch("pt_tt_dm10TESDown", &pt_tt_dm10TESDown);
    T->Branch("m_vis_dm10TESDown", &m_vis_dm10TESDown);
    T->Branch("pfmet_type1_dm10TESDown_Pt", &pfmet_type1_dm10TESDown_Pt);
    T->Branch("P_chi_pf_dm10TESDown", &P_chi_pf_dm10TESDown);
    T->Branch("LPT_dm10TESDown", &LPT_dm10TESDown);
    T->Branch("DeltaR_leg1_leg2", &DeltaR_leg1_leg2);
    T->Branch("cos_DeltaPhi_leg1_leg2", &cos_DeltaPhi_leg1_leg2);
    T->Branch("cos_DeltaPhi_PFMET_Higgs_dm10TESDown", &cos_DeltaPhi_PFMET_Higgs_dm10TESDown);

    T->Branch("mvaVar_mt_MZP600A0400_dm10TESUp", &mvaVar_mt_MZP600A0400_dm10TESUp);
    T->Branch("mvaVar_mt_MZP800A0400_dm10TESUp", &mvaVar_mt_MZP800A0400_dm10TESUp);
    T->Branch("mvaVar_mt_MZP1000A0400_dm10TESUp", &mvaVar_mt_MZP1000A0400_dm10TESUp);
    T->Branch("mvaVar_mt_MZP1200A0400_dm10TESUp", &mvaVar_mt_MZP1200A0400_dm10TESUp);

    T->Branch("mvaVar_et_MZP600A0400_dm10TESUp", &mvaVar_et_MZP600A0400_dm10TESUp);
    T->Branch("mvaVar_et_MZP800A0400_dm10TESUp", &mvaVar_et_MZP800A0400_dm10TESUp); 
    T->Branch("mvaVar_et_MZP1000A0400_dm10TESUp", &mvaVar_et_MZP1000A0400_dm10TESUp);
    T->Branch("mvaVar_et_MZP1200A0400_dm10TESUp", &mvaVar_et_MZP1200A0400_dm10TESUp);
        
    T->Branch("mvaVar_tt_MZP600A0400_dm10TESUp", &mvaVar_tt_MZP600A0400_dm10TESUp);
    T->Branch("mvaVar_tt_MZP800A0400_dm10TESUp", &mvaVar_tt_MZP800A0400_dm10TESUp);
    T->Branch("mvaVar_tt_MZP1000A0400_dm10TESUp", &mvaVar_tt_MZP1000A0400_dm10TESUp);
    T->Branch("mvaVar_tt_MZP1200A0400_dm10TESUp", &mvaVar_tt_MZP1200A0400_dm10TESUp);

    T->Branch("mvaVar_mt_MZP600A0400_dm10TESDown", &mvaVar_mt_MZP600A0400_dm10TESDown);
    T->Branch("mvaVar_mt_MZP800A0400_dm10TESDown", &mvaVar_mt_MZP800A0400_dm10TESDown);
    T->Branch("mvaVar_mt_MZP1000A0400_dm10TESDown", &mvaVar_mt_MZP1000A0400_dm10TESDown);
    T->Branch("mvaVar_mt_MZP1200A0400_dm10TESDown", &mvaVar_mt_MZP1200A0400_dm10TESDown);

    T->Branch("mvaVar_et_MZP600A0400_dm10TESDown", &mvaVar_et_MZP600A0400_dm10TESDown);
    T->Branch("mvaVar_et_MZP800A0400_dm10TESDown", &mvaVar_et_MZP800A0400_dm10TESDown); 
    T->Branch("mvaVar_et_MZP1000A0400_dm10TESDown", &mvaVar_et_MZP1000A0400_dm10TESDown);
    T->Branch("mvaVar_et_MZP1200A0400_dm10TESDown", &mvaVar_et_MZP1200A0400_dm10TESDown);
        
    T->Branch("mvaVar_tt_MZP600A0400_dm10TESDown", &mvaVar_tt_MZP600A0400_dm10TESDown);
    T->Branch("mvaVar_tt_MZP800A0400_dm10TESDown", &mvaVar_tt_MZP800A0400_dm10TESDown);
    T->Branch("mvaVar_tt_MZP1000A0400_dm10TESDown", &mvaVar_tt_MZP1000A0400_dm10TESDown);
    T->Branch("mvaVar_tt_MZP1200A0400_dm10TESDown", &mvaVar_tt_MZP1200A0400_dm10TESDown);

}

double generateH2TauSyncTree::mtTotCalc(TLorentzVector l1, TLorentzVector l2, TLorentzVector met )
{
        double value = 0.;

        double part1 = 2 * l1.Pt() * met.Pt() * (1-cos(l1.DeltaPhi(met)));
        double part2 = 2 * l2.Pt() * met.Pt() * (1-cos(l2.DeltaPhi(met)));
        double part3 = 2 * l1.Pt() * l2.Pt()  * (1-cos(l1.DeltaPhi(l2)));

        if( part1 + part2 + part3 > 0) value = sqrt( part1 + part2 + part3 );

        return value;
}

double generateH2TauSyncTree::pzetaVisCalc(TLorentzVector e, TLorentzVector mu)
{
	double value = 0;

	double den = sqrt(pow((cos(e.Phi())+cos(mu.Phi())),2) + pow((sin(e.Phi())+sin(mu.Phi())),2));
	double num = (e.Px()+mu.Px())*(cos(e.Phi())+cos(mu.Phi()));
	num += (e.Py()+mu.Py())*(sin(e.Phi())+cos(mu.Phi()));

	if(den!=0) value = num/den;


	return value;
}

double generateH2TauSyncTree::pzetaMissCalc(TLorentzVector e, TLorentzVector mu, TLorentzVector met)
{

	double value = 0;

	double den = sqrt(pow((cos(e.Phi())+cos(mu.Phi())),2) + pow((sin(e.Phi())+sin(mu.Phi())),2));
	double num = 	(met.Px())*(cos(e.Phi())+cos(mu.Phi()));
	num += (met.Py())*(sin(e.Phi())+sin(mu.Phi()));

	if(den!=0) value = num/den;
	return value;

}

void generateH2TauSyncTree::reset()
{

	/* in general set :

		-- weights to 1.0
		-- floats and doubles to -999.0  ** be careful when cutting to exclude -999.0 for example iso < 0.1 will pass!!! **
		-- ints to -999  ** be careful when cutting to exclude -999.0 for example pT < 100 will pass!!! **
		-- strings to "NULL"
		-- bool to 0
		-- and use .clear() on vectors

	*/


	final_weight = 1.0;
    
    sf_IDISO = 1.0;
    sf_TRIG = 1.0;
    sf_TRACK = 1.0;
    sf_ALD = 1.0;
    
	nominalCrossSection_Weight = 1.0;
	puWeight_Weight = 1.0;
	TopQuarkPtWeight_Weight = 1.0;
	ZReWeight_Weight = 1.0;
    ZReWeight_WeightUp = 1.0;
    ZReWeight_WeightDown = 1.0;
    KReWeight_Weight = 1.0;
    KReWeight_WeightUp = 1.0;
    KReWeight_WeightDown = 1.0;
    ZZReWeight_Weight = 1.0;           /* indluded in final_weight for ZZ*/
    ZZReWeight_WeightUp = 1.0;
    ZZReWeight_WeightDown = 1.0;
    WWReWeight_WeightUp = 1.0;
    JTF_WeightUp = 1.0;
    JTF_WeightDown = 1.0;
	NLOReWeight_Weight = 1.0;
	ScaleFactorsForPair_Weight = 1.0;
	ScaleFactorsForPair_WeightUp = 1.0;
	ScaleFactorsForPair_WeightDown = 1.0;
	QCDWeightForEleMuChannel_Weight = 1.0;
	QCDWeightForEleMuChannel_WeightUp = 1.0;
	QCDWeightForEleMuChannel_WeightDown = 1.0;
	QCDWeightForEleMuChannelNoPZetaCut_Weight = 1.0;
	QCDWeightForEleMuChannelNoPZetaCut_WeightUp = 1.0;
	QCDWeightForEleMuChannelNoPZetaCut_WeightDown = 1.0;
	highPtTauEff_WeightUp = 1.0;
	highPtTauEff_WeightDown = 1.0;

    flag_MVAEventType = -999;
    randNum = -999;
    
    mvaVar_mt_MZP600A0400= -999.;
    mvaVar_mt_MZP800A0400= -999.;
    mvaVar_mt_MZP1000A0400= -999.;
    mvaVar_mt_MZP1200A0400= -999.;
    
    mvaVar_et_MZP600A0400= -999.;
    mvaVar_et_MZP800A0400= -999.;
    mvaVar_et_MZP1000A0400= -999.;
    mvaVar_et_MZP1200A0400= -999.;
    
    mvaVar_tt_MZP600A0400= -999.;
    mvaVar_tt_MZP800A0400= -999.;
    mvaVar_tt_MZP1000A0400= -999.;
    mvaVar_tt_MZP1200A0400= -999.;
    
    mvaVar_mt_MZP600A0400_UESUp= -999.;
    mvaVar_mt_MZP800A0400_UESUp= -999.;
    mvaVar_mt_MZP1000A0400_UESUp= -999.;
    mvaVar_mt_MZP1200A0400_UESUp= -999.;
    
    mvaVar_et_MZP600A0400_UESUp= -999.;
    mvaVar_et_MZP800A0400_UESUp= -999.;
    mvaVar_et_MZP1000A0400_UESUp= -999.;
    mvaVar_et_MZP1200A0400_UESUp= -999.;
    
    mvaVar_tt_MZP600A0400_UESUp= -999.;
    mvaVar_tt_MZP800A0400_UESUp= -999.;
    mvaVar_tt_MZP1000A0400_UESUp= -999.;
    mvaVar_tt_MZP1200A0400_UESUp= -999.;
    
    mvaVar_mt_MZP600A0400_UESDown= -999.;
    mvaVar_mt_MZP800A0400_UESDown= -999.;
    mvaVar_mt_MZP1000A0400_UESDown= -999.;
    mvaVar_mt_MZP1200A0400_UESDown= -999.;
    
    mvaVar_et_MZP600A0400_UESDown= -999.;
    mvaVar_et_MZP800A0400_UESDown= -999.;
    mvaVar_et_MZP1000A0400_UESDown= -999.;
    mvaVar_et_MZP1200A0400_UESDown= -999.;
    
    mvaVar_tt_MZP600A0400_UESDown= -999.;
    mvaVar_tt_MZP800A0400_UESDown= -999.;
    mvaVar_tt_MZP1000A0400_UESDown= -999.;
    mvaVar_tt_MZP1200A0400_UESDown= -999.;
    
    rpt_1 = -999.;
    rpt_2 = -999.;
    rpfmt_1 = -999.;
    rmt_tot = -999.;
    rpt_tt = -999.;
    rm_vis = -999.;
    rmet = -999.;
    rP_chi_pf = -999.;
    rLPT = -999.;
    rDeltaR_leg1_leg2 = -999.;
    rcos_DeltaPhi_leg1_leg2 = -999.;
    rcos_DeltaPhi_PFMET_Higgs = -999.;

    rratio_weight = -999.;
    rfinal_weight = -999.;
    rnpu = -999.;
    revent = -999.;
    rrandNum = -999.;
    rDataCardInt = -999.;
    rIsZTT = -999.;
    rIsZJ = -999.;
    rIsZL = -999.;
    rIsTTT = -999.;

    originalXWGTUP = 1.0;			/* always init a weight to 1.0 */
    theory_scale_factors.clear();	/* std::vectors are reset using clear() */

	pairRank = 999;
	isOsPair = -999;
    isBoostedChannelPair = 0;
    DataCardInt = 0;
    
    genBosonTotal_pt = -999.;
    genBosonTotal_eta = -999.;
    genBosonTotal_phi = -999.;
    genBosonTotal_M = -999.;

	genBosonVisible_pt = -999.;
	genBosonVisible_eta = -999.;
	genBosonVisible_phi = -999.;
	genBosonVisible_M = -999.;
    
    genBosonTotal_Wpt = -999.;
    
    leg1_GENMOTHERpdgId = -999;
    leg2_GENMOTHERpdgId = -999;

	genTopPt1 = -999.;
	genTopPt2 = -999.;

	run = 0;
	event = 0;
	evt = 0;
	lumi = 0;

	npv =  0;
	npu =  -999.0;
	rho =  -999.0;
	puweight = 1.0;

	pt_1 = -999.0;
	phi_1 = -999.0;
	eta_1 = -999.0;
	m_1 = -999.0;
    pt_1_flat = -999.0;
	phi_1_flat = -999.0;
	eta_1_flat = -999.0;
	m_1_flat = -999.0;
	iso_1 = -999.0;
	dZ_1 = -999.0;
	dzTauVertex_1 = -999.0;
	d0_1 = -999.0;
	q_1 = -999;
	id_e_mva_nt_loose_1 = -999.0;
	tau_decay_mode_1 = -999;
	ZimpactTau_1 = -999.0;
	mt_1 = -999.0;
	pfmt_1 = -999.0;
    pfmt_1_UESUp = -999.0;
    pfmt_1_UESDown = -999.0;
    pfmt_1_JEnUp = -999.0;
    pfmt_1_JEnDown = -999.0;
	puppimt_1 = -999.0;
	mt_uncorr_1 = -999.0;
	responseUP_MTmvaMET_1 = -999.0;
	responseDOWN_MTmvaMET_1 = -999.0;
	resolutionUP_MTmvaMET_1 = -999.0;
	resolutionDOWN_MTmvaMET_1 = -999.0;
	gen_match_1 = -999;
	genMCmatch_pt_1 = -999.0;
	genMCmatch_eta_1 = -999.0;
	genMCmatch_phi_1 = -999.0;
	genMCmatch_M_1 = -999.0;
	MCMatchPdgId_1 = -999;
	byIsolationMVArun2v1DBdR03oldDMwLTraw_1 = -999.0;
	byTightIsolationMVArun2v1DBdR03oldDMwLT_1 = -999.0;
	byVTightIsolationMVArun2v1DBdR03oldDMwLT_1 = -999.0;
	byLooseIsolationMVArun2v1DBdR03oldDMwLT_1 = -999.0;
	byMediumIsolationMVArun2v1DBdR03oldDMwLT_1 = -999.0;
	byVLooseIsolationMVArun2v1DBdR03oldDMwLT_1 = -999.0;
	byVVTightIsolationMVArun2v1DBdR03oldDMwLT_1 = -999.0;
    byLooseIsolationMVArun2v1DBoldDMwLT_1 = -999.0;
    byTightIsolationMVArun2v1DBoldDMwLT_1 = -999.0;
	againstElectronVLooseMVA6_1 = -999.0;
	againstMuonTight3_1 = -999.0;
	againstElectronTightMVA6_1 = -999.0;
	againstMuonLoose3_1 = -999.0;
	decayModeFinding_1 = -999.0;
	//byIsolationMVA3oldDMwLTraw_1 = -999.0;
	byCombinedIsolationDeltaBetaCorrRaw3Hits_1 = -999.0;
	byIsolationMVArun2v1DBnewDMwLTraw_1 = -999.0;
	decayModeFindingNewDMs_1 = -999.0;

	pt_2 = -999.0;
	phi_2 = -999.0;
	eta_2 = -999.0;
	m_2 = -999.0;
    pt_2_flat = -999.0;
	phi_2_flat = -999.0;
	eta_2_flat = -999.0;
	m_2_flat = -999.0;
	iso_2 = -999.0;
	dZ_2 = -999.0;
	dzTauVertex_2 = -999.0;
	d0_2 = -999.0;
	q_2 = -999;
	id_e_mva_nt_loose_2 = -999.0;
	tau_decay_mode_2 = -999;
	ZimpactTau_2 = -999.0;
	mt_2 = -999.0;
	pfmt_2 = -999.0;
	puppimt_2 = -999.0;
	mt_uncorr_2 = -999.0;
	responseUP_MTmvaMET_2 = -999.0;
	responseDOWN_MTmvaMET_2 = -999.0;
	resolutionUP_MTmvaMET_2 = -999.0;
	resolutionDOWN_MTmvaMET_2 = -999.0;
	gen_match_2 = -999;
	genMCmatch_pt_2 = -999.0;
	genMCmatch_eta_2 = -999.0;
	genMCmatch_phi_2 = -999.0;
	genMCmatch_M_2 = -999.0;
	MCMatchPdgId_2 = -999;
	byIsolationMVArun2v1DBdR03oldDMwLTraw_2 = -999.0;
	byTightIsolationMVArun2v1DBdR03oldDMwLT_2 = -999.0;
	byVTightIsolationMVArun2v1DBdR03oldDMwLT_2 = -999.0;
	byLooseIsolationMVArun2v1DBdR03oldDMwLT_2 = -999.0;
	byMediumIsolationMVArun2v1DBdR03oldDMwLT_2 = -999.0;
	byVLooseIsolationMVArun2v1DBdR03oldDMwLT_2 = -999.0;
	byVVTightIsolationMVArun2v1DBdR03oldDMwLT_2 = -999.0;
    byLooseIsolationMVArun2v1DBoldDMwLT_2 = -999.0;
    byTightIsolationMVArun2v1DBoldDMwLT_2 = -999.0;
	againstElectronVLooseMVA6_2 = -999.0;
	againstMuonTight3_2 = -999.0;
	againstElectronTightMVA6_2 = -999.0;
	againstMuonLoose3_2 = -999.0;
	decayModeFinding_2 = -999.0;
	//byIsolationMVA3oldDMwLTraw_2 = -999.0;
	byCombinedIsolationDeltaBetaCorrRaw3Hits_2 = -999.0;
	byIsolationMVArun2v1DBnewDMwLTraw_2 = -999.0;
	decayModeFindingNewDMs_2 = -999.0;

    pt_1_TESUp = -999.0;
    mt_1_TESUp = -999.0;
    pfmt_1_TESUp = -999.0;

    pt_2_TESUp = -999.0;
    mt_2_TESUp = -999.0;
    pfmt_2_TESUp = -999.0;

    mt_tot_TESUp = -999.0;
    pt_tt_TESUp = -999.0;
    m_vis_TESUp = -999.0;
    pfmet_type1_TESUp_Pt = -999.0;
    P_chi_pf_TESUp = -999.0;
    LPT_TESUp = -999.0;
    DeltaR_leg1_leg2 = -999.0;
    cos_DeltaPhi_leg1_leg2 = -999.0;
    cos_DeltaPhi_PFMET_Higgs_TESUp = -999.0;

    pt_1_TESDown = -999.0;
    mt_1_TESDown = -999.0;
    pfmt_1_TESDown = -999.0;

    pt_2_TESDown = -999.0;
    mt_2_TESDown = -999.0;
    pfmt_2_TESDown = -999.0;

    mt_tot_TESDown = -999.0;
    pt_tt_TESDown = -999.0;
    m_vis_TESDown = -999.0;
    pfmet_type1_TESDown_Pt = -999.0;
    P_chi_pf_TESDown = -999.0;
    LPT_TESDown = -999.0;
    DeltaR_leg1_leg2 = -999.0;
    cos_DeltaPhi_leg1_leg2 = -999.0;
    cos_DeltaPhi_PFMET_Higgs_TESDown = -999.0;
    
    pt_1_L1TESUp = -999.0;
    mt_1_L1TESUp = -999.0;
    pfmt_1_L1TESUp = -999.0;

    pt_2_L1TESUp = -999.0;
    mt_2_L1TESUp = -999.0;
    pfmt_2_L1TESUp = -999.0;

    mt_tot_L1TESUp = -999.0;
    pt_tt_L1TESUp = -999.0;
    m_vis_L1TESUp = -999.0;
    pfmet_type1_L1TESUp_Pt = -999.0;
    P_chi_pf_L1TESUp = -999.0;
    LPT_L1TESUp = -999.0;
    DeltaR_leg1_leg2 = -999.0;
    cos_DeltaPhi_leg1_leg2 = -999.0;
    cos_DeltaPhi_PFMET_Higgs_L1TESUp = -999.0;

    pt_1_L1TESDown = -999.0;
    mt_1_L1TESDown = -999.0;
    pfmt_1_L1TESDown = -999.0;

    pt_2_L1TESDown = -999.0;
    mt_2_L1TESDown = -999.0;
    pfmt_2_L1TESDown = -999.0;

    mt_tot_L1TESDown = -999.0;
    pt_tt_L1TESDown = -999.0;
    m_vis_L1TESDown = -999.0;
    pfmet_type1_L1TESDown_Pt = -999.0;
    P_chi_pf_L1TESDown = -999.0;
    LPT_L1TESDown = -999.0;
    DeltaR_leg1_leg2 = -999.0;
    cos_DeltaPhi_leg1_leg2 = -999.0;
    cos_DeltaPhi_PFMET_Higgs_L1TESDown = -999.0;
    
    pt_1_L2TESUp = -999.0;
    mt_1_L2TESUp = -999.0;
    pfmt_1_L2TESUp = -999.0;

    pt_2_L2TESUp = -999.0;
    mt_2_L2TESUp = -999.0;
    pfmt_2_L2TESUp = -999.0;

    mt_tot_L2TESUp = -999.0;
    pt_tt_L2TESUp = -999.0;
    m_vis_L2TESUp = -999.0;
    pfmet_type1_L2TESUp_Pt = -999.0;
    P_chi_pf_L2TESUp = -999.0;
    LPT_L2TESUp = -999.0;
    DeltaR_leg1_leg2 = -999.0;
    cos_DeltaPhi_leg1_leg2 = -999.0;
    cos_DeltaPhi_PFMET_Higgs_L2TESUp = -999.0;

    pt_1_L2TESDown = -999.0;
    mt_1_L2TESDown = -999.0;
    pfmt_1_L2TESDown = -999.0;

    pt_2_L2TESDown = -999.0;
    mt_2_L2TESDown = -999.0;
    pfmt_2_L2TESDown = -999.0;

    mt_tot_L2TESDown = -999.0;
    pt_tt_L2TESDown = -999.0;
    m_vis_L2TESDown = -999.0;
    pfmet_type1_L2TESDown_Pt = -999.0;
    P_chi_pf_L2TESDown = -999.0;
    LPT_L2TESDown = -999.0;
    DeltaR_leg1_leg2 = -999.0;
    cos_DeltaPhi_leg1_leg2 = -999.0;
    cos_DeltaPhi_PFMET_Higgs_L2TESDown = -999.0;

    mvaVar_mt_MZP600A0400_TESUp = -999.0;
    mvaVar_mt_MZP800A0400_TESUp = -999.0;
    mvaVar_mt_MZP1000A0400_TESUp = -999.0;
    mvaVar_mt_MZP1200A0400_TESUp = -999.0;

    mvaVar_et_MZP600A0400_TESUp = -999.0;
    mvaVar_et_MZP800A0400_TESUp = -999.0; 
    mvaVar_et_MZP1000A0400_TESUp = -999.0;
    mvaVar_et_MZP1200A0400_TESUp = -999.0;
        
    mvaVar_tt_MZP600A0400_TESUp = -999.0;
    mvaVar_tt_MZP800A0400_TESUp = -999.0;
    mvaVar_tt_MZP1000A0400_TESUp = -999.0;
    mvaVar_tt_MZP1200A0400_TESUp = -999.0;

    mvaVar_mt_MZP600A0400_TESDown = -999.0;
    mvaVar_mt_MZP800A0400_TESDown = -999.0;
    mvaVar_mt_MZP1000A0400_TESDown = -999.0;
    mvaVar_mt_MZP1200A0400_TESDown = -999.0;

    mvaVar_et_MZP600A0400_TESDown = -999.0;
    mvaVar_et_MZP800A0400_TESDown = -999.0; 
    mvaVar_et_MZP1000A0400_TESDown = -999.0;
    mvaVar_et_MZP1200A0400_TESDown = -999.0;
        
    mvaVar_tt_MZP600A0400_TESDown = -999.0;
    mvaVar_tt_MZP800A0400_TESDown = -999.0;
    mvaVar_tt_MZP1000A0400_TESDown = -999.0;
    mvaVar_tt_MZP1200A0400_TESDown = -999.0;
    
    pt_1_dm0TESUp = -999.0;
    mt_1_dm0TESUp = -999.0;
    pfmt_1_dm0TESUp = -999.0;

    pt_2_dm0TESUp = -999.0;
    mt_2_dm0TESUp = -999.0;
    pfmt_2_dm0TESUp = -999.0;

    mt_tot_dm0TESUp = -999.0;
    pt_tt_dm0TESUp = -999.0;
    m_vis_dm0TESUp = -999.0;
    pfmet_type1_dm0TESUp_Pt = -999.0;
    P_chi_pf_dm0TESUp = -999.0;
    LPT_dm0TESUp = -999.0;
    DeltaR_leg1_leg2 = -999.0;
    cos_DeltaPhi_leg1_leg2 = -999.0;
    cos_DeltaPhi_PFMET_Higgs_dm0TESUp = -999.0;

    pt_1_dm0TESDown = -999.0;
    mt_1_dm0TESDown = -999.0;
    pfmt_1_dm0TESDown = -999.0;

    pt_2_dm0TESDown = -999.0;
    mt_2_dm0TESDown = -999.0;
    pfmt_2_dm0TESDown = -999.0;

    mt_tot_dm0TESDown = -999.0;
    pt_tt_dm0TESDown = -999.0;
    m_vis_dm0TESDown = -999.0;
    pfmet_type1_dm0TESDown_Pt = -999.0;
    P_chi_pf_dm0TESDown = -999.0;
    LPT_dm0TESDown = -999.0;
    DeltaR_leg1_leg2 = -999.0;
    cos_DeltaPhi_leg1_leg2 = -999.0;
    cos_DeltaPhi_PFMET_Higgs_dm0TESDown = -999.0;

    mvaVar_mt_MZP600A0400_dm0TESUp = -999.0;
    mvaVar_mt_MZP800A0400_dm0TESUp = -999.0;
    mvaVar_mt_MZP1000A0400_dm0TESUp = -999.0;
    mvaVar_mt_MZP1200A0400_dm0TESUp = -999.0;

    mvaVar_et_MZP600A0400_dm0TESUp = -999.0;
    mvaVar_et_MZP800A0400_dm0TESUp = -999.0; 
    mvaVar_et_MZP1000A0400_dm0TESUp = -999.0;
    mvaVar_et_MZP1200A0400_dm0TESUp = -999.0;
        
    mvaVar_tt_MZP600A0400_dm0TESUp = -999.0;
    mvaVar_tt_MZP800A0400_dm0TESUp = -999.0;
    mvaVar_tt_MZP1000A0400_dm0TESUp = -999.0;
    mvaVar_tt_MZP1200A0400_dm0TESUp = -999.0;

    mvaVar_mt_MZP600A0400_dm0TESDown = -999.0;
    mvaVar_mt_MZP800A0400_dm0TESDown = -999.0;
    mvaVar_mt_MZP1000A0400_dm0TESDown = -999.0;
    mvaVar_mt_MZP1200A0400_dm0TESDown = -999.0;

    mvaVar_et_MZP600A0400_dm0TESDown = -999.0;
    mvaVar_et_MZP800A0400_dm0TESDown = -999.0; 
    mvaVar_et_MZP1000A0400_dm0TESDown = -999.0;
    mvaVar_et_MZP1200A0400_dm0TESDown = -999.0;
        
    mvaVar_tt_MZP600A0400_dm0TESDown = -999.0;
    mvaVar_tt_MZP800A0400_dm0TESDown = -999.0;
    mvaVar_tt_MZP1000A0400_dm0TESDown = -999.0;
    mvaVar_tt_MZP1200A0400_dm0TESDown = -999.0;
    
    pt_1_dm1TESUp = -999.0;
    mt_1_dm1TESUp = -999.0;
    pfmt_1_dm1TESUp = -999.0;

    pt_2_dm1TESUp = -999.0;
    mt_2_dm1TESUp = -999.0;
    pfmt_2_dm1TESUp = -999.0;

    mt_tot_dm1TESUp = -999.0;
    pt_tt_dm1TESUp = -999.0;
    m_vis_dm1TESUp = -999.0;
    pfmet_type1_dm1TESUp_Pt = -999.0;
    P_chi_pf_dm1TESUp = -999.0;
    LPT_dm1TESUp = -999.0;
    DeltaR_leg1_leg2 = -999.0;
    cos_DeltaPhi_leg1_leg2 = -999.0;
    cos_DeltaPhi_PFMET_Higgs_dm1TESUp = -999.0;

    pt_1_dm1TESDown = -999.0;
    mt_1_dm1TESDown = -999.0;
    pfmt_1_dm1TESDown = -999.0;

    pt_2_dm1TESDown = -999.0;
    mt_2_dm1TESDown = -999.0;
    pfmt_2_dm1TESDown = -999.0;

    mt_tot_dm1TESDown = -999.0;
    pt_tt_dm1TESDown = -999.0;
    m_vis_dm1TESDown = -999.0;
    pfmet_type1_dm1TESDown_Pt = -999.0;
    P_chi_pf_dm1TESDown = -999.0;
    LPT_dm1TESDown = -999.0;
    DeltaR_leg1_leg2 = -999.0;
    cos_DeltaPhi_leg1_leg2 = -999.0;
    cos_DeltaPhi_PFMET_Higgs_dm1TESDown = -999.0;

    mvaVar_mt_MZP600A0400_dm1TESUp = -999.0;
    mvaVar_mt_MZP800A0400_dm1TESUp = -999.0;
    mvaVar_mt_MZP1000A0400_dm1TESUp = -999.0;
    mvaVar_mt_MZP1200A0400_dm1TESUp = -999.0;

    mvaVar_et_MZP600A0400_dm1TESUp = -999.0;
    mvaVar_et_MZP800A0400_dm1TESUp = -999.0; 
    mvaVar_et_MZP1000A0400_dm1TESUp = -999.0;
    mvaVar_et_MZP1200A0400_dm1TESUp = -999.0;
        
    mvaVar_tt_MZP600A0400_dm1TESUp = -999.0;
    mvaVar_tt_MZP800A0400_dm1TESUp = -999.0;
    mvaVar_tt_MZP1000A0400_dm1TESUp = -999.0;
    mvaVar_tt_MZP1200A0400_dm1TESUp = -999.0;

    mvaVar_mt_MZP600A0400_dm1TESDown = -999.0;
    mvaVar_mt_MZP800A0400_dm1TESDown = -999.0;
    mvaVar_mt_MZP1000A0400_dm1TESDown = -999.0;
    mvaVar_mt_MZP1200A0400_dm1TESDown = -999.0;

    mvaVar_et_MZP600A0400_dm1TESDown = -999.0;
    mvaVar_et_MZP800A0400_dm1TESDown = -999.0; 
    mvaVar_et_MZP1000A0400_dm1TESDown = -999.0;
    mvaVar_et_MZP1200A0400_dm1TESDown = -999.0;
        
    mvaVar_tt_MZP600A0400_dm1TESDown = -999.0;
    mvaVar_tt_MZP800A0400_dm1TESDown = -999.0;
    mvaVar_tt_MZP1000A0400_dm1TESDown = -999.0;
    mvaVar_tt_MZP1200A0400_dm1TESDown = -999.0;
    
    pt_1_dm10TESUp = -999.0;
    mt_1_dm10TESUp = -999.0;
    pfmt_1_dm10TESUp = -999.0;

    pt_2_dm10TESUp = -999.0;
    mt_2_dm10TESUp = -999.0;
    pfmt_2_dm10TESUp = -999.0;

    mt_tot_dm10TESUp = -999.0;
    pt_tt_dm10TESUp = -999.0;
    m_vis_dm10TESUp = -999.0;
    pfmet_type1_dm10TESUp_Pt = -999.0;
    P_chi_pf_dm10TESUp = -999.0;
    LPT_dm10TESUp = -999.0;
    DeltaR_leg1_leg2 = -999.0;
    cos_DeltaPhi_leg1_leg2 = -999.0;
    cos_DeltaPhi_PFMET_Higgs_dm10TESUp = -999.0;

    pt_1_dm10TESDown = -999.0;
    mt_1_dm10TESDown = -999.0;
    pfmt_1_dm10TESDown = -999.0;

    pt_2_dm10TESDown = -999.0;
    mt_2_dm10TESDown = -999.0;
    pfmt_2_dm10TESDown = -999.0;

    mt_tot_dm10TESDown = -999.0;
    pt_tt_dm10TESDown = -999.0;
    m_vis_dm10TESDown = -999.0;
    pfmet_type1_dm10TESDown_Pt = -999.0;
    P_chi_pf_dm10TESDown = -999.0;
    LPT_dm10TESDown = -999.0;
    DeltaR_leg1_leg2 = -999.0;
    cos_DeltaPhi_leg1_leg2 = -999.0;
    cos_DeltaPhi_PFMET_Higgs_dm10TESDown = -999.0;

    mvaVar_mt_MZP600A0400_dm10TESUp = -999.0;
    mvaVar_mt_MZP800A0400_dm10TESUp = -999.0;
    mvaVar_mt_MZP1000A0400_dm10TESUp = -999.0;
    mvaVar_mt_MZP1200A0400_dm10TESUp = -999.0;

    mvaVar_et_MZP600A0400_dm10TESUp = -999.0;
    mvaVar_et_MZP800A0400_dm10TESUp = -999.0; 
    mvaVar_et_MZP1000A0400_dm10TESUp = -999.0;
    mvaVar_et_MZP1200A0400_dm10TESUp = -999.0;
        
    mvaVar_tt_MZP600A0400_dm10TESUp = -999.0;
    mvaVar_tt_MZP800A0400_dm10TESUp = -999.0;
    mvaVar_tt_MZP1000A0400_dm10TESUp = -999.0;
    mvaVar_tt_MZP1200A0400_dm10TESUp = -999.0;

    mvaVar_mt_MZP600A0400_dm10TESDown = -999.0;
    mvaVar_mt_MZP800A0400_dm10TESDown = -999.0;
    mvaVar_mt_MZP1000A0400_dm10TESDown = -999.0;
    mvaVar_mt_MZP1200A0400_dm10TESDown = -999.0;

    mvaVar_et_MZP600A0400_dm10TESDown = -999.0;
    mvaVar_et_MZP800A0400_dm10TESDown = -999.0; 
    mvaVar_et_MZP1000A0400_dm10TESDown = -999.0;
    mvaVar_et_MZP1200A0400_dm10TESDown = -999.0;
        
    mvaVar_tt_MZP600A0400_dm10TESDown = -999.0;
    mvaVar_tt_MZP800A0400_dm10TESDown = -999.0;
    mvaVar_tt_MZP1000A0400_dm10TESDown = -999.0;
    mvaVar_tt_MZP1200A0400_dm10TESDown = -999.0;

	pt_tt = -999.0;
	DeltaR_leg1_leg2 = -999.0;
    cos_DeltaPhi_leg1_leg2 = -999.;
    cos_DeltaPhi_PFMET_Higgs = -999.;
    cos_DeltaPhi_PFMET_Higgs_UESUp = -999.;
    cos_DeltaPhi_PFMET_Higgs_UESDown = -999.;
    cos_DeltaPhi_PFMET_Higgs_JEnUp = -999.;
    cos_DeltaPhi_PFMET_Higgs_JEnDown = -999.;
    mt_tot = -999.0;
    weight_ttPtUp = -999.;
    weight_ttPtDown = -999.;
	mt_tot_UESUp = -999.0;
    mt_tot_UESDown = -999.0;
    mt_tot_JEnUp = -999.0;
    mt_tot_JEnDown = -999.0;
	m_vis = -999.0;
    m_vis_flat = -999.0;
	m_sv = -999.0;
	mt_sv = -999.0;
	SVFit_mvaMET_diTau_pt = -999.0;
	SVFit_mvaMET_diTau_eta = -999.0;
	SVFit_mvaMET_diTau_phi = -999.0;
	SVFit_mvaMET_FittedMET = -999.0;
	SVFit_mvaMET_FittedMETphi = -999.0;
	mvamet = -999.0;
	mvametphi = -999.0;
	met = -999.0;
	metphi = -999.0;
	puppimet = -999.0;
	puppimetphi = -999.0;
	uncorr_mvamet = -999.0;
	uncorr_mvametphi = -999.0;
	responseUP_mvaMET = -999.0;
	responseUP_mvaMETphi = -999.0;
	responseDOWN_mvaMET = -999.0;
	responseDOWN_mvaMETphi = -999.0;
	resolutionUP_mvaMET = -999.0;
	resolutionUP_mvaMETphi = -999.0;
	resolutionDOWN_mvaMET = -999.0;
	resolutionDOWN_mvaMETphi = -999.0;
	mvacov00 = -999.0;
	mvacov01 = -999.0;
	mvacov10 = -999.0;
	mvacov11 = -999.0;
	metcov00 = -999.0;
	metcov01 = -999.0;
	metcov10 = -999.0;
	metcov11 = -999.0;
    
	genMET = -999.0;
	genMETphi = -999.0;
	genMETeta = -999.0;
	genMETmass = -999.0;

	pfmet_raw_Pt = -999.0;
	pfmet_raw_Phi = -999.0;
	pfmet_raw_MT1 = -999.0;
	pfmet_raw_MT2 = -999.0;

	pfmet_type1_Pt = -999.0;
	pfmet_type1_Phi = -999.0;
	pfmet_type1_MT1 = -999.0;
	pfmet_type1_MT2 = -999.0;

	pfmet_type1_JetResUp_Pt = -999.0;
	pfmet_type1_JetResUp_Phi = -999.0;
	pfmet_type1_JetResUp_MT1 = -999.0;
	pfmet_type1_JetResUp_MT2 = -999.0;
	pfmet_type1_JetResDown_Pt = -999.0;
	pfmet_type1_JetResDown_Phi = -999.0;
	pfmet_type1_JetResDown_MT1 = -999.0;
	pfmet_type1_JetResDown_MT2 = -999.0;
	pfmet_type1_JetEnUp_Pt = -999.0;
	pfmet_type1_JetEnUp_Phi = -999.0;
	pfmet_type1_JetEnUp_MT1 = -999.0;
	pfmet_type1_JetEnUp_MT2 = -999.0;
	pfmet_type1_JetEnDown_Pt = -999.0;
	pfmet_type1_JetEnDown_Phi = -999.0;
	pfmet_type1_JetEnDown_MT1 = -999.0;
	pfmet_type1_JetEnDown_MT2 = -999.0;
	pfmet_type1_MuonEnUp_Pt = -999.0;
	pfmet_type1_MuonEnUp_Phi = -999.0;
	pfmet_type1_MuonEnUp_MT1 = -999.0;
	pfmet_type1_MuonEnUp_MT2 = -999.0;
	pfmet_type1_MuonEnDown_Pt = -999.0;
	pfmet_type1_MuonEnDown_Phi = -999.0;
	pfmet_type1_MuonEnDown_MT1 = -999.0;
	pfmet_type1_MuonEnDown_MT2 = -999.0;
	pfmet_type1_ElectronEnUp_Pt = -999.0;
	pfmet_type1_ElectronEnUp_Phi = -999.0;
	pfmet_type1_ElectronEnUp_MT1 = -999.0;
	pfmet_type1_ElectronEnUp_MT2 = -999.0;
	pfmet_type1_ElectronEnDown_Pt = -999.0;
	pfmet_type1_ElectronEnDown_Phi = -999.0;
	pfmet_type1_ElectronEnDown_MT1 = -999.0;
	pfmet_type1_ElectronEnDown_MT2 = -999.0;
	pfmet_type1_UnclusteredEnUp_Pt = -999.0;
	pfmet_type1_UnclusteredEnUp_Phi = -999.0;
	pfmet_type1_UnclusteredEnUp_MT1 = -999.0;
	pfmet_type1_UnclusteredEnUp_MT2 = -999.0;
	pfmet_type1_UnclusteredEnDown_Pt = -999.0;
	pfmet_type1_UnclusteredEnDown_Phi = -999.0;
	pfmet_type1_UnclusteredEnDown_MT1 = -999.0;
	pfmet_type1_UnclusteredEnDown_MT2 = -999.0;
	pfmet_type1_PhotonEnUp_Pt = -999.0;
	pfmet_type1_PhotonEnUp_Phi = -999.0;
	pfmet_type1_PhotonEnUp_MT1 = -999.0;
	pfmet_type1_PhotonEnUp_MT2 = -999.0;
	pfmet_type1_PhotonEnDown_Pt = -999.0;
	pfmet_type1_PhotonEnDown_Phi = -999.0;
	pfmet_type1_PhotonEnDown_MT1 = -999.0;
	pfmet_type1_PhotonEnDown_MT2 = -999.0;

	pzetavis = -999.0;
	pzetamiss = -999.0;
	pfpzetamiss = -999.0;
	puppipzetamiss = -999.0;
	pzetamiss_responseUP = -999.0;
	pzetamiss_responseDOWN = -999.0;
	pzetamiss_resolutionUP = -999.0;
	pzetamiss_resolutionDOWN = -999.0;
	pzetamiss_uncorr = -999.0;

	njets = -999;
	njetspt20 = -999;
	mjj = -999.0;
	jdeta = -999.0;
	njetingap = -999.0;
	njetingap20 = -999.0;
	jdphi = -999.0;
	jpt_1 = -999.0;
	jeta_1 = -999.0;
	jphi_1 = -999.0;
	jm_1 = -999.0;
	jmva_1 = -999.0;
	jpt_2 = -999.0;
	jeta_2 = -999.0;
	jphi_2 = -999.0;
	jm_2 = -999.0;
	jmva_2 = -999.0;
	nbtag = -999;
	nbtag_oneSigmaUp = -999;
	nbtag_oneSigmaDown = -999;
	bpt_1 = -999.0;
	beta_1 = -999.0;
	bphi_1 = -999.0;
	bm_1 = -999.0;
	bmva_1 = -999.0;
	bcsv_1 = -999.0;
	bpt_2 = -999.0;
	beta_2 = -999.0;
	bphi_2 = -999.0;
	bm_2 = -999.0;
	bmva_2 = -999.0;
	bcsv_2 = -999.0;
	nbtag_LooseWp = -999;
	nbtag_LooseWp_oneSigmaUp = -999;
	nbtag_LooseWp_oneSigmaDown = -999;
	bpt_1_LooseWp = -999.0;
	beta_1_LooseWp = -999.0;
	bphi_1_LooseWp = -999.0;
	bm_1_LooseWp = -999.0;
	bmva_1_LooseWp = -999.0;
	bcsv_1_LooseWp = -999.0;
	bpt_2_LooseWp = -999.0;
	beta_2_LooseWp = -999.0;
	bphi_2_LooseWp = -999.0;
	bm_2_LooseWp = -999.0;
	bmva_2_LooseWp = -999.0;
	bcsv_2_LooseWp = -999.0;
	nbtag_TightWp = -999;
	nbtag_TightWp_oneSigmaUp = -999;
	nbtag_TightWp_oneSigmaDown = -999;
	bpt_1_TightWp = -999.0;
	beta_1_TightWp = -999.0;
	bphi_1_TightWp = -999.0;
	bm_1_TightWp = -999.0;
	bmva_1_TightWp = -999.0;
	bcsv_1_TightWp = -999.0;
	bpt_2_TightWp = -999.0;
	beta_2_TightWp = -999.0;
	bphi_2_TightWp = -999.0;
	bm_2_TightWp = -999.0;
	bmva_2_TightWp = -999.0;
	bcsv_2_TightWp = -999.0;
	njets_JECshiftedUp = -999;
	njetspt20_JECshiftedUp = -999;
	mjj_JECshiftedUp = -999.0;
	jdeta_JECshiftedUp = -999.0;
	njetingap_JECshiftedUp = -999.0;
	njetingap20_JECshiftedUp = -999.0;
	jdphi_JECshiftedUp = -999.0;
	jpt_1_JECshiftedUp = -999.0;
	jeta_1_JECshiftedUp = -999.0;
	jphi_1_JECshiftedUp = -999.0;
	jm_1_JECshiftedUp = -999.0;
	jmva_1_JECshiftedUp = -999.0;
	jpt_2_JECshiftedUp = -999.0;
	jeta_2_JECshiftedUp = -999.0;
	jphi_2_JECshiftedUp = -999.0;
	jm_2_JECshiftedUp = -999.0;
	jmva_2_JECshiftedUp = -999.0;
	nbtag_JECshiftedUp = -999;
	bpt_1_JECshiftedUp = -999.0;
	beta_1_JECshiftedUp = -999.0;
	bphi_1_JECshiftedUp = -999.0;
	bm_1_JECshiftedUp = -999.0;
	bmva_1_JECshiftedUp = -999.0;
	bcsv_1_JECshiftedUp = -999.0;
	bpt_2_JECshiftedUp = -999.0;
	beta_2_JECshiftedUp = -999.0;
	bphi_2_JECshiftedUp = -999.0;
	bm_2_JECshiftedUp = -999.0;
	bmva_2_JECshiftedUp = -999.0;
	bcsv_2_JECshiftedUp = -999.0;
	nbtag_LooseWp_JECshiftedUp = -999;
	bpt_1_LooseWp_JECshiftedUp = -999.0;
	beta_1_LooseWp_JECshiftedUp = -999.0;
	bphi_1_LooseWp_JECshiftedUp = -999.0;
	bm_1_LooseWp_JECshiftedUp = -999.0;
	bmva_1_LooseWp_JECshiftedUp = -999.0;
	bcsv_1_LooseWp_JECshiftedUp = -999.0;
	bpt_2_LooseWp_JECshiftedUp = -999.0;
	beta_2_LooseWp_JECshiftedUp = -999.0;
	bphi_2_LooseWp_JECshiftedUp = -999.0;
	bm_2_LooseWp_JECshiftedUp = -999.0;
	bmva_2_LooseWp_JECshiftedUp = -999.0;
	bcsv_2_LooseWp_JECshiftedUp = -999.0;
	nbtag_TightWp_JECshiftedUp = -999;
	bpt_1_TightWp_JECshiftedUp = -999.0;
	beta_1_TightWp_JECshiftedUp = -999.0;
	bphi_1_TightWp_JECshiftedUp = -999.0;
	bm_1_TightWp_JECshiftedUp = -999.0;
	bmva_1_TightWp_JECshiftedUp = -999.0;
	bcsv_1_TightWp_JECshiftedUp = -999.0;
	bpt_2_TightWp_JECshiftedUp = -999.0;
	beta_2_TightWp_JECshiftedUp = -999.0;
	bphi_2_TightWp_JECshiftedUp = -999.0;
	bm_2_TightWp_JECshiftedUp = -999.0;
	bmva_2_TightWp_JECshiftedUp = -999.0;
	bcsv_2_TightWp_JECshiftedUp = -999.0;
	njets_JECshiftedDown = -999;
	njetspt20_JECshiftedDown = -999;
	mjj_JECshiftedDown = -999.0;
	jdeta_JECshiftedDown = -999.0;
	njetingap_JECshiftedDown = -999.0;
	njetingap20_JECshiftedDown = -999.0;
	jdphi_JECshiftedDown = -999.0;
	jpt_1_JECshiftedDown = -999.0;
	jeta_1_JECshiftedDown = -999.0;
	jphi_1_JECshiftedDown = -999.0;
	jm_1_JECshiftedDown = -999.0;
	jmva_1_JECshiftedDown = -999.0;
	jpt_2_JECshiftedDown = -999.0;
	jeta_2_JECshiftedDown = -999.0;
	jphi_2_JECshiftedDown = -999.0;
	jm_2_JECshiftedDown = -999.0;
	jmva_2_JECshiftedDown = -999.0;
	nbtag_JECshiftedDown = -999;
	bpt_1_JECshiftedDown = -999.0;
	beta_1_JECshiftedDown = -999.0;
	bphi_1_JECshiftedDown = -999.0;
	bm_1_JECshiftedDown = -999.0;
	bmva_1_JECshiftedDown = -999.0;
	bcsv_1_JECshiftedDown = -999.0;
	bpt_2_JECshiftedDown = -999.0;
	beta_2_JECshiftedDown = -999.0;
	bphi_2_JECshiftedDown = -999.0;
	bm_2_JECshiftedDown = -999.0;
	bmva_2_JECshiftedDown = -999.0;
	bcsv_2_JECshiftedDown = -999.0;
	nbtag_LooseWp_JECshiftedDown = -999;
	bpt_1_LooseWp_JECshiftedDown = -999.0;
	beta_1_LooseWp_JECshiftedDown = -999.0;
	bphi_1_LooseWp_JECshiftedDown = -999.0;
	bm_1_LooseWp_JECshiftedDown = -999.0;
	bmva_1_LooseWp_JECshiftedDown = -999.0;
	bcsv_1_LooseWp_JECshiftedDown = -999.0;
	bpt_2_LooseWp_JECshiftedDown = -999.0;
	beta_2_LooseWp_JECshiftedDown = -999.0;
	bphi_2_LooseWp_JECshiftedDown = -999.0;
	bm_2_LooseWp_JECshiftedDown = -999.0;
	bmva_2_LooseWp_JECshiftedDown = -999.0;
	bcsv_2_LooseWp_JECshiftedDown = -999.0;
	nbtag_TightWp_JECshiftedDown = -999;
	bpt_1_TightWp_JECshiftedDown = -999.0;
	beta_1_TightWp_JECshiftedDown = -999.0;
	bphi_1_TightWp_JECshiftedDown = -999.0;
	bm_1_TightWp_JECshiftedDown = -999.0;
	bmva_1_TightWp_JECshiftedDown = -999.0;
	bcsv_1_TightWp_JECshiftedDown = -999.0;
	bpt_2_TightWp_JECshiftedDown = -999.0;
	beta_2_TightWp_JECshiftedDown = -999.0;
	bphi_2_TightWp_JECshiftedDown = -999.0;
	bm_2_TightWp_JECshiftedDown = -999.0;
	bmva_2_TightWp_JECshiftedDown = -999.0;
	bcsv_2_TightWp_JECshiftedDown = -999.0;
	njets_JERup = -999;
	njetspt20_JERup = -999;
	mjj_JERup = -999.0;
	jdeta_JERup = -999.0;
	njetingap_JERup = -999.0;
	njetingap20_JERup = -999.0;
	jdphi_JERup = -999.0;
	jpt_1_JERup = -999.0;
	jeta_1_JERup = -999.0;
	jphi_1_JERup = -999.0;
	jm_1_JERup = -999.0;
	jmva_1_JERup = -999.0;
	jpt_2_JERup = -999.0;
	jeta_2_JERup = -999.0;
	jphi_2_JERup = -999.0;
	jm_2_JERup = -999.0;
	jmva_2_JERup = -999.0;
	nbtag_JERup = -999;
	bpt_1_JERup = -999.0;
	beta_1_JERup = -999.0;
	bphi_1_JERup = -999.0;
	bm_1_JERup = -999.0;
	bmva_1_JERup = -999.0;
	bcsv_1_JERup = -999.0;
	bpt_2_JERup = -999.0;
	beta_2_JERup = -999.0;
	bphi_2_JERup = -999.0;
	bm_2_JERup = -999.0;
	bmva_2_JERup = -999.0;
	bcsv_2_JERup = -999.0;
	nbtag_LooseWp_JERup = -999;
	bpt_1_LooseWp_JERup = -999.0;
	beta_1_LooseWp_JERup = -999.0;
	bphi_1_LooseWp_JERup = -999.0;
	bm_1_LooseWp_JERup = -999.0;
	bmva_1_LooseWp_JERup = -999.0;
	bcsv_1_LooseWp_JERup = -999.0;
	bpt_2_LooseWp_JERup = -999.0;
	beta_2_LooseWp_JERup = -999.0;
	bphi_2_LooseWp_JERup = -999.0;
	bm_2_LooseWp_JERup = -999.0;
	bmva_2_LooseWp_JERup = -999.0;
	bcsv_2_LooseWp_JERup = -999.0;
	nbtag_TightWp_JERup = -999;
	bpt_1_TightWp_JERup = -999.0;
	beta_1_TightWp_JERup = -999.0;
	bphi_1_TightWp_JERup = -999.0;
	bm_1_TightWp_JERup = -999.0;
	bmva_1_TightWp_JERup = -999.0;
	bcsv_1_TightWp_JERup = -999.0;
	bpt_2_TightWp_JERup = -999.0;
	beta_2_TightWp_JERup = -999.0;
	bphi_2_TightWp_JERup = -999.0;
	bm_2_TightWp_JERup = -999.0;
	bmva_2_TightWp_JERup = -999.0;
	bcsv_2_TightWp_JERup = -999.0;
	njets_JERdown = -999;
	njetspt20_JERdown = -999;
	mjj_JERdown = -999.0;
	jdeta_JERdown = -999.0;
	njetingap_JERdown = -999.0;
	njetingap20_JERdown = -999.0;
	jdphi_JERdown = -999.0;
	jpt_1_JERdown = -999.0;
	jeta_1_JERdown = -999.0;
	jphi_1_JERdown = -999.0;
	jm_1_JERdown = -999.0;
	jmva_1_JERdown = -999.0;
	jpt_2_JERdown = -999.0;
	jeta_2_JERdown = -999.0;
	jphi_2_JERdown = -999.0;
	jm_2_JERdown = -999.0;
	jmva_2_JERdown = -999.0;
	nbtag_JERdown = -999;
	bpt_1_JERdown = -999.0;
	beta_1_JERdown = -999.0;
	bphi_1_JERdown = -999.0;
	bm_1_JERdown = -999.0;
	bmva_1_JERdown = -999.0;
	bcsv_1_JERdown = -999.0;
	bpt_2_JERdown = -999.0;
	beta_2_JERdown = -999.0;
	bphi_2_JERdown = -999.0;
	bm_2_JERdown = -999.0;
	bmva_2_JERdown = -999.0;
	bcsv_2_JERdown = -999.0;
	nbtag_LooseWp_JERdown = -999;
	bpt_1_LooseWp_JERdown = -999.0;
	beta_1_LooseWp_JERdown = -999.0;
	bphi_1_LooseWp_JERdown = -999.0;
	bm_1_LooseWp_JERdown = -999.0;
	bmva_1_LooseWp_JERdown = -999.0;
	bcsv_1_LooseWp_JERdown = -999.0;
	bpt_2_LooseWp_JERdown = -999.0;
	beta_2_LooseWp_JERdown = -999.0;
	bphi_2_LooseWp_JERdown = -999.0;
	bm_2_LooseWp_JERdown = -999.0;
	bmva_2_LooseWp_JERdown = -999.0;
	bcsv_2_LooseWp_JERdown = -999.0;
	nbtag_TightWp_JERdown = -999;
	bpt_1_TightWp_JERdown = -999.0;
	beta_1_TightWp_JERdown = -999.0;
	bphi_1_TightWp_JERdown = -999.0;
	bm_1_TightWp_JERdown = -999.0;
	bmva_1_TightWp_JERdown = -999.0;
	bcsv_1_TightWp_JERdown = -999.0;
	bpt_2_TightWp_JERdown = -999.0;
	beta_2_TightWp_JERdown = -999.0;
	bphi_2_TightWp_JERdown = -999.0;
	bm_2_TightWp_JERdown = -999.0;
	bmva_2_TightWp_JERdown = -999.0;
	bcsv_2_TightWp_JERdown = -999.0;

    /* init the btag product SF to 1.0 */
    BtagEventSFproduct_looseWpDown = 1.0;
    BtagEventSFproduct_looseWpCentral = 1.0;
    BtagEventSFproduct_looseWpUp = 1.0;
    BtagEventSFproduct_mediumWpDown = 1.0;
    BtagEventSFproduct_mediumWpCentral = 1.0;
    BtagEventSFproduct_mediumWpUp = 1.0;
    BtagEventSFproduct_tightWpDown = 1.0;
    BtagEventSFproduct_tightWpCentral = 1.0;
    BtagEventSFproduct_tightWpUp = 1.0;
    BtagEventSFproduct_looseWpCentral_JECshiftedUp = 1.0;
    BtagEventSFproduct_mediumWpCentral_JECshiftedUp = 1.0;
    BtagEventSFproduct_tightWpCentral_JECshiftedUp = 1.0;
    BtagEventSFproduct_looseWpCentral_JECshiftedDown = 1.0;
    BtagEventSFproduct_mediumWpCentral_JECshiftedDown = 1.0;
    BtagEventSFproduct_tightWpCentral_JECshiftedDown = 1.0;
    BtagEventSFproduct_looseWpCentral_JERup = 1.0;
    BtagEventSFproduct_mediumWpCentral_JERup = 1.0;
    BtagEventSFproduct_tightWpCentral_JERup = 1.0;
    BtagEventSFproduct_looseWpCentral_JERdown = 1.0;
    BtagEventSFproduct_mediumWpCentral_JERdown = 1.0;
    BtagEventSFproduct_tightWpCentral_JERdown = 1.0;
    
    BtagEventSFproduct_or_DataTag_Central = 0.0;
    BtagEventSFproduct_or_DataTag_Up = 0.0;
    BtagEventSFproduct_or_DataTag_Down = 0.0;
    BtagEventSFproduct_or_DataTag_Central_JECshiftedUp = 0.0;
    BtagEventSFproduct_or_DataTag_Central_JECshiftedDown = 0.0;
    
	dilepton_veto = -999.0;
	extraelec_veto = -999.0;
	extramuon_veto = -999.0;

	HBHENoiseFilter = 0;
	HBHENoiseIsoFilter = 0;
	CSCTightHalo2015Filter = 0;
	EcalDeadCellTriggerPrimitiveFilter = 0;
	goodVerticesFilter = 0;
	eeBadScFilter = 0;
	chargedHadronTrackResolutionFilter = 0;
	muonBadTrackFilter = 0;
	globalTightHalo2016Filter = 0;
    BadChargedCandidateFilter = 0;
    BadPFMuonFilter = 0;

    BadMuonTaggedMoriond17 = 0;
    DuplicateMuonTaggedMoriond17 = 0;

	NUP = -999;
	weight = 1.0;
	generatorEventWeight = 1.0;
	lheHT = -999.0;
	lheOutGoingPartons = -999;
	lheZmass = -999.0;

	IsZTT = -999;
	IsZL = -999;
	IsZJ = -999;
	IsZLL = -999;
    
    IsTTT = 0;

	DataSet= "NULL";
	EventTotal = -999;
	EventType= "NULL";
	KeyName= "NULL";
	DataCard= "NULL";
	CrossSection = -999.0;
	FilterEff = -999.0;
	isSmallTree = 0;
	TauEsNumberSigmasShifted = -999.0;
	ElectronEsNumberSigmasShifted = -999.0;

	leg1_HLT_DoubleMediumIsoPFTau35_Trk1_eta2p1_Reg = 0.;
    leg1_HLT_DoubleMediumCombinedIsoPFTau35_Trk1_eta2p1_Reg = 0.;
	leg1_HLT_Ele25_eta2p1_WPTight_Gsf = 0.;
    leg1_HLT_Ele45_WPLoose_Gsf_L1JetTauSeeded = 0.;
	leg1_HLT_IsoMu22 = 0.;
	leg1_HLT_IsoTkMu22 = 0.;
    leg1_HLT_IsoMu22_eta2p1 = 0.;
    leg1_HLT_IsoTkMu22_eta2p1 = 0.;
    leg1_HLT_IsoMu24 = 0.;
    leg1_HLT_IsoTkMu24 = 0.;
    
	leg2_HLT_DoubleMediumIsoPFTau35_Trk1_eta2p1_Reg = 0.;
    leg2_HLT_DoubleMediumCombinedIsoPFTau35_Trk1_eta2p1_Reg = 0.;
	leg2_HLT_Ele25_eta2p1_WPTight_Gsf = 0.;
    leg2_HLT_Ele45_WPLoose_Gsf_L1JetTauSeeded = 0.;
	leg2_HLT_IsoMu22 = 0.;
	leg2_HLT_IsoTkMu22 = 0.;
    leg2_HLT_IsoMu22_eta2p1 = 0.;
    leg2_HLT_IsoTkMu22_eta2p1 = 0.;
    leg2_HLT_IsoMu24 = 0.;
    leg2_HLT_IsoTkMu24 = 0.;
    
    pairGoodForTrigger = 0;

	veto_leptonType.clear();
	veto_pt.clear();
	veto_eta.clear();
	veto_phi.clear();
	veto_M.clear();
	veto_charge.clear();
	veto_dxy.clear();
	veto_dz.clear();
	veto_RelIso.clear();
    veto_passesLooseMuonId.clear();
    veto_passesMediumMuonId_ICHEP16.clear();
    veto_passesMediumMuonId_Moriond17.clear();
	veto_passesMediumMuonId.clear();
    veto_passesTightMuonId.clear();
	veto_passElectronMVA80.clear();
	veto_passElectronMVA90.clear();
	veto_passVetoElectronCutBased.clear();
    veto_passTightElectronCutBased.clear();
	veto_isTrackerGlobalPFMuon.clear();
	veto_numberOfMissingInnerHits.clear();
	veto_numberOfMissingOuterHits.clear();
	veto_passConversionVeto.clear();
	veto_LeptonPassesThirdElectronVetoCuts.clear();
	veto_LeptonPassesThirdMuonVetoCuts.clear();
	veto_LeptonPassesDiElectronVetoCuts.clear();
	veto_LeptonPassesDiMuonVetoCuts.clear();

    LPT = -999.;

	//  expanded for sys. + alternate mets

    P_chi = -999.;
    M_min = -999.;
    
	P_chi_pf = -999.0;
    P_chi_pf_UESUp = -999.0;
    P_chi_pf_UESDown = -999.0;
    P_chi_pf_JEnUp = -999.0;
    P_chi_pf_JEnDown = -999.0;
	M_min_pf = -999.0;
	P_chi_puppi = -999.0;
	M_min_puppi = -999.0;
	P_chi_uncorr = -999.0;
	M_min_uncorr = -999.0;
	P_chi_responseUP = -999.0;
	M_min_responseUP = -999.0;
	P_chi_responseDOWN = -999.0;
	M_min_responseDOWN = -999.0;
	P_chi_resolutionUP = -999.0;
	M_min_resolutionUP = -999.0;
	P_chi_resolutionDOWN = -999.0;
	M_min_resolutionDOWN = -999.0;

}

std::vector<double> generateH2TauSyncTree::GetLeg1Leg2McTriggerWeights(TLorentzVector l1,TLorentzVector l2,
			int candType, int sysShift)
{
	std::vector<double> returnVec; /* index 0 is leg1 sf, index 1 is leg2 sf */
	returnVec.clear();
	return returnVec;

}

double generateH2TauSyncTree::CBeff(double pt, std::array <double, 5> par)
{
	/* note : arg must be ordered in this way :
	 m0    = par[0]
	 sigma = par[1] 
	 alpha = par[2]
	 n     = par[3]
	 norm =  par[4]
	*/

	return CBeff(pt, par[0], par[1], par[2] , par[3] , par[4]);

}



double generateH2TauSyncTree::CBeff(double pt, double m0, double sigma, double alpha, double n, double norm)
{
	/* from python in https://indico.cern.ch/event/489915/contributions/
	1167924/attachments/1227657/1798150/htt_12_2_2016.pdf#page=6*/
	
	double sqrtPiOver2 = std::sqrt(TMath::PiOver2());
	double sqrt2 = std::sqrt(2.);
	double sig = fabs(sigma);
	double t = (pt - m0)/sig * alpha /fabs(alpha);
	double absAlpha = fabs(alpha/sig); 

	double a = TMath::Power(n/absAlpha, n) * TMath::Exp(-0.5*absAlpha*absAlpha); 
	double b = absAlpha - n/absAlpha;
	double arg = absAlpha/sqrt2; 

	double ApproxErf = 1.;

	if (arg > 5.)  		  ApproxErf = 1.; 
	else if (arg < -5.)   ApproxErf = - 1.;
	else  				  ApproxErf = TMath::Erf(arg);

	double leftArea    = ( 1. +  ApproxErf) * sqrtPiOver2; 
	double rightArea   = ( a *1./TMath::Power(absAlpha-b, n-1) ) / (n-1) ;
	double area        = leftArea + rightArea; 


	if (t <= absAlpha )
	{
		arg = t / sqrt2;
		if (arg > 5.)        ApproxErf = 1.;
		else if (arg < -5.)  ApproxErf = -1.;
		else  			     ApproxErf = TMath::Erf(arg);
	
		return norm * (1+ ApproxErf) * sqrtPiOver2 / area; 
	}


	else
	{
		double num = norm * (leftArea + a * (1./TMath::Power(t-b,n-1) - 1./TMath::Power(absAlpha - b,n-1)) / (1-n));
		return (num / area);
	}
}



void generateH2TauSyncTree::initScaleFactorParametersRunII()
{

	// Run II efficiencies for HLT_DoubleMediumIsoPFTau35_Trk1_eta2p1_Reg

    Run2_TauTau_legTriggerEff_DataReal = {3.81919e+01 , 5.38746e+00 , 4.44730e+00 , 7.39646e+00 , 9.33402e-01};
    Run2_TauTau_legTriggerEff_DataFake = {3.90677e+01 , 7.03152e+00 , 1.11690e+01 , 1.29314e+00 , 9.99999e-01};

    Run2_TauTau_legTriggerEff_DataReal_dm0 = {39.04720616800577 ,7.264937515910671 ,7.481036846472791 ,2.2186432555766045 ,0.9979282164307409 };
    Run2_TauTau_legTriggerEff_DataReal_dm1 = {36.84199904252395 ,4.753143697071493 ,5.531763963597042 ,1.6898901594761204 ,0.9999999995735358 };
    Run2_TauTau_legTriggerEff_DataReal_dm10 = {40.75183952916798 ,5.3504479220811 ,2.313487662004634 ,9.97183433262574 ,0.9999999999882239 };
    Run2_TauTau_legTriggerEff_DataFake_dm0 = {37.12956370824129 ,6.983025697293738 ,5.363894244046358 ,2.4407568585097326 ,0.8911595667452936 };
    Run2_TauTau_legTriggerEff_DataFake_dm1 = {36.97927866559607 ,6.644502687462261 ,9.841922721381033 ,1.3676315442177307 ,0.9389735758585416 };
    Run2_TauTau_legTriggerEff_DataFake_dm10 = {41.43164872255737 ,6.705919352986819 ,4.1052875425000055 ,5.778073936840931 ,0.9338411883350912 };
    
    Run2_TauTau_legTriggerEff_MCReal_dm0 = {38.15007577090958 ,6.311019635208034 ,5.494801948596283 ,2.63419150959934 ,0.937851738589565 };
    Run2_TauTau_legTriggerEff_MCReal_dm1 = {36.3666829038776 ,4.5131284656300705 ,4.5285997287953785 ,1.6296625819252648 ,0.9999999997812011 };
    Run2_TauTau_legTriggerEff_MCReal_dm10 = {39.50277114020348 ,4.853831493248512 ,1.4594666703565793 ,117.98597869527406 ,0.9851724068805524 };
    Run2_TauTau_legTriggerEff_MCFake_dm0 = {35.32703569423916 ,5.691102365767336 ,3.9167310355092373 ,1.5068105014095727 ,0.9999999531102266 };
    Run2_TauTau_legTriggerEff_MCFake_dm1 = {36.250161105296414 ,4.8515088360634815 ,4.712192817046861 ,1.5331124616660115 ,0.9999999999586987 };
    Run2_TauTau_legTriggerEff_MCFake_dm10 = {41.85875828593461 ,6.291301353413331 ,4.38982391097697 ,2.4534222687442755 ,0.9999998258581441 };

	Run2_TauTau_legTriggerEff_Data = {3.45412e+01 , 5.63353e+00 , 2.49242e+00 , 3.35896e+00 , 1.00000e+00};
	Run2_TauTau_legTriggerEff_DataUP = {3.31713e+01 , 5.66551e+00 , 1.87175e+00 , 8.07790e+00 , 1.00000e+00};
	Run2_TauTau_legTriggerEff_DataDOWN = {3.56264e+01 , 5.30711e+00 , 2.81591e+00 , 2.40649e+00 , 9.99958e-01};
	Run2_TauTau_legTriggerEff_Mc = {3.60274e+01 , 5.89434e+00 , 5.82870e+00 , 1.83737e+00 , 9.58000e-01};
	Run2_TauTau_legTriggerEff_McUP = {3.56012e+01 , 5.97209e+00 , 6.09604e+00 , 1.68740e+00 , 9.87653e-01};
	Run2_TauTau_legTriggerEff_McDOWN = {3.62436e+01 , 5.58461e+00 , 5.12924e+00 , 2.05921e+00 , 9.32305e-01};
}


double generateH2TauSyncTree::GetTransverseMass(TLorentzVector L, TLorentzVector T)
{
  double pTxMET = sqrt(L.X()*L.X()+L.Y()*L.Y())*sqrt(T.X()*T.X()+T.Y()*T.Y());
  double CosDphi = cos(L.DeltaPhi(T));
  double MtSq = (2 * pTxMET*(1-CosDphi));
  return sqrt(MtSq);

}

void generateH2TauSyncTree::fillJetBranches(bool eventHasNominalLeptonEnergyScales, bool eventIsNotSmallTree,
											std::string variantString, jetDescription &jetDesc)
{

	////////////////////////////////////////////////////////////////////////
	// based on variantString, figure out which jet collection to access  //

	assert( variantString == "" ||\
	   variantString == "_JECshiftedUp" ||\
	   variantString == "_JECshiftedDown" ||\
	   variantString == "_JERup" ||\
	   variantString == "_JERdown" );

  	/////////////////////////////////////////////////////
  	//  JETS + btagging  + alternate btags      //
  	/////////////////////////////////////////////////////


	/* these jet variables are accessible in FlatTuple only 
	   if dealing with nominal jets OR nominal Tau/Electron Energy Scale */

	std::vector<double> jets_pt; 						jets_pt.clear();
	std::vector<double> jets_eta; 						jets_eta.clear();
	std::vector<double> jets_phi; 						jets_phi.clear();
	std::vector<double> jets_M; 						jets_M.clear();
	std::vector<double> jets_IsBTagged_LooseWpCentral; 	jets_IsBTagged_LooseWpCentral.clear();
	std::vector<double> jets_IsBTagged_MediumWpCentral; jets_IsBTagged_MediumWpCentral.clear();
	std::vector<double> jets_IsBTagged_TightWpCentral; 	jets_IsBTagged_TightWpCentral.clear();


	if( eventHasNominalLeptonEnergyScales || variantString=="" )
	{

		jets_pt 						=  R.getVD("jets_pt"+variantString);
		jets_eta 						=  R.getVD("jets_eta"+variantString);
		jets_phi 						=  R.getVD("jets_phi"+variantString);
		jets_M 							=  R.getVD("jets_M"+variantString);
		jets_IsBTagged_LooseWpCentral 	=  R.getVD("jets_IsBTagged_LooseWpCentral"+variantString);
		jets_IsBTagged_MediumWpCentral	=  R.getVD("jets_IsBTagged_MediumWpCentral"+variantString);
		jets_IsBTagged_TightWpCentral 	=  R.getVD("jets_IsBTagged_TightWpCentral"+variantString);

	}


	/* these jet variables are accessible in FlatTuple only if dealing 
	   with nominal jets AND nominal Tau/Electron Energy Scale */


	std::vector<double>	jets_IsBTagged_LooseWpUp;	 jets_IsBTagged_LooseWpUp.clear();
	std::vector<double> jets_IsBTagged_LooseWpDown;	 jets_IsBTagged_LooseWpDown.clear();
	std::vector<double>	jets_IsBTagged_MediumWpUp;	 jets_IsBTagged_MediumWpUp.clear();
	std::vector<double>	jets_IsBTagged_MediumWpDown; jets_IsBTagged_MediumWpDown.clear();
	std::vector<double>	jets_IsBTagged_TightWpUp;	 jets_IsBTagged_TightWpUp.clear();
	std::vector<double>	jets_IsBTagged_TightWpDown;	 jets_IsBTagged_TightWpDown.clear();


	/* only kept b-tag shifts in nominal case */
	if( eventHasNominalLeptonEnergyScales && variantString=="" )
	{
		jets_IsBTagged_LooseWpUp  	= R.getVD("jets_IsBTagged_LooseWpUp"+variantString);
		jets_IsBTagged_LooseWpDown  = R.getVD("jets_IsBTagged_LooseWpDown"+variantString);
		jets_IsBTagged_MediumWpUp  	= R.getVD("jets_IsBTagged_MediumWpUp"+variantString);
		jets_IsBTagged_MediumWpDown = R.getVD("jets_IsBTagged_MediumWpDown"+variantString);
		jets_IsBTagged_TightWpUp  	= R.getVD("jets_IsBTagged_TightWpUp"+variantString);
		jets_IsBTagged_TightWpDown  = R.getVD("jets_IsBTagged_TightWpDown"+variantString);
	}

	/* these jet variables are accessible in FlatTuple only if dealing 
	   with nominal jets */


	std::vector<float> jets_defaultBtagAlgorithm_RawScore; 		jets_defaultBtagAlgorithm_RawScore.clear();
	std::vector<double> jets_PU_jetIdRaw; 						jets_PU_jetIdRaw.clear();

	if( variantString=="" )
	{
		jets_defaultBtagAlgorithm_RawScore 	=  R.getVF("jets_defaultBtagAlgorithm_RawScore"+variantString);
		jets_PU_jetIdRaw 					=  R.getVD("jets_PU_jetIdRaw"+variantString);
	}

	
	/* counters for nominal JEC and JER jets */
	/* reset the b-tag counters to zero since we are counting them now */

	jetDesc.m_njets				= 0;
	jetDesc.m_njetspt20			= 0;

	jetDesc.m_njetingap 		= 0;	// Only filled if njetspt20>=2 
	jetDesc.m_njetingap20  		= 0; 	// Only filled if njetspt20>=2 

	/* b-tag counters */

   	jetDesc.m_nbtag 				= 0; // (medium WP central)
    jetDesc.m_nbtag_oneSigmaUp 		= 0; // (medium WP UP : one sigma + shift in btag sys)
    jetDesc.m_nbtag_oneSigmaDown	= 0; // (medium WP DOWN : one sigma - shift in btag sys)

	jetDesc.m_nbtag_LooseWp					= 0;  // (loose WP central)
	jetDesc.m_nbtag_LooseWp_oneSigmaUp		= 0;  // (loose WP UP : one sigma + shift in btag sys)
	jetDesc.m_nbtag_LooseWp_oneSigmaDown	= 0;  // (loose WP DOWN : one sigma - shift in btag sys)

	jetDesc.m_nbtag_TightWp					= 0;  // (tight WP central)
	jetDesc.m_nbtag_TightWp_oneSigmaUp 		= 0;  // (tight WP UP : one sigma + shift in btag sys)
	jetDesc.m_nbtag_TightWp_oneSigmaDown	= 0;  // (tight WP DOWN : one sigma - shift in btag sys)

	/* vectors to hold B-tagged jet pair index */
	
	std::vector < std::size_t > Bjet_pair;			// (medium Wp central)
	Bjet_pair.clear();

	std::vector < std::size_t > Bjet_pairUP;		// (medium Wp central + 1 sigma shift UP on btag sys.)
	Bjet_pairUP.clear();

	std::vector < std::size_t > Bjet_pairDOWN;		// (medium Wp central - 1 sigma shift DOWN on btag sys.)
	Bjet_pairDOWN.clear();

	std::vector < std::size_t > Bjet_pair_looseWp;	// (loose Wp central)
	Bjet_pair_looseWp.clear();

	std::vector < std::size_t > Bjet_pair_looseWpUP;  // (loose Wp central shift UP)
	Bjet_pair_looseWpUP.clear();

	std::vector < std::size_t > Bjet_pair_looseWpDOWN; // (loose Wp central shift DOWN)
	Bjet_pair_looseWpDOWN.clear();

	std::vector < std::size_t > Bjet_pair_tightWp;		// (tight Wp central)
	Bjet_pair_tightWp.clear();

	std::vector < std::size_t > Bjet_pair_tightWpUP;	// (tight Wp central shift UP)
	Bjet_pair_tightWpUP.clear();

	std::vector < std::size_t > Bjet_pair_tightWpDOWN;	// (tight Wp central shift DOWN)
	Bjet_pair_tightWpDOWN.clear();


	/* begin a loop on the nominal jets to figure out b-tagging */

	for(std::size_t i = 0; i<jets_pt.size();++i)
 	{ 		

 		if( fabs(jets_eta[i]) < 4.7 && jets_pt[i] > 30)  jetDesc.m_njets++;		
		if( fabs(jets_eta[i]) < 4.7 && jets_pt[i] > 20)  jetDesc.m_njetspt20++;

        // bjets at pt 30 for monohiggs
 		if( fabs(jets_eta[i]) < 2.4 && jets_pt[i] > 30)
 		{	
 			if( eventHasNominalLeptonEnergyScales || variantString=="" )
 			{
	 			if(jets_IsBTagged_LooseWpCentral[i] > 0.5)	Bjet_pair_looseWp.push_back(i);
	 			if(jets_IsBTagged_MediumWpCentral[i] > 0.5)	Bjet_pair.push_back(i);
	 			if(jets_IsBTagged_TightWpCentral[i] > 0.5)	Bjet_pair_tightWp.push_back(i);
	 		}

			if( eventHasNominalLeptonEnergyScales  && variantString=="" )
			{
				if(jets_IsBTagged_LooseWpUp[i] > 0.5)		Bjet_pair_looseWpUP.push_back(i);  	
				if(jets_IsBTagged_LooseWpDown[i] > 0.5)		Bjet_pair_looseWpDOWN.push_back(i);  
				if(jets_IsBTagged_MediumWpUp[i] > 0.5)		Bjet_pairUP.push_back(i);  	
				if(jets_IsBTagged_MediumWpDown[i] > 0.5)	Bjet_pairDOWN.push_back(i); 
				if(jets_IsBTagged_TightWpUp[i] > 0.5)		Bjet_pair_tightWpUP.push_back(i);  	
				if(jets_IsBTagged_TightWpDown[i] > 0.5)		Bjet_pair_tightWpDOWN.push_back(i);  
			}
		}
 	}


 	/* begin computation of branches for nomina jets & b-tags */


	//////////////////////////////////////////////////////////////////////////////
	// some 4-vectors we will need

	TLorentzVector   j1_20(0.,0.,0.,0.);
	TLorentzVector   j2_20(0.,0.,0.,0.);

	if(jets_pt.size() >= 1) 
	{
		j1_20.SetPtEtaPhiM(jets_pt[0],jets_eta[0],jets_phi[0],jets_M[0]);
		jetDesc.m_jpt_1 = jets_pt[0];
		jetDesc.m_jeta_1 = jets_eta[0];
		jetDesc.m_jphi_1 = jets_phi[0];
		jetDesc.m_jm_1 = jets_M[0];
		if( variantString=="" ) jetDesc.m_jmva_1 = jets_PU_jetIdRaw[0];
	}


	if(jets_pt.size() >= 2)
	{	
		j2_20.SetPtEtaPhiM(jets_pt[1],jets_eta[1],jets_phi[1],jets_M[1]);
		jetDesc.m_jpt_2 = jets_pt[1];
		jetDesc.m_jeta_2 = jets_eta[1];
		jetDesc.m_jphi_2 = jets_phi[1];
		jetDesc.m_jm_2 = jets_M[1];
		if( variantString=="" ) jetDesc.m_jmva_2 = jets_PU_jetIdRaw[1];


		jetDesc.m_mjj = (j1_20+j2_20).M();
	 	jetDesc.m_jdeta = (j1_20.Eta()-j2_20.Eta());
	 	jetDesc.m_jdphi = (j1_20.Phi()-j2_20.Phi());



		 if(jets_eta.size()>2)
		 {
		 	for(std::size_t i = 2; i<jets_eta.size();++i)
		 	{
		 		double current_eta = jets_eta[i];
		 		double low_bound = std::min(jets_eta[0],jets_eta[1]);
		 		double high_bound = std::max(jets_eta[0],jets_eta[1]);

		 		if(current_eta>low_bound && high_bound>current_eta)
		 		{
		 			if(jets_pt[i]>20) jetDesc.m_njetingap20++;
		 			if(jets_pt[i]>30) jetDesc.m_njetingap++;
		 		}
		 	}
		 }
	} // at least 2 jets

	/* handle the b-tag counts, since b-tag counts are small and +, ignore std::size_t vs int difference here */

   	jetDesc.m_nbtag 						= Bjet_pair.size();
    jetDesc.m_nbtag_oneSigmaUp 				= Bjet_pairUP.size();
    jetDesc.m_nbtag_oneSigmaDown			= Bjet_pairDOWN.size();

	jetDesc.m_nbtag_LooseWp					= Bjet_pair_looseWp.size();		
	jetDesc.m_nbtag_LooseWp_oneSigmaUp		= Bjet_pair_looseWpUP.size();
	jetDesc.m_nbtag_LooseWp_oneSigmaDown	= Bjet_pair_looseWpDOWN.size();

	jetDesc.m_nbtag_TightWp					= Bjet_pair_tightWp.size();
	jetDesc.m_nbtag_TightWp_oneSigmaUp 		= Bjet_pair_tightWpUP.size();
	jetDesc.m_nbtag_TightWp_oneSigmaDown	= Bjet_pair_tightWpDOWN.size();


	/* handle b-jet 1 and 2 quantities */

	if(Bjet_pair_looseWp.size() > 0)
	{
		jetDesc.m_bpt_1_LooseWp 		= jets_pt[Bjet_pair_looseWp[0]];
		jetDesc.m_beta_1_LooseWp 		= jets_eta[Bjet_pair_looseWp[0]];
		jetDesc.m_bphi_1_LooseWp 		= jets_phi[Bjet_pair_looseWp[0]];
		jetDesc.m_bm_1_LooseWp 			= jets_M[Bjet_pair_looseWp[0]];
		

		if( variantString=="" ) 
		{
			jetDesc.m_bmva_1_LooseWp 		= jets_PU_jetIdRaw[Bjet_pair_looseWp[0]];
			jetDesc.m_bcsv_1_LooseWp 		= jets_defaultBtagAlgorithm_RawScore[Bjet_pair_looseWp[0]];
		}	

		if(Bjet_pair_looseWp.size() > 1)
		{
			jetDesc.m_bpt_2_LooseWp 		= jets_pt[Bjet_pair_looseWp[1]];
			jetDesc.m_beta_2_LooseWp 		= jets_eta[Bjet_pair_looseWp[1]];
			jetDesc.m_bphi_2_LooseWp 		= jets_phi[Bjet_pair_looseWp[1]];
			jetDesc.m_bm_2_LooseWp 			= jets_M[Bjet_pair_looseWp[1]];
			if( variantString=="" ) 
			{
				jetDesc.m_bmva_2_LooseWp 		= jets_PU_jetIdRaw[Bjet_pair_looseWp[1]];
				jetDesc.m_bcsv_2_LooseWp 		= jets_defaultBtagAlgorithm_RawScore[Bjet_pair_looseWp[1]];
			}
		} 
	}
	if(Bjet_pair.size() > 0)
	{
		jetDesc.m_bpt_1 		= jets_pt[Bjet_pair[0]];
		jetDesc.m_beta_1 		= jets_eta[Bjet_pair[0]];
		jetDesc.m_bphi_1 		= jets_phi[Bjet_pair[0]];
		jetDesc.m_bm_1 			= jets_M[Bjet_pair[0]];
		if( variantString=="" ) 			{
			jetDesc.m_bmva_1 		= jets_PU_jetIdRaw[Bjet_pair[0]];
			jetDesc.m_bcsv_1 		= jets_defaultBtagAlgorithm_RawScore[Bjet_pair[0]];
		}	

		if(Bjet_pair.size() > 1)
		{

			jetDesc.m_bpt_2 		= jets_pt[Bjet_pair[1]];
			jetDesc.m_beta_2 		= jets_eta[Bjet_pair[1]];
			jetDesc.m_bphi_2 		= jets_phi[Bjet_pair[1]];
			jetDesc.m_bm_2 			= jets_M[Bjet_pair[1]];
			if( variantString=="" ) 
			{
				jetDesc.m_bmva_2 		= jets_PU_jetIdRaw[Bjet_pair[1]];
				jetDesc.m_bcsv_2 		= jets_defaultBtagAlgorithm_RawScore[Bjet_pair[1]];
			}
		} 
	}
	if(Bjet_pair_tightWp.size() > 0)
	{

		jetDesc.m_bpt_1_TightWp 		= jets_pt[Bjet_pair_tightWp[0]];
		jetDesc.m_beta_1_TightWp 		= jets_eta[Bjet_pair_tightWp[0]];
		jetDesc.m_bphi_1_TightWp 		= jets_phi[Bjet_pair_tightWp[0]];
		jetDesc.m_bm_1_TightWp 			= jets_M[Bjet_pair_tightWp[0]];
		if( variantString=="" ) 
		{
			jetDesc.m_bmva_1_TightWp 		= jets_PU_jetIdRaw[Bjet_pair_tightWp[0]];
			jetDesc.m_bcsv_1_TightWp 		= jets_defaultBtagAlgorithm_RawScore[Bjet_pair_tightWp[0]];
		}

		if(Bjet_pair_tightWp.size() > 1)
		{

			jetDesc.m_bpt_2_TightWp 		= jets_pt[Bjet_pair_tightWp[1]];
			jetDesc.m_beta_2_TightWp 		= jets_eta[Bjet_pair_tightWp[1]];
			jetDesc.m_bphi_2_TightWp 		= jets_phi[Bjet_pair_tightWp[1]];
			jetDesc.m_bm_2_TightWp 			= jets_M[Bjet_pair_tightWp[1]];
			if( variantString=="" ) 
			{
				jetDesc.m_bmva_2_TightWp 		= jets_PU_jetIdRaw[Bjet_pair_tightWp[1]];
				jetDesc.m_bcsv_2_TightWp 		= jets_defaultBtagAlgorithm_RawScore[Bjet_pair_tightWp[1]];
			}
		} 
	}	
	
}


double generateH2TauSyncTree::getFinalWeight(bool verbose_, TLorentzVector l1, TLorentzVector l2)
{

	// init weight

	double returnWeight_ = 1.0;

	/* include the nominal cross section / events weight */

	returnWeight_ *= getNominalWeight(verbose_);

	/* for MC include the puWeight */

	if(R.getB("isRealData") == 0)
	{
		returnWeight_ *= R.getD("puWeight");
	}

	/* include top pt rewight (this returns 1.0 for non ttbar sample ) */

	//returnWeight_ *=  getTopQuarkPtWeight(verbose_);

	/* include V pt reweights for W,DY samples. Monojet k factors for W, NLO weight for Z*/

    returnWeight_ *= getKFactor(verbose_);
    
    returnWeight_ *= getZReWeight(verbose_);
    
    returnWeight_ *= getWWFactor(0);
    
    returnWeight_ *= getZZFactor(0);
    
	/* include susy ggH NLO weight (1.0 for non valid samples ) 
		susy signal samples use 10 as tan beta --- need to double check on this
	*/
    //only for SUSY gluglu samples
	//returnWeight_ *= getNLOReWeight(verbose_, 10 );

    /* include anti-lepton disc factor for fake taus */
    
    returnWeight_ *= getALDScaleFactors(verbose_);

	/* include trigger x id x iso scale factor */
    ///useMuonCentral ON
	returnWeight_ *= getFinalScaleFactorsForPair(0,0,1,0,0,0,l1,l2);

	return returnWeight_;
}


/* function: getNominalWeight 
		-- nomial weight = (1000 * cross-section) / (crab_job_eff * EventTotal * FilterEff) in pb
		-- valid for all samples except data, W+jets, DY+jets, (returns to 1.0, for data)
		-- for DY and W+jets returns stiching weights		
		-- for SUSY Higgs & mono-Higgs the cross-sections are set to 1.0 			
*/

double generateH2TauSyncTree::getNominalWeight(bool verbose_)
{

	/* return 1.0 for data */

	if(R.getB("isRealData") == 1) 
	{
		if(verbose_) std::cout<<" real data, nominal weight returned as 1.0 \n";
		return 1.0;
	}

    bool stitch = 0;
    double returnWeight_ = 1.;
    
    /* include generator event weight for AMCNLO */
    
    returnWeight_ *= R.getD("generatorEventWeight");


    //including LOtoNLO for DY
    if(R.getS("KeyName") == "DYJetsToLL_M-50" ||\
    R.getS("KeyName") == "DYJetsToLL_M-50ext1-v2" ||\
    R.getS("KeyName") == "DY1JetsToLL_M-50" ||\
    R.getS("KeyName") == "DY2JetsToLL_M-50" ||\
    R.getS("KeyName") == "DY3JetsToLL_M-50"  ||\
    R.getS("KeyName") == "DY4JetsToLL_M-50")
    {
        stitch = 1;
        int outgoingJets_ = R.getI("lheOutGoingPartons");
        
        if (outgoingJets_ == 0 ) returnWeight_ = (0.039542/(R.getD("FilterEff")));
        if (outgoingJets_ == 1 ) returnWeight_ = (0.012748/(R.getD("FilterEff")));
        if (outgoingJets_ == 2 ) returnWeight_ = (0.013013/(R.getD("FilterEff")));
        if (outgoingJets_ == 3 ) returnWeight_ = (0.013380/(R.getD("FilterEff")));
        if (outgoingJets_ == 4 ) returnWeight_ = (0.010970/(R.getD("FilterEff")));
    }
    
    else if (R.getS("KeyName") == "WJetsToLNu" ||\
    R.getS("KeyName") == "WJetsToLNuext2-v1" ||\
    R.getS("KeyName") == "WJetsToLNu_HT-100To200orig" ||\
    R.getS("KeyName") == "WJetsToLNu_HT-100To200" ||\
    R.getS("KeyName") == "WJetsToLNu_HT-100To200ext2-v1" ||\
    R.getS("KeyName") == "WJetsToLNu_HT-200To400orig" ||\
    R.getS("KeyName") == "WJetsToLNu_HT-200To400ext1-v1" ||\
    R.getS("KeyName") == "WJetsToLNu_HT-200To400" ||\
    R.getS("KeyName") == "WJetsToLNu_HT-400To600orig" ||\
    R.getS("KeyName") == "WJetsToLNu_HT-400To600"  ||\
    R.getS("KeyName") == "WJetsToLNu_HT-600To800orig" ||\
    R.getS("KeyName") == "WJetsToLNu_HT-600To800" ||\
    R.getS("KeyName") == "WJetsToLNu_HT-800To1200orig" ||\
    R.getS("KeyName") == "WJetsToLNu_HT-800To1200" ||\
    R.getS("KeyName") == "WJetsToLNu_HT-1200To2500orig" ||\
    R.getS("KeyName") == "WJetsToLNu_HT-1200To2500" ||\
    R.getS("KeyName") == "WJetsToLNu_HT-2500ToInforig" ||\
    R.getS("KeyName") == "WJetsToLNu_HT-2500ToInf")
    {
        stitch = 1.0;
        if (R.getD("lheHT") < 100.0) {returnWeight_ = (0.5808909/(R.getD("FilterEff")));}
        else if (R.getD("lheHT") > 100.0 && R.getD("lheHT") < 200.0) {returnWeight_ = (0.0164683/(R.getD("FilterEff")));}
        else if (R.getD("lheHT") > 200.0 && R.getD("lheHT") < 400.0) {returnWeight_ = (0.0089255/(R.getD("FilterEff")));}
        else if (R.getD("lheHT") > 400.0 && R.getD("lheHT") < 600.0) {returnWeight_ = (0.0062354/(R.getD("FilterEff")));}
        else if (R.getD("lheHT") > 600.0 && R.getD("lheHT") < 800.0) {returnWeight_ = (0.0006453/(R.getD("FilterEff")));}
        else if (R.getD("lheHT") > 800.0 && R.getD("lheHT") < 1200.0) {returnWeight_ = (0.0007118/(R.getD("FilterEff")));}
        else if (R.getD("lheHT") > 1200.0 && R.getD("lheHT") < 2500.0) {returnWeight_ = (0.0001934/(R.getD("FilterEff")));}
        else if (R.getD("lheHT") > 2500.0)  {returnWeight_ = (0.0000122/(R.getD("FilterEff")));}
    }
    
    else if (R.getS("KeyName") == "EWKWMinus2Jets_WToLNu_M-50" ||\
    R.getS("KeyName") == "EWKWMinus2Jets_WToLNu_M-50--ext1-v1" ||\
    R.getS("KeyName") == "EWKWMinus2Jets_WToLNu_M-50--ext2-v1")
    {
        stitch = 1.0;
        returnWeight_ = 0.00468/(R.getD("FilterEff"));
    }
    
    else if (R.getS("KeyName") == "EWKWPlus2Jets_WToLNu_M-50" ||\
    R.getS("KeyName") == "EWKWPlus2Jets_WToLNu_M-50--ext1-v1" ||\
    R.getS("KeyName") == "EWKWPlus2Jets_WToLNu_M-50--ext2-v1")
    {
        stitch = 1.0;
        returnWeight_ = 0.00519/(R.getD("FilterEff"));
    }
    
    else if (R.getS("KeyName") == "EWKZ2Jets_ZToLL_M-50" ||\
    R.getS("KeyName") == "EWKZ2Jets_ZToLL_M-50--ext1-v1" ||\
    R.getS("KeyName") == "EWKZ2Jets_ZToLL_M-50--ext2-v1")
    {
        stitch = 1.0;
        returnWeight_ = 0.00399/(R.getD("FilterEff"));
    }
    
    else if (R.getS("KeyName") == "EWKZ2Jets_ZToNuNu" ||\
    R.getS("KeyName") == "EWKZ2Jets_ZToNuNu--ext1-v1")
    {
        stitch = 1.0;
        returnWeight_ = 0.01677/(R.getD("FilterEff"));
    }

    /* Stitch with negative weight sums */
    else if (R.getS("KeyName") == "WWTo1L1Nu2Q")
    {
        stitch = 1.0;
        returnWeight_ = (1000.0 * R.getD("CrossSection") * R.getD("generatorEventWeight")) / (435302995.93 * R.getD("FilterEff"));
    }
    else if (R.getS("KeyName") == "WZTo1L1Nu2Q")
    {
        stitch = 1.0;
        returnWeight_ = (1000.0 * R.getD("CrossSection") * R.getD("generatorEventWeight")) / (420511499.74 * R.getD("FilterEff"));
    }
    else if (R.getS("KeyName") == "WZTo1L3Nu")
    {
        stitch = 1.0;
        returnWeight_ = (1000.0 * R.getD("CrossSection") * R.getD("generatorEventWeight")) / (9357734.34 * R.getD("FilterEff"));
    }
    else if (R.getS("KeyName") == "WZTo2L2Q")
    {
        stitch = 1.0;
        returnWeight_ = (1000.0 * R.getD("CrossSection") * R.getD("generatorEventWeight")) / (234019002.42 * R.getD("FilterEff"));
    }
    else if (R.getS("KeyName") == "ZZTo2L2Q")
    {
        stitch = 1.0;
        returnWeight_ = (1000.0 * R.getD("CrossSection") * R.getD("generatorEventWeight")) / (77855498.35 * R.getD("FilterEff"));
    }
    else if (R.getS("KeyName") == "ZZTo2Q2Nu")
    {
        stitch = 1.0;
        returnWeight_ = (1000.0 * R.getD("CrossSection") * R.getD("generatorEventWeight")) / (196347735.89 * R.getD("FilterEff"));
    }
    else if (R.getS("KeyName") == "ZZTo4L--ext1-v1")
    {
        stitch = 1.0;
        returnWeight_ = (1000.0 * R.getD("CrossSection") * R.getD("generatorEventWeight")) / (20476437.37 * R.getD("FilterEff"));
    }
    else if (R.getS("KeyName") == "WWW_4F")
    {
        stitch = 1.0;
        returnWeight_ = (1000.0 * R.getD("CrossSection") * R.getD("generatorEventWeight")) / (50012.91 * R.getD("FilterEff"));
    }
    else if (R.getS("KeyName") == "WWZ")
    {
        stitch = 1.0;
        returnWeight_ = (1000.0 * R.getD("CrossSection") * R.getD("generatorEventWeight")) / (41171.81 * R.getD("FilterEff"));
    }
    else if (R.getS("KeyName") == "WZZ")
    {
        stitch = 1.0;
        returnWeight_ = (1000.0 * R.getD("CrossSection") * R.getD("generatorEventWeight")) / (13736.41 * R.getD("FilterEff"));
    }
    else if (R.getS("KeyName") == "ZZZ")
    {
        stitch = 1.0;
        returnWeight_ = (1000.0 * R.getD("CrossSection") * R.getD("generatorEventWeight")) / (3499.44 * R.getD("FilterEff"));
    }
    
    //new
    else if (R.getS("KeyName") == "TTZToLLNuNu_M-10ext1" ||\
    R.getS("KeyName") == "TTZToLLNuNu_M-10ext2" ||\
    R.getS("KeyName") == "TTZToLLNuNu_M-10ext3")
    {
        stitch = 1.0;
        returnWeight_ = (1000.0 * R.getD("CrossSection") * R.getD("generatorEventWeight")) / ( 3511199.52 * R.getD("FilterEff"));
    }
    else if (R.getS("KeyName") == "TTWJetsToLNu")
    {
        stitch = 1.0;
        returnWeight_ = (1000.0 * R.getD("CrossSection") * R.getD("generatorEventWeight")) / ( 731877.67 * R.getD("FilterEff"));
    }
    else if (R.getS("KeyName") == "ttHJetToNonbb")
    {
        stitch = 1.0;
        returnWeight_ = (1000.0 * R.getD("CrossSection") * R.getD("generatorEventWeight")) / ( 10129526.16 * R.getD("FilterEff"));
    }
    else if (R.getS("KeyName") == "WZTo2L2Q")
    {
        stitch = 1.0;
        returnWeight_ = (1000.0 * R.getD("CrossSection") * R.getD("generatorEventWeight")) / ( 234019002.42 * R.getD("FilterEff"));
    }
    
    /* fix wrong XS and stitch*/
    else if (R.getS("KeyName") == "WWToLNuQQ" ||\
    R.getS("KeyName") == "WWToLNuQQext1")
    {
        stitch = 1.0;
        returnWeight_ = 0.00484 / (R.getD("FilterEff"));
    }
    
	/* for other samples return the nominal weight */
	if (stitch == 0) returnWeight_ = ( 1000.0 * R.getD("CrossSection") ) / ( R.getI("EventTotal") * R.getD("FilterEff"));
    
	return returnWeight_;

}


/* function: getTopQuarkPtWeight 
		-- returns  top_wt = sqrt(exp(0.156-0.00137*top_1_pt)*exp(0.156-0.00137*top_2_pt)) (with top_x_pt capped at 400)
		   for ttbar samples
		-- returns 1.0 for all other samples
*/

double generateH2TauSyncTree::getTopQuarkPtWeight(bool verbose_)
{

	double returnWeight_ = 1.0;

	if(R.getS("KeyName") == "TT_TuneCUETP8M2T4")
	{

		double top_1_pt = std::min(R.getD("genTopPt1"), 400.0);
		double top_2_pt = std::min(R.getD("genTopPt2"), 400.0);

		returnWeight_ = sqrt(exp(0.0615-0.0005*top_1_pt)*exp(0.0615-0.0005*top_2_pt));

		if(verbose_) std::cout<<" in top sample with pt1, pt2 = "<<top_1_pt<<" , "<<top_2_pt<<" get top pt weight = "<<returnWeight_<<"\n";

	}

	return returnWeight_;

}          

double generateH2TauSyncTree::getZReWeight(bool verbose_)
{
	double returnWeight_ = 1.0;

	bool calc = 0;

    if(R.getS("KeyName") == "DYJetsToLL_M-50" ||\
    R.getS("KeyName") == "DYJetsToLL_M-50ext1-v2" ||\
    R.getS("KeyName") == "DY1JetsToLL_M-50" ||\
    R.getS("KeyName") == "DY2JetsToLL_M-50" ||\
    R.getS("KeyName") == "DY3JetsToLL_M-50"  ||\
    R.getS("KeyName") == "DY4JetsToLL_M-50"  ||\
    R.getS("KeyName") == "EWKZ2Jets_ZToLL_M-50"  ||\
    R.getS("KeyName") == "EWKZ2Jets_ZToNuNu"
    )	calc = 1;

	if(calc == 0)
	{
		if(verbose_) std::cout<<" non dy sample, return z reweight as 1.0 \n";
		return 1.0;
	}

	double genpT 	= R.getD("genBosonTotal_pt");
	double genMass 	= R.getD("genBosonTotal_M");
	returnWeight_ = zweightHist->GetBinContent(zweightHist->GetXaxis()->FindBin(genMass),
											   zweightHist->GetYaxis()->FindBin(genpT));

	if(verbose_) 
	{
		std::cout<<"  for dy sample, return z reweight wt( "<<genMass<<" , "<<	genpT <<" ) = "<<returnWeight_<<"\n";

	}

	return returnWeight_;

}        

/* function: getHighPtTauUncertainty(bool)
		-- returns  vector with element [0] = 1 + 0.2 * (gen_tauPt)/1000.0 and element [1] = 1 - 0.2 * (gen_tauPt)/1000.0
		-- for events with 2 hadronically decaying taus, uses the formula eff(leg1) + eff(leg2) - eff(leg1)*eff(leg2)
*/

std::vector<double> generateH2TauSyncTree::getHighPtTauUncertainty(bool verbose_)
{

    // used to be 20% per fraction TeV per leg, now 5% up and 35% down
	std::vector<double> returnVec_; 
	returnVec_.clear();

	double unc1 = 0.;
	double unc2 = 0.;
	double eff1 = 1.0;
	double eff2 = 1.0;

	/* leg1 */
	if(R.getI("leg1_leptonType") == 3 && R.getI("leg1_MCMatchType")==5)
	{
		unc1 = R.getD("leg1_genMCmatch_pt")/1000.0;
		if(verbose_) std::cout<<" applying high pt tau eff unc. to leg1 \n";
	}

	/* leg2 */
	if(R.getI("leg2_leptonType") == 3 && R.getI("leg2_MCMatchType")==5)
	{
		unc2 = R.getD("leg2_genMCmatch_pt")/1000.0;
		if(verbose_) std::cout<<" applying high pt tau eff unc. to leg2 \n";
	}

	/* if MC compute the uncertainty */
	if(R.getB("isRealData") == 0)
	{
		/* get the up shift */
		eff1 = 1 + 0.05*(unc1 + unc2);

		/* get the down shift */
		eff2 = 1 - 0.35*(unc1 + unc2);
	}

	returnVec_.clear();
	returnVec_.push_back(eff1);
	returnVec_.push_back(eff2);

	return returnVec_;

}   


/* function: getNLOReWeight 
		-- returns  NLO Higgs Pt based reweight for SUSY signal samples			   
		-- returns 1.0 for all other samples
*/

double generateH2TauSyncTree::getNLOReWeight(bool verbose_, int tanBeta_)
{
	bool eval_ = 0;
	int mass_ = 0;

		if( R.getS("KeyName") == "Fall15_SUSYggHTauTau_M80") {eval_ = 1; mass_ = 80;}
	else if( R.getS("KeyName") == "Fall15_SUSYggHTauTau_M90") {eval_ = 1; mass_ = 90;}
	else if( R.getS("KeyName") == "Fall15_SUSYggHTauTau_M100") {eval_ = 1; mass_ = 100;}
	else if( R.getS("KeyName") == "Fall15_SUSYggHTauTau_M110") {eval_ = 1; mass_ = 110;}
	else if( R.getS("KeyName") == "Fall15_SUSYggHTauTau_M120") {eval_ = 1; mass_ = 120;}
	else if( R.getS("KeyName") == "Fall15_SUSYggHTauTau_M130") {eval_ = 1; mass_ = 130;}
	else if( R.getS("KeyName") == "Fall15_SUSYggHTauTau_M140") {eval_ = 1; mass_ = 140;}
	else if( R.getS("KeyName") == "Fall15_SUSYggHTauTau_M160") {eval_ = 1; mass_ = 160;}
	else if( R.getS("KeyName") == "Fall15_SUSYggHTauTau_M180") {eval_ = 1; mass_ = 180;}
	else if( R.getS("KeyName") == "Fall15_SUSYggHTauTau_M200") {eval_ = 1; mass_ = 200;}
	else if( R.getS("KeyName") == "Fall15_SUSYggHTauTau_M250") {eval_ = 1; mass_ = 250;}
	else if( R.getS("KeyName") == "Fall15_SUSYggHTauTau_M300") {eval_ = 1; mass_ = 300;}
	else if( R.getS("KeyName") == "Fall15_SUSYggHTauTau_M350") {eval_ = 1; mass_ = 350;}
	else if( R.getS("KeyName") == "Fall15_SUSYggHTauTau_M400") {eval_ = 1; mass_ = 400;}
	else if( R.getS("KeyName") == "Fall15_SUSYggHTauTau_M450") {eval_ = 1; mass_ = 450;}
	else if( R.getS("KeyName") == "Fall15_SUSYggHTauTau_M500") {eval_ = 1; mass_ = 500;}
	else if( R.getS("KeyName") == "Fall15_SUSYggHTauTau_M600") {eval_ = 1; mass_ = 600;}
	else if( R.getS("KeyName") == "Fall15_SUSYggHTauTau_M700") {eval_ = 1; mass_ = 700;}
	else if( R.getS("KeyName") == "Fall15_SUSYggHTauTau_M800") {eval_ = 1; mass_ = 800;}
	else if( R.getS("KeyName") == "Fall15_SUSYggHTauTau_M900") {eval_ = 1; mass_ = 900;}
	else if( R.getS("KeyName") == "Fall15_SUSYggHTauTau_M1000") {eval_ = 1; mass_ = 1000;}
	else if( R.getS("KeyName") == "Fall15_SUSYggHTauTau_M1200") {eval_ = 1; mass_ = 1200;}
	else if( R.getS("KeyName") == "Fall15_SUSYggHTauTau_M1400") {eval_ = 1; mass_ = 1400;}
	else if( R.getS("KeyName") == "Fall15_SUSYggHTauTau_M1500") {eval_ = 1; mass_ = 1500;}
	else if( R.getS("KeyName") == "Fall15_SUSYggHTauTau_M1600") {eval_ = 1; mass_ = 1600;}
	else if( R.getS("KeyName") == "Fall15_SUSYggHTauTau_M1800") {eval_ = 1; mass_ = 1800;}
	else if( R.getS("KeyName") == "Fall15_SUSYggHTauTau_M2000") {eval_ = 1; mass_ = 2000;}
	else if( R.getS("KeyName") == "Fall15_SUSYggHTauTau_M2300") {eval_ = 1; mass_ = 2300;}
	else if( R.getS("KeyName") == "Fall15_SUSYggHTauTau_M2600") {eval_ = 1; mass_ = 2600;}
	else if( R.getS("KeyName") == "Fall15_SUSYggHTauTau_M2900") {eval_ = 1; mass_ = 2900;}
	else if( R.getS("KeyName") == "Fall15_SUSYggHTauTau_M3200") {eval_ = 1; mass_ = 3200;}


	/* return 1.0 for all non SUSY GGH samples */
	if(!eval_) return 1.0;


	/* else evalualte the NLO reweight */

	return NLO_returnNLOweight(mass_, tanBeta_, R.getD("genBosonTotal_pt"));

}




/* adapted from https://github.com/CMS-HTT/NLOReweightingTool */
void  generateH2TauSyncTree::NLO_ReadFile()
{

  TFile *file = new TFile("Reweight.root");

  const int num_of_tb = 60;  

  int imass = 0;
  for(auto mass: NLO_marray){
    for(int tanb=0; tanb < num_of_tb; tanb++){

      TString wname = "weight_MSSM_";
      wname += mass;
      wname += "_tanb_";
      wname += tanb + 1;

      //      std::cout << wname << " " << NLO_func[imass][tanb] << " " << NLO_func[imass][tanb]->GetBinContent(1) << std::endl;
      // NLO_func[imass][tanb] = (TGraphErrors*) gROOT->FindObject(wname); /* original code, replaced by Get */

	  NLO_func[imass][tanb] = (TGraphErrors*) file->Get(wname);

    }    
    imass++;
  }

}


/* adapted from https://github.com/CMS-HTT/NLOReweightingTool */
float generateH2TauSyncTree::NLO_returnNLOweight(Int_t mass, Int_t tanb, Double_t pt)
{


  if(pt > 800){
    //    std::cout << "[INFO] pT = " << pt << " exceeds the range --> set it to 800." << std::endl;    
    pt = 800;
  }

  auto iter = std::find(NLO_marray.begin(), NLO_marray.end(), mass);
  size_t index = std::distance(NLO_marray.begin(), iter);

  if(index == NLO_marray.size()){
    std::cout << "[WARNING] Invalid mass point ... " << mass << " -> return weight 1" << std::endl;    
    return 1;
  }
  
  if(tanb <1 || tanb > 60){
    std::cout << "[WARNING] Invalid tan(beta) point ... " << tanb << " -> return weight 1" << std::endl;
    return 1;
  }

  return NLO_func[index][tanb-1]->Eval(pt) ;

}

/* 
    function: getQCDWeightForEleMuChannel(bool)
	-- returns a size 6 double vector with elements :

	double qcdweight @ element 0
	double qcdweightup @ element 1 
	double qcdweightdown @ element 2 
	double qcdweightNoDZeta @ element 3 
	double qcdweightupNoDZeta @ element 4
	double qcdweightdownNoDZeta @ element 5 

	all set to 1.0 in case not e+mu channel
*/

std::vector<double> generateH2TauSyncTree::getQCDWeightForEleMuChannel(bool verbose, TLorentzVector l1, TLorentzVector l2)
{


	std::vector <double> returnVec;
	returnVec.clear();

	// accessing OS/SS extrapolation factor as a function
	// of pt(e), pt(mu), and dR(e,mu)

	if(R.getI("CandidateEventType") != 2) 
	{
		if(verbose) std::cout<<" qcd weights for ele mu channel called for non ele-mu channel returning <1.0> \n";
		for(int k=0; k<6; ++k) returnVec.push_back(1.0);
		return returnVec;
	}

	if(verbose) std::cout<<" qcd weights for ele mu channel called for ele-mu channel returning weights vector \n";
		
	double pt_e = l1.Pt();
	double pt_m = l2.Pt();
	double dR = R.getD("DeltaR_leg1_leg2");

	/* should always have leg1 = e and leg2 = mu, but just in case order changes in the future */
	if(R.getI("leg1_leptonType")!=1 && R.getI("leg2_leptonType") == 1)  pt_e = l2.Pt();
	if(R.getI("leg2_leptonType")!=2 && R.getI("leg1_leptonType") == 2)  pt_m = l1.Pt();


	//double qcdweight = qcdWeights->getWeight(pt_e,pt_m,dR);
    double qcdweight = 1.0;

	// accessing OS/SS extrapolation factor corresponding
	// to +1sigma systematic variation of QCD background shape

	//double qcdweightup = qcdWeights->getWeightUp(pt_e,pt_m,dR);
    double qcdweightup = 1.0;

	// We suggest to compute OS/SS extrapolation factor
	// corresponding to the -1sigma systematic variation
	// of the QCD background shape as

	double qcdweightdown = 1; 
	if(qcdweightup != 0) qcdweightdown = qcdweight * qcdweight / qcdweightup;

	//double qcdweightNoDZeta = qcdWeightsNoDZeta->getWeight(pt_e,pt_m,dR);
    double qcdweightNoDZeta = 1.0;

	// accessing OS/SS extrapolation factor corresponding
	// to +1sigma systematic variation of QCD background shape

	//double qcdweightupNoDZeta = qcdWeightsNoDZeta->getWeightUp(pt_e,pt_m,dR);
    double qcdweightupNoDZeta = 1.0;

	// We suggest to compute OS/SS extrapolation factor
	// corresponding to the -1sigma systematic variation
	// of the QCD background shape as

	double qcdweightdownNoDZeta = 1; 
	if(qcdweightupNoDZeta != 0) qcdweightdownNoDZeta = qcdweightNoDZeta * qcdweightNoDZeta / qcdweightupNoDZeta;

	returnVec.push_back(qcdweight);
	returnVec.push_back(qcdweightup);
	returnVec.push_back(qcdweightdown);
	returnVec.push_back(qcdweightNoDZeta);
	returnVec.push_back(qcdweightupNoDZeta);
	returnVec.push_back(qcdweightdownNoDZeta);


	return returnVec;
}

double generateH2TauSyncTree::getCentralMuonFactor(Double_t eta, Double_t pt, bool trig)
{
    float returnWeight_ = 1.0;
    float periodBCDEFweight = 0.557065;
    float periodGHweight = 0.442935;
    
    if(trig)
    {
        float TRIGsfBCDEF = sfHisto_Muon_Trigger_BCDEF->GetBinContent(sfHisto_Muon_Trigger_BCDEF->GetXaxis()->FindBin(eta),sfHisto_Muon_Trigger_BCDEF->GetYaxis()->FindBin(std::min(pt,499.)));
        
        float TRIGsfGH = sfHisto_Muon_Trigger_GH->GetBinContent(sfHisto_Muon_Trigger_GH->GetXaxis()->FindBin(eta),sfHisto_Muon_Trigger_GH->GetYaxis()->FindBin(std::min(pt,499.)));
        
        returnWeight_ *= (periodBCDEFweight*TRIGsfBCDEF + periodGHweight*TRIGsfGH);
    }
    else
    {
        float IDsfBCDEF = sfHisto_Muon_TightID2016_BCDEF->GetBinContent(sfHisto_Muon_TightID2016_BCDEF->GetXaxis()->FindBin(eta),sfHisto_Muon_TightID2016_BCDEF->GetYaxis()->FindBin(std::min(pt,119.)));
        float ISOsfBCDEF = sfHisto_Muon_TightIso_BCDEF->GetBinContent(sfHisto_Muon_TightIso_BCDEF->GetXaxis()->FindBin(eta),sfHisto_Muon_TightIso_BCDEF->GetYaxis()->FindBin(std::min(pt,119.)));
        
        float IDsfGH = sfHisto_Muon_TightID2016_GH->GetBinContent(sfHisto_Muon_TightID2016_GH->GetXaxis()->FindBin(eta),sfHisto_Muon_TightID2016_GH->GetYaxis()->FindBin(std::min(pt,119.)));
        float ISOsfGH = sfHisto_Muon_TightIso_GH->GetBinContent(sfHisto_Muon_TightIso_GH->GetXaxis()->FindBin(eta),sfHisto_Muon_TightIso_GH->GetYaxis()->FindBin(std::min(pt,119.)));
        
        returnWeight_ *= ((periodBCDEFweight*IDsfBCDEF + periodGHweight*IDsfGH) * (periodBCDEFweight*ISOsfBCDEF + periodGHweight*ISOsfGH));
    }
    return returnWeight_;
}

double generateH2TauSyncTree::getKFactor(bool verbose)
{
    double returnWeight_ = 1.0;
    
    // Add W jets key names only

    if(R.getS("KeyName") == "WJetsToLNu" ||\
	   R.getS("KeyName") == "WJetsToLNuext2-v1" ||\
	   R.getS("KeyName") == "WJetsToLNu_HT-70To100" ||\
	   R.getS("KeyName") == "WJetsToLNu_HT-100To200orig" ||\
	   R.getS("KeyName") == "WJetsToLNu_HT-100To200" ||\
       R.getS("KeyName") == "WJetsToLNu_HT-100To200ext2-v1" ||\
       R.getS("KeyName") == "WJetsToLNu_HT-200To400orig" ||\
       R.getS("KeyName") == "WJetsToLNu_HT-200To400ext1-v1" ||\
       R.getS("KeyName") == "WJetsToLNu_HT-200To400" ||\
       R.getS("KeyName") == "WJetsToLNu_HT-400To600orig" ||\
       R.getS("KeyName") == "WJetsToLNu_HT-400To600" ||\
       R.getS("KeyName") == "WJetsToLNu_HT-600To800orig" ||\
       R.getS("KeyName") == "WJetsToLNu_HT-600To800" ||\
       R.getS("KeyName") == "WJetsToLNu_HT-800To1200orig" ||\
       R.getS("KeyName") == "WJetsToLNu_HT-800To1200" ||\
       R.getS("KeyName") == "WJetsToLNu_HT-1200To2500orig" ||\
       R.getS("KeyName") == "WJetsToLNu_HT-1200To2500" ||\
       R.getS("KeyName") == "WJetsToLNu_HT-2500ToInforig" ||\
       R.getS("KeyName") == "WJetsToLNu_HT-2500ToInf")
       
    {
        double genBosonTotal_Wpt = R.getD("MaxPtGenBoson_WisconinStyle_pt");
        if (genBosonTotal_Wpt < 150.) genBosonTotal_Wpt = 151.;
        double k = EWK_Wcorr->GetBinContent(EWK_Wcorr->GetXaxis()->FindBin(genBosonTotal_Wpt));
        std::cout << " k factor: " << k << std::endl;
        returnWeight_ *= k;
        if (verbose)std::cout<<" W k-factor: " << returnWeight_ << " \n";
    }
    
    /*
    else if(R.getS("KeyName") == "DYJetsToLL_M-50" ||\
	   R.getS("KeyName") == "DYJetsToLL_M-50ext1-v2" ||\
	   R.getS("KeyName") == "DY1JetsToLL_M-50" ||\
	   R.getS("KeyName") == "DY2JetsToLL_M-50" ||\
	   R.getS("KeyName") == "DY3JetsToLL_M-50"  ||\
	   R.getS("KeyName") == "DY4JetsToLL_M-50")
	{
        // turn off k factor for DY now, using Z reweight
        returnWeight_ = 1.;
        //if(genBosonMass_ > 5.0) {returnWeight_ *= EWK_Zcorr->GetBinContent(EWK_Zcorr->FindBin(genBosonTotal_pt));}
        //else {returnWeight_ *= EWK_Gcorr->GetBinContent(EWK_Gcorr->FindBin(genBosonTotal_pt));}
        if (verbose) {std::cout<<" Z/G k-factor: " << returnWeight_ << " \n";}
    }
    */
    
    return returnWeight_;
    
}

double generateH2TauSyncTree::getKFactorSyst(bool verbose, bool down)
{
    double returnWeight_ = 1.0;
    // Add W jets key names only

    if(R.getS("KeyName") == "WJetsToLNu" ||\
	   R.getS("KeyName") == "WJetsToLNuext2-v1" ||\
	   R.getS("KeyName") == "WJetsToLNu_HT-70To100" ||\
	   R.getS("KeyName") == "WJetsToLNu_HT-100To200orig" ||\
	   R.getS("KeyName") == "WJetsToLNu_HT-100To200" ||\
       R.getS("KeyName") == "WJetsToLNu_HT-100To200ext2-v1" ||\
       R.getS("KeyName") == "WJetsToLNu_HT-200To400orig" ||\
       R.getS("KeyName") == "WJetsToLNu_HT-200To400ext1-v1" ||\
       R.getS("KeyName") == "WJetsToLNu_HT-200To400" ||\
       R.getS("KeyName") == "WJetsToLNu_HT-400To600orig" ||\
       R.getS("KeyName") == "WJetsToLNu_HT-400To600" ||\
       R.getS("KeyName") == "WJetsToLNu_HT-600To800orig" ||\
       R.getS("KeyName") == "WJetsToLNu_HT-600To800" ||\
       R.getS("KeyName") == "WJetsToLNu_HT-800To1200orig" ||\
       R.getS("KeyName") == "WJetsToLNu_HT-800To1200" ||\
       R.getS("KeyName") == "WJetsToLNu_HT-1200To2500orig" ||\
       R.getS("KeyName") == "WJetsToLNu_HT-1200To2500" ||\
       R.getS("KeyName") == "WJetsToLNu_HT-2500ToInforig" ||\
       R.getS("KeyName") == "WJetsToLNu_HT-2500ToInf")
       
    {
        double genBosonTotal_Wpt = R.getD("MaxPtGenBoson_WisconinStyle_pt");
        if (genBosonTotal_Wpt < 150) genBosonTotal_Wpt = 151;
        double k = EWK_Wcorr_dNLO->GetBinContent(EWK_Wcorr_dNLO->GetXaxis()->FindBin(genBosonTotal_Wpt));
        returnWeight_ *= k;
        if (down) {returnWeight_ = 2.0-k;}
        if (verbose)std::cout<<" W k-factor Systematic (EWK factor only) : " << returnWeight_ << " \n";
    }
    
    /*
    else if(R.getS("KeyName") == "DYJetsToLL_M-50" ||\
	   R.getS("KeyName") == "DYJetsToLL_M-50ext1-v2" ||\
	   R.getS("KeyName") == "DY1JetsToLL_M-50" ||\
	   R.getS("KeyName") == "DY2JetsToLL_M-50" ||\
	   R.getS("KeyName") == "DY3JetsToLL_M-50"  ||\
	   R.getS("KeyName") == "DY4JetsToLL_M-50")
	{
        // turn off k factor for DY now, using Z reweight
        returnWeight_ = 1.;
        //if(genBosonMass_ > 5.0) {returnWeight_ *= EWK_Zcorr->GetBinContent(EWK_Zcorr->FindBin(genBosonTotal_pt));}
        //else {returnWeight_ *= EWK_Gcorr->GetBinContent(EWK_Gcorr->FindBin(genBosonTotal_pt));}
        if (verbose) {std::cout<<" Z/G k-factor: " << returnWeight_ << " \n";}
    }
    */
    return returnWeight_;
}

double generateH2TauSyncTree::getWWFactor(bool verbose_)
{
    double returnWeight_ = 1.0;
    
    if(R.getS("KeyName") == "WWTo2L2Nu" || R.getS("KeyName") == "WWJJToLNuLNu_EWK_noTop" || R.getS("KeyName") == "WWTo2L2Nu_DoubleScattering" || R.getS("KeyName") == "WpWpJJ_EWK" || R.getS("KeyName") == "WpWpJJ_QCD" || R.getS("KeyName") == "GluGluWWTo2L2Nu")
    {
        double gen_met = R.getD("genMET");
        if (gen_met < 20.) {gen_met = 21.;}
        double xval = WWcorr->GetXaxis()->FindBin(gen_met);
        double k = WWcorr->GetBinContent(xval);
        if (verbose_) {std::cout << "WW up weight factor: " << k << std::endl;}
        returnWeight_ *= k;
    }
    return returnWeight_;
}

double generateH2TauSyncTree::getZZFactor(int var)
{
    double returnWeight_ = 1.0;
    
    if(R.getS("KeyName") == "ZZTo2L2Nu")
    {
        double genBosonTotal_Zpt = R.getD("MaxPtGenBoson_WisconinStyle_pt");
        float xval = ZZcorr->GetXaxis()->FindBin(genBosonTotal_Zpt);
    
        if (var==-1)
        {
            double j = ZZcorr->GetBinContent(xval,4);
            double k = ZZcorr->GetBinContent(xval,2);
            returnWeight_ = j/k;
        }
        else if (var==0)
        {
            double j = ZZcorr->GetBinContent(xval,2);
            double k = ZZcorr->GetBinContent(xval,1);
            returnWeight_ = j/k;
        }
        else if (var==1)
        {
            double j = ZZcorr->GetBinContent(xval,3);
            double k = ZZcorr->GetBinContent(xval,2);
            returnWeight_ = j/k;
        }
    }
    
    return returnWeight_;
}

double generateH2TauSyncTree::getJetTauFakeFactor(bool verbose, int variant, TLorentzVector l1, TLorentzVector l2)
{
    // remember to normalize tau-tau histogram
    double returnWeight_ = 1.0;
    double pt1 = l1.Pt();
    double pt2 = l2.Pt();
    if (pt1 > 200.0){pt1 = 200.0;}
    if (pt2 > 200.0){pt2 = 200.0;}

    if (R.getI("CandidateEventType")==5 || R.getI("CandidateEventType")==3)
    {
        if (variant == 1 && R.getI("leg2_MCMatchType") == 6)
        {
            returnWeight_ *= 1.09103-0.0039767*pt2+0.0000083736*pt2*pt2;
        }
        else if (variant == -1 && R.getI("leg2_MCMatchType") == 6)
        {
            returnWeight_ *= 1/(1.09103-0.0039767*pt2+0.0000083736*pt2*pt2);
        }
    }
    else if (R.getI("CandidateEventType")==6)
    {
        if (variant == 1)
        {
            if(R.getI("leg1_MCMatchType") == 6 && R.getI("leg2_MCMatchType") == 6)
            {
                returnWeight_ *= (1.09103-0.0039767*pt1+0.0000083736*pt1*pt1)*(1.09103-0.0039767*pt2+0.0000083736*pt2*pt2);
            }
            else if(R.getI("leg1_MCMatchType") == 6 && R.getI("leg2_MCMatchType") != 6)
            {
                returnWeight_ *= 1.09103-0.0039767*pt1+0.0000083736*pt1*pt1;
            }
            else if(R.getI("leg1_MCMatchType") != 6 && R.getI("leg2_MCMatchType") == 6)
            {
                returnWeight_ *= 1.09103-0.0039767*pt2+0.0000083736*pt2*pt2;
            }
        }
        else if (variant == -1)
        {
            if(R.getI("leg1_MCMatchType") == 6 && R.getI("leg2_MCMatchType") == 6)
            {
                returnWeight_ *= 1/((1.09103-0.0039767*pt1+0.0000083736*pt1*pt1)*(1.09103-0.0039767*pt2+0.0000083736*pt2*pt2));
            }
            else if(R.getI("leg1_MCMatchType") == 6 && R.getI("leg2_MCMatchType") != 6)
            {
                returnWeight_ *= 1/(1.09103-0.0039767*pt1+0.0000083736*pt1*pt1);
            }
            else if(R.getI("leg1_MCMatchType") != 6 && R.getI("leg2_MCMatchType") == 6)
            {
                returnWeight_ *= 1/(1.09103-0.0039767*pt2+0.0000083736*pt2*pt2);
            }
        }
    }
    return returnWeight_;
}

double generateH2TauSyncTree::getALDScaleFactors(bool verbose)
{

    double returnSF_ = 1.0;
    
    if (R.getI("CandidateEventType")==5)
    {
        if (R.getI("leg2_MCMatchType") == 1 || R.getI("leg2_MCMatchType") == 3)
        {
            if (fabs(R.getD("leg2_eta")) < 1.460) returnSF_ *= 1.213;
            else if (fabs(R.getD("leg2_eta")) > 1.558) returnSF_ *= 1.375;
        }
        else if (R.getI("leg2_MCMatchType") == 2 || R.getI("leg2_MCMatchType") == 4)
        {
            if (fabs(R.getD("leg2_eta")) < 0.4) returnSF_ *= 1.263;
            else if (fabs(R.getD("leg2_eta")) > 0.4 && fabs(R.getD("leg2_eta")) < 0.8) returnSF_ *= 1.364;
            else if (fabs(R.getD("leg2_eta")) > 0.8 && fabs(R.getD("leg2_eta")) < 1.2) returnSF_ *= 0.854;
            else if (fabs(R.getD("leg2_eta")) > 1.2 && fabs(R.getD("leg2_eta")) < 1.7) returnSF_ *= 1.712;
            else if (fabs(R.getD("leg2_eta")) > 1.7 && fabs(R.getD("leg2_eta")) < 2.3) returnSF_ *= 2.324;
        }
    }
    else if (R.getI("CandidateEventType")==3)
    {
        if (R.getI("leg2_MCMatchType") == 1 || R.getI("leg2_MCMatchType") == 3)
        {
            if (fabs(R.getD("leg2_eta")) < 1.460) returnSF_ *= 1.402;
            else if (fabs(R.getD("leg2_eta")) > 1.558) returnSF_ *= 1.900;
        }
        else if (R.getI("leg2_MCMatchType") == 2 || R.getI("leg2_MCMatchType") == 4)
        {
            if (fabs(R.getD("leg2_eta")) < 0.4) returnSF_ *= 1.010;
            else if (fabs(R.getD("leg2_eta")) > 0.4 && fabs(R.getD("leg2_eta")) < 0.8) returnSF_ *= 1.007;
            else if (fabs(R.getD("leg2_eta")) > 0.8 && fabs(R.getD("leg2_eta")) < 1.2) returnSF_ *= 0.870;
            else if (fabs(R.getD("leg2_eta")) > 1.2 && fabs(R.getD("leg2_eta")) < 1.7) returnSF_ *= 1.154;
            else if (fabs(R.getD("leg2_eta")) > 1.7 && fabs(R.getD("leg2_eta")) < 2.3) returnSF_ *= 2.281;
        }
    }
    else if (R.getI("CandidateEventType")==6)
    {
        if (R.getI("leg1_MCMatchType") == 1 || R.getI("leg1_MCMatchType") == 3)
        {
            if (fabs(R.getD("leg1_eta")) < 1.460) returnSF_ *= 1.213;
            else if (fabs(R.getD("leg1_eta")) > 1.558) returnSF_ *= 1.375;
        }
        else if (R.getI("leg1_MCMatchType") == 2 || R.getI("leg1_MCMatchType") == 4)
        {
            if (fabs(R.getD("leg1_eta")) < 0.4) returnSF_ *= 1.010;
            else if (fabs(R.getD("leg1_eta")) > 0.4 && fabs(R.getD("leg1_eta")) < 0.8) returnSF_ *= 1.007;
            else if (fabs(R.getD("leg1_eta")) > 0.8 && fabs(R.getD("leg1_eta")) < 1.2) returnSF_ *= 0.870;
            else if (fabs(R.getD("leg1_eta")) > 1.2 && fabs(R.getD("leg1_eta")) < 1.7) returnSF_ *= 1.154;
            else if (fabs(R.getD("leg1_eta")) > 1.7 && fabs(R.getD("leg1_eta")) < 2.3) returnSF_ *= 2.281;
        }
        
        if (R.getI("leg2_MCMatchType") == 1 || R.getI("leg2_MCMatchType") == 3)
        {
            if (fabs(R.getD("leg2_eta")) < 1.460) returnSF_ *= 1.213;
            else if (fabs(R.getD("leg2_eta")) > 1.558) returnSF_ *= 1.375;
        }
        else if (R.getI("leg2_MCMatchType") == 2 || R.getI("leg2_MCMatchType") == 4)
        {
            if (fabs(R.getD("leg2_eta")) < 0.4) returnSF_ *= 1.010;
            else if (fabs(R.getD("leg2_eta")) > 0.4 && fabs(R.getD("leg2_eta")) < 0.8) returnSF_ *= 1.007;
            else if (fabs(R.getD("leg2_eta")) > 0.8 && fabs(R.getD("leg2_eta")) < 1.2) returnSF_ *= 0.870;
            else if (fabs(R.getD("leg2_eta")) > 1.2 && fabs(R.getD("leg2_eta")) < 1.7) returnSF_ *= 1.154;
            else if (fabs(R.getD("leg2_eta")) > 1.7 && fabs(R.getD("leg2_eta")) < 2.3) returnSF_ *= 2.281;
        }
    }
    
    return returnSF_;

}

std::vector <double> generateH2TauSyncTree::getTauShift(int dm, int lt, int shift)
{
    std::vector<double> returnVector;
	returnVector.clear();

    double returnVal0_ = 1.0;
    double returnVal1_ = 1.0;
    
    double returnShift0_ = 1.0;
    double returnShift1_ = 1.0;
    
    /*
    double ESdm0 = .995;
    double ESdm1 = 1.011
    double ESdm10 = 1.006
    */
    
    double ESdm0 = 1.0;
    double ESdm1 = 1.0;
    double ESdm10 = 1.0;
    
    if (lt == 3)
    {
        if (shift == 1) {returnShift0_ = 1.012; returnShift1_ = 1.012;}
        else if (shift == -1) {returnShift0_ = 0.988; returnShift1_ = 0.988;}
        else {returnShift0_ = 1.0; returnShift1_ = 1.0;}
        
        if (dm==0) {returnShift1_ = 1.0; returnVal0_ = ESdm0 * returnShift0_; returnVal1_ = 1.0;}
        else if (dm==1) {returnVal0_ = ESdm1 * returnShift0_; returnVal1_ = ESdm1 * returnShift1_;}
        else if (dm==10) {returnVal0_ = ESdm10 * returnShift0_; returnVal1_ = ESdm10 * returnShift1_;}
    }
    
    if (R.getB("isRealData")==0)
    {
        returnVector.clear();
        returnVector.push_back(returnVal0_);
        returnVector.push_back(returnVal1_);
        returnVector.push_back(returnShift0_);
        returnVector.push_back(returnShift1_);
    }
    else
    {
        returnVector.clear();
        returnVector.push_back(1.0);
        returnVector.push_back(1.0);
        returnVector.push_back(1.0);
        returnVector.push_back(1.0);
    }
    return returnVector;
    
}

double generateH2TauSyncTree::getFinalScaleFactorsForPair(bool verbose, int sysShift, bool useMuonCentral, bool returnIDISO, bool returnTRIG, bool returnTRACK, TLorentzVector l1, TLorentzVector l2)
{

    double tauID1 = 1.0;  //current tau ID scale factor (just DMF)
    double tauID2 = 1.0;
    
	double returnSF = 1.0;
    
	/* muon + tau */
	if(R.getI("CandidateEventType") == 5)
	{
        if (R.getI("leg1_MCMatchType")==5) tauID1=0.95;
        if (R.getI("leg2_MCMatchType")==5) tauID2=0.95;
    
		if(verbose) std::cout<<" id x iso x trigger sf for Muon + Tau \n";

		double muonID = 1.0;
		double muonTrigger = 1.0;
        double muonTrack = 1.0;

		/* muon */
		double pt1 = 0.0;
		double eta1 = 0.0;
        
		/* should always have leg1 = mu and leg2 = tau, but just in case order changes in the future */
		if(R.getI("leg1_leptonType")==2 && R.getI("leg2_leptonType") == 3)  
		{
			pt1 = l1.Pt();
			eta1 = l1.Eta();
		}

		if(R.getI("leg2_leptonType")==2 && R.getI("leg1_leptonType") == 3)   
		{
			pt1 = l2.Pt();
			eta1 = l2.Eta();
		}
    
        tw->var("m_eta")->setVal(eta1);
        muonTrack = tw->function("m_trk_ratio")->getVal();
        
        if(useMuonCentral)
        {
            muonID = getCentralMuonFactor(fabs(eta1), pt1, 0);
            muonTrigger = getCentralMuonFactor(fabs(eta1), pt1, 1);
        }
        else
        {
            muonID = sfTool_Muon_IdIso0p15_eff->get_ScaleFactor(pt1,eta1);
            muonTrigger = sfTool_Muon_SingleMu_eff->get_ScaleFactor(pt1,eta1);
        }
        returnSF = muonID * muonTrigger * muonTrack * tauID2;
        if(returnIDISO) return muonID;
        else if(returnTRIG) return muonTrigger;
        else if(returnTRACK) return muonTrack;
        else return returnSF;
        
	/////////////////////////////////////
	}

	/* electron + tau */
	if(R.getI("CandidateEventType") == 3)
	{
    
        if (R.getI("leg1_MCMatchType")==5) tauID1=0.95;
        if (R.getI("leg2_MCMatchType")==5) tauID2=0.95;

		if(verbose) std::cout<<" id x iso x trigger sf for Electron + Tau \n";

		double electronID = 1.0;
		double electronTrigger = 1.0;
        double electronTrack = 1.0;

		/* electron */
		double pt1 = 0.0;
		double eta1 = 0.0;

		/* should always have leg1 = electron and leg2 = tau, but just in case order changes in the future */
		if(R.getI("leg1_leptonType")==1 && R.getI("leg2_leptonType") == 3)  
		{
			pt1 = l1.Pt();
			eta1 = l1.Eta();
		}

		if(R.getI("leg2_leptonType")==1 && R.getI("leg1_leptonType") == 3)   
		{
			pt1 = l2.Pt();
			eta1 = l2.Eta();
		}

        tw->var("e_eta")->setVal(eta1);
        electronTrack = tw->function("e_trk_ratio")->getVal();

		electronID = sfTool_Electron_IdIso0p10_eff->get_ScaleFactor(pt1,eta1);
		electronTrigger = sfTool_Electron_SingleEle_eff->get_ScaleFactor(pt1,eta1);

        returnSF = electronID * electronTrigger * electronTrack * tauID2;
        
        if(returnIDISO) return electronID;
        else if(returnTRIG) return electronTrigger;
        else if(returnTRACK) return electronTrack;
        else return returnSF;

	}

	/* electron + muon */
	if(R.getI("CandidateEventType") == 2)
	{

		if(verbose) std::cout<<" id x iso x trigger sf for Electron + Muon \n";

		double effData = 1.0;
		double effMC = 1.0;
		double electronID = 1.0;
		double muonID = 1.0;
		double SF = 1.0;

		/* electron */
		double pt1 = 0.0;
		double eta1 = 0.0;

		/* muon */
		double pt2 = 0.0;
		double eta2 = 0.0;

		/* should always have leg1 = electron and leg2 = muon, but just in case order changes in the future */
		if(R.getI("leg1_leptonType")==1 && R.getI("leg2_leptonType") == 2)  
		{
			pt1 = l1.Pt();
			eta1 = l1.Eta();

			pt2 = l2.Pt();
			eta2 = l2.Eta();
		}

		if(R.getI("leg2_leptonType")==1 && R.getI("leg1_leptonType") == 2)   
		{
			pt1 = l2.Pt();
			eta1 = l2.Eta();

			pt2 = l1.Pt();
			eta2 = l1.Eta();
		}

		electronID = sfTool_Electron_IdIso0p15_eff->get_ScaleFactor(pt1,eta1);
		muonID = sfTool_Muon_IdIso0p20_eff->get_ScaleFactor(pt2,eta2);

 		/* effData = eff_data(Mu17)*eff_data(Ele12)
 		+eff_data(Mu8)*eff_data(Ele17)-eff_data(Mu17)*eff_data(Ele17) */

		effData =  sfTool_Muon_Mu17_eff->get_EfficiencyData(pt2, eta2)*sfTool_Electron_Ele12_eff->get_EfficiencyData(pt1, eta1);
		effData += sfTool_Muon_Mu8_eff->get_EfficiencyData(pt2, eta2)*sfTool_Electron_Ele17_eff->get_EfficiencyData(pt1, eta1);
		effData -= sfTool_Muon_Mu17_eff->get_EfficiencyData(pt2, eta2) * sfTool_Electron_Ele17_eff->get_EfficiencyData(pt1, eta1);  		

 		/* effMC = eff_MC(Mu17)*eff_MC(Ele12)
 		+eff_MC(Mu8)*eff_MC(Ele17)-eff_MC(Mu17)*eff_MC(Ele17)  */

		effMC =  sfTool_Muon_Mu17_eff->get_EfficiencyMC(pt2, eta2)*sfTool_Electron_Ele12_eff->get_EfficiencyMC(pt1, eta1);
		effMC += sfTool_Muon_Mu8_eff->get_EfficiencyMC(pt2, eta2)*sfTool_Electron_Ele17_eff->get_EfficiencyMC(pt1, eta1);
		effMC -= sfTool_Muon_Mu17_eff->get_EfficiencyMC(pt2, eta2) * sfTool_Electron_Ele17_eff->get_EfficiencyMC(pt1, eta1);  		

		if(effMC!=0.) SF = effData/effMC;

		returnSF = electronID * muonID * SF;

        /////////////////////////////////////
		return returnSF;
	}

	/* tau + tau */
	if(R.getI("CandidateEventType") == 6)
	{
    
        if (R.getI("leg1_MCMatchType")==5) tauID1=0.94;
        if (R.getI("leg2_MCMatchType")==5) tauID2=0.94;

		if(verbose) std::cout<<"  trigger sf for Tau + Tau with systematic shifted "<<sysShift<<" sigmas \n";
        
		double effData1 = 1.0;
		double effMC1 = 1.0;
		double SF1 = 1.0;

		double effData2 = 1.0;
		double effMC2 = 1.0;
		double SF2 = 1.0;

		/* tau1 */
		double pt1 = l1.Pt();
        //double mt1 = TMath::Sqrt((R.getD("leg1_M"))*(R.getD("leg1_M"))+pt1*pt1);
        double dm1 = R.getI("leg1_decayMode");

		/* tau2 */
		double pt2 = l2.Pt();
        //double mt2 = TMath::Sqrt((R.getD("leg2_M"))*(R.getD("leg2_M"))+pt2*pt2);
        double dm2 = R.getI("leg2_decayMode");
        
		/* make sure sysShift is valid */
		assert (sysShift==0 || sysShift==1 || sysShift == -1);

		if(sysShift == 0)
		{
        
            w->var("t_pt")->setVal(pt1);
            w->var("t_dm")->setVal(dm1);
            SF1 = w->function("t_genuine_TightIso_tt_ratio")->getVal();
            
            //if using fake and real
            /*if(mt1<30.0)
            {
                SF1 = w->function("t_genuine_TightIso_tt_ratio")->getVal();
            }
            else
            {
                SF1 = w->function("t_fake_TightIso_tt_ratio")->getVal();
            }*/
            
            w->var("t_pt")->setVal(pt2);
            w->var("t_dm")->setVal(dm2);
            SF2 = w->function("t_genuine_TightIso_tt_ratio")->getVal();
            
            //if using fake and real
            /*if(mt2<30.0)
            {
                SF2 = w->function("t_genuine_TightIso_tt_ratio")->getVal();
            }
            else
            {
                SF2 = w->function("t_fake_TightIso_tt_ratio")->getVal();
            }*/
            
            /*
            if(mt1<30.0)
            {
                if (dm1==0) effData1 = CBeff( pt1, Run2_TauTau_legTriggerEff_DataReal_dm0);
                if (dm1==1) effData1 = CBeff( pt1, Run2_TauTau_legTriggerEff_DataReal_dm1);
                if (dm1==10) effData1 = CBeff( pt1, Run2_TauTau_legTriggerEff_DataReal_dm10);
            }
            else
            {
                if (dm1==0) effData1 = CBeff( pt1, Run2_TauTau_legTriggerEff_DataFake_dm0);
                if (dm1==1) effData1 = CBeff( pt1, Run2_TauTau_legTriggerEff_DataFake_dm1);
                if (dm1==10) effData1 = CBeff( pt1, Run2_TauTau_legTriggerEff_DataFake_dm10);
            }
            
            if(mt2<30.0)
            {
                if (dm2==0) effData2 = CBeff( pt2, Run2_TauTau_legTriggerEff_DataReal_dm0);
                if (dm2==1) effData2 = CBeff( pt2, Run2_TauTau_legTriggerEff_DataReal_dm1);
                if (dm2==10) effData2 = CBeff( pt2, Run2_TauTau_legTriggerEff_DataReal_dm10);
            }
            else
            {
                if (dm2==0) effData2 = CBeff( pt2, Run2_TauTau_legTriggerEff_DataFake_dm0);
                if (dm2==1) effData2 = CBeff( pt2, Run2_TauTau_legTriggerEff_DataFake_dm1);
                if (dm2==10) effData2 = CBeff( pt2, Run2_TauTau_legTriggerEff_DataFake_dm10);
            }
            
            
            if(mt1<30.0)
            {
                if (dm1==0) effMC1 = CBeff( pt1, Run2_TauTau_legTriggerEff_MCReal_dm0);
                if (dm1==1) effMC1 = CBeff( pt1, Run2_TauTau_legTriggerEff_MCReal_dm1);
                if (dm1==10) effMC1 = CBeff( pt1, Run2_TauTau_legTriggerEff_MCReal_dm10);
            }
            else
            {
                if (dm1==0) effMC1 = CBeff( pt1, Run2_TauTau_legTriggerEff_MCFake_dm0);
                if (dm1==1) effMC1 = CBeff( pt1, Run2_TauTau_legTriggerEff_MCFake_dm1);
                if (dm1==10) effMC1 = CBeff( pt1, Run2_TauTau_legTriggerEff_MCFake_dm10);
            }
            
            if(mt2<30.0)
            {
                if (dm2==0) effMC2 = CBeff( pt2, Run2_TauTau_legTriggerEff_MCReal_dm0);
                if (dm2==1) effMC2 = CBeff( pt2, Run2_TauTau_legTriggerEff_MCReal_dm1);
                if (dm2==10) effMC2 = CBeff( pt2, Run2_TauTau_legTriggerEff_MCReal_dm10);
            }
            else
            {
                if (dm2==0) effMC2 = CBeff( pt2, Run2_TauTau_legTriggerEff_MCFake_dm0);
                if (dm2==1) effMC2 = CBeff( pt2, Run2_TauTau_legTriggerEff_MCFake_dm1);
                if (dm2==10) effMC2 = CBeff( pt2, Run2_TauTau_legTriggerEff_MCFake_dm10);
            }
            */
            
			returnSF = SF1 * SF2 * tauID1 * tauID2;
		}

        //No Shifts defined in 2016, keeping in case needed later
		else if( sysShift == 1) /* defined as up data over down mc */
		{
			effData1 = CBeff( pt1, Run2_TauTau_legTriggerEff_DataUP);
			effData2 = CBeff( pt2, Run2_TauTau_legTriggerEff_DataUP);

			effMC1 = CBeff( pt1, Run2_TauTau_legTriggerEff_McDOWN);
			effMC2 = CBeff( pt2, Run2_TauTau_legTriggerEff_McDOWN);

			if(effMC1!=0) SF1 = effData1/effMC1;
			if(effMC2!=0)  SF2 = effData2/effMC2;

			returnSF = SF1 * SF2;
		}
		else if( sysShift == -1) /* defined as down data over up mc */
		{
			effData1 = CBeff( pt1, Run2_TauTau_legTriggerEff_DataDOWN);
			effData2 = CBeff( pt2, Run2_TauTau_legTriggerEff_DataDOWN);

			effMC1 = CBeff( pt1, Run2_TauTau_legTriggerEff_McUP);
			effMC2 = CBeff( pt2, Run2_TauTau_legTriggerEff_McUP);

			if(effMC1!=0) SF1 = effData1/effMC1;
			if(effMC2!=0)  SF2 = effData2/effMC2;

			returnSF = SF1 * SF2;
		}

        if(returnIDISO) return tauID1 * tauID2;
        else if(returnTRIG) return SF1 * SF2;
        else if(returnTRACK) return 1.0;
        else return returnSF;
	}


	/* if reaches here then print a warning */
	std::cout<<"WARNING ---- Trigger x ID x ISO sf not available for CandidateEventType "<<R.getI("CandidateEventType")<<" using 1.0 \n";
	return returnSF;
}

