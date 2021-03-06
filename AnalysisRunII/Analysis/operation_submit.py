#!/usr/local/bin/python
import subprocess
import os
import sys, getopt
import textwrap
import time

supportFileList = os.getcwd() + '/supportFiles/pDistPlots.root' + ',' + os.getcwd() + '/ZReweight/zpt_weights_2016_BtoH.root' + ',' + os.getcwd() + '/supportFiles/kfactors.root' + ',' + os.getcwd() + '/supportFiles/zzKfactorGenPt.root' + ',' + os.getcwd() + '/supportFiles/WWsigma.root'+ ',' + os.getcwd() + '/HTT-utilities/QCDModelingEMu/data/QCD_weight_emu.root'+ ',' + os.getcwd() + '/HTT-utilities/QCDModelingEMu/data/QCD_weight_emu_nodzeta.root' + ',' + os.getcwd() + '/HTT-utilities/LepEffInterface/data/Muon/Run2016BtoH/Muon_IdIso_IsoLt0p15_2016BtoH_eff.root' + ',' + os.getcwd() + '/HTT-utilities/LepEffInterface/data/Muon/Run2016BtoH/Muon_IsoMu24_OR_TkIsoMu24_2016BtoH_eff.root' + ',' + os.getcwd() + '/HTT-utilities/LepEffInterface/data/Muon/Run2016BtoH/Muon_IdIso_IsoLt0p2_2016BtoH_eff.root' + ',' + os.getcwd() + '/HTT-utilities/LepEffInterface/data/Electron/Run2016BtoH/Electron_IdIso_IsoLt0p1_eff.root' + ',' + os.getcwd() + '/HTT-utilities/LepEffInterface/data/Electron/Run2016BtoH/Electron_Ele25WPTight_eff.root' + ',' + os.getcwd() + '/HTT-utilities/LepEffInterface/data/Electron/Run2016BtoH/Electron_IdIso_IsoLt0p15_eff.root' + ',' + os.getcwd() + '/SUSYHiggsPtReweight/Reweight.root' + ',' + os.getcwd() + '/supportFiles/EfficienciesAndSF_ID_BCDEF.root' + ',' + os.getcwd() + '/supportFiles/EfficienciesAndSF_ISO_BCDEF.root' + ',' + os.getcwd() + '/supportFiles/EfficienciesAndSF_ID_GH.root' + ',' + os.getcwd() + '/supportFiles/EfficienciesAndSF_ISO_GH.root' + ',' + os.getcwd() + '/supportFiles/EfficienciesAndSF_TRIG_GH.root' + ',' + os.getcwd() + '/supportFiles/EfficienciesAndSF_TRIG_BCDEF.root' + ',' + os.getcwd() + '/CorrectionsWorkspace/htt_scalefactors_v16_4.root' + ',' + os.getcwd() + '/CorrectionsWorkspace/htt_scalefactors_sm_moriond_v2.root' + ',' + os.getcwd() + '/CorrectionsWorkspace/CrystalBallEfficiency.cxx' + ',' + os.getcwd() + '/CorrectionsWorkspace/CrystalBallEfficiency.h' + ',' + os.getcwd() + '/CorrectionsWorkspace/CrystalBallEfficiency_cxx.d' + ',' + os.getcwd() + '/CorrectionsWorkspace/CrystalBallEfficiency_cxx.so' + ',' + os.getcwd() + '/CorrectionsWorkspace/CrystalBallEfficiency_cxx_ACLiC_dict_rdict.pcm' + ',' + os.getcwd() + '/xmlM18/TMVAClassification_MLPBNN_MuTau_MZP600_MA0400_MDM100_nodes30_epoch1000.weights.xml' + ',' + os.getcwd() + '/xmlM18/TMVAClassification_MLPBNN_MuTau_MZP800_MA0400_MDM100_nodes30_epoch1000.weights.xml' + ',' + os.getcwd() + '/xmlM18/TMVAClassification_MLPBNN_MuTau_MZP1000_MA0400_MDM100_nodes30_epoch1000.weights.xml' + ',' + os.getcwd() + '/xmlM18/TMVAClassification_MLPBNN_MuTau_MZP1200_MA0400_MDM100_nodes30_epoch1000.weights.xml' + ',' + os.getcwd() + '/xmlM18/TMVAClassification_MLPBNN_EleTau_MZP600_MA0400_MDM100_nodes30_epoch1000.weights.xml' + ',' + os.getcwd() + '/xmlM18/TMVAClassification_MLPBNN_EleTau_MZP800_MA0400_MDM100_nodes30_epoch1000.weights.xml' + ',' + os.getcwd() + '/xmlM18/TMVAClassification_MLPBNN_EleTau_MZP1000_MA0400_MDM100_nodes30_epoch1000.weights.xml' + ',' + os.getcwd() + '/xmlM18/TMVAClassification_MLPBNN_EleTau_MZP1200_MA0400_MDM100_nodes30_epoch1000.weights.xml' + ',' + os.getcwd() + '/xmlM18/TMVAClassification_MLPBNN_TauTau_MZP600_MA0400_MDM100_nodes30_epoch1000.weights.xml' + ',' + os.getcwd() + '/xmlM18/TMVAClassification_MLPBNN_TauTau_MZP800_MA0400_MDM100_nodes30_epoch1000.weights.xml' + ',' + os.getcwd() + '/xmlM18/TMVAClassification_MLPBNN_TauTau_MZP1000_MA0400_MDM100_nodes30_epoch1000.weights.xml' + ',' + os.getcwd() + '/xmlM18/TMVAClassification_MLPBNN_TauTau_MZP1200_MA0400_MDM100_nodes30_epoch1000.weights.xml' + '\n'

execTitle = ""

def getOptions(argv):
  useCondor = 0
  tTree = ""
  inFile = ""
  operationList = []
  outputDir = ""
  try:
    opts, args = getopt.getopt(argv,"hct:i:o:d:",["tTree=","iFile=","opList=","outDir="])
  except getopt.GetoptError:
    print 'Usage: python operation_submit.py -c (to use condor) -i [TTree in FlatTuple] -i [input name and path to FlatTuple root file or txt file with list of root files] -o [operations to perform] -d [output path]'
    sys.exit(2)
  for opt, arg in opts:
    if opt == '-h':
      print 'Usage: python operation_submit.py -c (to use condor) -t [TTree in FlatTuple] -i [input name and path to FlatTuple root file or txt file with list of root files] -o [operations to perform] -d [output path]'
      sys.exit()
    elif opt == '-c':
      useCondor = 1
    elif opt in ("-t", "--tTree"):
      tTree = arg
    elif opt in ("-i", "--iFile"):
      inFile = arg
    elif opt in ("-o", "--opList"):
      operationList = arg.split()
    elif opt in ("-d", "--outDir"):
      outputDir = arg
  return [useCondor, tTree , inFile, operationList, outputDir]

def editCondorConfig(List):
  condorConfigFileTmp = open(os.getcwd() + '/run_condor.txt.tmp')
  condorConfigFileUse = open(os.getcwd() + '/run_condor.txt', 'w')
  for line in condorConfigFileTmp:
    if 'Place Arguments' in line:
      line = 'arguments = ' + List[1] + ' ' + List[2] + ' ' + List[3][0] + '\n'
    if 'Place Directory' in line:
      line = 'Initialdir = ' + List[4] + '\n'
    if 'Place Input' in line:
      if '.txt' in List[2]:
        line = 'transfer_input_files = ' + os.getcwd() + '/' + List[2] + ',' + supportFileList
      if '.root' in List[2]:
        line = 'transfer_input_files = ' + supportFileList
    condorConfigFileUse.write(line)
  condorConfigFileUse.close()
  condorConfigFileTmp.close()

allParams = getOptions(sys.argv[1:])
print allParams

if allParams[0] == 1:
  editCondorConfig(allParams)
  if os.path.isdir(allParams[4]) == 0:
    subprocess.call('mkdir -p ' + allParams[4], shell = True)
  subprocess.call('condor_submit run_condor.txt', shell = True)
else:
  subprocess.call('./run_analysis ' + allParams[1] + ' ' + allParams[2] + ' ' + allParams[3][0], shell = True)
