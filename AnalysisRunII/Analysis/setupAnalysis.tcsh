git clone https://github.com/CMS-HTT/LeptonEff-interface.git HTT-utilities
cd HTT-utilities/LepEffInterface
git clone https://github.com/CMS-HTT/LeptonEfficiencies.git data
cd ..
git clone https://github.com/CMS-HTT/QCDModelingEMu.git QCDModelingEMu
cp QCDModelingEMu/data/QCD_weight_emu.root $CMSSW_BASE/src
cp QCDModelingEMu/data/QCD_weight_emu_nodzeta.root $CMSSW_BASE/src
cd ..

git clone https://github.com/CMS-HTT/CorrectionsWorkspace.git
cd CorrectionsWorkspace
git submodule update --init
python makeCorrectionsWorkspace.py
cp CrystalBallEfficiency_cxx.so ../
cd -

