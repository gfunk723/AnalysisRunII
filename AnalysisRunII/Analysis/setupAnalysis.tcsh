git clone https://github.com/CMS-HTT/LeptonEff-interface.git HTT-utilities
cd HTT-utilities/LepEffInterface
git clone https://github.com/CMS-HTT/LeptonEfficiencies.git data
cd -

git clone https://github.com/CMS-HTT/CorrectionsWorkspace.git
cd CorrectionsWorkspace
git submodule update --init
python makeCorrectionsWorkspace.py
cd -

