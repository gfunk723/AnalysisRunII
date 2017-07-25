#include "makePlot.C"


void exec()
{
  makeLTauStack( "ET", "xtt_et.inputs-13TeV-mt.root", "et_inclusive", 3, "Total Mt", "GeV", 0, "#tau_{e}#tau_{h}", "Golden", 1, 0, 1, 200, 0);
  makeLTauStack( "MT", "xtt_mt.inputs-13TeV-mt.root", "mt_inclusive", 3, "Total Mt", "GeV", 0, "#tau_{m}#tau_{h}", "Golden", 1, 0, 1, 200, 0);
  makeLTauStack( "TT", "xtt_tt.inputs-13TeV-mt.root", "tt_inclusive", 3, "Total Mt", "GeV", 0, "#tau_{h}#tau_{h}", "Golden", 1, 0, 1, 200, 0);
}
