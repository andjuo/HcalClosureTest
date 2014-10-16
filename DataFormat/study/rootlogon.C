{
  gROOT->ProcessLine(".L dijet_PFNtuple.C+");
  if (1) {
    gROOT->ProcessLine(".L link_DijetRespCorrData.cc+");
    gROOT->ProcessLine(".L DijetRespCorrDataExtended.cc+");
  }
  if (1) {
    gROOT->ProcessLine(".L CPlot.cc+");
    gROOT->ProcessLine(".L Analyzer.cc+");
  }
}
