{

  gROOT->ProcessLine(".x MitStyleRemix.cc+");
  gROOT->ProcessLine(".x CPlot.cc+");

  if (0) { // to study dijet data. Old version (before Nov 20, 2014)
    gROOT->ProcessLine(".L link_DijetRespCorrData.cc+");
  }

  if (0) { // to skim dijet data
    gROOT->ProcessLine(".L link_PFDijetTree.h+");
    gROOT->ProcessLine(".x link_GammaJetFit.cc+");
  }

  if (1) { // study gamma+jet data
    gROOT->ProcessLine(".x pf_gammajettree.cc+");
    gROOT->ProcessLine(".x link_GammaJetFit.cc+");
    gROOT->ProcessLine(".x link_GammaJetFitAnalyzer.cc+");
  }

}
