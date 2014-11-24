#ifndef helper_HH
#define helper_HH

#include <iostream>
#include <vector>
#include <TCanvas.h>
#include <TString.h>
#include <TSystem.h>

// -----------------------------------------------

#ifndef study_helper_HH
inline
void HERE(const char *msg) { std::cout << msg << std::endl; }
#endif

// -----------------------------------------------

template<class T>
inline
void printVec(const char *msg, const std::vector<T> &v) {
  std::cout << msg << std::endl;
  std::cout << "printVec of size=" << v.size() << "\n";
  for (unsigned int i=0; i<v.size(); ++i) {
    std::cout << "i=" << i << " " << v[i] << "\n";
  }
  std::cout << std::endl;
  return ;
}

// -----------------------------------------------

#ifndef HistoCollector_H
inline
void SaveCanvas(TCanvas* canv, const TString &canvName, TString destDir="plots") {
  gSystem->mkdir(destDir,kTRUE);
  gSystem->mkdir(destDir+TString("/png"),kTRUE);
  gSystem->mkdir(destDir+TString("/pdf"),kTRUE);
  gSystem->mkdir(destDir+TString("/root"),kTRUE);

  TString saveName=destDir+TString("/png/");
  saveName+=canvName;
  saveName+=".png";
  canv->SaveAs(saveName);
  saveName.ReplaceAll("png","pdf");
  canv->SaveAs(saveName);
  saveName.ReplaceAll("pdf","root");
  canv->SaveAs(saveName);
  return;
}
#endif

// -----------------------------------------------

#endif
