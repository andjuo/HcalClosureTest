#ifndef HistoCollector_H
#define HistoCollector_H

#include <TROOT.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <vector>
#include <iostream>

// -----------------------------------------------------------

std::vector<TString> packMessages(int msgLineCount, ...);
void SaveCanvas(TCanvas* canv, const TString &canvName, TString destDir="plots");

// -----------------------------------------------------------
// -----------------------------------------------------------

class HistoCollector_t {
 protected:
  std::vector<TH1D*> fHistos1Dv;
  std::vector<TH2D*> fHistos2Dv;
  std::vector<TCanvas*> fCanvasV;
  std::vector<TString> fMsgV;
  TString fTag; // not used

 public:
 HistoCollector_t() :
  fHistos1Dv(), fHistos2Dv(), fCanvasV(), fMsgV(), fTag()
    {}

 HistoCollector_t(HistoCollector_t &hc) :
  fHistos1Dv(hc.fHistos1Dv), fHistos2Dv(hc.fHistos2Dv),
    fCanvasV(hc.fCanvasV), fMsgV(hc.fMsgV), fTag(hc.fTag)
    {}

  void Clear();

  TCanvas* GetCanvas(int idx) const { return fCanvasV[idx]; }

  void Add(TH1D *h) { fHistos1Dv.push_back(h); }
  void Add(TH2D *h) { fHistos2Dv.push_back(h); }
  void Add(TCanvas *c) { fCanvasV.push_back(c); }
  void Add(TH1D *h, TCanvas *c) { fHistos1Dv.push_back(h); fCanvasV.push_back(c); }
  void Add(TH1D *h1, TH1D *h2, TCanvas *c) { fHistos1Dv.push_back(h1); fHistos1Dv.push_back(h2); fCanvasV.push_back(c); }
  void Add(TH2D *h, TCanvas *c) { fHistos2Dv.push_back(h); fCanvasV.push_back(c); }

  void Add(TString msg) { fMsgV.push_back(msg); }
  void Add(const std::vector<TString> &vec);
  int  Add(int msgCount, ... );
  
  int SaveToFile(TString fName);
  int SaveCanvases(TString outFileNameTag) const;
};

// -----------------------------------------------------------

#endif
