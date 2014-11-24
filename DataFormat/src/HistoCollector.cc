#ifdef __localRun
#  include "../interface/HistoCollector.h"
#else
#  include "HcalClosureTest/DataFormat/interface/HistoCollector.h"
#endif

#include <TFile.h>
#include <TObjString.h>
#include <TSystem.h>

#ifdef __localRun
#ifdef __CINT__
#pragma link C++ class HistoCollector_t;
#endif
#endif

// ----------------------------------------------------------------

std::vector<TString> packMessages(int msgLineCount, ...)
{
  std::vector<TString> knownMessages;
  if (msgLineCount>0) {
    knownMessages.reserve(msgLineCount);
    va_list vl;
    va_start(vl,msgLineCount);
    for (int i=0; i<msgLineCount; ++i) {
      typedef const char *constCharPtr;
      TString temp= TString(va_arg(vl,constCharPtr));
      knownMessages.push_back(temp);
    }
    va_end(vl);
  }
  return knownMessages;
}

// ----------------------------------------------------------------

void SaveCanvas(TCanvas* canv, const TString &canvName, TString destDir)
{
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

// ----------------------------------------------------------------
// ----------------------------------------------------------------

void HistoCollector_t::Clear()
{
  fHistos1Dv.clear();
  fHistos2Dv.clear();
  fCanvasV.clear();
  fMsgV.clear();
  fTag.Clear();
}

// ----------------------------------------------------------------

void HistoCollector_t::Add(const std::vector<TString> &vec) {
  fMsgV.reserve(fMsgV.size() + vec.size());
  for (unsigned int i=0; i<vec.size(); ++i) {
    fMsgV.push_back(vec[i]);
  }
  return;
}

// ----------------------------------------------------------------

int HistoCollector_t::Add(int msgLineCount, ... )
{
  int count=0;
  if (msgLineCount>0) {
    fMsgV.reserve(fMsgV.size() + msgLineCount);
    va_list vl;
    va_start(vl,msgLineCount);
    for (int i=0; i<msgLineCount; ++i) {
      typedef const char *constCharPtr;
      TString temp= TString(va_arg(vl,constCharPtr));
      if (temp.Length()) {
	count++;
	fMsgV.push_back(temp);
      }
    }
    va_end(vl);
  }
  return count;
}

// ----------------------------------------------------------------

int HistoCollector_t::SaveToFile(TString fName) 
{
  TFile fout(fName,"recreate");
  if (!fout.IsOpen()) {
    std::cout << "SaveVectorsToFile(" << fName
	      << "): failed to create the file\n";
    return 0;
  }
  fout.cd();

  for (unsigned int i=0; i<fHistos1Dv.size(); ++i) fHistos1Dv[i]->Write();
  for (unsigned int i=0; i<fHistos2Dv.size(); ++i) fHistos2Dv[i]->Write();
  for (unsigned int i=0; i<fCanvasV.size(); ++i) fCanvasV[i]->Write();

  if (fMsgV.size()) {
    TObjString ostr;
    ostr.Write("info_fields_follow");

    for (unsigned int i=0; i<fMsgV.size(); ++i) {
      ostr.Write(fMsgV[i]);
    }
  }

  // put time stamp
  time_t ltime;
  ltime=time(NULL);
  TString str = TString(asctime(localtime(&ltime)));
  if (str[str.Length()-1]=='\n') str.Remove(str.Length()-1,1);
  TObjString date(str);
  date.Write(str.Data());

  fout.Close();
  std::cout << "File <" << fout.GetName() << "> created\n";
  return 1;

}

// ----------------------------------------------------------------

int HistoCollector_t::SaveCanvases(TString outFileNameTag) const
{
  TString destDir="plots_" + outFileNameTag;
  for (unsigned int i=0; i<fCanvasV.size(); ++i) {
    TString figName="fig-" + outFileNameTag + TString("-") +
      TString(fCanvasV[i]->GetName());
    SaveCanvas(fCanvasV[i],figName,destDir);
  }
  return int(fCanvasV.size());
}

// ----------------------------------------------------------------
