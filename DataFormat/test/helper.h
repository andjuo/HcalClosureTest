#ifndef helper_HH
#define helper_HH

#include <iostream>
#include <vector>
#include <TCanvas.h>
#include <TString.h>
#include <TSystem.h>

// -----------------------------------------------

inline void helper() { std::cout << "helper.h loaded\n"; }

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

// Assuming "x,y|z,t,w" format put x,y and z,t,w into two arrays
//
inline
int identifyTriggerIndices(const TString str,
			   std::vector<int> &arr1, std::vector<int> &arr2,
			   int debug=0) {
  arr1.clear(); arr2.clear();
  if (str.Length()<=1) {
    std::cout << "incorrect string <" << str << ">\n";
    return 0;
  }

  int divPos= str.Index('|');
  if (divPos==-1) {
    std::cout << "divPos could not be located\n";
    return 0;
  }

  const char *s= str.Data();
  if (s[0]!='|') arr1.push_back(atoi(s));
  for (int i=0; i<str.Length(); ++i) {
    if ((s[i]==',') || (s[i]=='|')) {
      if (i<divPos) arr1.push_back(atoi(s+i+1));
      else {
	if (i+1 < str.Length())
	  arr2.push_back(atoi(s+i+1));
      }
    }
  }

  if (debug) {
    std::cout << "identified " << arr1.size() << " indices I: ";
    for (unsigned int i=0; i<arr1.size(); ++i) std::cout << " " << arr1[i];
    std::cout << "\n";
    std::cout << "identified " << arr2.size() << " indices II: ";
    for (unsigned int i=0; i<arr2.size(); ++i) std::cout << " " << arr2[i];
    std::cout << "\n";
  }

  return 1;
}

// -----------------------------------------------

#endif
