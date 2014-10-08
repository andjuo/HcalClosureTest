#ifndef _HCALCLOSURETEST_DATAFORMAT_RUNPFJETCORR_H_
#define _HCALCLOSURETEST_DATAFORMAT_RUNPFJETCORR_H_

#include <vector>
#include <iostream>

#include "TFile.h"
#include "TChain.h"
#include "TH1D.h"
//#include "TSystem.h"
//#include "TROOT.h"

#ifndef __localRun
#  include "HcalClosureTest/DataFormat/src/DijetRespCorrData.cc"
#else
#  include "../src/DijetRespCorrData.cc"
#endif


#ifdef __localRun
#ifdef __CINT__
#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ class DijetRespCorrDatum;
#pragma link C++ class DijetRespCorrData;
#endif
#endif


#endif
