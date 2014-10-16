#include "link_DijetRespCorrData.h"
#include "../src/DijetRespCorrData.cc"

#ifdef __CINT__
#pragma link C++ class DijetRespCorrDatum;
#pragma link C++ class DijetRespCorrData;
#endif


// ----------------------------------------------------------------------

std::ostream& operator<<(std::ostream& out, const DijetRespCorrDatum &d) {
  out << "DijetRespCorrDatum:\n";
  out << " weight=" << d.GetWeight() << "\n";
  out << " tag jet: (eta,phi) = (" << d.GetTagEta() << ","
      << d.GetTagPhi() << ")\n";
  out << " tag EcalE=" << d.GetTagEcalE() << "\n";
  std::map<Int_t,Double_t> E;
  d.GetTagHcalE(E);
  out << " tag HcalE= ";
  for (std::map<Int_t,Double_t>::const_iterator it= E.begin();
       it != E.end(); ++it) {
    out << " (" << it->first << "," << it->second << ")";
  }
  out << "\n";

  out << " probe jet: (eta,phi) = (" << d.GetProbeEta() << ","
      << d.GetProbePhi() << ")\n";
  out << " probe EcalE=" << d.GetProbeEcalE() << "\n";

  d.GetProbeHcalE(E);
  out << " probe HcalE= ";
  for (std::map<Int_t,Double_t>::const_iterator it= E.begin();
       it != E.end(); ++it) {
    out << " (" << it->first << "," << it->second << ")";
  }
  out << "\n";

  out << " third jet (px,py)=(" << d.GetThirdJetPx() << ","
      << d.GetThirdJetPy() << ")\n";
  out << " candTrack info ignored (N=" << d.GetCandTrackN() << ")\n";
  return out;
}

// ----------------------------------------------------------------------
