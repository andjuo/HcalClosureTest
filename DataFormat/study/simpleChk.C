#include "dijet_PFNtuple.h"
#include <cmath>

void simpleChk(int ientry=0, int checkCase=-1) {
  dijet_PFNtuple inpData;
  Long64_t nentries = inpData.fChain->GetEntriesFast();
  if (nentries<=0) return;

  inpData.LoadTree(ientry);
  inpData.GetEntry(ientry);

  if (0 || (checkCase==1)) {
    // check pT of the jet

    double E= inpData.tpfjet_E;
    double p= inpData.tpfjet_p;
    double pT=inpData.tpfjet_pt;
    double eta=inpData.tpfjet_eta;
    double phi=inpData.tpfjet_phi;

    std::cout << "E=" << E << ", p=" << p << ", pT=" << pT
	      << ", eta=" << eta << ", phi=" << phi << "\n";
    std::cout << "E/cosh(eta) = " << (E/std::cosh(eta)) << "\n";
    std::cout << "p/cosh(eta) = " << (p/std::cosh(eta))
	      << " (expected pT=" << pT << ")\n";
    std::cout << "sqrt(E^2-p^2) = " << sqrt(E*E-p*p) << "\n";

  }
  else if (0 || (checkCase==2)) {
    // check energy additivity

    double E= inpData.tpfjet_E;
    double p= inpData.tpfjet_p;
    double pT=inpData.tpfjet_pt;
    double eta=inpData.tpfjet_eta;
    double phi=inpData.tpfjet_phi;
    double nonHadE= inpData.tpfjet_unkown_E + inpData.tpfjet_electron_E +
      inpData.tpfjet_muon_E + inpData.tpfjet_photon_E;
    double hadE=0;
    for (int i=0; i<inpData.tpfjet_had_n; ++i) {
      hadE += inpData.tpfjet_had_E->at(i);
    }

    std::cout << "E=" << E << ", p=" << p << ", pT=" << pT
	      << ", eta=" << eta << ", phi=" << phi << "\n";
    std::cout << "E/cosh(eta) = " << (E/std::cosh(eta)) << "\n";
    std::cout << "p/cosh(eta) = " << (p/std::cosh(eta)) << "\n";
    std::cout << "sqrt(E^2-p^2) = " << sqrt(E*E-p*p) << "\n";
    std::cout << " nonHadE+hadE= " << nonHadE << " + " << hadE
	      << " = " << (nonHadE+hadE) << " (expected E=" << E << ")\n";
  }
  else if (0 || (checkCase==3)) {
    // check pT components

    double E= inpData.tpfjet_E;
    double p= inpData.tpfjet_p;
    double pT=inpData.tpfjet_pt;
    double eta=inpData.tpfjet_eta;
    double phi=inpData.tpfjet_phi;
    double nonHadPx= inpData.tpfjet_unkown_px + inpData.tpfjet_electron_px +
      inpData.tpfjet_muon_px + inpData.tpfjet_photon_px;
    double nonHadPy= inpData.tpfjet_unkown_py + inpData.tpfjet_electron_py +
      inpData.tpfjet_muon_py + inpData.tpfjet_photon_py;
    double hadPx=0, hadPy=0;
    for (int i=0; i<inpData.tpfjet_had_n; ++i) {
      hadPx += inpData.tpfjet_had_px->at(i);
      hadPy += inpData.tpfjet_had_py->at(i);
    }

    std::cout << "E=" << E << ", p=" << p << ", pT=" << pT
	      << ", eta=" << eta << ", phi=" << phi << "\n";
    std::cout << "E/cosh(eta) = " << (E/std::cosh(eta)) << "\n";
    std::cout << "p/cosh(eta) = " << (p/std::cosh(eta)) << "\n";
    std::cout << "sqrt(E^2-p^2) = " << sqrt(E*E-p*p) << "\n";
    std::cout << " nonHadPx+hadPx= " << nonHadPx << " + " << hadPx
	      << " = " << (nonHadPx+hadPx) << "\n";
    std::cout << " nonHadPy+hadPy= " << nonHadPy << " + " << hadPy
	      << " = " << (nonHadPy+hadPy) << "\n";
    std::cout << " nonHadPt+hadPt= sqrt[("
	      << sqrt(pow(nonHadPx,2)+pow(nonHadPy,2)) << ")^2 + ("
	      << sqrt(pow(hadPx,2)+pow(hadPy,2)) << ")^2] = "
	    << sqrt(pow(nonHadPx,2)+pow(nonHadPy,2)+pow(hadPx,2)+pow(hadPy,2))
	      << "\n";
    std::cout << " nonHadPt+hadPt= sqrt[("
	      << (nonHadPx+hadPx) << ")^2 + ("
	      << (nonHadPy+hadPy) << ")^2] = "
	      << sqrt(pow(nonHadPx+hadPx,2)+pow(nonHadPy+hadPy,2))
	      << "\n";
  }
  else if (0 || (checkCase==4)) {
    // check the meaning of rawHcalE and EcalE
    double sum=0;
    for (int i=0; i<inpData.tpfjet_had_n; ++i) {
      if (0) {
	std::cout << "tpfjet_had_n=" << inpData.tpfjet_had_n << ", sizes="
		  << inpData.tpfjet_had_E->size() << " "
		  << inpData.tpfjet_had_EcalE->size() << " "
		  << inpData.tpfjet_had_rawHcalE->size() << " "
		  << inpData.tpfjet_had_emf->size() << "\n";
      }

      double had_E= inpData.tpfjet_had_E->at(i);
      double had_EcalE= inpData.tpfjet_had_EcalE->at(i);
      double had_rawHcalE= inpData.tpfjet_had_rawHcalE->at(i);
      double had_emf= inpData.tpfjet_had_emf->at(i);
      int had_Id= inpData.tpfjet_had_id->at(i);
      int trackIdx= inpData.tpfjet_had_candtrackind->at(i);
      int nTowers= inpData.tpfjet_had_ntwrs->at(i);
      std::cout << "i=" << i << ", had_Id=" <<  had_Id
		<< ", had_E=" << had_E
		<< ", had_EcalE=" << had_EcalE << ", had_rawHcalE="
		<< had_rawHcalE << ", had_emf=" << had_emf;
      if ((nTowers==0) && (trackIdx>-1)) std::cout << "; no recHits, has track";
      std::cout << "\n";
      sum+= had_E;
    }
    std::cout << "(sum had_E)=" << sum << "\n";
  }
  else if (0 || (checkCase==5)) {
    // check decomposition into towers
    std::cout << "tpfjet_ntwrs=" << inpData.tpfjet_ntwrs << "\n"; 
    int oldHadInd=-1;
    double totSumWithFrac=0;
    double sum=0, sumWithFrac=0;
    for (int i=0; i<inpData.tpfjet_ntwrs; ++i) {
      int hadIdx= inpData.tpfjet_twr_hadind->at(i);
      if (hadIdx != oldHadInd) {
	if (oldHadInd!=-1) {
	  std::cout << "oldHadInd=" << oldHadInd
		    << ", sum=" << sum
		    << ", sumWithFrac=" << sumWithFrac
		    << "\n";
	}
	oldHadInd= hadIdx;
	sum=0;
	sumWithFrac=0;
      }
      std::cout << "iTwr=" << i << ", hadIdx=" << hadIdx
		<< ", hade=" << inpData.tpfjet_twr_hade->at(i)
		<< ", frac=" << inpData.tpfjet_twr_frac->at(i)
		<< ", subdet=" << inpData.tpfjet_twr_subdet->at(i)
		<< "\n";
      double hadE= inpData.tpfjet_twr_hade->at(i);
      double hadEWithFrac= hadE * inpData.tpfjet_twr_frac->at(i);
      sum+=  hadE;
      sumWithFrac += hadEWithFrac;
      totSumWithFrac += hadEWithFrac;
    }
    std::cout << "oldHadInd=" << oldHadInd
	      << ", sum=" << sum
	      << ", sumWithFrac=" << sumWithFrac
	      << "\n";
    std::cout << " --- totSumWithFrac=" << totSumWithFrac << "\n";
  }
  else if (1 || (checkCase==6)) {
    // check the balancing equation
    std::cout << "tpfjet_ntwrs=" << inpData.tpfjet_ntwrs << "\n"; 
    double sumHcalE=0;
    double sumHadEcalE=0;
    double sumCandNoRecHits=0;
    // for hadronic energy -- sum over towers
    for (int i=0; i<inpData.tpfjet_ntwrs; ++i) {
      double hadE= inpData.tpfjet_twr_hade->at(i);
      if (hadE<0) continue;
      int clusterIdx= inpData.tpfjet_twr_clusterind->at(i);
      if (clusterIdx ||
	  inpData.tpfjet_cluster_dR->at(clusterIdx)
	  ) {
	sumHcalE+= hadE * inpData.tpfjet_twr_frac->at(i);
      }
    }

    //. for remaining energy -- sum over hadrons
    for (int i=0; i<inpData.tpfjet_had_n; ++i) {
      if (inpData.tpfjet_had_id->at(i) < 2) {
	// it's either charged or neutral hadron
	sumHadEcalE += inpData.tpfjet_had_EcalE->at(i);
      }
    }

    for (int i=0; i<inpData.tpfjet_had_n; ++i) {
      int trackIdx=inpData.tpfjet_had_candtrackind->at(i);
      if ((inpData.tpfjet_had_ntwrs->at(i) == 0) && // no recHits
	  (trackIdx > -1)) // charged hadron
	{
	  double trPx= inpData.tpfjet_candtrack_px->at(trackIdx);
	  double trPy= inpData.tpfjet_candtrack_py->at(trackIdx);
	  double trPz= inpData.tpfjet_candtrack_pz->at(trackIdx);
	  double trE = sqrt(pow(trPx,2) + pow(trPy,2) + pow(trPz,2));
	  double ecalE= inpData.tpfjet_had_EcalE->at(i);
	  std::cout << "i=" << i << ", trE=" << trE << ", ecalE=" << ecalE
		    << ", diff=" << (trE-ecalE) << "\n";
	  sumCandNoRecHits += (trE-ecalE);
      }
    }

    std::cout << "sumHcalE=" << sumHcalE << "\n";
    std::cout << "sumHadEcalE=" << sumHadEcalE << "\n";
    std::cout << "sumCandNoRecHits=" << sumCandNoRecHits << "\n";
    double hadronEn= sumHcalE + sumHadEcalE + sumCandNoRecHits;
    double otherEn = inpData.tpfjet_unkown_E + inpData.tpfjet_electron_E
      + inpData.tpfjet_muon_E + inpData.tpfjet_photon_E;
    std::cout << "totalEn= hadronEn + otherEn = "
	      << hadronEn << " + " << otherEn << " = "
	      << (hadronEn + otherEn)
	      << " (expected E=" << inpData.tpfjet_E << ")\n";
  }
}
