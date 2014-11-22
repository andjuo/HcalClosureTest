#ifndef study_helper_HH
#define study_helper_HH

#include <TH1D.h>
#include<iostream>
#include <vector>

// ---------------------------------------------------------

inline
void HERE(const char *msg) { std::cout << msg << std::endl; }

// ---------------------------------------------------------

template<class T>
inline
void HERE(const char *msg, const T& val) {
  std::cout << msg << std::endl;
  std::cout << val << std::endl;
}

// ---------------------------------------------------------

template<class T>
inline
void printVec(const char *msg, const std::vector<T> &vec, int printEol=0) {
  std::cout << msg << " [" << vec.size() << "]: ";
  if (printEol) std::cout << "\n";
  for (unsigned int i=0; i<vec.size(); ++i) {
    std::cout << " " << vec[i];
    if (printEol) std::cout << "\n";
  }
  if (!printEol) std::cout << "\n";
}

// ---------------------------------------------------------

#ifndef GammaJetFitData_H_
template <class type_t>
inline
std::ostream& operator<<(std::ostream& out, const std::vector<type_t> *vec) {
  out << " [" << vec->size() << "]: ";
  for (unsigned int i=0; i<vec->size(); ++i) {
    out << " " << (*vec)[i];
  }
  return out;
}
#endif

// ---------------------------------------------------------

//------------------------------------------------------------------------------------------------------------------------

inline
int printHisto(std::ostream& out, const TH1D* histo, int exponent=0,
	       int maxLines=-1) {
  if (!histo) {
    out << "printHisto: histo is null\n";
    return 0;
  }
  char buf[100];
  const char *format= (exponent) ?
    " %5.2f-%5.2f    %e    %e\n" :
    " %5.2f-%5.2f    %f    %f\n";

  out << "values of " << histo->GetName() << "\n";
  int imax=histo->GetNbinsX();
  int truncated=0;
  if ((maxLines>0) && (imax>maxLines)) { imax=maxLines; truncated=1; }
  for(int i=1; i<=imax; i++) {
    double x=histo->GetBinLowEdge(i);
    double w=histo->GetBinWidth(i);
    sprintf(buf,format,
	    x,x+w,histo->GetBinContent(i),histo->GetBinError(i));
    out << buf;
  }
  if (truncated) out << "... output truncated to " << maxLines << " lines\n";
  return 1;
}

//------------------------------------------------------------------------------------------------------------------------

inline int printHisto(const TH1D* histo, int exponent=0, int maxLines=-1)
{ return printHisto(std::cout, histo, exponent, maxLines); }

//------------------------------------------------------------------------------------------------------------------------

#endif
