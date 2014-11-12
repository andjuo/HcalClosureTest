#ifndef helper_HH
#define helper_HH

#include <iostream>
#include <vector>

// -----------------------------------------------

inline
void HERE(const char *msg) { std::cout << msg << std::endl; }

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

#endif
