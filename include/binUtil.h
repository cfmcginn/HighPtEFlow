#ifndef BINUTIL_H
#define BINUTIL_H

//vector
#include <vector>

//ROOT
#include "TMath.h"

inline void getLinBins(const Float_t lower, const Float_t higher, const Int_t nBins, Double_t bins[])
{
  bins[0] = lower;
  bins[nBins] = higher;

  Float_t interval = (bins[nBins] - bins[0])/nBins;

  for(Int_t iter = 1; iter < nBins; iter++){
    bins[iter] = bins[0] + iter*interval;
  }

  return;
}

inline void getLinBins(const Float_t lower, const Float_t higher, const Int_t nBins, std::vector<double>* bins)
{
  bins->clear();
  for(int bIX = 0; bIX < nBins+1; ++bIX){
    bins->push_back(0.0);
  }

  bins->at(0) = lower;
  bins->at(nBins) = higher;

  Float_t interval = (bins->at(nBins) - bins->at(0))/nBins;

  for(Int_t iter = 1; iter < nBins; iter++){
    bins->at(iter) = bins->at(0) + iter*interval;
  }

  return;
}

inline void getLogBins(const Float_t lower, const Float_t higher, const Int_t nBins, Double_t bins[])
{
  Float_t logBins[nBins+1];
  bins[0] = lower;
  bins[nBins] = higher;

  logBins[0] = TMath::Log10(lower);
  logBins[nBins] = TMath::Log10(higher);

  Float_t interval = (logBins[nBins] - logBins[0])/nBins;

  for(Int_t iter = 1; iter < nBins; iter++){
    logBins[iter] = logBins[0] + iter*interval;
    bins[iter] = TMath::Power(10, logBins[iter]);
  }

  return;
}

#endif
