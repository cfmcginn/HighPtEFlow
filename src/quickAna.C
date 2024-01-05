//Author C. McGinn 2024.01.05

//c and cpp
#include <algorithm>
#include <iostream>
#include <map>
#include <string>

//ROOT
#include "TDatime.h"
#include "TEnv.h"
#include "TFile.h"
#include "TH1D.h"
#include "TMath.h"
#include "TSystem.h"
#include "TTree.h"

//include
#include "include/binUtil.h"
#include "include/etaPhiFunc.h"
#include "include/stringUtil.h"

int quickAna(std::string configName)
{
  //Grab date stamp
  TDatime date;
  const std::string dateStr = std::to_string(date.GetDate());

  //Define some maximums for param checking
  const Int_t nCentBinsMax = 8;  

  //Get our config and params for the job
  TEnv* config_p = new TEnv(configName.c_str());
  std::string inFileName = config_p->GetValue("INFILENAME", "");
  std::string jetColl = config_p->GetValue("JETCOLL", "");
  const ULong64_t entriesPerPrint = config_p->GetValue("ENTRIESPERPRINT", 100);
  std::string outFileName = config_p->GetValue("OUTFILENAME", "out.root");
  Int_t nCentBins = config_p->GetValue("NCENTBINS", 100);//Default above max to throw return 1 if not defined  
  std::vector<int> centBinsLow = strToVectI(config_p->GetValue("CENTBINSLOW", ""));
  std::vector<int> centBinsHigh = strToVectI(config_p->GetValue("CENTBINSHIGH", ""));
  
  std::vector<std::string> histNames = {"leadingJtPt", "subleadingJtPt", "diJtDPhi", "diJtAJ"};
  std::vector<std::string> histLabels = {"Leading Jet p_{T} [GeV]", "Subleading Jet p_{T} [GeV]", "#Delta#phi_{1,2}", "A_{J}"};  
  std::vector<std::string> histNamesUpper;

  std::map<std::string, int> binsNMap;
  std::map<std::string, double> binsLowHighMap;  
  std::map<std::string, bool> binsBoolMap;

  for(auto const & histName : histNames){
    std::string histNameUpper = histName;
    std::transform(histNameUpper.begin(), histNameUpper.end(), histNameUpper.begin(), ::toupper);

    histNamesUpper.push_back(histNameUpper);
    binsNMap[histNameUpper + "NBINS"] = config_p->GetValue((histNameUpper + "NBINS").c_str(), -1);
    binsLowHighMap[histNameUpper + "LOW"] = config_p->GetValue((histNameUpper + "LOW").c_str(), -1.0);
    binsLowHighMap[histNameUpper + "HIGH"] = config_p->GetValue((histNameUpper + "HIGH").c_str(), -1.0);
    binsBoolMap[histNameUpper + "DOLOG"] = config_p->GetValue((histNameUpper + "DOLOG").c_str(), false);    
    binsBoolMap[histNameUpper + "DOLIN"] = config_p->GetValue((histNameUpper + "DOLIN").c_str(), false);    
  }

  //Dump the bins map
  /*
  for(auto const & histName : histNamesUpper){
    std::cout << histName << std::endl;
    std::cout << " " << histName << "NBINS: " << binsNMap[histName + "NBINS"] << std::endl;
    std::cout << " " << histName << "LOW: " << binsLowHighMap[histName + "LOW"] << std::endl;
    std::cout << " " << histName << "HIGH: " << binsLowHighMap[histName + "HIGH"] << std::endl;
    std::cout << " " << histName << "DOLOG: " << binsBoolMap[histName + "DOLOG"] << std::endl;
    std::cout << " " << histName << "DOLIN: " << binsBoolMap[histName + "DOLIN"] << std::endl;
  }
  */
  
  //Output handling
  std::string outDir = "";
  if(outFileName.find("/") == std::string::npos){//if no output dir, make one w/ timestamp
    outDir = "output/" + dateStr;
    gSystem->Exec(("mkdir -p " + outDir).c_str());

    outFileName = outDir + "/" + outFileName;
  }
  
  TFile* outFile_p = new TFile(outFileName.c_str(), "RECREATE");
  TH1D *leadingJtPt_h[nCentBinsMax], *subleadingJtPt_h[nCentBinsMax], *diJtDPhi_h[nCentBinsMax], *diJtAJ_h[nCentBinsMax];
  
  const Int_t nMaxBins = 200;
  Double_t bins[nMaxBins+1];
  
  //Create our histograms looping over cent
  for(Int_t cI = 0; cI < nCentBins; ++cI){
    std::string centBinStr = "CentBin" + std::to_string(cI);

    for(unsigned int hI = 0; hI < histNamesUpper.size(); ++hI){
      std::string name = histNames[hI] + "_" + centBinStr + "_h";
      std::string label = ";" + histLabels[hI] + ";Counts";
      const Int_t nBins = binsNMap[histNamesUpper[hI] + "NBINS"];
      const Float_t binsLow = binsLowHighMap[histNamesUpper[hI] + "LOW"];
      const Float_t binsHigh = binsLowHighMap[histNamesUpper[hI] + "HIGH"];
      const Bool_t binsDoLog = binsBoolMap[histNamesUpper[hI] + "DOLOG"];
      const Bool_t binsDoLin = binsBoolMap[histNamesUpper[hI] + "DOLIN"];

      if(binsDoLin){
	getLinBins(binsLow, binsHigh, nBins, bins);
      }
      else if(binsDoLog){
	getLogBins(binsLow, binsHigh, nBins, bins);
      }
      else{
	std::cout << "UHOH: DO LIN/LOG BINS NOT DEFINED. return 1" << std::endl;
	return 1;	
      }
      
      if(histNamesUpper[hI] == "LEADINGJTPT") leadingJtPt_h[cI] = new TH1D(name.c_str(), label.c_str(), nBins, bins);
      else if(histNamesUpper[hI] == "SUBLEADINGJTPT") subleadingJtPt_h[cI] = new TH1D(name.c_str(), label.c_str(), nBins, bins);
      else if(histNamesUpper[hI] == "DIJTDPHI") diJtDPhi_h[cI] = new TH1D(name.c_str(), label.c_str(), nBins, bins);
      else if(histNamesUpper[hI] == "DIJTAJ") diJtAJ_h[cI] = new TH1D(name.c_str(), label.c_str(), nBins, bins);
      else{
	std::cout << "UHOH: HISTNAMESUPPER: " << histNamesUpper[hI] << " is undefined, return 1" << std::endl;
	return 1;
      }
    }
  }

  
  //Grab file and ttrees
  TFile* inFile_p = new TFile(inFileName.c_str(), "READ");
  TTree* jetTree_p = (TTree*)inFile_p->Get(jetColl.c_str());
  TTree* hiTree_p = (TTree*)inFile_p->Get("hiEvtAnalyzer/HiTree");
  
  //define variables for reading ttrees
  const Int_t nMaxJets = 500;
  Int_t nref_;
  Float_t jtpt_[nMaxJets], jteta_[nMaxJets], jtphi_[nMaxJets];

  Int_t hiBin_;
  
  //Set status/addresses
  jetTree_p->SetBranchStatus("*", 0);
  jetTree_p->SetBranchStatus("nref", 1);
  jetTree_p->SetBranchStatus("jtpt", 1);
  jetTree_p->SetBranchStatus("jteta", 1);
  jetTree_p->SetBranchStatus("jtphi", 1);

  jetTree_p->SetBranchAddress("nref", &nref_);
  jetTree_p->SetBranchAddress("jtpt", jtpt_);
  jetTree_p->SetBranchAddress("jteta", jteta_);
  jetTree_p->SetBranchAddress("jtphi", jtphi_);

  hiTree_p->SetBranchStatus("*", 0);
  hiTree_p->SetBranchStatus("hiBin", 1);

  hiTree_p->SetBranchAddress("hiBin", &hiBin_);
  
  //Get nentries and start loopin
  const ULong64_t nEntries = jetTree_p->GetEntries();

  std::cout << "Processing " << nEntries << " entries..." << std::endl;
  for(ULong64_t entry = 0; entry < nEntries; ++entry){
    if(entry%entriesPerPrint == 0) std::cout << " Entry " << entry << "/" << nEntries << "..." << std::endl;

    jetTree_p->GetEntry(entry);
    hiTree_p->GetEntry(entry);

    int centBinPos = -1;
    for(int cI = 0; cI < nCentBins; ++cI){
      if(hiBin_ >= centBinsLow[cI]*2 && hiBin_ < centBinsHigh[cI]*2){
	centBinPos = cI;
	break;
      }
    }
    
    if(centBinPos < 0){
      //      std::cout << "No Pos for hiBin: " << hiBin_ << std::endl;
      continue;
    }

    Float_t leadPt = -999;
    Float_t subleadPt = -999;
    Int_t leadPos = -1;
    Int_t subleadPos = -1;

    for(Int_t jI = 0; jI < nref_; ++jI){
      if(jtpt_[jI] > leadPt){
	subleadPt = leadPt;
	subleadPos = leadPos;

	leadPt = jtpt_[jI];
	leadPos = jI;
      }
      else if(jtpt_[jI] > subleadPt){
	subleadPt = jtpt_[jI];
	subleadPos = jI;
      }
    }

    if(leadPt < binsLowHighMap["LEADINGJTPTLOW"]) continue;
    if(subleadPt < binsLowHighMap["SUBLEADINGJTPTLOW"]) continue;
    
    leadingJtPt_h[centBinPos]->Fill(leadPt);
    subleadingJtPt_h[centBinPos]->Fill(subleadPt);

    Float_t dPhi = TMath::Abs(getDPHI(jtphi_[leadPos], jtphi_[subleadPos]));

    diJtDPhi_h[centBinPos]->Fill(dPhi);

    Float_t aj = (leadPt - subleadPt)/(leadPt + subleadPt);
    diJtAJ_h[centBinPos]->Fill(aj);
  }
  
  
  //cleanup
  inFile_p->Close();
  delete inFile_p;

  outFile_p->cd();
  for(Int_t cI = 0; cI < nCentBins; ++cI){
    leadingJtPt_h[cI]->Write("", TObject::kOverwrite);
    subleadingJtPt_h[cI]->Write("", TObject::kOverwrite);
    diJtDPhi_h[cI]->Write("", TObject::kOverwrite);
    diJtAJ_h[cI]->Write("", TObject::kOverwrite);

    delete leadingJtPt_h[cI];
    delete subleadingJtPt_h[cI];
    delete diJtDPhi_h[cI];
    delete diJtAJ_h[cI];
  }
  outFile_p->Close();
  delete outFile_p;
  
  delete config_p;
    
  return 0;
}

int main(int argc,  char* argv[])
{
  if(argc != 2){
    std::cout << "Usage: ./bin/quickAna.exe <configName>. return 1" << std::endl;
    return 1;
  }
  
  int retVal = 0;
  retVal += quickAna(argv[1]);
  return retVal;
}
