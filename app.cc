#include "simple_t2k.hh"
#include "TMath.h"

#include "RooArgList.h"
#include "RooAbsReal.h"
#include "RooAbsCategory.h"
#include "RooRealVar.h"
#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "TString.h"
#include <TVectorD.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TFrame.h>

using namespace std;

int main(int argc, char**argv){

 RooFitResult* res;
 Sterile * rep = new Sterile ("_rep");
   char formula[10];

 std::cout<<"start to run "<<std::endl;

 // vecInput1 is the CV for pulls while vecInput2 is the unc. for pulls
 TH1D* vecInput1 = new TH1D("","",11,0,11);
 TH1D* vecInput2 = new TH1D("","",11,0,11);
 //vecInput1->Print();

 vecInput1 ->SetBinContent(1, 1);
 vecInput2 ->SetBinContent(1, 1);

 std::cout<<"'ve set some inputs "<<std::endl;

 // binned from 0.5-9 MeV with bin width of 0.25 MeV
 int nBins = 30;
 TH1D* binHist = new TH1D("","",nBins+1,0,nBins+1);
 for(Int_t i=0;i<nBins+1;i++){
   binHist->SetBinContent(i+1, 0.5 + 0.25*i);
 }

 //TString fileLocation = "/gpfs/projects/McGrewGroup/gyang/REACTOR-related/";
 TString fileLocation = "./";
 rep->setFileLocation(fileLocation);

 rep->SetBinning(binHist);

 rep->SetMatrixNameDC(fileLocation+"tbd");
 std::vector<TH1D*> tempPredList = rep->preparePrediction(rep->getPullList(), false);
 rep->prepareData(rep->preparePrediction(rep->getPullList(), false));

 rep->setPull(vecInput1); 
 rep->setPullUnc(vecInput2);

 std::cout<<"ended up with setting us basic stuff "<<std::endl;
 rep->FillEv(rep->getPullList());

 RooArgList list("list");
 list.add(*rep);
 sprintf(formula,"%s","@0");
 RooFormulaVar* fcn = new RooFormulaVar("fit","fit",formula,list);

 // ******************************** Important setup here *************************************
 // *******************************************************************************************
 rep->setSysts(false);
 rep->ifEShiftHist(true);
 rep->setEShiftHist("Escalefraction.root");
 // *******************************************************************************************
 // ******************************************************************************************* 

 rep->getParVar(0)->setConstant(false);

 std::cout<<"------------  Getting current spectra "<<std::endl;
 std::vector<TH1D*> outPrediction = tempPredList; // rep->GetCurrentPrediction();
 std::vector<TH1D*> outData = rep->GetCurrentData(outPrediction);
 std::cout<<"------------  Have got current spectra "<<std::endl;

 // save the un-oscillated standard prediction and data spectra
 TFile* outputFile = new TFile("outputFigs.root","RECREATE");
 for(Int_t i=0;i<outPrediction.size();i++)
 {
   outPrediction[i]->Write(Form("outPrediction[%d]",i));
   if(i<5) outData[i]->Write(Form("outData[%d]",i));
 }
 outputFile -> Close();
 std::cout<<"first thing saved "<<std::endl;

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

 RooMinuit m(*fcn);
 m.setStrategy(2);
 Double_t callsEDM[2] = {10500., 1.e-6};
 Int_t irf = 0;

 gMinuit->mnexcm("MIGRAD",callsEDM,2,irf);
 m.migrad();
 res = m.save();
 double bestFit = res->minNll();
 std::cout<<"fit status code is : "<< res->status()<<std::endl;	
 //status = 0    : OK
 //status = 1    : Covariance was mad  epos defined
 //status = 2    : Hesse is invalid
 //status = 3    : Edm is above max
 //status = 4    : Reached call limit
 //status = 5    : Any other failure
 std::cout<<"quiality code of covariance matrix is : "<< res->covQual()<<std::endl;
 //status = -1 :  not available (inversion failed or Hesse failed)
 //status =  0 : available but not positive defined
 //status =  1 : covariance only approximate
 //status =  2 : full matrix but forced pos def
 //status =  3 : full accurate matrix

 //////////////////////	
 //inline Int_t status() const 
 //inline Int_t covQual() const 
 //inline Int_t numInvalidNLL() const 
 //inline Double_t edm() const 
 //inline Double_t minNll() const 
 /////////////////////

 outPrediction = rep->GetCurrentPrediction();
 //outData = rep->GetCurrentData(outPrediction);

 //std::cout<<"list of reactor pulls : "<<std::endl;
 //for(Int_t i=18;i<41;i++){std::cout<<" "<<rep->getPar(i)<<std::endl;}

 std::cout<<"size of output prediction list "<<outPrediction.size()<<std::endl;
 std::cout<<"size of output data list "<<outData.size()<<std::endl;
 //std::cout<<atoi(argv[1])<<" "<<atoi(argv[2])<<" "<<bestFit<<std::endl;

 //std::cout<<"result list: "<<std::endl;
 //std::cout<<"chi2: "<<bestFit <<std::endl;

 double bb =  rep->getParVar(2)->getAsymErrorLo();
 double dd =  rep->getParVar(2)->getAsymErrorHi();

 //cout<<"errors are "<<bb<<" "<<dd<<endl;

}
