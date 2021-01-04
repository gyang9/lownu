/*
 *  low-nu simple combined fit
 *
 *  Author: Guang Yang
 */
#include "simple_fit.hh"
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
#include <TFile.h>

using namespace std;

  Sterile ::Sterile (const char* name) 
  : RooAbsReal(name,name)
{

// there will be: pull 0-9 totally 10 flux systs.

  _pulls     = new RooListProxy("_pulls","_pulls",this);
  RooRealVar* Par1 = new RooRealVar("flux0","par1",0,0,100);
  RooRealVar* Par2 = new RooRealVar("flux1","par2",0,0,100); // TMath::ASin(TMath::Sqrt(0.95))/2.,0,100);
  RooRealVar* Par3 = new RooRealVar("flux2","par3",0,0,100);

  Par1->setConstant(false);
  Par2->setConstant(false);
  Par3->setConstant(false);

  _parlist.add(*Par1);
  _parlist.add(*Par2);
  _parlist.add(*Par3);

  this->addServerList(*_pulls);

  dataDC = new TH1D("","dataDC",30, 0.5, 8.);
};

Sterile ::~Sterile ()
{;}

//=================================================================================================================================
TMatrixD* Sterile::prepareCovMatrix(Int_t nBins, TVectorD* fVec) const
{

  TFile fMatrixDC(fileNameDC);

  TMatrixD* outMat = new TMatrixD( nBins , nBins);

  if(inSyst){
      // 3% error, diagonal only
      double errlist[100];
      for (int temp=0;temp<10;temp++)
	errlist[temp] = fLownuErr;
      TVectorD* errList = new TVectorD(nBins);
      for(Int_t i = 0; i< nBins; i++)
      {
          (*errList)[i] = errlist[i];
      }
      //-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
      for(Int_t i = 0; i< nBins; i++)
      {
        for(Int_t j =i ;j<nBins; j++) 
        {
	  if(i == j)	
            (*outMat)(i,j) = (*errList)[i] * (*errList)[j] * (*fVec)[i] * (*fVec)[j];
        }
      }
  }
  for(Int_t i = 0; i< nBins ; i++)
  {
    if((*outMat)(i,i) == 0) (*outMat)(i,i) += 0.0000000001;
  }

return outMat ;
}

Double_t Sterile ::FillEv( RooListProxy* _pulls ) const 
{

   //std::cout<<"in FillEv() "<<std::endl;
   int nBins = _nBins;

   std::vector<TH1D*> tempPredList = this->preparePrediction(_pulls, true);
   //std::cout<<"filled in new pediction "<<std::endl;

   //std::vector<TH1D*> tempDataList = this->prepareData();
   TH1D* predDC = tempPredList[0];

   TVectorD* fVec = new TVectorD(nBins );
   TVectorD* fData = new TVectorD(nBins );

   TH1D* tempVec[5]; 
   TH1D* tempDat[5]; 
   for(Int_t i=0;i<5;i++){
     tempVec[i] = new TH1D("","",100,0,10);
     tempDat[i] = new TH1D("","",100,0,10);
   }

   for(Int_t i=0;i<nBins;i++){	 
         (*fVec)[i]               = predDC->GetBinContent(i+1) ; //* this->surv_Prob( (thisE) , _pulls, 400)  ;
         tempVec[0] -> SetBinContent(i+1, (*fVec)[i]);
   }
   for(Int_t i=0;i<nBins;i++){
   	(*fData)[i]             = dataDC->GetBinContent(i+1)   ;
        tempDat[0] -> SetBinContent(i+1, (*fData)[i]);
   }
   //std::cout<<"data ready also "<<std::endl;
   //--------------------------------------------------------------------------------------------------------------------------------Data------------------------------------

   // scale to same total rate, doing shape only analysis 
   double scaling1 = tempDat[0]->Integral() / tempVec[0]->Integral();

   for(Int_t i=0;i<nBins;i++){
      (*fData)[i]             = 0; //dataDC->GetBinContent(i+1)   - (*fVec)[i] * scaling1  ;
   }

   for(Int_t i=0;i<nBins;i++){ 
      (*fData)[i]             = TMath::Abs( dataDC->GetBinContent(i+1)   - (*fVec)[i] * scaling1  );
      //std::cout<<dataDC->GetBinContent(i+1)<<" "<<dataDYB->GetBinContent(i+1)<<" "<<dataRENO->GetBinContent(i+1)<<" "<<dataNEOS->GetBinContent(i+1)<<" "<<dataPROS->GetBinContent(i+1)<<std::endl;
   }

   TMatrixD* covMat = this->prepareCovMatrix(nBins , fVec);
   covMat->Invert();

   TVectorD mulVec(*fData);
   mulVec *= (*covMat);

   Double_t currentResult = TMath::Abs(mulVec*(*fData));
   //std::cout<<"DC12_chi2 sans pull "<<currentResult<<std::endl;
  
   return (Double_t) currentResult ; 
}

//================================================================================================================================================5. Fill the Ev, Prediction===============
Double_t Sterile ::evaluate() const
{ 

Double_t matPart = this->FillEv(_pulls);//original FillEv is matPart

Double_t extraPull = this -> ExtraPull (_pulls);//same variable extraPull
Double_t tot = matPart + extraPull; //If needed, add pull terms here.

return tot;

}

Double_t Sterile ::ExtraPull (RooListProxy* _pulls) const
{
Double_t pullAdd = 0;
for(Int_t i=0;i<11;i++){
 pullAdd += TMath::Power(( ((RooAbsReal*)_pulls->at(i))->getVal() - (*pullCV)[i] ),2) / TMath::Power( (*pullUnc)[i],2) ;
    }
 //std::cout<<"DC13_extra pull penalty: "<<pullAdd<<std::endl;
 return pullAdd;
}

std::vector<TH1D*> Sterile:: preparePrediction(RooListProxy* _pulls, bool Iosc) const
{

  ((RooAbsReal*)_pulls->at(i+12))->getVal()	

  std::vector<TH1D*> predictionList;
  predictionList.clear();
/*
  predictionList.push_back(predDC);
  predictionList.push_back(predDYB);
  predictionList.push_back(predRENO);
  predictionList.push_back(predNEOS);
  predictionList.push_back(fpredDC);
  predictionList.push_back(fpredDYB);
  predictionList.push_back(fpredRENO);
  predictionList.push_back(fpredNEOS);
  predictionList.push_back(fpredPROS);
*/
  float receNuE;
  TFile file("");
  TTree* tree = (TTree*)file->Get("tree");
  tree->SetBranchAddress("recoNeutrinoE", &receNuE);
  TH1F predNuE("","",16,0,8);
  for (int i = 0; i < tree->GetEntries(); i++) {
      predNuE.Fill(recoNuE + ((RooAbsReal*)_pulls->at(i+12))->getVal() * ND_numubar_RHC* recoNuE);
  }
  predictionList.push_back(predNuE);

  //for (std::vector<TH1D*>::iterator pred = predictionList.begin();
  //        pred != predictionList.end()
  //        pred++) {
  //    *pred = pred + ((RooAbsReal*)_pulls->at(i+12))->getVal()* 
  //}

  //std::cout<<"************* before folding "<<predPROS->Integral()<<std::endl;
  //TH1D* temp(predPROS);
  //TH1D* fpredPROS = this->folding(temp);
  // seems pushing back the prompt energy spectrum
  return predictionList;

}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////==============================7.5 print prediction with folding===


std::vector<TH1D*> Sterile:: prepareData(std::vector<TH1D*> tempPredList) const// 
{

  //dataDC = new TH1D("","",_nBins,binEdge[0],binEdge[_nBins]);
  //dataDYB = new TH1D("","",_nBins,binEdge[0],binEdge[_nBins]);
  //dataRENO = new TH1D("","",_nBins,binEdge[0],binEdge[_nBins]);
  //dataNEOS = new TH1D("","",_nBins,binEdge[0],binEdge[_nBins]);

  TGraph* gradataDC   = new TGraph(fileLocation+"/data/dataDC.txt", "%*lg %lg %lg", "");
  TGraph* gradataDYB  = new TGraph(fileLocation+"/data/dataDYB.txt", "%*lg %lg %lg", "");
  TGraph* gradataRENO = new TGraph(fileLocation+"/data/dataRENO.txt", "%*lg %lg %lg", "");
  TGraph* gradataNEOS = new TGraph(fileLocation+"/data/dataNEOS.txt", "%*lg %lg %lg", "");
  TGraph* gradataPROS = new TGraph(fileLocation+"/data/dataPROS.txt", "%*lg %lg %lg", "");

  TH1D* predDC = tempPredList[0];
  TH1D* predDYB = tempPredList[1];
  TH1D* predRENO = tempPredList[2];
  TH1D* predNEOS = tempPredList[3];
  TH1D* predPROS = tempPredList[4];

  for(Int_t i=0;i<predDC->GetNbinsX();i++)
  {
    dataDC -> SetBinContent(i+1, predDC->GetBinContent(i+1)  * gradataDC->Eval(predDC->GetBinCenter(i+1)) );
    dataDYB -> SetBinContent(i+1,  predDYB->GetBinContent(i+1)  * gradataDYB->Eval(predDYB->GetBinCenter(i+1)) );
    dataRENO -> SetBinContent(i+1,  predRENO->GetBinContent(i+1)  * gradataRENO->Eval(predRENO->GetBinCenter(i+1)) );
    dataNEOS -> SetBinContent(i+1,  predNEOS->GetBinContent(i+1)  * gradataNEOS->Eval(predNEOS->GetBinCenter(i+1)) );
    dataPROS -> SetBinContent(i+1,  predPROS->GetBinContent(i+1)  * gradataPROS->Eval(predPROS->GetBinCenter(i+1)) );
    //std::cout<<"DC23.5 "<<predPROS->GetBinContent(i+1)<<" "<<predPROS->GetBinContent(i+1)  * gradataPROS->Eval(predPROS->GetBinCenter(i+1))<<std::endl;
  }

  std::vector<TH1D*> dataList;
  dataList.push_back(dataDC);
  dataList.push_back(dataDYB);
  dataList.push_back(dataRENO);
  dataList.push_back(dataNEOS);
  dataList.push_back(dataPROS);

  return dataList;

}

//===========================================================================================================================================================8. prepareData ========================

Double_t Sterile ::getPar(int i) {
(((RooAbsReal*)_pulls->at(i))->getVal());
}

RooRealVar* Sterile ::getParVar(int i) {
return ((RooRealVar*)_pulls->at(i));
}


void Sterile :: setSyst(Double_t syst){
_syst = syst;
}

void Sterile :: setdm2CV(Double_t dm2CV){
_dm2CV = dm2CV;
}

void Sterile :: setdm2Unc(Double_t dm2Unc){
_dm2Unc = dm2Unc;
}

void Sterile :: addSK(Bool_t wSK){
withSK = wSK;
}

void Sterile :: setAtmBaseline(Double_t AtmBaseline){
_AtmBaseline = AtmBaseline;
}

void Sterile :: setDensity(Double_t Density){
_Density = Density;
}

void Sterile :: setNBins(Double_t Bins){
_Bins= Bins;
}

void Sterile :: setTime(Double_t time){
_time= time;
}

void Sterile :: setPull(TH1D* pullvecCV){
pullCV = new TVectorD(11);
for(Int_t i=0;i<11;i++){
(*pullCV)[i] =  pullvecCV->GetBinContent(i+1);
    }
}

void Sterile :: setPullUnc(TH1D* pullvecUnc){
pullUnc = new TVectorD(11);
for(Int_t i=0;i<11;i++){
(*pullUnc)[i] = pullvecUnc->GetBinContent(i+1);
    }
}

Double_t Sterile::getPullUnc(Int_t pN){
return (*pullUnc)[pN];
}

void Sterile::DataSwitch(Bool_t dataSwitch) const
{
Bool_t _dataSwitch = dataSwitch;
}

Bool_t Sterile::getDataSwitch() const
{
return _dataSwitch;
}

RooListProxy* Sterile::getPullList() const
{
return _pulls;
}

void Sterile::SetBinning(TH1D* binHist)
{
for(Int_t i=0;i< binHist->GetNbinsX(); i++)
{
binEdge[i] = binHist->GetBinContent(i+1);
}
_nBins = binHist->GetNbinsX()-1;
}

void Sterile::SetFissionFraction(TH1D* fissionHist)
{
for(Int_t i=0; i< fissionHist->GetNbinsX();i++)
fissionFraction[i] = fissionHist->GetBinContent(i+1);
}

void Sterile::SetMatrixNameDC(TString matrixName)
{
fileNameDC = matrixName;
}

void Sterile::SetMatrixNameDYB(TString matrixName)
{
fileNameDYB = matrixName;
}

void Sterile::SetMatrixNameNEOS(TString matrixName)
{
fileNameNEOS = matrixName;
}

void Sterile::SetMatrixNamePROS(TString matrixName)
{
fileNamePROS = matrixName;
}


void Sterile::SetMatrixNameRENO(TString matrixName)
{
fileNameRENO = matrixName;
}

void Sterile::SetModelList(std::vector<TString> mlist)
{
modelList = mlist;
}

std::vector<TH1D*> Sterile:: GetCurrentPrediction()
{
return this->preparePrediction(this->getPullList(), true);
}

std::vector<TH1D*> Sterile:: GetCurrentData(std::vector<TH1D*> predAList)
{
return this->prepareData(predAList);
}

void Sterile:: fitSingleExp(TString input)
{
singleExp = input;
}

void Sterile::setBaselineDC(Double_t bl)
{
baselineDC = bl;
}
void Sterile::setBaselineDYB(Double_t bl)
{ 
baselineDYB = bl;
}
void Sterile::setBaselineRENO(Double_t bl)
{ 
baselineRENO = bl;
}
void Sterile::setBaselineNEOS(Double_t bl)
{ 
baselineNEOS = bl;
}
void Sterile::setBaselinePROS(Double_t bl)
{
baselinePROS = bl;
}

void Sterile::ifEqualIso(bool iso)
{
equalIso = iso;
}

bool Sterile::GetEqualIso()
{
return equalIso;
}

void Sterile::setFileLocation(TString fileL)
{
fileLocation = fileL;
}

void Sterile::setSysts(bool syst)
{
inSyst = syst;
}

bool Sterile::GetSysts()
{
return inSyst;
}

void Sterile::ifEShiftHist(bool eshifthist)
{
ifEHist = eshifthist;
}

bool Sterile::getIfEShiftHist()
{
return ifEHist ;
}

void Sterile::setEShiftHist(TString file)
{
TFile eshift(file);
std::cout<<"taking escale file "<<file<<std::endl;
TH1D* histEscale = (TH1D*)eshift.Get("Escalefraction");
vecEscale = new TVectorD(histEscale->GetNbinsX());
for(int i=0;i<histEscale->GetNbinsX();i++)
  (*vecEscale)[i] = histEscale->GetBinContent(i+1);
}

TVectorD* Sterile::getTestVec(){
return testVec;
}

//////////////////////////////////////////////////////////////////////////////////////////////==========================================9. making short function in Sterile variable============

TMatrixD* Sterile:: ConversionMatrix(TString inputFile, TString inputTree)
{
TFile f(inputFile);
TTree* t = (TTree*)f.Get(inputTree);
//   double binEdge[100];
//   Int_t  _nBins;

fHist = new TH2D("","",_nBins,binEdge,_nBins,binEdge);
TH2D* cfHist = new TH2D("","",_nBins,0.5,0.5+_nBins*0.25,_nBins,0.5,0.5+_nBins*0.25);
//TH2D* fHist = new TH2D("","",34,0.5,9,34,0.5,9);

fMatrix = new TMatrixD(_nBins, _nBins);
TMatrixD* cfMatrix = new TMatrixD(_nBins, _nBins);
uMatrix = new TMatrixD(_nBins, _nBins);
unfoldingMatrix = new TMatrixD(_nBins, _nBins);

double nu, prompt;
t->SetBranchAddress("myNeutrinoEnergy_Th",&nu);
t->SetBranchAddress("myPromptEvisID",&prompt);

for(Int_t i=0;i<t->GetEntries();i++){
  t->GetEntry(i);
  fHist->Fill(nu, prompt);
  cfHist->Fill(nu, prompt);
  //std::cout<<" true vs. prompt "<<nu<<" "<<prompt<<std::endl;
}

for(Int_t i=0;i<fHist->GetNbinsX();i++){
  double summ = 0;
  for(Int_t j=0;j<fHist->GetNbinsY();j++){
     summ += fHist->GetBinContent(i+1,j+1);
  }
  for(Int_t j=0;j<fHist->GetNbinsY();j++){
     if(summ>0) fHist->SetBinContent(i+1,j+1,fHist->GetBinContent(i+1,j+1)/summ);
     (*fMatrix)(j,i) = fHist->GetBinContent(i+1,j+1);
     (*cfMatrix)(j,i) = fHist->GetBinContent(i+1,j+1);
     (*uMatrix)(j,i) = fHist->GetBinContent(i+1,j+1);
  }
}

for(Int_t i=0;i<fHist->GetNbinsY();i++){
  double summ = 0;
  for(Int_t j=0;j<fHist->GetNbinsX();j++){
     summ += fHist->GetBinContent(j+1,i+1);
  }
  for(Int_t j=0;j<fHist->GetNbinsX();j++){
     if(summ>0) fHist->SetBinContent(j+1,i+1,fHist->GetBinContent(j+1,i+1)/summ);
     (*unfoldingMatrix)(i,j) = fHist->GetBinContent(j+1,i+1);
  }
}


t->Delete();
f.Close();
return cfMatrix;
}

TH1D* Sterile:: folding(TH1D* input) const{
TH1D* output(input);
for(Int_t i=0;i<uMatrix->GetNrows();i++){	
  double sum = 0;
  for(Int_t j=0;j<input->GetNbinsX();j++){
    sum += (*uMatrix)(i,j)*input->GetBinContent(j+1);
  }	  
  output->SetBinContent(i+1,sum);
}
//std::cout<<"folded "<<std::endl;
return output;
}
