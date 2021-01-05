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

Lownu ::Lownu (const char* name) 
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
    _pulls->add(_parlist);

    this->addServerList(*_pulls);

    dataDC = new TH1D("","dataDC",30, 0.5, 8.);
};

Lownu ::~Lownu ()
{;}

//=================================================================================================================================
TMatrixD* Lownu::prepareCovMatrix(Int_t nBins, TVectorD* fVec) const
{

    //TFile fMatrixDC(fileNameDC);

    TMatrixD* outMat = new TMatrixD( nBins , nBins);

    if(inSyst){
        // 3% error, diagonal only
        double fLownuErr;
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
	(*outMat)(i,i) += (*fVec)[i] ;    
        if((*outMat)(i,i) == 0) (*outMat)(i,i) += 0.0000000001;
    }

    return outMat ;
}

Double_t Lownu ::FillEv( RooListProxy* _pulls ) const 
{

    //std::cout<<"in FillEv() "<<std::endl;
    int nBins = _nBins;

    std::vector<TH1D*> tempPredList = this->preparePrediction(_pulls, true);
    //std::vector<std::vector<float>> tempPredList = this->preparePrediction(_pulls, true);
    std::cout<<"filled in new pediction "<<std::endl;

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
        //(*fVec)[i]               = predDC->GetBinContent(i+1) ; //* this->surv_Prob( (thisE) , _pulls, 400)  ;
        tempVec[0] -> SetBinContent(i+1, (*fVec)[i]);
    }
    for(Int_t i=0;i<nBins;i++){
        (*fData)[i]             = dataDC->GetBinContent(i+1)   ;
        tempDat[0] -> SetBinContent(i+1, (*fData)[i]);
    }
    std::cout<<"data ready also "<<std::endl;
    //--------------------------------------------------------------------------------------------------------------------------------Data------------------------------------

    // scale to same total rate, doing shape only analysis 
    double scaling1 = tempDat[0]->Integral() / tempVec[0]->Integral();

    for(Int_t i=0;i<nBins;i++){ 
        (*fData)[i]             = TMath::Abs( dataDC->GetBinContent(i+1)   - (*fVec)[i]  );
        //std::cout<<dataDC->GetBinContent(i+1)<<" "<<dataDYB->GetBinContent(i+1)<<" "<<dataRENO->GetBinContent(i+1)<<" "<<dataNEOS->GetBinContent(i+1)<<" "<<dataPROS->GetBinContent(i+1)<<std::endl;
    }

    TMatrixD* covMat = this->prepareCovMatrix(nBins , fVec);
    covMat->Invert();

    TVectorD mulVec(*fData);
    mulVec *= (*covMat);

    Double_t currentResult = TMath::Abs(mulVec*(*fData));
    std::cout<<"DC12_chi2 sans pull "<<currentResult<<std::endl;

    return (Double_t) currentResult ; 
}

//================================================================================================================================================5. Fill the Ev, Prediction===============
Double_t Lownu ::evaluate() const
{ 

    Double_t matPart = this->FillEv(_pulls);//original FillEv is matPart

    Double_t extraPull = this -> ExtraPull (_pulls);//same variable extraPull
    Double_t tot = matPart + extraPull; //If needed, add pull terms here.

    return tot;

}

Double_t Lownu ::ExtraPull (RooListProxy* _pulls) const
{
    Double_t pullAdd = 0;
    for(Int_t i=0;i<1;i++){
        pullAdd += TMath::Power(( ((RooAbsReal*)_pulls->at(i))->getVal() - (*pullCV)[i] ),2) / TMath::Power( (*pullUnc)[i],2) ;
    }
    std::cout<<"extra pull penalty: "<<pullAdd<<std::endl;
    return pullAdd;
}

std::vector<TH1D*> Lownu:: preparePrediction(RooListProxy* _pulls, bool Iosc) const
//std::vector<std::vector<float>> Lownu:: preparePrediction(RooListProxy* _pulls, bool Iosc) const
{

    //((RooAbsReal*)_pulls->at(i+12))->getVal()	;
    std::cout<<"updating preparePrediction() .."<<std::endl;

    //std::vector<std::vector<float>> predictionList;
    //std::vector<float> predNuE;
    TFile file2(fileLocation + "flux_shifts.root");
    TH1D* ND_numubar_RHC = (TH1D*)file2.Get("syst0/ND_numubar_RHC");
    int numBins = ND_numubar_RHC->GetNbinsX();
    double minimum = ND_numubar_RHC->GetBinLowEdge(1);
    double maximum = ND_numubar_RHC->GetBinLowEdge(numBins) + ND_numubar_RHC->GetBinWidth(numBins);
    std::vector<TH1D*> predictionList;
    predictionList.clear();
    TH1D* predNuE = new TH1D("","",numBins,minimum,maximum);

    for (int i = 0; i < this->inputTree->GetEntries(); i++) {
        int temp = 0;
        this->inputTree->GetEntry(i);
        std::cout << this->recoNuE << std::endl;
        for (int ii = 1; ii < ND_numubar_RHC->GetNbinsX()+1; ii++)
        {
            if (ND_numubar_RHC->GetBinLowEdge(ii) + ND_numubar_RHC->GetBinWidth(ii) > trueNuE && ND_numubar_RHC->GetBinLowEdge(ii) < trueNuE) {
                temp = ii;
            }
        }
        //predNuE.push_back(recoNuE, 1 + ((RooAbsReal*)_pulls->at(0))->getVal() * ND_numubar_RHC->GetBinContent(temp+1));
        predNuE->Fill(recoNuE, 1 + ((RooAbsReal*)_pulls->at(0))->getVal() * ND_numubar_RHC->GetBinContent(temp+1));
    }
    predictionList.push_back(predNuE);

    std::cout<<"return preprePrediction()"<<std::endl;

    //std::cout<<"************* before folding "<<predPROS->Integral()<<std::endl;
    //TH1D* temp(predPROS);
    //TH1D* fpredPROS = this->folding(temp);
    // seems pushing back the prompt energy spectrum
    return predictionList;

}

void Lownu::SetInputTree(TString fileName)
{
    file = std::make_unique<TFile> (fileName);
    this->inputTree = (TTree*)file.get()->Get("tree");
    if (!inputTree)
        throw std::runtime_error("asdasdad");
    this->inputTree->SetBranchAddress("recoNeutrinoE", &recoNuE);
    this->inputTree->SetBranchAddress("trueNeutrinoE", &trueNuE);
    std::cout << __func__ << std::endl;
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////==============================7.5 print prediction with folding===


std::vector<TH1D*> Lownu:: prepareData(std::vector<TH1D*> tempPredList) const// 
{
    return tempPredList;
}

//===========================================================================================================================================================8. prepareData ========================

Double_t Lownu ::getPar(int i) {
    (((RooAbsReal*)_pulls->at(i))->getVal());
}

RooRealVar* Lownu ::getParVar(int i) {
    return ((RooRealVar*)_pulls->at(i));
}


void Lownu :: setSyst(Double_t syst){
    _syst = syst;
}

void Lownu :: setdm2CV(Double_t dm2CV){
    _dm2CV = dm2CV;
}

void Lownu :: setdm2Unc(Double_t dm2Unc){
    _dm2Unc = dm2Unc;
}

void Lownu :: addSK(Bool_t wSK){
    withSK = wSK;
}

void Lownu :: setAtmBaseline(Double_t AtmBaseline){
    _AtmBaseline = AtmBaseline;
}

void Lownu :: setDensity(Double_t Density){
    _Density = Density;
}

void Lownu :: setNBins(Double_t Bins){
    _Bins= Bins;
}

void Lownu :: setTime(Double_t time){
    _time= time;
}

void Lownu :: setPull(TH1D* pullvecCV){
    pullCV = new TVectorD(11);
    for(Int_t i=0;i<11;i++){
        (*pullCV)[i] =  pullvecCV->GetBinContent(i+1);
    }
}

void Lownu :: setPullUnc(TH1D* pullvecUnc){
    pullUnc = new TVectorD(11);
    for(Int_t i=0;i<11;i++){
        (*pullUnc)[i] = pullvecUnc->GetBinContent(i+1);
    }
}

Double_t Lownu::getPullUnc(Int_t pN){
    return (*pullUnc)[pN];
}

void Lownu::DataSwitch(Bool_t dataSwitch) const
{
    Bool_t _dataSwitch = dataSwitch;
}

Bool_t Lownu::getDataSwitch() const
{
    return _dataSwitch;
}

RooListProxy* Lownu::getPullList() const
{
    return _pulls;
}

void Lownu::SetBinning(TH1D* binHist)
{
    for(Int_t i=0;i< binHist->GetNbinsX(); i++)
    {
        binEdge[i] = binHist->GetBinContent(i+1);
    }
    _nBins = binHist->GetNbinsX()-1;
}

void Lownu::SetFissionFraction(TH1D* fissionHist)
{
    for(Int_t i=0; i< fissionHist->GetNbinsX();i++)
        fissionFraction[i] = fissionHist->GetBinContent(i+1);
}

void Lownu::SetMatrixNameDC(TString matrixName)
{
    fileNameDC = matrixName;
}

void Lownu::SetMatrixNameDYB(TString matrixName)
{
    fileNameDYB = matrixName;
}

void Lownu::SetMatrixNameNEOS(TString matrixName)
{
    fileNameNEOS = matrixName;
}

void Lownu::SetMatrixNamePROS(TString matrixName)
{
    fileNamePROS = matrixName;
}


void Lownu::SetMatrixNameRENO(TString matrixName)
{
    fileNameRENO = matrixName;
}

void Lownu::SetModelList(std::vector<TString> mlist)
{
    modelList = mlist;
}

//std::vector<std::vector<float>> Lownu:: GetCurrentPrediction()
std::vector<TH1D*> Lownu:: GetCurrentPrediction()
{
    return this->preparePrediction(this->getPullList(), true);
}

std::vector<TH1D*> Lownu:: GetCurrentData(std::vector<TH1D*> predAList)
{
    return this->prepareData(predAList);
}

void Lownu:: fitSingleExp(TString input)
{
    singleExp = input;
}

void Lownu::setBaselineDC(Double_t bl)
{
    baselineDC = bl;
}
void Lownu::setBaselineDYB(Double_t bl)
{ 
    baselineDYB = bl;
}
void Lownu::setBaselineRENO(Double_t bl)
{ 
    baselineRENO = bl;
}
void Lownu::setBaselineNEOS(Double_t bl)
{ 
    baselineNEOS = bl;
}
void Lownu::setBaselinePROS(Double_t bl)
{
    baselinePROS = bl;
}

void Lownu::ifEqualIso(bool iso)
{
    equalIso = iso;
}

bool Lownu::GetEqualIso()
{
    return equalIso;
}

void Lownu::setFileLocation(TString fileL)
{
    fileLocation = fileL;
}

void Lownu::setSysts(bool syst)
{
    inSyst = syst;
}

bool Lownu::GetSysts()
{
    return inSyst;
}

void Lownu::ifEShiftHist(bool eshifthist)
{
    ifEHist = eshifthist;
}

bool Lownu::getIfEShiftHist()
{
    return ifEHist ;
}

void Lownu::SetInputName(TString fileName)
{
    ffileName = fileName;
}

void Lownu::setEShiftHist(TString file)
{
    TFile eshift(file);
    std::cout<<"taking escale file "<<file<<std::endl;
    TH1D* histEscale = (TH1D*)eshift.Get("Escalefraction");
    vecEscale = new TVectorD(histEscale->GetNbinsX());
    for(int i=0;i<histEscale->GetNbinsX();i++)
        (*vecEscale)[i] = histEscale->GetBinContent(i+1);
}

TVectorD* Lownu::getTestVec(){
    return testVec;
}

//////////////////////////////////////////////////////////////////////////////////////////////==========================================9. making short function in Lownu variable============

TMatrixD* Lownu:: ConversionMatrix(TString inputFile, TString inputTree)
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

TH1D* Lownu:: folding(TH1D* input) const{
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
