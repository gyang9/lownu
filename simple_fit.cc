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

Lownu::Lownu (const char* name, int numPars) 
    : RooAbsReal(name,name)
    , mNumberOfParameters(numPars)
{
    _pulls = new RooListProxy("_pulls","_pulls",this);

    //Par[numPars+1] is for bkg
    RooRealVar* Par[numPars+1];
    for (int i = 0; i < numPars + 1; i++) {
        if (i == numPars)
            Par[i] = new RooRealVar("background", Form("par%d", i+1), 0, 0, 100);
        else
            Par[i] = new RooRealVar(Form("flux%d", i), Form("par%d", i+1), 0, 0, 100);
        Par[i]->setConstant(false);
        _parlist.add(*(Par[i]));
    }
    _pulls->add(_parlist);
    this->addServerList(*_pulls);
};

Lownu::~Lownu()
{;}

//=================================================================================================================================
TMatrixD* Lownu::prepareCovMatrix(Int_t nBins, TVectorD* pred) const
{
    std::cout << "now: " << __func__ << std::endl;
    //TFile fMatrixDC(fileNameDC);

    TMatrixD* outMat = new TMatrixD(nBins , nBins);

    if (inSyst) {
        // 3% error, diagonal only
        double fLownuErr = 0.03;
        double errlist[100];
        for (int temp = 0; temp < 10; temp++)
            errlist[temp] = fLownuErr;
        TVectorD* errList = new TVectorD(nBins);
        for (Int_t i = 0; i< nBins; i++) {
            (*errList)[i] = errlist[i];
        }
        //-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
        for(Int_t i = 0; i < nBins; i++) {
            for(Int_t j = i; j < nBins; j++) {
                if(i == j)	
                    (*outMat)(i, j) = (*errList)[i] * (*errList)[j] * (*pred)[i] * (*pred)[j];
            }
        }
    }

    for(Int_t i = 0; i < nBins ; i++) {
        (*outMat)(i, i) += (*pred)[i];    
        if((*outMat)(i, i) == 0) 
            (*outMat)(i, i) += 0.0000000001;
    }
    return outMat ;
}

Double_t Lownu::FillEv(RooListProxy* _pulls) const 
{
    //std::cout<<"in FillEv() "<<std::endl;
    std::cout << "now: " << __func__ << std::endl;
    std::vector<TH1D> tempPredList = this->preparePrediction(_pulls, true);
    //std::cout<<"filled in new pediction "<<std::endl;

    TVectorD* pred = new TVectorD(_nBins);
    std::vector<TVectorD> tempPrediction;
    TVectorD* fData = new TVectorD(_nBins);
    std::vector<TVectorD> tempData;
    TVectorD* difference = new TVectorD(_nBins);

    for (Int_t i = 0; i< _nBins; i++) {	 
        (*pred)[i] = tempPredList[0].GetBinContent(i+1) ; //* this->surv_Prob( (thisE) , _pulls, 400)  ;
    }
    for (Int_t i = 0; i < _nBins; i++) {
        (*fData)[i] = this->mData->GetBinContent(i+1)   ;
    }
    //data - prediciton
    for (Int_t i = 0; i < _nBins; i++) { 
        (*difference)[i] = TMath::Abs((*fData)[i] - (*pred)[i]);
    }

    /*
    TMatrixD* covMat = this->prepareCovMatrix(_nBins, pred);
    covMat->Invert();

    TVectorD mulVec(*difference);
    mulVec *= (*covMat);

    Double_t currentResult = TMath::Abs(mulVec*(*difference));
    */
    Double_t currentResult;
    for (Int_t i = 0; i < _nBins; i++) {
        if ((*pred)[i] == 0)
            continue;
        currentResult += TMath::Power((*difference)[i], 2)/(*pred)[i];
    }
    std::cout << "DC12_chi2 sans pull " << currentResult << std::endl;

    return (Double_t) currentResult; 
}

//================================================================================================================================================5. Fill the Ev, Prediction===============
Double_t Lownu::evaluate() const
{ 
    std::cout << "now: " << __func__ << std::endl;
    Double_t matPart = this->FillEv(_pulls);//original FillEv is matPart
    Double_t extraPull = this->ExtraPull(_pulls);//same variable extraPull
    Double_t tot = matPart + extraPull; //If needed, add pull terms here.

    return tot;
}

Double_t Lownu::ExtraPull(RooListProxy* _pulls) const
{
    std::cout << "now: " << __func__ << std::endl;
    Double_t pullAdd = 0;
    for(Int_t i = 0; i < this->GetNumberOfParameters() + 1; i++) {
        pullAdd += TMath::Power((((RooAbsReal*)_pulls->at(i))->getVal() - (*pullCV)[i]), 2)/TMath::Power((*pullUnc)[i], 2) ;
    }
    std::cout << "extra pull penalty: " << pullAdd << std::endl;
    return pullAdd;
}

std::vector<TH1D> Lownu::preparePrediction(RooListProxy* _pulls, bool Iosc) const
{
    std::cout << "now: " << __func__ << std::endl;
    //std::cout << "updating preparePrediction() .." << std::endl;

    std::vector<TH1D> predictionList;
    predictionList.clear();

    int numBins = this->syst[0].GetNbinsX();
    double minimum = this->syst[0].GetBinLowEdge(1);
    double maximum = this->syst[0].GetBinLowEdge(numBins) + this->syst[0].GetBinWidth(numBins);
    TH1D predNuE("", "", numBins, minimum, maximum);
    for (int event = 0; event < this->inputTree->GetEntries(); event++) {
        inputTree->GetEntry(event);
        int temp = 0;
        double weight = 1;
        //signal
        if (category == 1) {
            for (int tempPar = 0; tempPar < this->GetNumberOfParameters(); tempPar++) {

                for (int tempBin = 1; tempBin < this->syst[tempPar].GetNbinsX() + 1; tempBin++) {
                    if (this->syst[tempPar].GetBinLowEdge(tempBin) + this->syst[tempPar].GetBinWidth(tempBin) > trueNuE/1000. 
                        && this->syst[tempPar].GetBinLowEdge(tempBin) < trueNuE/1000.) {
                        temp = tempBin;
                    }
                }
                weight = weight + ((RooAbsReal*)_pulls->at(tempPar))->getVal() * this->syst[tempPar].GetBinContent(temp);
            }
            predNuE.Fill(recoNuE/1000., weight);
        }
        //bkg
        if (category == 3) {
            predNuE.Fill(recoNuE/1000., ((RooAbsReal*)_pulls->at(this->GetNumberOfParameters()))->getVal() * 1);
        }
    }
    predictionList.push_back(predNuE);

    //std::cout << "return preprePrediction()" << std::endl;

    //std::cout<<"************* before folding "<<predPROS->Integral()<<std::endl;
    //TH1D* temp(predPROS);
    //TH1D* fpredPROS = this->folding(temp);
    // seems pushing back the prompt energy spectrum
    return predictionList;
}

void Lownu::SetInputTree(TString fileName)
{
    std::cout << "now: " << __func__ << std::endl;
    file = std::make_unique<TFile> (fileName);
    this->inputTree = (TTree*)file.get()->Get("tree");
    if (!inputTree)
        throw std::runtime_error("asdasdad");
    this->inputTree->SetBranchAddress("recoNeutrinoE", &recoNuE);
    this->inputTree->SetBranchAddress("trueNeutrinoE", &trueNuE);
    this->inputTree->SetBranchAddress("category", &category);

    this->mData = new TH1D("mData", "mData", 16, 0, 8);

    for (int i = 0; i < inputTree->GetEntries(); i++) {
        inputTree->GetEntry(i);
        this->mData->Fill(recoNuE/1000.);
    }
    //std::cout << "now: " << __func__ << std::endl;
}

void Lownu::SetInputSyst(TString fileName)
{
    std::cout << "now: " << __func__ << std::endl;
    flux_shifts = std::make_unique<TFile> (fileName);
    TH1D temp_ND_numubar_RHC[10];
    for (int i = 0; i < this->GetNumberOfParameters(); i++)
    {
        TH1D* temp = (TH1D*)flux_shifts.get()->Get(Form("syst%d/ND_numubar_RHC",i));
        temp_ND_numubar_RHC[i] = *temp;
        this->syst.push_back(temp_ND_numubar_RHC[i]);
    }
    //std::cout << "now: " << __func__ << std::endl;
}

void Lownu::SetNumberOfParameters(int inNum)
{
    std::cout << "now: " << __func__ << std::endl;
    if (inNum < 1)
        throw std::runtime_error("SetNumberOfParameters(): invalid argument");
    this->mNumberOfParameters = inNum;
}

const int Lownu::GetNumberOfParameters() const
{
    return this->mNumberOfParameters;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////==============================7.5 print prediction with folding===


std::vector<TH1D> Lownu::prepareData(std::vector<TH1D> tempPredList) const// 
{
    std::cout << "now: " << __func__ << std::endl;
    return tempPredList;
}

//===========================================================================================================================================================8. prepareData ========================

Double_t Lownu::getPar(int i) {
    return (((RooAbsReal*)_pulls->at(i))->getVal());
}

RooRealVar* Lownu::getParVar(int i) {
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

void Lownu :: setPull(TH1D* pullvecCV) {
    this->pullCV = new TVectorD(this->GetNumberOfParameters() + 1);
    for(Int_t i = 0; i < this->GetNumberOfParameters() + 1; i++) {
        (*this->pullCV)[i] =  pullvecCV->GetBinContent(i+1);
    }
}

void Lownu :: setPullUnc(TH1D* pullvecUnc) {
    this->pullUnc = new TVectorD(this->GetNumberOfParameters() + 1);
    for(Int_t i = 0; i < this->GetNumberOfParameters() + 1; i++) {
        (*this->pullUnc)[i] = pullvecUnc->GetBinContent(i+1);
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
    for(Int_t i = 0; i < binHist->GetNbinsX(); i++)
    {
        binEdge[i] = binHist->GetBinContent(i + 1);
    }
    _nBins = binHist->GetNbinsX() - 1;
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

std::vector<TH1D> Lownu:: GetCurrentPrediction()
{
    return this->preparePrediction(this->getPullList(), true);
}

std::vector<TH1D> Lownu:: GetCurrentData(std::vector<TH1D> predAList)
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
    std::cout << "now: " << __func__ << std::endl;
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
    std::cout << "now: " << __func__ << std::endl;
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
