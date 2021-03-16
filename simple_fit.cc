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

//kSacling : anti neutrino CC 0pi event for 3DST per year, CDR
float kScaling = 2.4E+6 * 0.25; //1year
const float kCorelation = 0.00;
TVectorD energyScaleError(16);

Lownu::Lownu (const char* name, int numPars, double inError) 
    : RooAbsReal(name,name)
    , mNumberOfParameters(numPars)
    , mError(inError)
{
    _pulls = new RooListProxy("_pulls","_pulls",this);

    //Par[numPars+1] is for bkg
    //E' = a + b*E_reco
    //Par[numPars+2] : a
    //Par[numPars+3] : b
    //Par[numPars+ 4~20] : each energy bin
    
    RooRealVar* Par[numPars+3+16];
    for (int i = 0; i < numPars + 3; i++) {
        if (i == numPars)
            Par[i] = new RooRealVar("background", Form("par%d", i+1), 0, -100, 100);
        else if (i == numPars + 1)
            Par[i] = new RooRealVar("a", Form("par%d", i+2), 0, -100, 100);
        else if (i == numPars + 2)
            Par[i] = new RooRealVar("b", Form("par%d", i+3), 0, -100, 100);
        else
            Par[i] = new RooRealVar(Form("flux systematic %d", i), Form("par%d", i+1), 0, -100, 100);
        Par[i]->setConstant(false);
        _parlist.add(*(Par[i]));
    }
    for (int i = numPars + 3; i < numPars + 3 + 16; i++) {
        Par[i] = new RooRealVar(Form("energy bin %d", i - numPars - 2), Form("par%d", i), 0, -100, 100);
        Par[i]->setConstant(false);
        _parlist.add(*(Par[i]));
    }
    _pulls->add(_parlist);
    this->addServerList(*_pulls);
    mNuCut = 200;
};

Lownu::~Lownu()
{
    delete this->_pulls;
}

//=================================================================================================================================
TMatrixD* Lownu::prepareCovMatrix(Int_t nBins, TVectorD* pred) const
{
    //output covariant matrix
    TMatrixD* outMat = new TMatrixD(nBins , nBins);

    //cross section uncertainty
    //{
    TVectorD uncertainties(_nBins);
    for (int i = 0; i < _nBins; i++) {
        uncertainties[i] = mError;
    }
    TVectorD recoUncertainties(_nBins);
    recoUncertainties = *(this->mFolding)*uncertainties;
    //}

    //only fill diagonal element
    //(i, i) = statistic + cross section uncertainty^2
    //{
    for(Int_t i = 0; i < nBins ; i++) {
        (*outMat)(i, i) = (*pred)[i] * kScaling / this->mData->GetEntries();    
        if((*outMat)(i, i) == 0) 
            (*outMat)(i, i) += 0.0000000001;
    }
    for(Int_t i = 0; i < nBins ; i++) {
        for(Int_t j = 0; j < nBins ; j++) {
            if (i == j) {
                continue;
            }
            //(*outMat)(i, j) = kCorelation * std::pow((*outMat)(i,i), 0.5) * std::pow((*outMat)(j,j), 0.5);
        }
    }
    //}

    return outMat ;
}

TMatrixD* Lownu::prepareCovMatrix2(Int_t nBins) const
{
    //output covariant matrix
    TMatrixD* outMat = new TMatrixD(nBins , nBins);

    //only fill diagonal element
    //(i, i) = cross section uncertainty^2
    //{
    for(Int_t i = 0; i < nBins ; i++) {
        (*outMat)(i, i) = std::pow(mError, 2);
    }
    for(Int_t i = 0; i < nBins ; i++) {
        for(Int_t j = 0; j < nBins ; j++) {
            if (i == j) {
                continue;
            }
            (*outMat)(i, j) = kCorelation * std::pow((*outMat)(i,i), 0.5) * std::pow((*outMat)(j,j), 0.5);
        }
    }
    //}

    return outMat ;
}

Double_t Lownu::FillEv(RooListProxy* _pulls) const 
{
    std::vector<TH1D> tempPredList = this->preparePrediction(_pulls, true);

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
    //kSacling : anti neutrino CC 0pi event for 3DST per year, CDR
    for (Int_t i = 0; i < _nBins; i++) { 
        (*difference)[i] = kScaling/this->mData->GetEntries() * ((*pred)[i] - (*fData)[i]);
    }

    TMatrixD* covMat = this->prepareCovMatrix(_nBins, pred);
    covMat->Invert();

    TVectorD mulVec(*difference);
    mulVec *= (*covMat);

    Double_t currentResult = TMath::Abs(mulVec*(*difference));

    //std::cout << "DC12_chi2 sans pull " << currentResult << std::endl;

    return (Double_t) currentResult; 
}

Double_t Lownu::FillEv2(RooListProxy* _pulls) const 
{
    TVectorD* e_i = new TVectorD(_nBins);
    TVectorD* centralValue = new TVectorD(_nBins);
    TVectorD* difference = new TVectorD(_nBins);

    for (Int_t i = 0; i < _nBins; i++) {	 
        (*e_i)[i] = ((RooAbsReal*)_pulls->at(this->GetNumberOfParameters()+2+i+1))->getVal();
    }
    for (Int_t i = 0; i < _nBins; i++) {
        (*centralValue)[i] = 0;
    }

    //e_i - centralValue
    for (Int_t i = 0; i < _nBins; i++) { 
        (*difference)[i] = (*e_i)[i] - (*centralValue)[i];
    }

    TMatrixD* covMat = this->prepareCovMatrix2(_nBins);
    covMat->Invert();

    TVectorD mulVec(*difference);
    mulVec *= (*covMat);

    Double_t currentResult = TMath::Abs(mulVec*(*difference));

    //std::cout << "DC12_chi2 sans pull " << currentResult << std::endl;

    return (Double_t) currentResult; 
}

//================================================================================================================================================5. Fill the Ev, Prediction===============
Double_t Lownu::evaluate() const
{ 
    Double_t matPart = this->FillEv(_pulls);//original FillEv is matPart
    Double_t energyPart = this->FillEv2(_pulls);//(e_i - CV)^T * ( ) * (e_i - CV)
    Double_t extraPull = this->ExtraPull(_pulls);//same variable extraPull
    std::cout << "p - d side: " << matPart << std::endl;
    std::cout << "e - CV side: " << energyPart << std::endl;
    std::cout << "extraPull: " << extraPull << std::endl;
    Double_t tot = matPart + energyPart + extraPull; //If needed, add pull terms here.

    return tot;
}

Double_t Lownu::ExtraPull(RooListProxy* _pulls) const
{
    Double_t pullAdd = 0;
    for(Int_t i = 0; i < this->GetNumberOfParameters() + 1; i++) {
        pullAdd += TMath::Power((((RooAbsReal*)_pulls->at(i))->getVal() - (*pullCV)[i]), 2)/TMath::Power((*pullUnc)[i], 2) ;
    }
    pullAdd += TMath::Power((((RooAbsReal*)_pulls->at(this->GetNumberOfParameters() + 1))->getVal()), 2) / TMath::Power(1, 2);
    pullAdd += TMath::Power((((RooAbsReal*)_pulls->at(this->GetNumberOfParameters() + 2))->getVal()), 2) / TMath::Power(1, 2);
    for(Int_t i = this->GetNumberOfParameters() + 3; i < this->GetNumberOfParameters() + 3 + 16; i++) {
        pullAdd += TMath::Power(((RooAbsReal*)_pulls->at(i))->getVal() - 0, 2)/TMath::Power(1, 2);
    }
    return pullAdd;
}

std::vector<TH1D> Lownu::preparePrediction(RooListProxy* _pulls, bool Iosc) const
{
    std::vector<TH1D> predictionList;
    predictionList.clear();

    //using the same number of bins, flux systematic
    //{
    int numBins = this->syst[0].GetNbinsX();
    double minimum = this->syst[0].GetBinLowEdge(1);
    double maximum = this->syst[0].GetBinLowEdge(numBins) + this->syst[0].GetBinWidth(numBins);
    TH1D predNuE("", "", numBins, minimum, maximum);
    //}

    for (int event = 0; event < this->inputTree->GetEntries(); event++) {
        inputTree->GetEntry(event);
        int temp = 0;
        double weight = 1;
        double weightEnergyBin = 0;
        double a = 0;
        double b = 0;
        //signal
        //{
        if (trackNum == 1) {
            if ((category == 0 || category == 1) && recoNeutronKE < this->mNuCut) {// && isSecondary != 1) {
                for (int tempPar = 0; tempPar < this->GetNumberOfParameters(); tempPar++) {
                    for (int tempBin = 1; tempBin < this->syst[tempPar].GetNbinsX() + 1; tempBin++) {
                        //find corresponding bin, true neutrinoE
                        if (this->syst[tempPar].GetBinLowEdge(tempBin) + this->syst[tempPar].GetBinWidth(tempBin) > trueNuE/1000. 
                                && this->syst[tempPar].GetBinLowEdge(tempBin) < trueNuE/1000.) {
                            temp = tempBin;
                        }
                    }
                    weight *= 1 + ((RooAbsReal*)_pulls->at(tempPar))->getVal() * this->syst[tempPar].GetBinContent(temp);
                }
                //weight *= ((RooAbsReal*)_pulls->at(this->GetNumberOfParameters()+1))->getVal() * this->eScale->GetBinContent(temp);
                //weight *= ((RooAbsReal*)_pulls->at(this->GetNumberOfParameters()+2))->getVal() * this->smearing->GetBinContent(temp);
                weightEnergyBin = ((RooAbsReal*)_pulls->at(this->GetNumberOfParameters()+2+temp))->getVal();
                predNuE.Fill(recoNuE/1000., (weight) * (1 + weightEnergyBin));
            }
            //}

            //bkg
            //{
            if ((category == 2 || category == 3) && recoNeutronKE < this->mNuCut) {
                //use 100% (1) error for background
                predNuE.Fill(recoNuE/1000., (1 + ((RooAbsReal*)_pulls->at(this->GetNumberOfParameters()))->getVal()) * 1);
                //predNuE.Fill(recoNuE/1000);
            }
            //}
        }
    }
    predictionList.push_back(predNuE);

    return predictionList;
}

std::vector<TH1D> Lownu::TEST(RooListProxy* _pulls)
{
    std::vector<TH1D> predictionList;
    predictionList.clear();

    //using the same number of bins, flux systematic
    //{
    int numBins = this->syst[0].GetNbinsX();
    double minimum = this->syst[0].GetBinLowEdge(1);
    double maximum = this->syst[0].GetBinLowEdge(numBins) + this->syst[0].GetBinWidth(numBins);
    TH1D predNuE("", "", numBins, minimum, maximum);
    //}
    for (int bin = 0; bin < numBins; ++bin) {
        double tempBinValue = 0;
        double tempBinDinominator = 0;
        for (int ii = 0; ii < this->GetNumberOfParameters(); ++ii) {
            tempBinValue += std::pow((this->syst[ii].GetBinContent(bin+1) * ((RooRealVar*)_pulls->at(ii))->getError()), 2);
            tempBinDinominator += std::pow(this->syst[ii].GetBinContent(bin+1), 2);
        }
        predNuE.SetBinContent(bin+1, std::pow(tempBinValue, 0.5)/std::pow(tempBinDinominator, 0.5));
    }
    //predNuE.Add(&syst[0], ((RooRealVar*)_pulls->at(0))->getError());
    predictionList.push_back(predNuE);

    return predictionList;
}

void Lownu::SetInputTree(TString fileName)
{
    file = std::make_unique<TFile> (fileName);
    this->inputTree = (TTree*)file.get()->Get("tree");
    if (!inputTree)
        throw std::runtime_error("invalid input");
    this->inputTree->SetBranchAddress("recoNeutrinoE", &recoNuE);
    this->inputTree->SetBranchAddress("trueNeutrinoE", &trueNuE);
    this->inputTree->SetBranchAddress("category", &category);
    this->inputTree->SetBranchAddress("recoNeutronKE", &recoNeutronKE);
    this->inputTree->SetBranchAddress("trackNum", &trackNum);
    //this->inputTree->SetBranchAddress("isSecondary", &isSecondary);

    this->mData = new TH1D("mData", "mData", 16, 0, 8);
    this->mFolding = std::make_unique<TMatrixD> (_nBins, _nBins);
    this->eScale = new TH1D("","energy scale",16,0,8);
    eScale->SetBinContent(1,-0.1608301);
    eScale->SetBinContent(2,-0.1064572);
    eScale->SetBinContent(3,-0.08198026);
    eScale->SetBinContent(4,-0.06008165);
    eScale->SetBinContent(5,-0.02177312);
    eScale->SetBinContent(6,0.01526615);
    eScale->SetBinContent(7,0.08380262);
    eScale->SetBinContent(8,0.15233);
    eScale->SetBinContent(9,0.1639289);
    eScale->SetBinContent(10,0.213944);
    eScale->SetBinContent(11,0.1481013);
    eScale->SetBinContent(12,0.09745763);
    eScale->SetBinContent(13,0.1262887);
    eScale->SetBinContent(14,0.03533569);
    eScale->SetBinContent(15,0.015625);
    eScale->SetBinContent(16,-0.004651163);
    this->smearing = new TH1D("unnamed","smearing",16,0,8);
    smearing->SetBinContent(1,-0.05306868);
    smearing->SetBinContent(2,-0.01210618);
    smearing->SetBinContent(3,-0.004864977);
    smearing->SetBinContent(4,-0.00413987);
    smearing->SetBinContent(5,-0.002425308);
    smearing->SetBinContent(6,-0.002289431);
    smearing->SetBinContent(7,-0.004670003);
    smearing->SetBinContent(8,-0.008615162);
    smearing->SetBinContent(9,-0.01163406);
    smearing->SetBinContent(10,-0.03235897);
    smearing->SetBinContent(11,-0.05162117);
    smearing->SetBinContent(12,-0.1178922);
    smearing->SetBinContent(13,-0.1923525);
    smearing->SetBinContent(14,-0.1187281);
    smearing->SetBinContent(15,-0.2341523);
    smearing->SetBinContent(16,-0.06467317);

    //for each bin
    double numberOfTrueNuE[_nBins];
    double numberOfRecoNuE[_nBins][_nBins];
    for (int i = 0; i < _nBins; i++) {
        numberOfTrueNuE[i] = 0;
        for (int j = 0; j < _nBins; j++) {
            numberOfRecoNuE[i][j] = 0;
        }
    }

    for (int i = 0; i < inputTree->GetEntries(); i++) {
        inputTree->GetEntry(i);
        if (trackNum == 1) {
            if (recoNeutronKE < this->mNuCut) {
                if ((category == 0 || category == 1)) {
                    this->mData->Fill(recoNuE/1000.);
                }
                if (category == 2 || category == 3) {
                    this->mData->Fill(recoNuE/1000.);
                }

                for (int j = 0; j < _nBins; j++) {
                    if (0.5*j < trueNuE/1000. && trueNuE/1000. < 0.5*(j+1)) {
                        numberOfTrueNuE[j]++;
                        for (int k = 0; k < _nBins; k++) {
                            if (0.5*k < recoNuE/1000. && recoNuE/1000. < 0.5*(k+1)) {
                                numberOfRecoNuE[j][k]++;
                            }
                        }
                    }
                }
            }
        }
    }
    for(Int_t i = 0; i < _nBins; i++) {
        for(Int_t j = 0; j < _nBins; j++) {
            (*mFolding)(i, j) = (numberOfRecoNuE[j][i]/numberOfTrueNuE[j]);
        }
    }
}

void Lownu::SetInputSyst(TString fileName)
{
    fileFluxShifts = std::make_unique<TFile> (fileName);
    TH1D temp_ND_numubar_RHC[this->GetNumberOfParameters()];
    for (int i = 0; i < this->GetNumberOfParameters(); i++)
    {
        TH1D* temp = (TH1D*)fileFluxShifts.get()->Get(Form("syst%d/ND_numubar_RHC",i));
        temp_ND_numubar_RHC[i] = *temp;
        this->syst.push_back(temp_ND_numubar_RHC[i]);
    }
}

void Lownu::SetNumberOfParameters(int inNum)
{
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
