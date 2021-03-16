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
const int kNumberOfEnergyBin = 16;

Lownu::Lownu (const char* name, double inError) 
    : RooAbsReal(name,name)
    , mLownuCrossSectionError(inError)
{
    _pulls = new RooListProxy("_pulls","_pulls",this);

    RooRealVar* Par[kNumberOfEnergyBin + 19 + 1];
    for (int i = 0; i < kNumberOfEnergyBin; i++) {
        Par[i] = new RooRealVar(Form("energy bin %d", i), Form("par%d", i+1), 0, -100, 100);
        Par[i]->setConstant(false);
        _parlist.add(*(Par[i]));
    }
    for (int i = kNumberOfEnergyBin; i < kNumberOfEnergyBin + 19; i++) {
        Par[i] = new RooRealVar(Form("cov %d", i - 16), Form("par%d", i+1), 0, -100, 100);
        Par[i]->setConstant(false);
        _parlist.add(*(Par[i]));
    }
    for (int i = kNumberOfEnergyBin + 19; i < kNumberOfEnergyBin + 19 + 1; i++) {
        Par[i] = new RooRealVar("background", Form("par%d", i+1), 0, -100, 100);
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

std::vector<TH1D> Lownu::preparePrediction(RooListProxy* _pulls, bool Iosc) const
{
    std::vector<TH1D> predictionList;
    predictionList.clear();

    TH1D predNuE("", "", kNumberOfEnergyBin, 0, 8);

    for (int event = 0; event < this->inputTree->GetEntries(); event++) {
        inputTree->GetEntry(event);
        int temp = 0;
        int temp2 = 0;
        double weight = 0;
        double weight2 = 0;
        //signal
        //{
        if (trackNum == 1) {
            if ((category == 0 || category == 1) && recoNeutronKE < this->mNuCut) {// && isSecondary != 1) {
                for (int tempBin = 0; tempBin < kNumberOfEnergyBin; tempBin++) {
                    //find corresponding bin, true neutrinoE
                    if ((tempBin + 1) * 0.5 > trueNuE/1000. && tempBin * 0.5 < trueNuE/1000.) {
                        temp = tempBin;
                    }
                }
                //0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 7, 8, 12, 16, 20, 40, 116
                for (int tempBin = 0; tempBin < 11; tempBin++) {
                    //find corresponding bin, true neutrinoE
                    if ((tempBin + 1) * 0.5 > trueNuE/1000. && tempBin * 0.5 < trueNuE/1000.) {
                        temp2 = tempBin;
                    }
                }
                if (7 > trueNuE/1000. && 6 < trueNuE/1000.) {
                    temp2 = 13;
                }
                if (8 > trueNuE/1000. && 7 < trueNuE/1000.) {
                    temp2 = 14;
                }
                if (12 > trueNuE/1000. && 8 < trueNuE/1000.) {
                    temp2 = 15;
                }
                if (16 > trueNuE/1000. && 12 < trueNuE/1000.) {
                    temp2 = 16;
                }
                if (20 > trueNuE/1000. && 16 < trueNuE/1000.) {
                    temp2 = 17;
                }
                if (40 > trueNuE/1000. && 20 < trueNuE/1000.) {
                    temp2 = 18;
                }
                if (116 > trueNuE/1000. && 40 < trueNuE/1000.) {
                    temp2 = 19;
                }
                weight = ((RooAbsReal*)_pulls->at(temp))->getVal();
                weight2 = ((RooAbsReal*)_pulls->at(kNumberOfEnergyBin + temp2))->getVal();
                predNuE.Fill(recoNuE/1000., (1 + weight) * (1 + weight2));
            }
            //}

            //bkg
            //{
            if ((category == 2 || category == 3) && recoNeutronKE < this->mNuCut) {
                //use 100% (1) error for background
                predNuE.Fill(recoNuE/1000., (1 + ((RooAbsReal*)_pulls->at(kNumberOfEnergyBin + 19))->getVal()) * 1);
            }
            //}
        }
    }
    predictionList.push_back(predNuE);

    return predictionList;
}

Double_t Lownu::FillEv(RooListProxy* _pulls) const 
{
    std::vector<TH1D> tempPredList = this->preparePrediction(_pulls, true);

    TVectorD* pred = new TVectorD(kNumberOfEnergyBin);
    std::vector<TVectorD> tempPrediction;
    TVectorD* fData = new TVectorD(kNumberOfEnergyBin);
    std::vector<TVectorD> tempData;
    TVectorD* difference = new TVectorD(kNumberOfEnergyBin);

    for (Int_t i = 0; i < kNumberOfEnergyBin; i++) {	 
        (*pred)[i] = tempPredList[0].GetBinContent(i+1) ; //* this->surv_Prob( (thisE) , _pulls, 400)  ;
    }
    for (Int_t i = 0; i < kNumberOfEnergyBin; i++) {
        (*fData)[i] = this->mData->GetBinContent(i+1)   ;
    }

    //data - prediciton
    //kSacling : anti neutrino CC 0pi event for 3DST per year, CDR
    //only stat
    for (Int_t i = 0; i < kNumberOfEnergyBin; i++) { 
        (*difference)[i] = kScaling/this->mData->GetEntries() * ((*pred)[i] - (*fData)[i]);
    }

    TMatrixD* covMat = this->prepareCovMatrix(kNumberOfEnergyBin, pred);
    covMat->Invert();

    TVectorD mulVec(*difference);
    mulVec *= (*covMat);

    Double_t currentResult = TMath::Abs(mulVec*(*difference));

    //std::cout << "DC12_chi2 sans pull " << currentResult << std::endl;

    return (Double_t) currentResult; 
}
TMatrixD* Lownu::prepareCovMatrix(Int_t nBins, TVectorD* pred) const
{
    //output covariant matrix
    TMatrixD* outMat = new TMatrixD(nBins , nBins);

    //only fill diagonal element
    //(i, i) = statistic
    for(Int_t i = 0; i < nBins ; i++) {
        (*outMat)(i, i) = (*pred)[i] * kScaling / this->mData->GetEntries();    
        if((*outMat)(i, i) == 0) 
            (*outMat)(i, i) += 0.0000000001;
    }

    return outMat ;
}

Double_t Lownu::FillEv2(RooListProxy* _pulls) const 
{
    TVectorD* e_i = new TVectorD(kNumberOfEnergyBin);
    TVectorD* centralValue = new TVectorD(kNumberOfEnergyBin);
    TVectorD* difference = new TVectorD(kNumberOfEnergyBin);

    for (Int_t i = 0; i < kNumberOfEnergyBin; i++) {	 
        (*e_i)[i] = ((RooAbsReal*)_pulls->at(i))->getVal();
    }
    for (Int_t i = 0; i < kNumberOfEnergyBin; i++) {
        (*centralValue)[i] = 0;
    }

    //e_i - centralValue
    for (Int_t i = 0; i < kNumberOfEnergyBin; i++) { 
        (*difference)[i] = (*e_i)[i] - (*centralValue)[i];
    }

    TMatrixD* covMat = this->prepareCovMatrix2(kNumberOfEnergyBin);
    covMat->Invert();

    TVectorD mulVec(*difference);
    mulVec *= (*covMat);

    Double_t currentResult = TMath::Abs(mulVec*(*difference));

    return (Double_t) currentResult; 
}
TMatrixD* Lownu::prepareCovMatrix2(Int_t nBins) const
{
    //output covariant matrix
    TMatrixD* outMat = new TMatrixD(nBins , nBins);

    //(i, i) = cross section uncertainty^2
    //(i, j) = correlation
    for(Int_t i = 0; i < nBins ; i++) {
        (*outMat)(i, i) = std::pow(mLownuCrossSectionError, 2);
    }
    for(Int_t i = 0; i < nBins ; i++) {
        for(Int_t j = 0; j < nBins ; j++) {
            if (i == j) {
                continue;
            }
            (*outMat)(i, j) = kCorelation * std::pow((*outMat)(i,i), 0.5) * std::pow((*outMat)(j,j), 0.5);
        }
    }

    return outMat ;
}

Double_t Lownu::FillEv3(RooListProxy* _pulls) const 
{
    TVectorD* f_i = new TVectorD(19);
    TVectorD* centralValue = new TVectorD(19);
    TVectorD* difference = new TVectorD(19);

    for (Int_t i = 0; i < 19; i++) {	 
        (*f_i)[i] = ((RooAbsReal*)_pulls->at(kNumberOfEnergyBin + i))->getVal();
    }
    for (Int_t i = 0; i < 19; i++) {
        (*centralValue)[i] = 0;
    }

    //f_i - centralValue
    for (Int_t i = 0; i < 19; i++) { 
        (*difference)[i] = (*f_i)[i] - (*centralValue)[i];
    }

    TMatrixD* covMat = this->covMatrix;
    covMat->Invert();

    TVectorD mulVec(*difference);
    mulVec *= (*covMat);

    Double_t currentResult = TMath::Abs(mulVec*(*difference));

    return (Double_t) currentResult; 
}

Double_t Lownu::evaluate() const
{ 
    Double_t part1 = this->FillEv(_pulls);  //P-D
    Double_t part2 = this->FillEv2(_pulls); //e-CV
    Double_t part3 = this->FillEv3(_pulls); //f-CV
    std::cout << "p - d: " << part1 << std::endl;
    std::cout << "e - CV: " << part2 << std::endl;
    std::cout << "f - CV: " << part3 << std::endl;
    Double_t tot = part1 + part2 + part3;

    return tot;
}

void Lownu::SetInputTree(TString fileName)
{
    dataFile = std::make_unique<TFile> (fileName);  
    this->inputTree = (TTree*)dataFile.get()->Get("tree");
    if (!inputTree)
        throw std::runtime_error("invalid input");
    this->inputTree->SetBranchAddress("recoNeutrinoE", &recoNuE);
    this->inputTree->SetBranchAddress("trueNeutrinoE", &trueNuE);
    this->inputTree->SetBranchAddress("category", &category);
    this->inputTree->SetBranchAddress("recoNeutronKE", &recoNeutronKE);
    this->inputTree->SetBranchAddress("trackNum", &trackNum);

    this->mData = new TH1D("","data", kNumberOfEnergyBin, 0, 8);
    for (int i = 0; i < inputTree->GetEntries(); ++i) {
        inputTree->GetEntry(i);
        if (trackNum == 1 && recoNeutronKE < 200) {
            this->mData->Fill(recoNuE/1000.);
        }
    }
        
}
void Lownu::SetInputSyst(TString fileName)
{
    covRootFile = std::make_unique<TFile> (fileName);
    TMatrixD* RHC_cov_nom_mat = (TMatrixD*)covRootFile.get()->Get("RHC_cov_nom_mat");
    TMatrixD* RHC_cov_shp_mat = (TMatrixD*)covRootFile.get()->Get("RHC_cov_shp_mat");
    TMatrixD* input_RHC_chol_diag_nom_mat = (TMatrixD*)covRootFile.get()->Get("RHC_chol_diag_nom_mat");
    this->covMatrix = new TMatrixD(19, 19);
    this->RHC_chol_diag_nom_mat = new TMatrixD(19, 19);
    for (int i = 0; i < 19; ++i) {
        for (int j = 0; j < 19; ++j) {
            (*(this->covMatrix))(i, j) = (*RHC_cov_nom_mat)(i + 19, j + 19) + (*RHC_cov_shp_mat)(i + 19, j + 19);
            (*(this->RHC_chol_diag_nom_mat))(i,i) = (*input_RHC_chol_diag_nom_mat)(i + 19, i + 19);
        }
    }
}
