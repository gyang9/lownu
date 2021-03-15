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

std::vector<TH1D> histVec;
void loop(int numPar, double Error);

using namespace std;

int main(int argc, char**argv) {

    int numPar = 0;
    ////if (argc < 2) {
    ////    std::cout << "./app {number of parameters}" << std::endl;
    ////    return 0;
    ////}
    ////numPar = std::stoi(argv[1]);

    std::cout << "how many parameters? (1~256) : ";
    std::cin >> numPar;
    //while (numPar < 1 || numPar > 256) {
    //    std::cout << "invalied number of parameters, put again (1~256) : ";
    //    std::cin >> numPar;
    //}
    //double Error = 0;
    //std::cout << "how much error? : ";
    //std::cin >> Error;

    //RooFitResult* res;
    //std::unique_ptr<Lownu> rep = std::make_unique<Lownu> ("_rep", numPar, Error);
    //char formula[10];

    //std::cout << "start to run " << std::endl;

    //// vecInput1 is the CV for pulls while vecInput2 is the unc. for pulls
    //auto vecInput1 = std::make_unique<TH1D> ("", "", rep->GetNumberOfParameters() + 1, 0, rep->GetNumberOfParameters() + 1);
    //auto vecInput2 = std::make_unique<TH1D> ("", "", rep->GetNumberOfParameters() + 1, 0, rep->GetNumberOfParameters() + 1);

    //for (int i = 0; i < rep->GetNumberOfParameters() + 1; i++) {
    //    vecInput1->SetBinContent(i+1, 0);
    //    vecInput2->SetBinContent(i+1, 1);
    //}

    //std::cout << "you've set some inputs " << std::endl;

    //int nBins = 16;
    //std::unique_ptr<TH1D> binHist = std::make_unique<TH1D> ("", "", nBins+1 , 0, nBins + 1);
    //for(Int_t i = 0; i < nBins + 1; i++) {
    //    binHist->SetBinContent(i+1, 0.5 + 0.25 * i);
    //}

    //TString fileLocation = "./";
    //rep->setFileLocation(fileLocation);

    //rep->SetBinning(binHist.get());

    ////rep->SetInputName(fileLocation + "variableOutputAfterCut.root");
    //rep->SetInputTree(fileLocation + "variableOutputAfterCut.root");
    //rep->SetInputSyst(fileLocation + "flux_shifts.root");
    //gStyle->SetOptStat(false);
    //TCanvas can10;
    //rep->mFolding->Draw("colz text");
    //can10.SaveAs("folding.pdf");

    //std::vector<TH1D> tempPredList = rep->preparePrediction(rep->getPullList(), false);
    //std::cout << "tempPredList.size(): " << tempPredList.size() << std::endl;

    //TCanvas can4;
    //tempPredList[0].Draw();
    //rep->mData->SetLineColor(2);
    //rep->mData->Draw("same");
    //can4.SaveAs("beforeFit.png");

    //std::vector<TH1D> pred1 = rep->preparePrediction(rep->getPullList(), false);
    //TCanvas can8;
    //pred1[0].Draw();
    //rep->mData->Draw("same");
    //can8.SaveAs("1sigmashift.pdf");

    //rep->setPull(vecInput1.get()); 
    //rep->setPullUnc(vecInput2.get());

    //std::cout << "ended up with setting us basic stuff " << std::endl;
    //rep->FillEv(rep->getPullList());

    //RooArgList list("list");
    //list.add(*rep);
    //sprintf(formula, "%s", "@0");
    ////RooFormulaVar* fcn = new RooFormulaVar("fit","fit",formula,list);
    //std::unique_ptr<RooFormulaVar> fcn = std::make_unique<RooFormulaVar> ("fit", "fit", formula, list);

    //// ******************************** Important setup here *************************************
    //// *******************************************************************************************
    //rep->setSysts(false);
    ////rep->ifEShiftHist(true);
    ////rep->setEShiftHist("Escalefraction.root");
    //// *******************************************************************************************
    //// ******************************************************************************************* 

    //std::cout << "------------  Getting current spectra " << std::endl;
    //std::vector<TH1D> outPrediction = rep->GetCurrentPrediction(); //tempPredList; // rep->GetCurrentPrediction();
    //std::vector<TH1D> outData = rep->GetCurrentData(outPrediction);
    //std::cout << "------------  Have got current spectra " << std::endl;

    //// save the un-oscillated standard prediction and data spectra
    ////TFile* outputFile = new TFile("outputFigs.root","RECREATE");
    //std::unique_ptr<TFile> outputFile = std::make_unique<TFile> ("outputFigs.root","RECREATE");
    //std::cout << "outPrediction size: " << outPrediction.size() << std::endl;
    //for(Int_t i = 0; i < outPrediction.size(); i++) {
    //    outPrediction[i].Write(Form("outPrediction[%d]", i));
    //    if(i < 5) 
    //        outData[i].Write(Form("outData[%d]", i));
    //}
    //outputFile->Close();
    //std::cout << "first thing saved " << std::endl;

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    //RooMinuit m(*fcn);
    //m.setStrategy(2);
    //Double_t callsEDM[2] = {10500., 1.e-6};
    //Int_t irf = 0;

    //gMinuit->mnexcm("MIGRAD", callsEDM, 2, irf);
    //m.migrad();
    //res = m.save();
    //double bestFit = res->minNll();
    //std::cout << "fit status code is : " <<  res->status() << std::endl;	
    ////status = 0    : OK
    ////status = 1    : Covariance was mad  epos defined
    ////status = 2    : Hesse is invalid
    ////status = 3    : Edm is above max
    ////status = 4    : Reached call limit
    ////status = 5    : Any other failure
    //std::cout << "quiality code of covariance matrix is : " <<  res->covQual() << std::endl;
    ////status = -1 :  not available (inversion failed or Hesse failed)
    ////status =  0 : available but not positive defined
    ////status =  1 : covariance only approximate
    ////status =  2 : full matrix but forced pos def
    ////status =  3 : full accurate matrix
    ////////////////////////	
    ////inline Int_t status() const 
    ////inline Int_t covQual() const 
    ////inline Int_t numInvalidNLL() const 
    ////inline Double_t edm() const 
    ////inline Double_t minNll() const 
    ///////////////////////

    //outPrediction = rep->GetCurrentPrediction();
    ////outData = rep->GetCurrentData(outPrediction);

    ////std::cout << "list of reactor pulls : " << std::endl;
    ////for(Int_t i=18;i<41;i++){std::cout << " " << rep->getPar(i) << std::endl;}

    //std::cout << "size of output prediction list " << outPrediction.size() << std::endl;
    //std::cout << "size of output data list " << outData.size() << std::endl;
    ////std::cout << atoi(argv[1]) << " " << atoi(argv[2]) << " " << bestFit << std::endl;

    ////std::cout << "result list: " << std::endl;
    ////std::cout << "chi2: " << bestFit  << std::endl;

    //TCanvas can;
    //rep->mData->SetLineColor(2);
    //rep->mData->Draw();
    //can.SaveAs("data.png");

    //std::vector<TH1D> pred = rep->preparePrediction(rep->getPullList(), false);

    //TCanvas can1;
    //pred[0].Draw();
    //can1.SaveAs("pred.png");

    //TCanvas can2;
    //pred[0].Draw();
    //rep->mData->Draw("same");
    //can2.SaveAs("afterfit_together.png");
    //for (int i = 0; i < rep->GetNumberOfParameters() + 3; ++i) {
    //    std::cout << Form("((RooRealVar*)rep->_pulls->at(%d))->getError(): ",i) << ((RooRealVar*)rep->_pulls->at(i))->getError() << std::endl;
    //}

    //for(int i = 0; i < 3; ++i) {
        //loop(numPar, 0.20 * (i));
    //}

      gROOT->SetStyle("Plain");
  gStyle->SetTitleBorderSize(0);
  gStyle->SetOptStat("");

  gStyle->SetLabelFont(102,"");
  gStyle->SetLabelSize(0.06,"");
  gStyle->SetLabelFont(102,"xyz");
  gStyle->SetLabelSize(0.04,"xyz");
  gStyle->SetLabelOffset(0.001,"x");
  gStyle->SetLabelOffset(0.01,"y");

  gStyle->SetTitleFont(42,"xyz");
  gStyle->SetTitleFontSize(0.06);
  gStyle->SetTitleFont(42);
  gStyle->SetTitleSize(0.06,"xyz");
  gStyle->SetTitleOffset(0.9,"x");
  gStyle->SetTitleOffset(1.2,"y");

  gStyle->SetStripDecimals(kFALSE);

  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetPadBottomMargin(0.15);

  gStyle->SetStatW(0.35);
  gStyle->SetStatH(0.25);

  gStyle->SetPadTickX(kTRUE);
  gStyle->SetPadTickY(kTRUE);

  gStyle->SetPalette(1);
  gStyle->SetNumberContours(99);

  gStyle->SetHistLineWidth(2);
  gStyle->SetFrameLineWidth(2);
  gStyle->SetFuncWidth(2);

  gStyle->SetStatFont(42);
  gStyle->SetOptFit(1111);
  gStyle->SetOptStat(0);
  gStyle->SetOptStat(000000);

    loop(numPar, 0.);
    loop(numPar, 0.1);
    TLegend l;

    TCanvas can5;
    histVec.at(0).GetXaxis()->SetTitle("reco E_{#nu} (GeV)");
    histVec.at(0).GetYaxis()->SetTitle("post fit / pre fit");
    histVec.at(0).SetMinimum(0);
    histVec.at(0).SetMaximum(1);
    histVec.at(0).SetLineColor(1);
    histVec.at(0).Draw();
    l.AddEntry(&histVec.at(0), "0%");
    histVec.at(0).SetLineColor(2);
    histVec.at(0).Draw("same");
    l.AddEntry(&histVec.at(1), "10%");
    histVec.at(1).SetLineColor(6);
    histVec.at(1).Draw("same");
    //for(int i = 0; i < 2; ++i) {
        //l.AddEntry(&histVec.at(i), Form("%d",5*i));
        //histVec.at(i).SetLineColor(kViolet+i);
        //histVec.at(i).Draw("same");
    //}
    l.Draw();
    can5.SaveAs("3_lownu_uncertainty5.pdf");
    can5.SaveAs("3_lownu_uncertainty5.C");
}

void loop(int numPar, double Error)
{
    RooFitResult* res;
    std::unique_ptr<Lownu> rep = std::make_unique<Lownu> ("_rep", numPar, Error);
    char formula[10];

    // vecInput1 is the CV for pulls while vecInput2 is the unc. for pulls
    auto vecInput1 = std::make_unique<TH1D> ("", "", rep->GetNumberOfParameters() + 1, 0, rep->GetNumberOfParameters() + 1);
    auto vecInput2 = std::make_unique<TH1D> ("", "", rep->GetNumberOfParameters() + 1, 0, rep->GetNumberOfParameters() + 1);

    for (int i = 0; i < rep->GetNumberOfParameters() + 1; i++) {
        vecInput1->SetBinContent(i+1, 0);
        vecInput2->SetBinContent(i+1, 1);
    }

    int nBins = 16;
    std::unique_ptr<TH1D> binHist = std::make_unique<TH1D> ("", "", nBins+1 , 0, nBins + 1);
    for(Int_t i = 0; i < nBins + 1; i++) {
        binHist->SetBinContent(i+1, 0.5 + 0.25 * i);
    }

    TString fileLocation = "./";
    rep->setFileLocation(fileLocation);

    rep->SetBinning(binHist.get());

    rep->SetInputName(fileLocation + "variableOutputAfterCut.root");
    rep->SetInputTree(fileLocation + "variableOutputAfterCut.root");
    rep->SetInputSyst(fileLocation + "flux_shifts.root");

    std::vector<TH1D> tempPredList = rep->preparePrediction(rep->getPullList(), false);
    std::cout << "tempPredList.size(): " << tempPredList.size() << std::endl;

    rep->setPull(vecInput1.get()); 
    rep->setPullUnc(vecInput2.get());
    rep->FillEv(rep->getPullList());

    RooArgList list("list");
    list.add(*rep);
    sprintf(formula, "%s", "@0");
    std::unique_ptr<RooFormulaVar> fcn = std::make_unique<RooFormulaVar> ("fit", "fit", formula, list);
    rep->setSysts(false);
    std::unique_ptr<TFile> outputFile = std::make_unique<TFile> ("outputFigs.root","RECREATE");

    RooMinuit m(*fcn);
    m.setStrategy(2);
    Double_t callsEDM[2] = {10500., 1.e-6};
    Int_t irf = 0;

    gMinuit->mnexcm("MIGRAD", callsEDM, 2, irf);
    m.migrad();
    res = m.save();
    double bestFit = res->minNll();

    std::vector<TH1D> asdf = rep->TEST(rep->getPullList());
    asdf[0].SetStats(false);
    histVec.push_back(asdf[0]);

}
