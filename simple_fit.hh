/*
 *  Lownu for lownu fit header file.
 *
 *  Author: Guang Yang
 */

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <iomanip>
#include <ctype.h>
#include <string.h>
#include <math.h>
#include <sstream>
#include <TList.h>

#include <TROOT.h>
#include <TCanvas.h>
#include <TApplication.h>
#include <TH1.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TFile.h>
#include <TRint.h>
#include <TH2.h>
#include <TFormula.h>
#include <TF1.h>
#include <TF2.h>
#include <TMath.h>
#include <Math/DistFunc.h>
#include <TLine.h>
#include <TTree.h>
#include <TGraph.h>
#include <TRandom.h>
#include <TRandom3.h>
#include <TGraphErrors.h>
#include <TVirtualFFT.h>
#include <TFoamIntegrand.h>
#include <TMatrixD.h>
#include <TVectorT.h>
#include <TDecompChol.h>

#include <RooFit.h>
#include "RooGlobalFunc.h"
#include "RooAbsReal.h"
#include "RooAbsCategory.h"
#include "RooRealVar.h"
#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooDataSet.h"
#include "RooArgList.h"
#include "RooArgSet.h"
#include "RooGaussian.h"
#include "RooProdPdf.h"
#include "RooWorkspace.h"
#include "RooMinuit.h"
#include "RooNLLVar.h"
#include "RooFitResult.h"
#include "RooPlot.h"
#include "RooDataSet.h"
#include "RooExtendPdf.h"
#include "RooChi2Var.h"
#include "RooMinuit.h"
#include "RooDataHist.h"
#include "RooHistPdf.h"
#include "RooRandom.h"
#include <RooMsgService.h>
#include <RooHist.h>
#include <RooTrace.h>
#include <RooCategory.h>
#include "RooConstVar.h"
#include "RooBinning.h"

#include "TStopwatch.h"
#include "TFile.h"
#include "TMinuit.h"

#include "RooFit.h"
#include "RooMinuit.h"
#include "RooFitResult.h"
#include "TMinuit.h"
#include <RooRealVar.h>

using namespace RooFit;


  class Lownu : public RooAbsReal {

  public:

    Lownu (const char* name);

    Lownu (const Lownu & other, const char* name = 0): RooAbsReal(other,name) {};
    virtual TObject* clone(const char* newname) const {return new Lownu (*this, newname);};
    virtual ~Lownu () ;

    Lownu (const Lownu & Lownu );

    void randomGo(int newIsotope, RooListProxy* _pulls);

    Lownu & operator=(const Lownu & rhs);

    RooFormulaVar* Chi2() ;

    Double_t FillEv(RooListProxy* _pulls) const;

    Double_t ExtraPull(RooListProxy* _pulls) const;

    TMatrixD* prepareMatrix(TH2D* conv) const;
 
    void setSyst(Double_t syst) ;
    void setAtmBaseline(Double_t AtmBaseline) ;
    void setDensity(Double_t Density) ;

    void setBaselineDYB(Double_t bl) ;
    void setBaselineDC(Double_t bl) ;
    void setBaselineRENO(Double_t bl) ;
    void setBaselineNEOS(Double_t bl) ;
    void setBaselinePROS(Double_t bl) ;

    void setNBins(Double_t Bins) ;
    void setTime(Double_t time) ;
    void setCorr1(Double_t Corr1) ;
    void setCorr2(Double_t Corr2) ;
    void setCorr3(Double_t Corr3) ;
    void setCorr4(Double_t Corr4) ;

    void setdm2CV(Double_t dm2CV) ;
    void setdm2Unc(Double_t dm2Unc) ;

    void setData(RooListProxy* _pulls) const;
    void setSysts(bool inSyst);

    void setPull(TH1D* pullvecCV) ;
    void setPullUnc(TH1D* pullvecUnc) ;
    Double_t getPullUnc(Int_t pN) ;
    void addSK(Bool_t wSK) ;

    void SetBinning(TH1D* binHist) ;
    void SetFissionFraction(TH1D* fissionHist) ;
    void SetMatrixNameDYB(TString matrixNameDYB);
    void SetMatrixNameDC(TString matrixNameDC);
    void SetMatrixNameRENO(TString matrixNameRENO);
    void SetMatrixNameNEOS(TString matrixNameNEOS);
    void SetMatrixNamePROS(TString matrixNamePROS);
    void SetInputName(TString fileName);

    void SetModelList(std::vector<TString>);

    TMatrixD* prepareCovMatrix(TH2D* conv, TVectorD* vec, Double_t syst) const;

    std::vector<TH1D*> prepareData(std::vector<TH1D*> tempPredLis) const;
    //std::vector<TH1D*> preparePrediction(RooListProxy* _pulls) const;
    //std::vector<std::vector<float>> preparePrediction(RooListProxy* _pulls, bool Iosc ) const;
    std::vector<TH1D*> preparePrediction(RooListProxy* _pulls, bool Iosc ) const;
    std::vector<double> wonseokPreparePrediction;//(RooListProxy* _pulls, bool Iosc ) const;//May,20,2020
    //std::vector<TH1D*> preparePrediction(RooListProxy* _pulls, std::vector<TString> modelList) const;

    //std::vector<double> histForScale;//May,20,2020
    //double scalingForNEOS;//May,20,2020

    TMatrixD* prepareT2kCovMatrix(TMatrixD* covM_t2k, TVectorD* fVec_t2k) const;
    TMatrixD* prepareT2kCovMatrix(TMatrixD* covM_t2k, TVectorD* fVec_t2k, Int_t nBins) const;
    TMatrixD* prepareT2kCovMatrixND(TMatrixD* covM_t2k, TVectorD* fVec_t2k, TVectorD* fVec_t2k_nd, Int_t nBins) const;
    TMatrixD* prepareSkCovMatrix(TMatrixD* covM_sk, TVectorD* fVec_sk) const;
    TMatrixD* prepareJunoCovMatrix(Int_t nBins, TMatrixD* covM_JUNO, TVectorD* fVec_JUNO) const;

    TMatrixD* prepareCovMatrix(Int_t nBins, TVectorD* fVec) const;

    RooRealVar Par1 ;
    RooRealVar Par2 ;
    RooRealVar Par3 ;
    RooRealVar Par4 ;
    RooRealVar Par5 ;
    RooRealVar Par6 ;
    RooRealVar Par7 ;
    RooRealVar Par8 ;
    RooRealVar Par9 ;
    RooRealVar Par10 ;
    RooRealVar Par11 ;
    RooRealVar Par12 ;
    RooRealVar Par13 ;
    RooRealVar Par14 ;
    RooRealVar Par15 ;
    RooRealVar Par16 ;
    RooRealVar Par17 ;
    RooRealVar Par18 ;
    RooRealVar Par19 ;
    RooRealVar Par20 ;
    RooRealVar Par21 ;
    RooRealVar Par22 ;
    RooRealVar Par23 ;
    RooRealVar Par24 ;
    RooRealVar Par25 ;
    RooRealVar Par26 ;
    RooRealVar Par27 ;
    RooRealVar Par28 ;
    RooRealVar Par29 ;
    RooRealVar Par30 ;
    RooRealVar Par31 ;
    RooRealVar Par32 ;
    RooRealVar Par33 ;
    RooRealVar Par34 ;
    RooRealVar Par35 ;
    RooRealVar Par36 ;
    RooRealVar Par37 ;
    RooRealVar Par38 ;
    RooRealVar Par39 ;
    RooRealVar Par40 ;
    RooRealVar Par41 ;
    RooRealVar Par42 ;
    RooRealVar Par43 ;
    RooRealVar Par44 ;
    RooRealVar Par45 ;
    RooRealVar Par46 ;
    RooRealVar Par47 ;
    RooRealVar Par48 ;
    RooRealVar Par49 ;
    RooRealVar Par50 ;
    RooRealVar Par51 ;
    RooRealVar Par52 ;
    RooRealVar Par53 ;
    RooRealVar Par54 ;
    RooRealVar Par55 ;
    RooRealVar Par56 ;
    RooRealVar Par57 ;
    RooRealVar Par58 ;
    RooRealVar Par59 ;
    RooRealVar Par60 ;

    Double_t _par1;
    Double_t _par2;
    Double_t _par3;
    Double_t _par4;
    Double_t _par5;
    Double_t _par6;
    Double_t _par7;
    Double_t _par8;
    Double_t _par9;
    Double_t _par10;
    Double_t _par11;
    Double_t _par12;
    Double_t _par13;
    Double_t _par14;
    Double_t _par15;
    Double_t _par16;
    Double_t _par17;
    Double_t _par18;
    Double_t _par19;
    Double_t _par20;
    Double_t _par21;
    Double_t _par22;
    Double_t _par23;
    Double_t _par24;
    Double_t _par25;
    Double_t _par26;
    Double_t _par27;
    Double_t _par28;
    Double_t _par29;
    Double_t _par30;
    Double_t _par31;
    Double_t _par32;
    Double_t _par33;
    Double_t _par34;
    Double_t _par35;
    Double_t _par36;
    Double_t _par37;
    Double_t _par38;
    Double_t _par39;
    Double_t _par40;
    Double_t _par41;
    Double_t _par42;
    Double_t _par43;
    Double_t _par44;
    Double_t _par45;
    Double_t _par46;
    Double_t _par47;
    Double_t _par48;
    Double_t _par49;
    Double_t _par50;
    Double_t _par51;
    Double_t _par52;
    Double_t _par53;
    Double_t _par54;
    Double_t _par55;
    Double_t _par56;
    Double_t _par57;
    Double_t _par58;
    Double_t _par59;
    Double_t _par60;

    Double_t getPar(int i) ;
    void DataSwitch(Bool_t dataSwitch=true) const;
    Bool_t getDataSwitch() const;

    Double_t surv_t2k(Double_t L, Double_t E, RooListProxy* _pulls) const ;
    Double_t app_t2k(Double_t L, Double_t E, Double_t density, RooListProxy* _pulls) const ;
    Double_t survNue_sk(Double_t L, Double_t E, Double_t density, RooListProxy* _pulls) const ;
    Double_t survNumu_sk(Double_t L, Double_t E, RooListProxy* _pulls) const ;

    Double_t surv_t2k(Double_t L, Double_t E, TVectorD* parVec) const ;
    Double_t app_t2k(Double_t L, Double_t E, Double_t density, TVectorD* parVec) const ;
    Double_t survNue_sk(Double_t L, Double_t E, Double_t density, TVectorD* parVec) const ;
    Double_t survNumu_sk(Double_t L, Double_t E, TVectorD* parVec) const ;

    Double_t surv_JUNO(Double_t E, RooListProxy* _pulls) const ;
    Double_t surv_JUNO(Double_t E, TVectorD* parVec) const ;

    Double_t surv_JUNO(Double_t E, RooListProxy* _pulls, Double_t L) const ;
    Double_t surv_JUNO(Double_t E, TVectorD* parVec, Double_t L) const ;

    TF1* fitting_function;//June
    Double_t surv_Prob(Double_t E, RooListProxy* _pulls, Double_t L, TF1* func ) const ;
    //Double_t surv_Prob(Double_t E, RooListProxy* _pulls, Double_t L ) const ;
    //TF1* getFittingFunction();
    //void setFittingFunction(TF1*);

    RooRealVar* getParVar(int i) ;
    RooListProxy* getParVar() ;
    RooListProxy* getPullList() const;

    TVectorD* getSkVec() ;
    TVectorD* getT2kVec() ;
    TVectorD* getJunoVec() ;
    TVectorD* getJunoData() ;

    TH1D* getJunoHist() ;
    TH1D* getJunoDataHist() ;

    RooArgList _parlist;
    RooListProxy* _pulls;

    //std::vector<std::vector<float>> GetCurrentPrediction(); //혹시 이것때문에 출력은 절대 아닐 텐데...
    std::vector<TH1D*> GetCurrentPrediction(); //혹시 이것때문에 출력은 절대 아닐 텐데...
    std::vector<TH1D*> GetCurrentData(std::vector<TH1D*> pred);
    
    TF1* GetIBDXsecFormula() const;
    TGraph* GetIBDXsecPoints() const;
    std::vector<TH1D*> GetFluxPrediction(RooListProxy* _pulls, bool Iosc) ;

    //std::vector<TMatrixD*> ConversionMatrix(TString inputFile, TString inputTree);
    TMatrixD* ConversionMatrix(TString inputFile, TString inputTree);
    //TH2D* GetConversionHist() const {return fHist[0];}//need
    TH2D* GetConversionHist() const {return fHist;}//need
    //std::vector<TMatrixD*> GetConversionMatrix() const {return fMatrix;}//need
    TMatrixD* GetConversionMatrix() const {return fMatrix;}//need
    //TH1D* folding(TH1D* input, int something) const;
    TH1D* folding(TH1D* input) const;
    TVectorD* getTestVec();//May,20,2020
    Double_t getScaling4();//June,3,2020
    void fitSingleExp(TString input);

    void ifEqualIso(bool iso = false);
    bool GetEqualIso();
    bool GetSysts();

    void setFileLocation(TString fileL);

    void ifEShiftHist(bool eshifthist);
    void setEShiftHist(TString file);
    bool getIfEShiftHist();

    //TH2D* fHist[5];
    TH2D* fHist;
    //TH2D* fHist_DC;
    //TH2D* fHist_DYB;
    //TH2D* fHist_RENO;
    //TH2D* fHist_NEOS;
    //TH2D* fHist_PROS;

    //std::vector<TMatrixD*> fMatrix;
    TMatrixD* fMatrix;
    //TMatrixD* fMatrix_DC;
    //TMatrixD* fMatrix_DYB;
    //TMatrixD* fMatrix_RENO;
    //TMatrixD* fMatrix_NEOS;
    //TMatrixD* fMatrix_PROS;

    //std::vector<TMatrixD*> uMatrix;
    TMatrixD* uMatrix;
    //TMatrixD* uMatrix_DC;
    //TMatrixD* uMatrix_DYB;
    //TMatrixD* uMatrix_RENO;
    //TMatrixD* uMatrix_NEOS;
    //TMatrixD* uMatrix_PROS;

    TMatrixD* unfoldingMatrix;

    TH2D* conv;
    TH1D* nueAfter;
    TH1D* nueBefore;
    TH1D* numuTrue;

    Double_t _syst;
    Double_t _dm2CV;
    Double_t _dm2Unc;

    Double_t _Density;
    Double_t _AtmBaseline;
    Double_t _time;
    Double_t _Bins;
    Double_t _Corr1;
    Double_t _Corr2;
    Double_t _Corr3;
    Double_t _Corr4;
    Double_t dataSet;
    Bool_t _dataSwitch;

    TVectorD* pullCV; 
    TVectorD* pullUnc;
    TVectorD* fData_t2k;
    TVectorD* fData_sk;
    TVectorD* fData_JUNO;

    TVectorD* fVec_t2k;
    TVectorD* fVec_t2k_nd;
    TVectorD* fVec_sk;
    TVectorD* fVec_JUNO;
    TVectorD* fVecShadow_JUNO;
 
    TVectorD* MCDelegate;
    TVectorD* dataDelegate;

    TH1D* hData_JUNO;
    TH1D* hShadow_JUNO;
    
    Bool_t withSK ;

    TString fileName;
    TString fileNameRENO;
    TString fileNameNEOS;
    TString fileNameDC;
    TString fileNameDYB;
    TString fileNamePROS;
    TString ffileName;
    double fissionFraction[100];
    double binEdge[100];
    Int_t  _nBins;
    int newBin = 60; //new
    int newBin_DCDYB = 24; //new
    int newBin_RENO = 26; //new
    int newBin_NEOS = 61; //new
    int newBin_PROS = 29; //new
    int newBin_total = newBin_DCDYB + newBin_DCDYB + newBin_RENO + newBin_NEOS + newBin_PROS;//new
    double binWidth = 0.125; //new
    /*
    Double_t prob[10000];//May,30,2020 : 1000 oscillation in surv_Prob function
    double scaling44=123.456;//June,3,2020 : Print out Scaling4
    double predNEOS2[100];//new
    Double_t thisE_NEOS;//new
    */
    TString singleExp;
    double baselineDC;
    double baselineDYB;
    double baselineRENO;
    double baselineNEOS;
    double baselinePROS;

    bool inSyst;
    bool equalIso;
    TString fileLocation;

    bool ifEHist ;
    TString EHistFile ;
    TVectorD* vecEscale;    

    TH1D* dataDC;
    TH1D* dataDYB;
    TH1D* dataRENO;
    TH1D* dataNEOS;
    TH1D* dataPROS;
    //TH1D* predDC;
    //TH1D* predDYB;
    //TH1D* predRENO;
    //TH1D* predNEOS;

    std::vector<TString> modelList;

    std::vector<TH1D*> dataList;
    TVectorD* testVec = new TVectorD(100);//May,21,2020
    //double scaling4=123.456;//June,3,2020
    //  private:
    //
    float recoNuE;
    float trueNuE;
    TTree* inputTree = nullptr;
    std::unique_ptr<TFile> file;
    void SetInputTree(TString fileName);
   
   virtual  Double_t evaluate() const ;

   protected:


  };


