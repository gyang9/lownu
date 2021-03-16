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

        Lownu (const char* name, double inError);

        Lownu (const Lownu & other, const char* name = 0)
            : RooAbsReal(other,name) {};
        virtual TObject* clone(const char* newname) const {return new Lownu (*this, newname);};
        virtual ~Lownu () ;
        virtual  Double_t evaluate() const ;

        Lownu (const Lownu& Lownu);

        Lownu & operator=(const Lownu & rhs);
        void SetInputTree(TString fileName);
        void SetInputSyst(TString fileName);

        Double_t FillEv(RooListProxy* _pulls) const;
        Double_t FillEv2(RooListProxy* _pulls) const;
        Double_t FillEv3(RooListProxy* _pulls) const;
        TMatrixD* prepareCovMatrix(Int_t nBins, TVectorD* fVec) const;
        TMatrixD* prepareCovMatrix2(Int_t nBins) const;
        TMatrixD* prepareCovMatrix3() const;

        std::vector<TH1D> preparePrediction(RooListProxy* _pulls, bool Iosc ) const;

        RooArgList _parlist;
        RooListProxy* _pulls;

        float mNuCut;
        float recoNuE;
        float trueNuE;
        float recoNeutronKE;
        int category;
        int trackNum;

        TTree* inputTree = nullptr;
        std::unique_ptr<TFile> dataFile;
        std::unique_ptr<TFile> covRootFile;

        TH1D* mData;
        double mLownuCrossSectionError = 0;

        TMatrixD* covMatrix;
        TMatrixD* RHC_chol_diag_nom_mat;
};

