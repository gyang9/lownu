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

const double kLownuUncertainty = 0.1;

int main(int argc, char**argv) {

    RooFitResult* res;
    std::unique_ptr<Lownu> rep = std::make_unique<Lownu> ("_rep", kLownuUncertainty);
    char formula[10];

    rep->SetInputTree("variableOutputAfterCut.root");
    rep->SetInputSyst("flux_covariance_matrices.root");

    RooArgList list("list");
    list.add(*rep);
    sprintf(formula, "%s", "@0");
    std::unique_ptr<RooFormulaVar> fcn = std::make_unique<RooFormulaVar> ("fit", "fit", formula, list);

    RooMinuit m(*fcn);
    m.setStrategy(2);
    Double_t callsEDM[2] = {10500., 1.e-6};
    Int_t irf = 0;

    gMinuit->mnexcm("MIGRAD", callsEDM, 2, irf);
    m.migrad();
    res = m.save();
    double bestFit = res->minNll();

    TH1D* test = new TH1D("","fi/RHC_chol_diag_nom_mat",19,0,19);
    for (int i = 0; i < 19; ++i) {
        test->SetBinContent(i+1, ((RooRealVar*)rep->_pulls->at(16 + i))->getError() / (*(rep->RHC_chol_diag_nom_mat))(i,i));
    }

    TCanvas can;
    test->Draw();
    can.SaveAs("test.pdf");
}
