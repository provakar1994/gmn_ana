/*
  Script to learn-try-apply RDataFrame stuff!
*/

#include <iostream>

#include "TChain.h"

void rdf() {

  // TChain C("Tout");
  // C.Add("siout/qelas_ana_simu_sbs14_sbs70p_model2.root");

  // ROOT::RDataFrame df(C);

  ROOT::RDataFrame df("Tout", "siout/qelas_ana_simu_sbs14_sbs70p_model2*");
  
  // auto df_var = df
  //   .Define("W2","");

  // auto colNames = df.GetColumnNames();
  // for (auto &&colName : colNames) {
  //   std::cout << colName << std::endl;
  // }

  auto df_filtered = df.Filter("fiduCut");

  // TCanvas *c1 = new TCanvas("c1","c1",600,400);
  // c1->cd();
  auto h = df.Histo1D({"histName", "histTitle", 200, -3.5, 2.5}, "dx", "weight");
  h->DrawClone();

  auto h2 = df_filtered.Histo1D({"histName2", "histTitle", 200, -3.5, 2.5}, "dx", "weight");
  h2->SetLineColor(kRed);
  h2->DrawClone("same");

  auto d = df.Display("dx");
  d->Print();

}



