/* 
   We will try to extract fEvtHdr info.
*/

#include <iostream>
#include "TChain.h"

#include "../../src/SetROOTVar.cpp"

void fEvtHdr (const char *rootfile_w_path) {
  
  TChain *T = new TChain("T");
  T->Add(rootfile_w_path);

  UInt_t runnum = 0, evnum = 0;
  // T->SetMakeClass(1); // This is important for reasons unknown!!
  // T->SetBranchAddress("fEvtHdr.fRun", &fEvtHdr_fRun); 
  setrootvar::setbranch(T, "fEvtHdr", "fRun", &runnum);
  setrootvar::setbranch(T, "fEvtHdr", "fEvtNum", &evnum);

  T->GetEntry(10);

  std::cout << " Run # " << evnum << std::endl;
}
