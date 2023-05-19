#ifndef SIMU_JOB_H
#define SIMU_JOB_H

typedef struct SimuJob {

  std::string sfname;    // g4sbs output filename with path
  std::string rfname;    // Replayed ROOT filename with path
  std::string generator; // generator (SIMC/g4sbs)
  double ngenreq;        // # of simulated events requested
  double ntried;         // # of simulated events tried to get the requested
  double genvol;         // generation volume (sr for g4sbs, MeV*sr2 for SIMC)
  double lumi;           // luminosity (ub^-1 for SIMC & Hz/cm2 for g4sbs) 
  double charge;         // C, total charge per simulation job (rel. for simc)
  double ebeam;          // GeV, uncorrected beam energy (rel. for g4sbs)
  double ibeam;          // A, beam current (rel. for g4sbs)

  // constructor 
SimuJob(): 
  sfname("NONE"),rfname("NONE"),generator("NONE"),ngenreq(0),ntried(0),
    genvol(0),lumi(0),charge(0),ebeam(0),ibeam(0)
  {}

  // define an ostream operator to print to screen conveniently
  friend ostream& operator <<(ostream &out, const SimuJob& sjob) {
    out << " ------------" << std::endl;
    out << " Simu. file name    : " << sjob.sfname << std::endl;
    out << " Replayed file name : " << sjob.rfname << std::endl;
    out << " Generator          : " << sjob.generator << std::endl;
    out << " # events requested : " << sjob.ngenreq << std::endl;
    out << " # events tried     : " << sjob.ntried << std::endl;
    out << " Generation vol.    : " << sjob.genvol << std::endl;
    out << " Luminosity         : " << sjob.lumi << std::endl;
    if (sjob.generator.compare("simc")==0)
      out << " Charge (C)         : " << sjob.charge << std::endl;
    if (sjob.generator.compare("g4sbs")==0) {
      out << " Ebeam (GeV)        : " << sjob.ebeam << std::endl;
      out << " Ibeam (A)          : " << sjob.ibeam << std::endl;
    }
    out << " ------------" << std::endl << std::endl;
    return out;
  }

  // sets data by reading MC job summary (exclusively for util::ReadSimuJobSummary function)
  void SetDataSimuJob(std::vector<std::string> data) {
    sfname    = data[0];
    rfname    = data[1];
    generator = data[2];
    ngenreq   = stod(data[3]);
    ntried    = stod(data[4]);
    genvol    = stod(data[5]);
    lumi      = stod(data[6]);
    if (generator.compare("simc")==0)
      charge  = stod(data[7])/1000.;
    if (generator.compare("g4sbs")==0) {
      ebeam   = stod(data[7]);      
      ibeam   = stod(data[8])*1e-6;      
    }
  }

} SimuJob_t;  

#endif
