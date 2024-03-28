#ifndef SIMU_JOB_H
#define SIMU_JOB_H

typedef struct SimuJob {

  std::string sfname{0};    // g4sbs output filename with path
  std::string rfname{0};    // Replayed ROOT filename with path
  std::string generator{0}; // generator (simc,g4sbs,etc.)
  std::string process{0};   // (heep,deep,deen,deeN,inel,etc.)
  double ngenreq{0};        // # of simulated events requested
  double ntried{0};         // # of simulated events tried to get the requested
  double genvol{0};         // generation volume (sr for g4sbs, MeV*sr2 for SIMC)
  double lumi{0};           // luminosity (ub^-1 for SIMC & Hz/cm2 for g4sbs) 
  double ebeam{0};          // GeV, uncorrected beam energy (rel. for g4sbs)
  double charge{0};         // C, total charge per simulation job (rel. for simc)
  double ibeam{0};          // A, beam current (rel. for g4sbs)
  bool usingRS{0};          // flag indicating status of rejection sampling (RS) (rel. for simc)
  double maxwtRS{0};        // ubMeV^-1sr^-2, maximum weight for RS (rel. for simc)

  // constructor 
  SimuJob () {}

  // sets data by reading MC job summary (exclusively for util::ReadSimuJobSummary function)
  void SetDataSimuJob(std::vector<std::string> const &data) {
    sfname    = data[0];
    rfname    = data[1];
    generator = data[2];
    process   = data[3];
    ngenreq   = stod(data[4]);
    ntried    = stod(data[5]);
    genvol    = stod(data[6]);
    lumi      = stod(data[7]);
    ebeam     = stod(data[8]);      
    if (generator.compare("simc")==0) {
      charge  = stod(data[9])/1000.;
      if (data.size()>10) {
	usingRS = stoi(data[10]);
	maxwtRS = stod(data[11]);
      }
    }
    if (generator.compare("g4sbs")==0) {
      ibeam   = stod(data[9])*1e-6;      
    }
  }

  // define an ostream operator to print to screen conveniently
  friend std::ostream& operator <<(std::ostream &out, const SimuJob& sjob) {
    out << " ------------" << std::endl;
    out << " Simu. file name    : " << sjob.sfname << std::endl;
    out << " Replayed file name : " << sjob.rfname << std::endl;
    out << " Generator          : " << sjob.generator << std::endl;
    out << " Process            : " << sjob.process << std::endl;
    out << " # events requested : " << sjob.ngenreq << std::endl;
    out << " # events tried     : " << sjob.ntried << std::endl;
    out << " Generation vol.    : " << sjob.genvol << std::endl;
    out << " Luminosity         : " << sjob.lumi << std::endl;
    out << " Ebeam (GeV)        : " << sjob.ebeam << std::endl;
    if (sjob.generator.compare("simc")==0)
      out << " Charge (C)         : " << sjob.charge << std::endl;
      out << " Using RS           : " << sjob.usingRS << std::endl;
      out << " Max Weight RS      : " << sjob.maxwtRS << std::endl;
    if (sjob.generator.compare("g4sbs")==0) 
      out << " Ibeam (A)          : " << sjob.ibeam << std::endl;
    out << " ------------" << std::endl << std::endl;
    return out;
  }


} SimuJob_t;  

#endif
