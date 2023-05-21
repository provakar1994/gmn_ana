#ifndef CODA_RUN_H
#define CODA_RUN_H

typedef struct CodaRun {

  int runnum;
  int sbsconf;
  std::string target;
  int sbsmag;           // SBS magnet current (A)
  int bbmag;            // BB magnet current (A)
  double ebeam;         // GeV, avg. over entire run 
  double ebeam_std;     // GeV, std. over entire run 
  double charge;        // C, total charge collected by the run
  double DAQltime;      // %

  // constructor 
CodaRun(): 
  runnum(0),sbsconf(0),target("NONE"),sbsmag(0),bbmag(0),ebeam(0),ebeam_std(0),charge(0),DAQltime(0)
  {}

  // define an ostream operator to print to screen conveniently
  friend std::ostream& operator <<(std::ostream &out, const CodaRun& crun) {
    out << " ------------" << std::endl;
    out << " Run number        : " << crun.runnum << std::endl;
    out << " SBS config        : " << crun.sbsconf << std::endl;
    out << " Target            : " << crun.target << std::endl;
    out << " SBS mag. cur. (A) : " << crun.sbsmag << std::endl;
    out << " BB mag. cur. (A)  : " << crun.bbmag << std::endl;
    out << " Avg. ebeam (GeV)  : " << crun.ebeam << std::endl;
    out << " Ebeam std. (GeV)  : " << crun.ebeam_std << std::endl;
    out << " Tot. charge. (C)  : " << crun.charge << std::endl;
    out << " ------------" << std::endl << std::endl;
    return out;
  }

  // sets data by reading runsheet (exclusively for util::ReadRunList functions)
  void SetDataRunSheet(std::vector<std::string> data) {
    if (data.size()==8) {
        sbsconf   = stoi(data[0]);
        runnum    = stoi(data[1]);
        target    = data[2];
        sbsmag    = stoi(data[3]);
        bbmag     = stoi(data[4]);
        ebeam     = stod(data[5]);
        ebeam_std = stod(data[6]);
        charge    = stod(data[7]);
     } else 
      throw std::runtime_error("Potential NaN column entry in the run sheet!");
    }

} CodaRun_t;  

#endif
