#ifndef CODA_RUN_H
#define CODA_RUN_H

typedef struct CodaRun {

  int runnum;
  int sbsconf;
  std::string target;
  int sbsmag;           // SBS magnet current (A)
  int bbmag;            // BB magnet current (A)
  double ebeam;         // GeV, avg. over entire run 
  double charge;        // C, total charge collected by the run
  double DAQltime;      // %

  // constructor 
CodaRun(): 
  runnum(0),sbsconf(0),target("NONE"),sbsmag(0),bbmag(0),ebeam(0),charge(0),DAQltime(0)
  {}

  // sets data by reading runsheet (exclusively for util::ReadRunList functions)
  void SetDataRunSheet(std::vector<std::string> data) {
    sbsconf = stoi(data[0]);
    runnum = stoi(data[1]);
    target = data[2];
    sbsmag = stoi(data[3]);
    bbmag = stoi(data[4]);
    ebeam = stod(data[5]);
  }

} CodaRun_t;  

#endif
