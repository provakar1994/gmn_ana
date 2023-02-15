#include <vector>
#include <iostream>
#include <fstream>
#include <ostream>
#include <sstream>

using namespace std;

int closest(std::vector<int> const& vec, int value) {
  auto const it = std::lower_bound(vec.begin(), vec.end(), value);
  if (it == vec.end()) { return -1; }

  return *it;
}

struct CodaRun {
  int runnum;
  int evnum;
  double ebeam;
};

int main() 
{
  vector<CodaRun> codarVec;
  
  string fst = "test1";
  string lst = ".csv";
  //string file = to_string("test1") + to_string(".csv");
  string file = fst + lst;

  ifstream adcGain_data;
  adcGain_data.open(file);
  string readline;
  int elemID=0;
  if(adcGain_data.is_open()){
    //cout << " Reading run info from: "<< adcGain_rfile << endl;
    string skip_header;
    getline(adcGain_data, skip_header);
    while(getline(adcGain_data,readline)){
      istringstream tokenStream(readline);
      string token;
      char delimiter = ',';
      vector<string> temp;
      while(getline(tokenStream,token,delimiter)){
	string temptoken=token;
	//adcGain[elemID] = temptoken.Atof();
	temp.push_back(temptoken);
	//cout << stoi(temptoken) << " ";
      }
      //cout << stod(temp[2]) << std::endl;
      CodaRun temp_cr;
      if (stoi(temp[0]) == 12313) {// && stoi(temp[1]) == 625) {
      	temp_cr.runnum = stoi(temp[0]);
      	temp_cr.evnum = stoi(temp[1]);
      	temp_cr.ebeam = stod(temp[2]);

	codarVec.push_back(temp_cr);
      }

      temp.clear();
      //}
    }
  }else{
    cerr << " **!** No file : " << file  << endl;
    throw;
  }
  adcGain_data.close();

  cout << codarVec.size() << " " 
       << codarVec[0].runnum << " " 
       << codarVec[0].evnum << " " 
       << codarVec[0].ebeam << endl;


  std::string target = "LD2";
  std::string fst1 = "../DB/good_runList_GMn_nTPE_"; 
  std::string mid = "_pass_";
  std::string lst1 = ".csv";
  int replay_pass = 1;
  std::string run_spreadsheet = fst1 + target + mid + std::to_string(replay_pass) + lst1;
  cout << run_spreadsheet << endl;

  return 0;
}
