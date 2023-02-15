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
  CodaRun codar;
  codar.runnum = 5;
  codar.evnum = 5;
  codar.ebeam = 5.5;
  vector<codar> codarVec;

  ifstream adcGain_data;
  adcGain_data.open("test1.csv");
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
       //if (stoi(token) == 12313) {
	while(getline(tokenStream,token,delimiter)){
	  string temptoken=token;
	  //adcGain[elemID] = temptoken.Atof();
	  cout << stoi(temptoken) << " ";
	  elemID++;
	}
	cout << std::endl;
	//}
    }
  }else{
    cerr << " **!** No file : "  << endl;
    throw;
  }
  adcGain_data.close();

  return 0;
}
