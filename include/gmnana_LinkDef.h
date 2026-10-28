#if defined(__CINT__) || defined(__CLING__)

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

// ROOT >= 6.40 rootcling needs at least one positive selection rule
#pragma link C++ namespace constant;

// general utility namespace
//#pragma link C++ namespace gmnana;

//#pragma link C++ function setrootvar::setbranch(TTree*, std::string, std::string, void*);
//#pragma link C++ function setrootvar::setbranch(TTree*, std::string, std::vector<std::string>, std::vector<void*>);
//#pragma link C++ function setrootvar::setbranch(TTree*, std::string, std::vector<std::string>, std::vector<void*>, int);
//#pragma link C++ function setrootvar::setbranch(TTree*, std::string, std::vector<std::string>, std::vector<void*>, std::vector<int>);

#endif
