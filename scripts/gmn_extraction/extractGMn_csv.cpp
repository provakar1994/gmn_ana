/* Extracts GMn for all kinematic points of a results CSV file (e.g. scripts/lta/results.csv),
   for all EMFF models of the gmnx module, w/ full error propagation, and writes a CSV table.
   Input CSV (one row per kinematic point, columns found by name in the header):
     - required: q2 (GeV2), eps (acceptance-averaged Q2 & eps), R (n/p Born CS ratio), R_stat
     - systematics: every column after R_stat (absolute errors on R; all uncorrelated).
       Columns ending in "_sc" are scale errors, all others point-to-point. Empty -> 0 (warning).
     - labels: every other column before R_stat (e.g. conf, field, etheta) is copied to the output.
   Output CSV: one row per (kinematic point x model), all errors on GMn/(mu_n*GD):
     <labels>,Q2,eps,etheta_eff,R,R_stat,Model,TPE,Val,StatErr,<p2p sources>,pCS,GEn,
     <scale sources>,TotSys_p2p,TotSys_sc,StatP2P,TotErr
     A point-to-point source X is written as "X_p2p" if a scale source "X_sc" exists.
   Physics, models and error propagation: see README.md in this directory.
   Usage  : (from this directory)
            root -l -b -q ../../../install_gmn_ana/runenv/rootlogon.C \
              'extractGMn_csv.cpp("../lta/results.csv", "gmn_results.csv")'
*/
#include <algorithm>
#include <fstream>
#include <iostream>
#include <sstream>

#include "GMnExtraction.h"

namespace {
  // Columns of the input CSV w/ a fixed meaning (all others are labels or systematics)
  std::string const kColQ2 = "q2", kColEps = "eps", kColR = "R", kColRstat = "R_stat";
  std::string const kScaleSuffix = "_sc", kP2PSuffix = "_p2p";

  // Splits one CSV line at commas, keeping empty fields; strips a trailing '\r' (CRLF files)
  std::vector<std::string> SplitCSVLine(std::string line) {
    if (!line.empty() && line.back() == '\r') line.pop_back();
    std::vector<std::string> fields;
    std::stringstream ss(line);
    std::string field;
    while (std::getline(ss, field, ',')) fields.push_back(field);
    if (!line.empty() && line.back() == ',') fields.push_back(""); // trailing empty field
    return fields;
  }

  // Index of column name in header, -1 if not found
  int ColumnIndex(std::vector<std::string> const& header, std::string const& name) {
    auto it = std::find(header.begin(), header.end(), name);
    return it == header.end() ? -1 : static_cast<int>(it - header.begin());
  }

  bool EndsWith(std::string const& s, std::string const& suffix) {
    return s.size() >= suffix.size() && s.compare(s.size() - suffix.size(), suffix.size(), suffix) == 0;
  }
}

int extractGMn_csv(std::string const& inFile, std::string const& outFile) {

  // -- Read the header --
  std::ifstream in(inFile);
  if (!in.is_open()) { std::cerr << "ERROR: cannot open " << inFile << "\n"; return 1; }
  std::string line;
  std::getline(in, line);
  std::string const bom = "\xEF\xBB\xBF"; // UTF-8 byte order mark, if present
  if (line.compare(0, bom.size(), bom) == 0) line.erase(0, bom.size());
  std::vector<std::string> const header = SplitCSVLine(line);

  int const iQ2 = ColumnIndex(header, kColQ2), iEps = ColumnIndex(header, kColEps);
  int const iR = ColumnIndex(header, kColR), iRstat = ColumnIndex(header, kColRstat);
  if (iQ2 < 0 || iEps < 0 || iR < 0 || iRstat < 0) {
    std::cerr << "ERROR: " << inFile << " needs the columns " << kColQ2 << ", " << kColEps << ", "
	      << kColR << ", " << kColRstat << "\n";
    return 1;
  }
  // labels: columns before R_stat w/o a fixed meaning; systematics: all columns after R_stat
  std::vector<int> iLabels, iSys;
  std::vector<std::string> sysNames;
  for (int i = 0; i < static_cast<int>(header.size()); ++i) {
    if (i == iQ2 || i == iEps || i == iR || i == iRstat) continue;
    if (i < iRstat) iLabels.push_back(i);
    else { iSys.push_back(i); sysNames.push_back(header[i]); }
  }

  // -- Output columns --
  // p2p sources first (as "X_p2p" if "X_sc" exists), then pCS & GEn, then the scale sources
  std::vector<size_t> p2pOrder, scaleOrder; // indices into sysNames
  for (size_t j = 0; j < sysNames.size(); ++j)
    (EndsWith(sysNames[j], kScaleSuffix) ? scaleOrder : p2pOrder).push_back(j);
  auto outName = [&](size_t j) {
    std::string const& name = sysNames[j];
    bool const hasScalePartner = std::find(sysNames.begin(), sysNames.end(), name + kScaleSuffix) != sysNames.end();
    return hasScalePartner ? name + kP2PSuffix : name;
  };
  std::ofstream out(outFile);
  if (!out.is_open()) { std::cerr << "ERROR: cannot write " << outFile << "\n"; return 1; }
  for (int i : iLabels) out << header[i] << ",";
  out << "Q2,eps,etheta_eff,R,R_stat,Model,TPE,Val,StatErr";
  for (size_t j : p2pOrder) out << "," << outName(j);
  out << ",pCS,GEn";
  for (size_t j : scaleOrder) out << "," << sysNames[j];
  out << ",TotSys_p2p,TotSys_sc,StatP2P,TotErr\n";

  // -- Process every kinematic point --
  std::cout << "\nPrimary model per point, GMn/(mu_n*GD) +/- stat +/- p2p sys +/- scale sys (total):\n";
  int nPoints = 0, nRows = 0, lineNo = 1;
  while (std::getline(in, line)) {
    ++lineNo;
    if (line.empty() || line == "\r") continue;
    std::vector<std::string> f = SplitCSVLine(line);
    f.resize(header.size()); // missing trailing fields -> empty
    if (f[iQ2].empty() || f[iEps].empty() || f[iR].empty() || f[iRstat].empty()) {
      std::cerr << "ERROR: line " << lineNo << ": empty " << kColQ2 << "/" << kColEps << "/" << kColR
		<< "/" << kColRstat << ", skipping\n";
      continue;
    }
    std::vector<double> sysVals;
    for (size_t j = 0; j < iSys.size(); ++j) {
      if (f[iSys[j]].empty())
	std::cout << "WARNING: line " << lineNo << ": " << sysNames[j] << " is empty, using 0\n";
      sysVals.push_back(f[iSys[j]].empty() ? 0. : std::stod(f[iSys[j]]));
    }

    gmnx::Input const input = gmnx::MakeInput(std::stod(f[iQ2]), std::stod(f[iEps]), std::stod(f[iR]),
					      std::stod(f[iRstat]), sysNames, sysVals);
    std::vector<gmnx::ModelResult> const results = gmnx::ExtractAll(input);
    ++nPoints;

    // one output row per model
    for (auto const& r : results) {
      for (int i : iLabels) out << f[i] << ",";
      out << Form("%.4f,%.4f,%.4f,%.4f,%.4f,%s,%s,%.4f,%.4f",input.Q2,input.eps,input.etheta,input.R,input.dR_stat,
		  r.model.name.c_str(),gmnx::TPEGroupName(r.model.tpe).c_str(),r.val,r.stat);
      for (size_t j : p2pOrder) out << Form(",%.4f",r.sys[j]);
      out << Form(",%.4f,%.4f",r.pCS,r.GEnSys);
      for (size_t j : scaleOrder) out << Form(",%.4f",r.sys[j]);
      out << Form(",%.4f,%.4f,%.4f,%.4f\n",r.sysP2P,r.sysScale,r.statP2P,r.tot);
      ++nRows;
    }

    // screen summary: primary model
    gmnx::ModelResult const& p = results[gmnx::kPrimaryModel];
    std::string labels;
    for (int i : iLabels) labels += header[i] + "=" + f[i] + " ";
    std::cout << Form("%sQ2=%.3f eps=%.3f R=%.4f | %s: %.4f +/- %.4f +/- %.4f +/- %.4f (%.4f)\n",labels.c_str(),
		      input.Q2,input.eps,input.R,p.model.name.c_str(),p.val,p.stat,p.sysP2P,p.sysScale,p.tot);
  }

  std::cout << Form("\nWrote %d rows (%d points x %zu models) to %s\n\n",nRows,nPoints,
		    gmnx::ModelDefs(1.).size(),outFile.c_str());
  return 0;
}
