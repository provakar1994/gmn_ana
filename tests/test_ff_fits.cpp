/*
  Unit test for the EMFF fits in EMFFFits.h (Ye, Kelly, Galster).
  ----
  1. Ye2017 (w/ TPE) vs Ye's official lookup tables (DB/FF_fits/Ye2018/ye2018_wTPE_*.csv),
     at the table nodes with 0.5 <= Q2 <= 15 GeV2:
     - the parametrized FF value agrees w/ the table value to a small fraction of dG
       (the tables are printed w/ 7 significant digits);
     - GetFFwErr returns the table error dG;
     - GetFFwParErr returns the table's parametrized error dG_Par.
     The tables give G/GD (GEp, GEn) and G/(mu*GD) (GMp, GMn).
  2. Ye2017 no-TPE (placeholder, "noExtraTPE" table): returns the table GEp, GMp & errors.
  3. Ye2017 outside the table range: error value -1000.
  4. Kelly2004 and Galster1971 vs the formulas of the papers, written out independently:
     - Kelly, PRC 70, 068202 (2004), Eq. 2 & Table I; tau = Q2/(4 Mp^2) for all FFs:
         GEp, GMp/mu_p, GMn/mu_n = (1 + a1*tau) / (1 + b1*tau + b2*tau^2 + b3*tau^3)
         GEn = A*tau/(1 + B*tau) * GD,  A = 1.70, B = 3.30
     - Galster (as quoted in Perdrisat et al., Prog. Part. Nucl. Phys. 59, 694, Eq. 16);
       tau = Q2/(4 Mn^2):  GEn = -mu_n*tau/(1 + 5.6*tau) * GD
*/
#include "Constants.h"
#include "EMFFFits.h"
#include "LookUpTableReader.h"
#include "TestUtils.h"
#include <cmath>

namespace {
  double const kRoundOffTol = 1e-12; // same formula / same table entry: only round-off
  double const kValueTolInErr = 5e-3; // |Ye param - table value| in units of dG (table precision)
  double const kParErrTol = 5e-4;    // Ye parametrized error vs dG_Par (table precision)
  double const kQ2Min = 0.5, kQ2Max = 15.; // GeV2, Q2 range checked against the Ye tables
  double const kLambda2 = 0.71;     // GeV2, dipole mass parameter

  double GD(double Q2) { return 1./std::pow(1. + Q2/kLambda2, 2); }

  // Kelly's rational form (1 + a1*tau) / (1 + b1*tau + b2*tau^2 + b3*tau^3)
  double KellyRatio(double tau, double a1, double b1, double b2, double b3) {
    return (1. + a1*tau) / (1. + b1*tau + b2*tau*tau + b3*tau*tau*tau);
  }

  // One FF of the Ye tables: which table, column names, and the G(0) factor (1 or mu)
  struct YeTableFF {
    G_t kG;
    std::string file, valCol, errCol, parErrCol;
    double G0;
  };
}

int main() {
  // -- 1. Ye2017 (w/ TPE) vs the official lookup tables --
  Ye2017 ye;
  YeTableFF const yeFFs[] = {
    {G_t::kGEp, "ye2018_wTPE_proton_lookup.csv",  "GEp/GD",   "dGEp/GD",   "dGEp_Par/GD",   1.},
    {G_t::kGMp, "ye2018_wTPE_proton_lookup.csv",  "GMp/muGD", "dGMp/muGD", "dGMp_Par/muGD", constant::mup},
    {G_t::kGEn, "ye2018_wTPE_neutron_lookup.csv", "GEn/GD",   "dGEn/GD",   "dGEn_Par/GD",   1.},
    {G_t::kGMn, "ye2018_wTPE_neutron_lookup.csv", "GMn/muGD", "dGMn/muGD", "dGMn_Par/muGD", constant::mun}};
  for (auto const& ff : yeFFs) {
    LookUpTableReader table;
    table.readCSV(EMFFFits::GetDBPath("FF_fits/Ye2018/" + ff.file));
    int nNodes = 0;
    for (int row = 0; row < table.getRowCount(); ++row) {
      double const Q2 = table.GetValueByRowAndColumn(row, 0);
      if (Q2 < kQ2Min || Q2 > kQ2Max) continue;
      ++nNodes;
      double const scale = GD(Q2)*std::abs(ff.G0); // table entry -> FF (or error)
      double const tabVal = table.GetValueByRowAndColumnName(row, ff.valCol)*GD(Q2)*ff.G0;
      double const tabErr = table.GetValueByRowAndColumnName(row, ff.errCol)*scale;
      double const tabParErr = table.GetValueByRowAndColumnName(row, ff.parErrCol)*scale;
      std::string const at = " " + ff.valCol + " at Q2 = " + std::to_string(Q2);

      std::vector<double> const wErr = ye.GetFFwErr(ff.kG, Q2);
      std::vector<double> const wParErr = ye.GetFFwParErr(ff.kG, Q2);
      CHECK_ABS(wErr[0], tabVal, kValueTolInErr*tabErr, "Ye value vs table" + at);
      CHECK_REL(wErr[1], tabErr, kRoundOffTol, "Ye GetFFwErr = table dG" + at);
      CHECK_REL(wParErr[0], wErr[0], kRoundOffTol, "Ye GetFFwParErr value = GetFFwErr value" + at);
      CHECK_REL(wParErr[1], tabParErr, kParErrTol, "Ye GetFFwParErr = table dG_Par" + at);
    }
    CHECK_TRUE(nNodes > 100, "Ye table has nodes in the checked Q2 range: " + ff.valCol);
  }

  // -- 2. Ye2017 no-TPE (placeholder): the "noExtraTPE" table, pure FF values --
  LookUpTableReader noTPETable;
  noTPETable.readCSV(EMFFFits::GetDBPath("FF_fits/Ye2018/ye2018_noExtraTPE_proton_lookup.csv"));
  for (int row = 0; row < noTPETable.getRowCount(); ++row) {
    double const Q2 = noTPETable.GetValueByRowAndColumn(row, 0);
    if (Q2 < kQ2Min || Q2 > kQ2Max) continue;
    std::string const at = " at Q2 = " + std::to_string(Q2);
    std::vector<double> const GEp = ye.GetFFwErrNoTPE(G_t::kGEp, Q2);
    std::vector<double> const GMp = ye.GetFFwErrNoTPE(G_t::kGMp, Q2);
    CHECK_REL(GEp[0], noTPETable.GetValueByRowAndColumnName(row, "GEp"), kRoundOffTol, "Ye no-TPE GEp" + at);
    CHECK_REL(GEp[1], noTPETable.GetValueByRowAndColumnName(row, "dGEp"), kRoundOffTol, "Ye no-TPE dGEp" + at);
    CHECK_REL(GMp[0], noTPETable.GetValueByRowAndColumnName(row, "GMp"), kRoundOffTol, "Ye no-TPE GMp" + at);
    CHECK_REL(GMp[1], noTPETable.GetValueByRowAndColumnName(row, "dGMp"), kRoundOffTol, "Ye no-TPE dGMp" + at);
  }

  // -- 3. Ye2017 outside the table range (Q2 up to ~9.9e3 GeV2) --
  double const kQ2OutOfRange = 2e4;
  CHECK_TRUE(ye.GetFFwErr(G_t::kGMp, kQ2OutOfRange)[1] == -1000, "Ye GetFFwErr outside the table range");

  // -- 4. Kelly2004 and Galster1971 vs the paper formulas --
  Kelly2004 kelly;
  Galster1971 galster;
  for (double Q2 = kQ2Min; Q2 <= kQ2Max; Q2 += 0.5) {
    std::string const at = " at Q2 = " + std::to_string(Q2);
    double const tau_p = Q2/(4.*constant::Mp*constant::Mp);
    double const tau_n = Q2/(4.*constant::Mn*constant::Mn);
    CHECK_REL(kelly.GetFF(G_t::kGEp, Q2), KellyRatio(tau_p, -0.24, 10.98, 12.82, 21.97), kRoundOffTol, "Kelly GEp" + at);
    CHECK_REL(kelly.GetFF(G_t::kGMp, Q2), constant::mup*KellyRatio(tau_p, 0.12, 10.97, 18.86, 6.55), kRoundOffTol, "Kelly GMp" + at);
    CHECK_REL(kelly.GetFF(G_t::kGMn, Q2), constant::mun*KellyRatio(tau_p, 2.33, 14.72, 24.20, 84.1), kRoundOffTol, "Kelly GMn" + at);
    CHECK_REL(kelly.GetFF(G_t::kGEn, Q2), 1.70*tau_p/(1. + 3.30*tau_p)*GD(Q2), kRoundOffTol, "Kelly GEn" + at);
    CHECK_REL(galster.GetFF(G_t::kGEn, Q2), -constant::mun*tau_n/(1. + 5.6*tau_n)*GD(Q2), kRoundOffTol, "Galster GEn" + at);
  }

  return test::Summary("test_ff_fits");
}
