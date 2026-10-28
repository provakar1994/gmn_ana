/*
  Unit test for LookUpTableReader (used to read the EMFF lookup tables in DB/).
  ----
  Uses the small fixture tests/data/log_grid.csv:
    Q2   : 0.1, 1, 10, 100 (a log grid; the header line starts with a UTF-8 BOM)
    log10Q2 = log10(Q2)  -> linear in log(Q2), so log interpolation must be exact everywhere
    twoQ2   = 2*Q2       -> checks the interpolation weights between two rows
  Checks: BOM stripping, row count, exact values at the nodes, log-linear interpolation
  in between, the upper table edge, and the -99999 error value (out of range, bad column).
*/
#include "LookUpTableReader.h"
#include "TestUtils.h"
#include <cmath>

int main() {
  double const kTol = 1e-12;     // pure floating point round-off
  double const kErrValue = -99999; // value returned by LookUpTableReader on failure

  LookUpTableReader table;
  table.readCSV(std::string(GMNANA_TEST_DATA_DIR) + "/log_grid.csv");

  // -- header and size --
  CHECK_TRUE(table.getColumnIndex("Q2") == 0, "BOM is stripped from the 1st column name");
  CHECK_TRUE(table.getColumnIndex("twoQ2") == 2, "column index by name");
  CHECK_TRUE(table.getColumnIndex("missing") == -1, "unknown column gives -1");
  CHECK_TRUE(table.getRowCount() == 4, "number of rows read");

  // -- values at the nodes must be returned exactly --
  for (double q2 : {0.1, 1.0, 10.0, 100.0}) {
    CHECK_REL(table.GetLogInterpValueByColumnName(q2, "twoQ2"), 2.*q2, kTol, "exact at node Q2 = " + std::to_string(q2));
  }

  // -- in between: log10(Q2) is linear in log(Q2), so interpolation is exact --
  for (double q2 : {0.15, 0.7, 3.0, 42.0, 99.9}) {
    CHECK_REL(table.GetLogInterpValueByColumnName(q2, "log10Q2"), std::log10(q2), kTol,
	      "log-linear interpolation at Q2 = " + std::to_string(q2));
  }

  // -- interpolation weights: half way between Q2 = 1 and 10 in log(Q2) is sqrt(10) --
  // twoQ2 there is the mean of the two rows, (2 + 20)/2 = 11 (not 2*sqrt(10))
  CHECK_REL(table.GetLogInterpValueByColumnName(std::sqrt(10.), "twoQ2"), 11., kTol,
	    "log-linear weights half way between two rows");

  // -- error handling --
  CHECK_TRUE(table.GetLogInterpValueByColumnName(0.05, "twoQ2") == kErrValue, "below the table range");
  CHECK_TRUE(table.GetLogInterpValueByColumnName(200., "twoQ2") == kErrValue, "above the table range");
  CHECK_TRUE(table.GetLogInterpValueByColumnName(5., "missing") == kErrValue, "unknown column");

  return test::Summary("test_lookup_table");
}
