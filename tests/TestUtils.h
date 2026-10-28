/*
  Minimal helpers for the gmnana unit tests (no external test framework needed).
  ----
  Usage in a test:
    CHECK_REL(actual, expected, relTol, "what is checked");
    CHECK_ABS(actual, expected, absTol, "what is checked");
    CHECK_TRUE(condition, "what is checked");
    return test::Summary("test_name");   // 0 if all checks passed, 1 otherwise
  Every failed check is printed with its file and line. ctest treats a non-zero
  return value as a failed test.
*/
#ifndef GMNANA_TEST_UTILS_H
#define GMNANA_TEST_UTILS_H

#include <cmath>
#include <cstdio>
#include <string>

namespace test {

  // Number of checks run and failed so far (shared by all checks in one test executable)
  inline int& NChecks() { static int n = 0; return n; }
  inline int& NFailed() { static int n = 0; return n; }

  // Records one check; prints the details if it failed
  inline void Record(bool passed, std::string const& what, std::string const& detail,
		     char const* file, int line) {
    ++NChecks();
    if (passed) return;
    ++NFailed();
    std::printf("FAILED: %s\n   %s\n   at %s:%d\n", what.c_str(), detail.c_str(), file, line);
  }

  // Passes if |actual - expected| <= relTol*|expected|
  inline void CheckRel(double actual, double expected, double relTol, std::string const& what,
		       char const* file, int line) {
    double const diff = std::abs(actual - expected);
    bool const passed = diff <= relTol*std::abs(expected);
    char detail[200];
    std::snprintf(detail, sizeof(detail), "actual = %.12g, expected = %.12g, rel. diff = %.3g (tol %.1g)",
		  actual, expected, expected != 0 ? diff/std::abs(expected) : diff, relTol);
    Record(passed, what, detail, file, line);
  }

  // Passes if |actual - expected| <= absTol (use when the expected value can be ~0)
  inline void CheckAbs(double actual, double expected, double absTol, std::string const& what,
		       char const* file, int line) {
    double const diff = std::abs(actual - expected);
    char detail[200];
    std::snprintf(detail, sizeof(detail), "actual = %.12g, expected = %.12g, abs. diff = %.3g (tol %.1g)",
		  actual, expected, diff, absTol);
    Record(diff <= absTol, what, detail, file, line);
  }

  // Prints the final tally and returns the exit code for main()
  inline int Summary(std::string const& testName) {
    std::printf("%s: %d/%d checks passed\n", testName.c_str(), NChecks() - NFailed(), NChecks());
    return NFailed() == 0 ? 0 : 1;
  }

} //::test

#define CHECK_REL(actual, expected, relTol, what) test::CheckRel((actual), (expected), (relTol), (what), __FILE__, __LINE__)
#define CHECK_ABS(actual, expected, absTol, what) test::CheckAbs((actual), (expected), (absTol), (what), __FILE__, __LINE__)
#define CHECK_TRUE(cond, what) test::Record((cond), (what), "condition is false", __FILE__, __LINE__)

#endif
