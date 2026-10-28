# gmn_ana
SBS GMn analysis framework.

## Unit tests
Unit tests for the library live in `tests/` (one `test_*.cpp` executable per topic,
small inputs in `tests/data/`). They are built together with the library (disable with
`-DGMNANA_BUILD_TESTS=OFF`) and run from the build directory with:
```
ctest --output-on-failure
```
