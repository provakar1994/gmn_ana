# gmn_ana
Analysis framework for the SBS GMn / nTPE experiment (JLab Hall A, E12-09-019): extracts the
neutron magnetic form factor GMn from quasi-elastic d(e,e'n)/d(e,e'p) yield ratios, using
data/MC fits of the HCAL `dx` distribution. C++/ROOT shared library plus analysis macros.

**Start here: [docs/ANALYSIS_GUIDE.md](docs/ANALYSIS_GUIDE.md)** — analysis chain, kinematic
settings, conventions, and an index of every script with a how-to (run, inputs, outputs) per
analysis stage in [docs/howto/](docs/howto/).

## Quick start
```
mkdir build_gmn_ana && cd build_gmn_ana        # next to the source tree, not inside it
cmake -DCMAKE_INSTALL_PREFIX=$PWD/../install_gmn_ana ../gmn_ana
make -j4 && make install
ctest --output-on-failure          # unit tests (tests/); disable with -DGMNANA_BUILD_TESTS=OFF
```
Macros are run from `scripts/` (data, data/MC) or `scripts/simulation/` (MC) with a JSON config,
e.g. `root -l -b -q 'w2_offset_check.cpp("configs/conf_w2_offset_check.json")'`.
Details: [docs/howto/01_setup.md](docs/howto/01_setup.md).

## Layout
| Path | Content |
|---|---|
| `include/`, `src/` | shared library (`gmnana`) |
| `scripts/` | real-data analysis and data/MC comparisons; configs in `scripts/configs/` |
| `scripts/simulation/` | MC analysis; configs in `scripts/simulation/configs/` |
| `tests/` | unit tests (`test_*.cpp`, inputs in `tests/data/`) |
| `DB/` | run lists, form-factor tables |
| `docs/` | analysis guide and how-tos |

Every change is logged in [CHANGELOG.md](CHANGELOG.md).
