## This directory will hold scripts to deal with EPICS (E) and scaler (TSsbs) trees

#### Instructions to carry out various procedures 

- Calculate average ebam (and std) per run:
    1. Modify `conf_get_avgHALLA_p_prun.json` config file with desired values
    2. Execute `get_avgHALLA_p_prun.cpp` with the above config file
    3. Execute `get_avgHALLA_p_prun.ipynb` with the above config file
    4. Find the output CSV file in `../../DB/<target>` directory
    5. Use `../../DB/join_n_merge_rLists.ipynb` to concatenate and/or merge the new CSV file with the existing ones