#import numpy as np
#import pandas as pd

def readfile(macfile):
    '''Reads a .mac file and returns a list
    containing all the lines in the file'''
    # Defining as list to hold mac file data
    macfile_data = []
    # Opening and reading the mac file
    with open(macfile) as f:
        #macfile_data = f.read();
        for line in f:
            macfile_data.append(line)
    return macfile_data

def update_flag_value(macfile_data, 
                      flag, 
                      oldvalue, 
                      newvalue):
    '''Updates the value of a specific g4sbs flag'''
    # Defining an empty list to hold updated mac file data
    updated_macfile_data = []
    # Update the desired flag
    for item in macfile_data:
        if (flag in item) and (item[0] != '#'):
            # Raise exception if given flag value doesn't exist
            if oldvalue not in item:
                raise Exception("Default flag value provided is wrong!")
            # Update flag with desired value
            item = item.replace(oldvalue, newvalue)
        updated_macfile_data.append(item)
    return updated_macfile_data

def writefile(filename, 
              updated_macfile_data):
    '''Writes a mac file with updated info'''
    with open(filename, 'w') as f:
        for line in updated_macfile_data:
            f.write(line)


def main():

    # User inputs
    infilename = 'sbs14_sbs70p_QE.mac'
    outfilebase = 'sbs14_sbs70p_QE'
    oldvalue = '.82'
    list_of_newvalue = ['.83', '.84']

    # Read all the lines in the mac file
    macfile_data = readfile(infilename)

    for item in list_of_newvalue:
        # Update the desired flag value
        updated_macfile_data = update_flag_value(macfile_data, 
                                                 '/g4sbs/sbsmagfield',
                                                 oldvalue,
                                                 item)
        # Write a new mac file with updated values
        outfilename = outfilebase + '_' + item.replace('.', '') + 'T.mac'
        writefile(outfilename, updated_macfile_data)


if __name__=="__main__":
    main()
