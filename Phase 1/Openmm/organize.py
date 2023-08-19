# Import necessary libraries. These provide functionalities for various tasks such as:
# subprocess: Running shell commands from Python
# time: Time-related functions
# datetime: Date-related functions
# shutil: High-level file operations like copying and removal
# os: Interacting with the operating system
# math and random: Mathematical functions and random number generation
# numpy: Scientific computing and array manipulation
#import subprocess
#import time
import datetime
import shutil
import os
from math import *
from random import *
import numpy as np
import pandas as pd
import glob
import matplotlib.pyplot as plt
import re

parameters="multi_run.dat"
#num_run is first number, says how many runs are we doing?
#start is the second number, says where we are in that list
#marker is are we at the start or end of a task. default is S for starting


# define file names and directories
dest_file="1LC.txt"
dest_file_debug="debug_file.txt"
dest_file_OLC="OneLetterCode.dat"
dest_file_TLC="ThreeLetterCode.dat"
dest_required_file="stats_module.dat"
dest_required_file_csv="data_multi.csv"
cur_dir=os.getcwd()
dir_comp=f"{cur_dir}/completed_run"

        
# This function checks if a specified parameter file exists. If not, it creates the file and writes default values into it.
def check_multi(parameters,marker='S'):
    # If the specified parameter file does not exist, create it.
    if not os.path.exists(parameters):
        f = open(parameters, 'x')
        # Write default values into the newly created parameter file.
        f.write('1 #Number of proteins you care to run')
        f.write('\n1 #Start value (if you are just starting the process should be on 1 by default)')
        f.write(f'\n{marker} #Put S or E based on starting or end process (don\'t really need to touch), if its I then set the other 2 numbers and run the file again.') 
        f.close()# ^ the exception to S or E is I for initializing, just makes sure it's all there.

# This function checks if a specified destination file exists. If not, it creates the file.
def check_dest_file(dest_file='1LC.txt'):
    # If the specified destination file does not exist, create it.
    if not os.path.exists(dest_file):
        f = open(dest_file, 'x')
        f.close()

# This function checks if a specified debug file exists. If not, it creates the file.
def check_dest_file_debug(dest_file_debug="debug_file.txt"):
    # If the specified debug file does not exist, create it.
    if not os.path.exists(dest_file_debug):
        f = open(dest_file_debug, 'x')
        f.close()
        
# This function reads a multi-run parameter file and extracts the parameters from it.
def read_multi_run(parameters='multi_run.dat'):
    # Open the specified parameter file.
    with open(parameters,'r') as para:
        lines = []
        # Read the file line by line and store each line in a list.
        for line in para:
            lines.append(line)
    para.close()
    num_run_raw=lines[0]
    start_raw=lines[1]
    marker_raw=lines[2]
    
    # Extract the number of proteins to run.
    if '#Number of proteins you care to run' in num_run_raw:
        num_run=int(num_run_raw.replace('#Number of proteins you care to run','').strip())
    else:
        num_run=int(num_run_raw)

    # Extract the start value.
    if '#Start value (if you are just starting the process should be on 1 by default)' in start_raw:
        start=int(start_raw.replace('#Start value (if you are just starting the process should be on 1 by default)','').strip())
    else:
        start=int(start_raw)

    # Extract the marker.
    if '#Put S or E based on starting or end process (don\'t really need to touch), if its I then set the other 2 numbers and run the file again.' in marker_raw:
        marker=marker_raw.replace('#Put S or E based on starting or end process (don\'t really need to touch), if its I then set the other 2 numbers and run the file again.','')
    else:
        marker=marker_raw
    marker=marker.strip()

    # Remove the lines that have been processed.
    lines.remove(num_run_raw)
    lines.remove(start_raw)
    lines.remove(marker_raw)
  
    clean_list=[]
  
    # Go through the remaining lines and clean them up.
    for line in lines:
        line=line.strip()
        # If the line is not empty,
        if line.strip():
            # If the line contains a newline character, remove it.
                if '\n' in line:
                    # Remove newline character and any leading or trailing white spaces
                    line=line.replace('\n','').strip()
                    clean_list.append(line)
                else:
                    # Add the line into the clean list after removing leading or trailing white spaces
                    clean_list.append(line.strip())
                
    # Return the number of runs, start value, marker, and the cleaned list of lines
    return num_run,start,marker,clean_list

def edit_multi_run(parameters,num_run,start,marker,list_files):
    # Open the parameter file for writing
    with open(parameters,'w') as para:
        # Check marker status and edit start and marker accordingly
        if marker=='S':
            marker='E'
        elif marker=='E':
            start+=1
            marker='S'
            
        # Write the number of runs, start value, and marker to the file
        para.write(f'{num_run} #Number of proteins you care to run')
        para.write(f'\n{start} #Start value (if you are just starting the process should be on 1 by default)')
        para.write(f'\n{marker} #Put S or E based on starting or end process (don\'t really need to touch), if its I then set the other 2 numbers and run the file again.')

        # Write all the remaining files from the list_files to the parameter file
        for filing in list_files:
            para.write(f'\n{filing}')
    # Close the file
    para.close()
    
def cap_original_files(parameters='multi_run.dat'):
    # Get the current directory
    cur_dir=os.getcwd()
    
    # Get a list of all files in the current directory
    original_files = [f for f in os.listdir(cur_dir) if os.path.isfile(os.path.join(cur_dir, f))]
    
    # Append all original files to the parameter file
    with open(parameters,'a') as para:
        for o_file in original_files:
            para.write(f'\n{o_file}')
    
    # Close the file
    para.close()

def check_def_dir():
    # Get the current directory
    cur_dir=os.getcwd()
    
    # Define the paths to the run and completed directories
    dir_comp=f'{cur_dir}/completed_run'
  
    # Check if the directories exist, if not, create them
    if not os.path.exists(dir_comp):
        os.mkdir(dir_comp)
                     
def organized_end_dir_creator(start):
    # Get current directory
    Protein_Name=row_specific_info_csv(start,'Name')
    cur_dir=os.getcwd()
    # Define directory for completed run
    dir_comp=f'{cur_dir}/completed_run'
    # Define the end directory name
    end_dir_name=f'{dir_comp}/{Protein_Name}/'
    # Check if the directory exists
    if not os.path.exists(end_dir_name):
        # If not, create it
        os.mkdir(end_dir_name)
    else:
        end_dir_name=f'{dir_comp}/{Protein_Name}'
        i=2
        name_check=end_dir_name
        while True:
            name_check=f'{end_dir_name}_{i}/' 
            if not os.path.exists(name_check):#checking to make sure does not OverWrite dir info
                os.mkdir(name_check)
                end_dir_name=name_check
                break
            else:
                i+=1
    # Return the end directory name
    return end_dir_name

def change_e_to_s(parameters,num_run,start,marker,list_files,dest_file='1LC.txt'):
    # Create end directory and get its path
    current_end_dir=organized_end_dir_creator(start)
    # Get current directory
    cur_dir=os.getcwd()
    # Define the location of the one letter file
    location_One_letter=f"{cur_dir}/{dest_file}"
    # Define the end directory file for one letter
    end_dir_file_One_letter=f"{current_end_dir}/{dest_file}"
    # Copy the one letter file to the end directory
    shutil.copy2(location_One_letter,end_dir_file_One_letter)
    # Get a list of all files in the current directory
    current_list_of_files= [f for f in os.listdir(cur_dir) if os.path.isfile(os.path.join(cur_dir, f))]
    # Create a copy of the list
    files_in_need_for_transportation=current_list_of_files[:]
    # Edit multi run parameters
    edit_multi_run(parameters,num_run,start,marker,list_files)
    # Loop over all files in current directory
    for c_file in current_list_of_files:
        # And for all files in the list files
        for o_file in list_files:
            # If file from current directory matches a file from the list, remove it from the files needing transportation
            if c_file==o_file:
                files_in_need_for_transportation.remove(o_file)
    # Loop over all non default files
    for not_default_file in files_in_need_for_transportation:
        # Define the location of the non default file
        not_default_file_location=f'{cur_dir}/{not_default_file}'
        # Move the non default file to the end directory
        shutil.move(not_default_file_location,current_end_dir)

# This function checks whether the required directories and files exist and if not, creates them. 
# It then reads in the parameters from the specified file and runs the source file checker function.
# The source file name is also created based on the start value. The function returns the parameters read in.
def check_all(parameters):
    check_def_dir()  # Check if the default directories exist, create them if not.
    num_run,start,marker,list=read_multi_run(parameters)  # Read in parameters from the parameters file.
    marker=marker.strip()  # Remove any leading/trailing spaces from the marker.
    return num_run,start,marker,list

# I forgot comments on definitions have this format not the previous
def find_stokes_num_from_hoomd_file():
    """
    Scans through the current directory for files starting with 'hoomd', strips the string 'hoomd-' and '.out' from the file name and converts the resultant string to an integer.
    Returns this integer representing the hoomd number.
    """
    for file in os.listdir(os.getcwd()):  # Iterate over all files in the current directory
        if file.startswith("hoomd"):  # If file name starts with 'hoomd'
            file = file.replace('hoomd-', '').replace('.out', '').strip()  # Remove 'hoomd-' from the file name  # Remove '.out' from the file name
            hoomd_number = int(file)  # Convert file name to integer
    return hoomd_number

def remove_hoomd_file():
    """
    Iterates over all files in the current directory and removes the ones starting with 'hoomd'
    """
    for file in os.listdir(os.getcwd()):  # Iterate over all files in the current directory
        if file.startswith("hoomd"):  # If file name starts with 'hoomd'
            os.remove(file)  # Remove the file

def add_hoomd_file_multi(list):
    """
    Iterates over all files in the current directory and adds the one starting with 'hoomd' to a list
    """
    for item in list:
        if item.startswith("hoomd"):
            return list
    for file in os.listdir(os.getcwd()):  # Iterate over all files in the current directory
        if file.startswith("hoomd"):  # If file name starts with 'hoomd'
            list.append(file)  # Remove the file        
    return list

# Special use-case for plotting current inside
def clean_up():
    """
    Removes specific files in the directory based on the file name patterns, and moves 'debug' file to another directory.
    """
    num_run, start, marker, list_files = read_multi_run()  # Reading details from multiple run
    if start > num_run:  # Only proceed if the starting number is greater than the total number of runs
        read_dat_files_data_merger_create_csv()
        remove_hoomd_file()
        
def read_entire_csv_return_dict(dest_required_file_csv = "data_multi.csv"):
    '''Define function to read entire CSV and return a dictionary.'''
    # Using pandas to read the CSV file located at the destination provided.
    df = pd.read_csv(dest_required_file_csv)

    # Converting the dataframe into a list of dictionaries where each dictionary represents a row of data.
    data = df.to_dict('records')

    # Return the list of dictionaries
    return data

def grab_info_from_csv_dict_row(start):
    '''Define function to grab information from a specific row in the CSV.'''
    # Call the function to read the CSV file and store the data in csv_dictionary
    csv_dictionary=read_entire_csv_return_dict()

    # Adjusting the row number to align with zero-indexed lists.
    csv_row_num=start-1 

    # Getting the information from the specific row in the CSV.
    current_row_info=csv_dictionary[csv_row_num]

    # Return the row information
    return current_row_info

def row_specific_info_csv(start,column_name):
    '''Define function to grab specific column information from a specific row in the CSV.'''
    # Call the function to get the information of the specific row
    current_row_info=grab_info_from_csv_dict_row(start)

    # Get the information of the specific column from the row.
    specific_csv_row_info=current_row_info[column_name]

    # Return the specific column information
    return specific_csv_row_info

def search_csv_sequence_return_name(file_location_protein_sequence,dest_required_file_csv = "data_multi.csv"):
    '''Define function to grab specific column information from a specific row in the CSV.'''
    
    file = open(file_location_protein_sequence, 'r')
    read_sequence=file.readline()
    file.close()
    csv_dictionary=read_entire_csv_return_dict()
    i=0
    matching_rows=[]
    for row in csv_dictionary:
        if row['Sequence']==read_sequence:
            protein_name=row['Name']
            matching_rows.append(i)
        i+=1
    return protein_name,matching_rows

def write_csv_sequence_to_1LC(dest_file="1LC.txt"):
    num_run,start,marker,clean_list=read_multi_run()
    Protein_Sequence=row_specific_info_csv(start,'Sequence')
    with open(dest_file,'w') as file:
        # Check marker status and edit start and marker accordingly
        file.write(f'{Protein_Sequence}')
    # Close the file
    file.close()

def get_csv_headers(csv_file):
    """
    This function reads a CSV file and returns all the headers as a list.

    :param csv_file: Location of the CSV file
    :return: List of headers
    """
    try:
        df = pd.read_csv(csv_file)
        headers = df.columns.tolist()
        return headers
    except FileNotFoundError:
        print(f"File not found: {csv_file}")
        return []

def get_unique_and_matching_headers_relative_to_first_csv_file(csv_file1, csv_file2):
    """
    This function takes two CSV files, retrieves the headers from each file,
    and returns a list of headers that are in the first CSV file but not in the second CSV file.

    :param csv_file1: Location of the first CSV file
    :param csv_file2: Location of the second CSV file
    :return: List of unique headers in the first CSV file
    """
    headers1 = get_csv_headers(csv_file1)
    headers2 = get_csv_headers(csv_file2)

    unique_headers = [header for header in headers1 if header not in headers2]
    matching_headers=[header for header in headers1 if header in headers2]
    return unique_headers,matching_headers    
 
def read_dat_files_data_merger_create_csv(csv_file="data_multi.csv", dir_comp=dir_comp):
    """
    This function reads .dat files from the given directory, merges them,
    and then merges the resulting DataFrame with data from a CSV file.

    :param csv_file: Location of the CSV file to merge data with
    :param dir_comp: The directory where the .dat files are located
    """

    # Use glob to get all the .dat files in subdirectories
    config_stat_files = glob.glob(f'{dir_comp}/**/config_stat.dat', recursive=True)
    csv_data = pd.read_csv(csv_file)

    # Loop through the files and read them in with pandas
    dataframes = []  # a list to hold all the individual pandas DataFrames
    print(config_stat_files)
    for filename in config_stat_files:
        print(filename)
        df = pd.read_csv(filename)  # assumes .dat files are CSV-formatted
        seq_file_loc = f'{dir_comp}/{os.path.basename(os.path.dirname(filename))}/{dest_file}'
        print(seq_file_loc)
        protein_name, matching_rows = search_csv_sequence_return_name(seq_file_loc)

        dir_name = os.path.basename(os.path.dirname(filename))
        dir_name_without_protein = dir_name.replace(protein_name, '').strip()

        # Extract the instance number from the directory name
        instance_number = None
        match = re.search(r'_(\d+)$', dir_name_without_protein)
        if match:
            instance_number = int(match.group(1))
            matching_array_value = instance_number - 1
        else:
            matching_array_value = 0
        print(f"Matching Rows: {matching_rows}")
        print(f"Matching Array Value: {matching_array_value}")
        print(f"Matching Row given Matching Array Value: {matching_rows[matching_array_value]}")

        non_matching_headers = set(df.columns) - set(csv_data.columns)
        df_subset = df[non_matching_headers].copy()  # Select only the non-matching columns from df
        df_subset.index = [matching_rows[matching_array_value]]  # Wrap the value in a list to create a single-item collection

        # Format the values as floats with a maximum of 3 decimal places
        df_subset = df_subset.applymap(lambda x: int(x) if round(float(x), 3).is_integer() else round(float(x), 3))
        dataframes.append(df_subset)  # Append the subset DataFrame to the list of dataframes

    # Concatenate all the dataframes together
    merged_data = pd.concat(dataframes)

    # Merge with the original csv_data
    merged_data = pd.merge(csv_data, merged_data, left_index=True, right_index=True, how='left')
    merged_data['# ParticleN'] = merged_data['# ParticleN'].fillna(0).astype(int)

    # Write it out
    merged_data.to_csv('merged_config_stat_results.csv', index=False)
         
# Check if the parameters file exists
if not os.path.exists(parameters):
    # If the parameters file does not exist, initialize it with default values, and also initialize the destination and debug files.
    check_multi(parameters,'I')
    check_dest_file()
    check_dest_file_debug()
    
else:
    # If the parameters file exists, call the check_all function to set up the required directories/files and read in the parameters.
    num_run,start,marker,list_files=check_all(parameters)
    # If the start value is not greater than the number of runs, proceed
    if not start>num_run:
        # If the marker indicates initialization, write the parameters into the parameters file.
        if marker =='I':
            marker='S'
            with open(parameters,'w') as para:
                para.write(f'{num_run} #Number of proteins you care to run')
                para.write(f'\n{start} #Start value (if you are just starting the process should be on 1 by default)')
                para.write(f'\n{marker} #Put S or E based on starting or end process (don\'t really need to touch), if its I then set the other 2 numbers and run the file again.')
                para.close()   
            cap_original_files()  # Capture all original files
        # If the marker is 'S', then proceed with the generation of the three different codes and edit the parameters file.
        elif marker=='S':
            write_csv_sequence_to_1LC()  # Copy letter text
            list_files=add_hoomd_file_multi(list_files)
            edit_multi_run(parameters,num_run,start,marker,list_files)  # Edit the parameters file
        # If the marker is 'E', then move files to completed directory and edit the parameters file.
        elif marker=='E':
            change_e_to_s(parameters,num_run,start,marker,list_files)
            num_run,start,marker,list_files=read_multi_run()
            if start > num_run:
                clean_up()
        else:
            print('invalid marker')  # If the marker is not any of the valid values, print error message.
