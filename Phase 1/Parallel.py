import pandas as pd
import random
import os
import shutil
import subprocess
import time
import sys
from collections import defaultdict
import datetime
import glob

max_safe_group=10
over_guess_percent=30 

cur_dir=os.getcwd()
parent_dir = os.path.dirname(cur_dir)
check_dir_cal=f"{cur_dir}/Calculations"
check_dir_archive=f"{parent_dir}/Data_Archive"



file_path_data = "data.csv"
file_path_m_data = "data_multi.csv"



if os.path.exists(check_dir_cal) and len(os.listdir(check_dir_cal))!=0:
    in_process = True
else:
    in_process = False
    
if not os.path.exists(file_path_data):#preference is given to data.csv if both exist
    if not os.path.exists(file_path_m_data):
        print(f"{file_path_data} and {file_path_m_data} don't exist!\n You need at least one of thses files.")
        exit()
    else:
        file_path=file_path_m_data
else:
    file_path=file_path_data

if not os.path.exists(check_dir_archive):
    os.makedirs(check_dir_archive)
    

def execute_shell_command(command, output_file="debug_file.txt"):
    try:
        # Execute the shell command as a subprocess
        process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)

        # Open the output file in append mode
        with open(output_file, 'a') as file:
            # Read the output line by line
            for line in process.stdout:
                # Write the output to the file
                file.write(line)
                file.flush()  # Ensure immediate write to the file

        # Wait for the process to finish
        process.wait()

        # Check the return code
        if process.returncode != 0:
            print(f"Command '{command}' failed with return code {process.returncode}")
            # You can raise an exception here or handle the error condition as needed

    except Exception as e:
        print(f"An error occurred while executing the command: {str(e)}")
        # Handle the exception as needed

def execute_shell_command_start(command, output_file="debug_file.txt"):
    try:
        # Execute the shell command as a subprocess
        for _ in range(20):
            process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)

            # Open the output file in append mode
            with open(output_file, 'a') as file:
                # Read the output line by line
                for line in process.stdout:
                    # Write the output to the file
                    file.write(line)
                    file.flush()  # Ensure immediate write to the file

                    # Check if there is any response
                    if line.strip():
                        #print("Command received a non-empty response. Stopping further execution.")
                        return

            # Wait for the process to finish
            process.wait()

            # Check the return code
            if process.returncode != 0:
                print(f"Command '{command}' failed with return code {process.returncode}")
                # You can raise an exception here or handle the error condition as needed

    except Exception as e:
        print(f"An error occurred while executing the command: {str(e)}")
        # Handle the exception as needed
     
def est_time_equation(number_of_letters):
    estimated_time = ((2.8912195151139 * 10**-7) * number_of_letters**4.14262) + 10
    return estimated_time

def est_time_equation_open(number_of_letters):
    estimated_time = ((2.8912195151139 * 10**-7) * number_of_letters**4.14262) + 10
    estimated_time= estimated_time*4
    return estimated_time

def process_csv(file_path):
    # Read the CSV file
    df = pd.read_csv(file_path)

    # Create a temporary column with the number of letters in "Sequence"
    df['letters_count'] = df['Sequence'].apply(len)
    if file_path==file_path_data:
        df['estimated_time'] = df['letters_count'].apply(est_time_equation)
    elif file_path==file_path_m_data:
        df['estimated_time'] = df['letters_count'].apply(est_time_equation_open)
    

    # Calculate the total estimated time
    total_estimated_time = df['estimated_time'].sum()

    # Calculate the average estimated time
    average_estimated_time = df['estimated_time'].mean()

    # Calculate the maximum estimated time
    max_estimated_time = df['estimated_time'].max()

    # Calculate the number of rows
    number_of_rows = len(df)

    # Sort by the 'estimated_time' column
    df_sorted = df.sort_values(by='estimated_time', ascending=False)

    # Optionally, drop the temporary column if not needed
    df_sorted = df_sorted.drop(columns=['letters_count'])

    return df_sorted, total_estimated_time, average_estimated_time, max_estimated_time, number_of_rows

def create_optimized_groups(df_sorted, max_estimated_time, max_safe_group):
    # Sort the DataFrame by estimated time in descending order
    df_sorted = df_sorted.sort_values(by='estimated_time', ascending=False).reset_index(drop=True)

    # Initialize groups and group times
    groups = defaultdict(list)
    group_times = defaultdict(float)

    # Distribute rows across groups, trying to balance the estimated times
    for index, row in df_sorted.iterrows():
        # Find the group with the smallest current total estimated time
        min_group_id = min(group_times, key=group_times.get, default=0)

        # If adding the row would exceed max_estimated_time, create a new group (if allowed)
        if group_times[min_group_id] + row['estimated_time'] > max_estimated_time and len(groups) < max_safe_group:
            min_group_id = max(groups) + 1

        # Add the row to the group
        groups[min_group_id].append(index)
        group_times[min_group_id] += row['estimated_time']

    # Assign Comp_Group_ID based on the groups
    df_sorted['Comp_Group_ID'] = -1
    for group_id, group_indices in groups.items():
        for index in group_indices:
            df_sorted.at[index, 'Comp_Group_ID'] = group_id

    return df_sorted

def get_unique_group_ids(df_with_groups):
    unique_group_ids = df_with_groups['Comp_Group_ID'].unique().tolist()
    return unique_group_ids

def copy_files_and_create_csv_by_group(df_with_groups, cur_dir, unique_group_ids, file_path):
    
        
    # Create the "Calculations" directory inside the current directory
    calculations_directory = os.path.join(cur_dir, 'Calculations')
    os.makedirs(calculations_directory, exist_ok=True)

    if file_path==file_path_data:
        # List all files in the  directory (excluding CSV files)
        source_directory = os.path.join(cur_dir, 'Hoomd')
        filenames = [f for f in os.listdir(source_directory) if os.path.isfile(os.path.join(source_directory, f)) and not f.endswith('.csv')]
    else:
        source_directory = os.path.join(cur_dir, 'Openmm')
        filenames = [f for f in os.listdir(source_directory) if os.path.isfile(os.path.join(source_directory, f)) and not f.endswith('.csv')]

    for group_id in unique_group_ids:
        # Create a subdirectory for the group inside the "Calculations" directory
        group_directory = os.path.join(calculations_directory, str(group_id))
        os.makedirs(group_directory, exist_ok=True)

        # Get the rows corresponding to the group
        group_rows = df_with_groups[df_with_groups['Comp_Group_ID'] == group_id]

        # Copy the files to the group subdirectory
        for filename in filenames:
            source_path = os.path.join(source_directory, filename)
            destination_path = os.path.join(group_directory, filename)
            shutil.copy2(source_path, destination_path)

        # Create a new CSV file for the group inside the group subdirectory
        group_csv_path = os.path.join(group_directory, f'{file_path}')
        group_rows.to_csv(group_csv_path, index=False)

def update_run_sh(group_ids, cur_dir, over_guess_percent, df_with_groups):
    for group_id in group_ids:
        # Get the rows corresponding to the group
        group_rows = df_with_groups[df_with_groups['Comp_Group_ID'] == group_id]

        # Calculate the estimated group time in minutes using your equation
        estimated_group_time_minutes = group_rows['estimated_time'].sum()

        # Add the specified percentage extra to the time
        total_minutes = estimated_group_time_minutes * (1 + float(over_guess_percent) / 100)
        hours = int(total_minutes // 60)
        minutes = int(total_minutes % 60)

        # Prepare the new time string (seconds will be set to 00)
        new_time_string = f"{hours:02}:{minutes:02}:00"

        # Construct the path to the run.sh file for this group
        file_path = os.path.join(cur_dir, 'Calculations', str(group_id), 'run.sh')

        # Read the original file
        with open(file_path, 'r') as file:
            lines = file.readlines()
        file.close()

        # Replace Windows line endings with Unix line endings
        lines = [line.replace('\r\n', '\n') for line in lines] # Cleans up code that was copied in Windows OS

        # Update the time line
        for i, line in enumerate(lines):
            if line.startswith("#SBATCH --time"):
                lines[i] = f"#SBATCH --time={new_time_string}\n"
                break

        # Write the updated lines back to the file
        with open(file_path, 'w') as file:
            file.writelines(lines)
        file.close()

def jump_start(group_ids, cur_dir):
    # Iterate through each group ID directory
    for group_id in group_ids:
        group_directory = os.path.join(cur_dir, 'Calculations', str(group_id))

        # Change the working directory to the group directory
        os.chdir(group_directory)

        # Execute the specified shell commands
        print(f"Starting group {group_id}...")
        execute_shell_command("conda deactivate")
        execute_shell_command("module load anaconda/anaconda3 && source activate hoomd && conda activate uber_env && module load anaconda/anaconda3 && python3 Initialize.py")
        print(f"Group {group_id} initialization complete.")

    # Return to the original working directory
    os.chdir(cur_dir)
    job_id_gather(group_ids, cur_dir)
    
def job_id_gather(group_ids, cur_dir):
    time.sleep(150)
    job_list=[]
    # Iterate through each group ID directory
    
    for group_id in group_ids:
        group_directory = os.path.join(cur_dir, 'Calculations', str(group_id))

        # Change the working directory to the group directory
        os.chdir(group_directory)

        # Execute the specified shell commands
        print(f"Gathering Job Number of Group ID: {group_id}...")
        temp=find_stokes_num_from_hoomd_file()
        job_list.append(temp)
    # Return to the original working directory
    os.chdir(cur_dir)
    file_from_list(job_list,'Job_list.txt')

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
    
def find_stokes_num_from_job_file():
    """
    Scans through the current directory for files starting with 'job', strips the string 'job-' and '.out' from the file name and converts the resultant string to an integer.
    Returns this integer representing the job number.
    """
    for file in os.listdir(os.getcwd()):  # Iterate over all files in the current directory
        if file.startswith("job"):  # If file name starts with 'job'
            file = file.replace('job-', '').replace('.out', '').strip()  # Remove 'job-' from the file name  # Remove '.out' from the file name
            job_number = int(file)  # Convert file name to integer
    return job_number    
    
def file_from_list(temp_list,name_of_file):
    with open(name_of_file, 'w') as file:
        # Iterate through the list and write each element to the file
        for item in temp_list:
            temp=f"{item}\n"
            file.write(temp)
    file.close()

def display_group_estimated_times(df_with_groups, unique_group_ids):
    print("Estimated Total Time for Each Group:")
    for group_id in unique_group_ids:
        # Get the rows corresponding to the group
        group_rows = df_with_groups[df_with_groups['Comp_Group_ID'] == group_id]
        
        # Calculate the estimated total time for the group
        estimated_total_time_minutes = group_rows['estimated_time'].sum()

        # Convert to hours and minutes
        hours = int(estimated_total_time_minutes // 60)
        minutes = int(estimated_total_time_minutes % 60)
        
        print(f"Group ID: {group_id} - Estimated Time: {hours} hours, {minutes} minutes")

def data_gather(file_path, check_dir_cal):
    """
    Gathers data from the specified directory, looks for files matching the file_path variable,
    and creates a large dataset by combining all the data with similar headers.
    Then, it creates a new CSV file in the original directory with the name of the file_path variable.
    """
    # List all the directories in check_dir_cal
    sub_dirs = [os.path.join(check_dir_cal, d) for d in os.listdir(check_dir_cal) if os.path.isdir(os.path.join(check_dir_cal, d))]
    
    # Initialize an empty list to hold all the dataframes
    dataframes = []

    # Loop through each subdirectory
    for sub_dir in sub_dirs:
        # Look for files named 'merged_config_stat_results.csv' in each subdirectory
        file_pattern = os.path.join(sub_dir, 'merged_config_stat_results.csv')
        files = glob.glob(file_pattern)
            
        # Process only the first file if any files were found
        if files:
            df = pd.read_csv(files[0])
            dataframes.append(df)
                
    # Concatenate all the dataframes together
    merged_data = pd.concat(dataframes, ignore_index=True)
    
    # Write it out to a new CSV file in the original directory
    output_file = os.path.join(cur_dir, file_path)
    merged_data.to_csv(output_file, index=False)
    print(f"Data gathered and saved as '{output_file}'")
    
    # Copy the CSV file to the check_dir_cal directory with the same name
    shutil.copy(output_file, check_dir_cal)
    print(f"Copied '{output_file}' to '{check_dir_cal}'")
    
    # Copy the 'Job_list.txt' file to the check_dir_cal directory
    job_list_file = os.path.join(cur_dir, 'Job_list.txt')
    shutil.copy(job_list_file, check_dir_cal)
    print(f"Copied '{job_list_file}' to '{check_dir_cal}'")

def packing_up():
    """
    Packages up the 'completed_run' directory into a zip file. The zip file is named with the Job number (hoomd number) and the date-time when the packing is done. The hoomd file is removed before packing.
    """
    cur_dt = datetime.datetime.now()  # Get the current date and time
    formatted_dt = cur_dt.strftime("%Y-%m-%d_%H-%M")  # Formatting date and time to the accepted naming convention
    stokes_num = find_stokes_num_from_job_file()  # Get Stokes Job Number 
    zip_file_name = f'completed_run_Job({stokes_num})_TimeFinished({formatted_dt})'  # Define the zip file name
    zip_file_path = os.path.join(cur_dir, zip_file_name)  # Define the zip file path
    shutil.make_archive(zip_file_path, 'zip', check_dir_cal)  # Create the zip archive
    print(f"Zipped directory '{check_dir_cal}' and saved as '{zip_file_path}'")  # Inform the user of the completed task

    # Delete the check_dir_cal directory
    shutil.rmtree(check_dir_cal)
    print(f"Deleted directory '{check_dir_cal}'")

    # Delete any file in cur_dir that ends with .out
    out_files = glob.glob(os.path.join(cur_dir, '*.out'))
    for file in out_files:
        os.remove(file)
        print(f"Deleted file '{file}'")
    return stokes_num

def create_directory(directory, name):
    """
    Creates a new directory with the specified name in the specified directory if it doesn't already exist.

    :param directory: The directory where the new directory should be created
    :param name: The name of the new directory
    """
    # Join the directory and name to form the name of the new directory
    new_dir = os.path.join(directory, name)
    
    # Check if the new directory already exists
    if not os.path.exists(new_dir):
        # If it doesn't exist, create it
        os.makedirs(new_dir)
        print(f"Created directory '{new_dir}'")
    else:
        print(f"Directory '{new_dir}' already exists")
    return new_dir

def migrate(stokes_num,quick_data_dir):
    # Create a temporary path for the renamed file in the same directory as old_path
    cur_dt = datetime.datetime.now()  # Get the current date and time
    formatted_dt = cur_dt.strftime("%Y-%m-%d_%H-%M")  # Formatting date and time to the accepted naming convention
    new_name=f'{stokes_num}_TimeFinished({formatted_dt}).csv'
    new_file_path = os.path.join(quick_data_dir, new_name)
    old_file_path = os.path.join(cur_dir, file_path)
    os.rename(old_file_path, new_file_path)
    
def decide(file_path,in_process):
    if in_process:
        data_gather(file_path,check_dir_cal)
        stokes_num=packing_up()
        zip_files = glob.glob(os.path.join(cur_dir, '*.zip'))
        if not zip_files:
            print("No zip file found!")
            exit()
        zip_file_path = zip_files[0]
    
    
        if file_path==file_path_data:
            # If running Hoomd
            print("Running Hoomd")
            hoomd_jobs_dir=create_directory(check_dir_archive,'Hoomd_Jobs')
            shutil.move(zip_file_path, hoomd_jobs_dir)
            print(f"Moved '{zip_file_path}' to '{hoomd_jobs_dir}'")
            quick_data_dir=create_directory(hoomd_jobs_dir,'Quick_Data')
            migrate(stokes_num,quick_data_dir)
            
        elif file_path==file_path_m_data:
            # If running OpenMM
            print("Running OpenMM")
            openmm_jobs_dir=create_directory(check_dir_archive,'OpenMM_Jobs')
            shutil.move(zip_file_path, openmm_jobs_dir)
            print(f"Moved '{zip_file_path}' to '{openmm_jobs_dir}'")
            quick_data_dir=create_directory(openmm_jobs_dir,'Quick_Data')
            migrate(stokes_num,quick_data_dir)
            
        else:
            print("No File Found!")
            exit()
    else:
        df_sorted, total_time, average_time, max_time, num_rows = process_csv(file_path)
        print(f"Total Estimated Time: {total_time}")
        print(f"Average Estimated Time: {average_time}")
        print(f"Number of Rows: {num_rows}")
        print(df_sorted.head())
        df_grouped=create_optimized_groups(df_sorted, max_time, max_safe_group)
        group_ids=get_unique_group_ids(df_grouped)
        copy_files_and_create_csv_by_group(df_grouped,cur_dir, group_ids,file_path)
        update_run_sh(group_ids, cur_dir, over_guess_percent, df_grouped)
        display_group_estimated_times(df_grouped, group_ids)
        jump_start(group_ids, cur_dir)




decide(file_path,in_process)