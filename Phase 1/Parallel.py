import pandas as pd
import random
import os
import shutil
import subprocess
import time
import sys
from collections import defaultdict


max_safe_group=6
over_guess_percent=20 

cur_dir=os.getcwd()
dir_cal=f"{cur_dir}/Calculations"

file_path_data = "data.csv"
file_path_m_data = "data_multi.csv"


if not os.path.exists(file_path_data):#preference is given to data.csv if both exist
    if not os.path.exists(file_path_m_data):
        print(f"{file_path_data} and {file_path_m_data} don't exist!\n You need at least one of thses files.")
        exit()
    else:
        file_path=file_path_m_data
else:
    file_path=file_path_data
    
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

def process_csv(file_path):
    # Read the CSV file
    df = pd.read_csv(file_path)

    # Create a temporary column with the number of letters in "Sequence"
    df['letters_count'] = df['Sequence'].apply(len)

    # Apply your equation to the temporary column to calculate the estimated time
    df['estimated_time'] = df['letters_count'].apply(est_time_equation)

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
    time.sleep(50)
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