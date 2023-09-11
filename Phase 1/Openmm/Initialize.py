import subprocess
import time
import shutil
import os
import sys
import pandas as pd


parameters='multi_run.dat'
cur_dir=os.getcwd()

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
    
def read_entire_csv_return_dict_and_num_rows(dest_required_file_csv = "data_multi.csv"):
    # Read the CSV data
    df = pd.read_csv(dest_required_file_csv)

    # Get the number of rows
    num_rows = df.shape[0]

    # Convert the DataFrame to a list of dictionaries
    data = df.to_dict('records')

    return data, num_rows

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

def check_dest_file_debug(dest_file_debug="debug_file.txt"):
    if not os.path.exists(dest_file_debug):
        f = open(dest_file_debug, 'x')
        f.close()
                  
#Start of what actually is being ran
execute_shell_command_start("conda deactivate")
execute_shell_command("module load anaconda/anaconda3 && source activate hoomd && source $HOME/my-envs/uber_env/bin/activate && python3 organize.py")
print("Initializing process started\n")


# In the furture maybe add the ability to just certain rows ect to run (that would go here)


num_run,start,marker,clean_list=read_multi_run()
data,num_rows=read_entire_csv_return_dict_and_num_rows()
num_run=num_rows
edit_multi_run(parameters,num_run,start,marker,clean_list)


execute_shell_command_start("conda deactivate")
execute_shell_command("module load anaconda/anaconda3 && source activate hoomd && source $HOME/my-envs/uber_env/bin/activate && python3 organize.py")
execute_shell_command_start("conda deactivate")
execute_shell_command("module purge && module load anaconda/anaconda3 && job_id=$(sbatch run.sh | awk '{print $4}') && output_file=\"hoomd-$job_id.out\" && touch \"$output_file\"")

