import subprocess
import time
import shutil
import os 

parameters='multi_run.dat'
#num_run is first number, says how many runs are we doing?
#start is the second number, says where we are in that list
#marker is are we at the start or end of a task. default is S for starting

dest_file='1LC.txt'
dest_file_debug="debug_file.txt"

cur_dir=os.getcwd()
dir_to_run=f"{cur_dir}/to_run"
dir_comp=f"{cur_dir}/completed_run"

def read_multi_run(parameters='multi_run.dat'):
    with open(parameters,'r') as para:
        lines = []
        for line in para:
            lines.append(line)
    para.close()
    
    num_run_raw=lines[0]
    start_raw=lines[1]
    marker_raw=lines[2]
    
    if '#Number of proteins you care to run' in num_run_raw:
        num_run=int(num_run_raw.replace('#Number of proteins you care to run','').strip())
    else:
        num_run=int(num_run_raw)
    if '#Start value (if you are just starting the process should be on 1 by default)' in start_raw:
        start=int(start_raw.replace('#Start value (if you are just starting the process should be on 1 by default)','').strip())
    else:
        start=int(start_raw)
    if '#Put S or E based on starting or end process (don\'t really need to touch), if its I then set the other 2 numbers and run the file again.' in marker_raw:
        marker=marker_raw.replace('#Put S or E based on starting or end process (don\'t really need to touch), if its I then set the other 2 numbers and run the file again.','')
    else:
        marker=marker_raw
    marker=marker.strip()

    lines.remove(num_run_raw)
    lines.remove(start_raw)
    lines.remove(marker_raw)
  
    clean_list=[]
  
    for line in lines:
        line=line.strip()
        if line.strip():
            if '\n' in line:
                line=line.replace('\n','').strip()
                clean_list.append(line)
            else:
                clean_list.append(line.strip())
                
    return num_run,start,marker,clean_list

def check_dest_file_debug(dest_file_debug="debug_file.txt"):
    if not os.path.exists(dest_file_debug):
        f = open(dest_file_debug, 'x')
        f.close()

check_dest_file_debug()

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

stokes_num=find_stokes_num_from_hoomd_file()

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

def execute_shell_command_no_communication(command):
    try:
        # Execute the shell command as a subprocess
        process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)

        # Read the output line by line, but don't do anything with it
        for line in process.stdout:
            pass  # No operation

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


# Define the commands
commands = [
    'module load anaconda/anaconda3 && source activate hoomd && source /home/sbrandon/my-envs/uber_env/bin/activate && python3 organize.py',
    'module load anaconda/anaconda3 && source activate hoomd && source /home/sbrandon/my-envs/uber_env/bin/activate && python3 ProteinBulk.py',
    'module load anaconda/anaconda3 && source activate hoomd && source /home/sbrandon/my-envs/uber_env/bin/activate && python3 long_config_analysis.py',
    'module load anaconda/anaconda3 && source activate hoomd && source /home/sbrandon/my-envs/uber_env/bin/activate && python3 cm_fix.py',
    'module load anaconda/anaconda3 && source activate hoomd && source /home/sbrandon/my-envs/uber_env/bin/activate && python3 organize.py'
]

num_run, start, marker, list_files = read_multi_run()

while start <= num_run:
    num_run, start, marker, list_files = read_multi_run()
    start_time = time.time()
    if num_run==1:
        for command in commands:
            execute_shell_command(command, output_file="debug_file.txt")
        break
    else:    
        for command in commands:
            execute_shell_command(command, output_file="debug_file.txt")
        
        if start==num_run:#goes into this after finishing commands so this tells us if we finished and leaves beofre attempting to read again
            break
execute_shell_command_no_communication(f'scancel {stokes_num}')