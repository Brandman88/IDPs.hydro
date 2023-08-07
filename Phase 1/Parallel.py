import pandas as pd
import random
import os
import shutil


max_safe_group=100
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
        
def your_equation(number_of_letters):
    # Your specific equation
    estimated_time = ((2.8912195151139 * 10**-7) * number_of_letters**4.14262) + 10
    return estimated_time

def process_csv(file_path):
    # Read the CSV file
    df = pd.read_csv(file_path)

    # Create a temporary column with the number of letters in "Sequence"
    df['letters_count'] = df['Sequence'].apply(len)

    # Apply your equation to the temporary column to calculate the estimated time
    df['estimated_time'] = df['letters_count'].apply(your_equation)

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

def create_optimized_groups(df_sorted, max_estimated_time, max_safe_group, iterations=10):
    best_groups = None
    best_variance = float('inf')

    for _ in range(iterations):
        # Shuffle the DataFrame
        df_shuffled = df_sorted.sample(frac=1).reset_index(drop=True)

        # Initialize groups and group times
        groups = []
        group_times = []

        # Greedy algorithm with random sampling
        for index, row in df_shuffled.iterrows():
            added_to_group = False

            # Check existing groups to see if the row can be added
            for i, group_time in enumerate(group_times):
                if group_time + row['estimated_time'] <= max_estimated_time:
                    groups[i].append(index)
                    group_times[i] += row['estimated_time']
                    added_to_group = True
                    break

            # If the row was not added to an existing group, create a new group (if allowed)
            if not added_to_group and len(groups) < max_safe_group:
                groups.append([index])
                group_times.append(row['estimated_time'])

        # Check if the current grouping has a better balance (lower variance)
        current_variance = pd.Series(group_times).var()
        if current_variance < best_variance:
            best_variance = current_variance
            best_groups = groups

    # Assign Comp_Group_ID based on the best groups
    df_sorted['Comp_Group_ID'] = 0
    for group_id, group_indices in enumerate(best_groups):
        for index in group_indices:
            df_sorted.at[index, 'Comp_Group_ID'] = group_id

    return df_sorted

def get_unique_group_ids(df_with_groups):
    unique_group_ids = df_with_groups['Comp_Group_ID'].unique().tolist()
    return unique_group_ids

def copy_files_and_create_csv_by_group(df_with_groups, source_directory, destination_directory, unique_group_ids,file_path):
    # Create the "Calculations" directory inside the destination
    calculations_directory = os.path.join(destination_directory, 'Calculations')
    os.makedirs(calculations_directory, exist_ok=True)

    for group_id in unique_group_ids:
        # Create a subdirectory for the group inside the "Calculations" directory
        group_directory = os.path.join(calculations_directory, str(group_id))
        os.makedirs(group_directory, exist_ok=True)

        # Get the rows and filenames corresponding to the group
        group_rows = df_with_groups[df_with_groups['Comp_Group_ID'] == group_id]
        filenames = group_rows['filename'].tolist()

        # Copy the files to the group subdirectory (excluding CSV files)
        for filename in filenames:
            if not filename.endswith('.csv'):
                source_path = os.path.join(source_directory, filename)
                destination_path = os.path.join(group_directory, filename)
                shutil.copy2(source_path, destination_path)

        # Create a new CSV file for the group inside the group subdirectory
        group_csv_path = os.path.join(group_directory, f'{file_path}')
        group_rows.to_csv(group_csv_path, index=False)

def update_run_sh(group_ids, estimated_group_times, cur_dir,over_guess_percent):
    for group_id, estimated_group_time_minutes in zip(group_ids, estimated_group_times):
        # Add 10% to the time
        total_minutes = estimated_group_time_minutes * (1+float(float(over_guess_percent)/100))
        hours = int(total_minutes // 60)
        minutes = int(total_minutes % 60)

        # Prepare the new time string (seconds will be set to 00)
        new_time_string = f"{hours:02}:{minutes:02}:00"

        # Construct the path to the run.sh file for this group
        file_path = os.path.join(cur_dir, 'Calculations', str(group_id), 'run.sh')

        # Read the original file
        with open(file_path, 'r') as file:
            lines = file.readlines()

        # Replace Windows line endings with Unix line endings
        lines = [line.replace('\r\n', '\n') for line in lines] #Cleans up code that was copied in Windows OS

        # Update the time line
        for i, line in enumerate(lines):
            if line.startswith("#SBATCH --time"):
                lines[i] = f"#SBATCH --time={new_time_string}\n"
                break

        # Write the updated lines back to the file
        with open(file_path, 'w') as file:
            file.writelines(lines)






file_path = "path/to/yourfile.csv"
df_sorted, total_time, average_time, num_rows = process_csv(file_path)
print(f"Total Estimated Time: {total_time}")
print(f"Average Estimated Time: {average_time}")
print(f"Number of Rows: {num_rows}")
print(df_sorted.head())