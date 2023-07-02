import datetime
import shutil
import os
import numpy as np
import pandas as pd
import glob
import matplotlib.pyplot as plt
import re
import statistics
import ast
from matplotlib.lines import Line2D
from math import *
from random import *
from scipy.stats import pearsonr
import inspect
from importlib.util import spec_from_file_location, module_from_spec



def get_definitions():
    definitions = {}
    current_module = inspect.currentframe().f_globals
    for name, obj in current_module.items():
        if inspect.isfunction(obj):
            definitions[name] = obj
    return definitions

def get_definitions_useful():
    definitions = {}
    current_module = inspect.currentframe().f_globals
    for name, obj in current_module.items():
        if inspect.isfunction(obj):
            definitions[name] = obj
    del definitions['pearsonr']
    del definitions['spec_from_file_location']
    del definitions['module_from_spec']
    del definitions['get_definitions']
    del definitions['extract_number']
    del definitions['ask_user_actions']
    del definitions['get_definitions_useful']
    del definitions['show_dictionary_keys_from_csv_choose']
    del definitions['show_dictionary_keys_from_csv']
    del definitions['list_csv_files_in_directory_choose']
    del definitions['list_py_files_in_directory_choose']
    del definitions['list_dat_files_in_directory_choose']
    del definitions['get_common_dictionary_keys_choose']
    del definitions['get_common_dictionary_keys']
    return definitions

def get_common_dictionary_keys_choose(file1, file2, name):
    # Read the data from the CSV files
    data1 = pd.read_csv(file1)
    data2 = pd.read_csv(file2)

    # Get the dictionary keys from both files
    keys1 = set(data1.columns)
    keys2 = set(data2.columns)

    # Find the common keys present in both files
    common_keys = list(keys1.intersection(keys2))

    if not common_keys:
        print("No common dictionary keys found in the CSV files.")
        return None

    print(f"Possbile {name}:")
    for index, key in enumerate(common_keys, start=1):
        print(f"{index}. {key}")

    # Prompt the user to choose a key
    while True:
        key_choice = input(f"Enter the number corresponding to the {name} you want to choose: ")

        if not key_choice.isdigit() or int(key_choice) not in range(1, len(common_keys) + 1):
            print("Invalid choice. Please try again.")
        else:
            chosen_key = common_keys[int(key_choice) - 1]
            print(f"Chosen key: {chosen_key}")
            break

    return chosen_key

def get_common_dictionary_keys(file1, file2):
    # Read the data from the CSV files
    data1 = pd.read_csv(file1)
    data2 = pd.read_csv(file2)

    # Get the dictionary keys from both files
    keys1 = set(data1.columns)
    keys2 = set(data2.columns)

    # Find the common keys present in both files
    common_keys = list(keys1.intersection(keys2))

    if not common_keys:
        print("No common dictionary keys found in the CSV files.")
        return None

    return common_keys

def get_definitions_other_file():
    name=list_py_files_in_directory_choose('getting list of python file definitions')
    spec = spec_from_file_location('module_name', name)
    module = module_from_spec(spec)
    spec.loader.exec_module(module)

    definitions = {}
    for name, obj in inspect.getmembers(module):
        if inspect.isfunction(obj):
            definitions[name] = obj
            
    print("Definitions:")
    for name in definitions:
        print(name, flush=True)
    return definitions

def show_dictionary_keys_from_csv_choose(file_path,name):
    # Read the data from the CSV file
    data_df = pd.read_csv(file_path)

    # Get the column names as dictionary keys
    keys = list(data_df.columns)

    if not keys:
        print("No dictionary keys found in the CSV file.")
        return

    for index, key in enumerate(keys, start=1):
        print(f"{index}. {key}")

    # Prompt the user to choose a key
    while True:
        key_choice = input(f"Choose for a number for {name} input: ")

        if not key_choice.isdigit() or int(key_choice) not in range(1, len(keys) + 1):
            print("Invalid choice. Please try again.")
        else:
            chosen_key = keys[int(key_choice) - 1]
            print(f"Chosen key: {chosen_key}")
            break

    return chosen_key

def show_dictionary_keys_from_csv(file_path):
    # Read the data from the CSV file
    data_df = pd.read_csv(file_path)

    # Get the column names as dictionary keys
    keys = list(data_df.columns)

    if not keys:
        print("No dictionary keys found in the CSV file.")
        return None

    return keys

def list_csv_files_in_directory_choose(name):
    # Get the current directory
    current_dir = os.getcwd()

    # Get all files in the directory ending with ".csv"
    csv_files = [file for file in os.listdir(current_dir) if file.endswith(".csv")]

    if not csv_files:
        print("No CSV files found in the current directory.")
        return None

    print("CSV files:")
    for index, file_name in enumerate(csv_files, start=1):
        print(f"{index}. {file_name}")

    # Prompt the user to choose a file
    while True:
        file_choice = input(f"Enter the number corresponding to the file you want to choose for {name}: ")

        if not file_choice.isdigit() or int(file_choice) not in range(1, len(csv_files) + 1):
            print("Invalid choice. Please try again.")
        else:
            chosen_file = csv_files[int(file_choice) - 1]
            print(f"Chosen file: {chosen_file}")
            break

    return chosen_file

def list_py_files_in_directory_choose(name):
    # Get the current directory
    current_dir = os.getcwd()

    # Get all files in the directory ending with ".py"
    py_files = [file for file in os.listdir(current_dir) if file.endswith(".py")]

    if not py_files:
        print("No Python files found in the current directory.")
        return None

    print("Python files:")
    for index, file_name in enumerate(py_files, start=1):
        print(f"{index}. {file_name}")

    # Prompt the user to choose a file
    while True:
        file_choice = input(f"Enter the number corresponding to the file you want to choose for {name}: ")

        if not file_choice.isdigit() or int(file_choice) not in range(1, len(py_files) + 1):
            print("Invalid choice. Please try again.")
        else:
            chosen_file = py_files[int(file_choice) - 1]
            print(f"Chosen file: {chosen_file}")
            break

    return chosen_file

def list_dat_files_in_directory_choose():
    # Get the current directory
    current_dir = os.getcwd()

    # Get all files in the directory ending with ".dat"
    dat_files = [file for file in os.listdir(current_dir) if file.endswith(".dat")]

    if not dat_files:
        print("No .dat files found in the current directory.")
        return None

    print(".dat files:")
    for index, file_name in enumerate(dat_files, start=1):
        print(f"{index}. {file_name}")

    # Prompt the user to choose a file
    while True:
        file_choice = input(f"Enter the number corresponding to the file you want to choose: ")

        if not file_choice.isdigit() or int(file_choice) not in range(1, len(dat_files) + 1):
            print("Invalid choice. Please try again.")
        else:
            chosen_file = dat_files[int(file_choice) - 1]
            print(f"Chosen file: {chosen_file}")
            break

    return chosen_file

def extract_number(s):
    # Find the first number in the string
    numbers = re.findall('\d+', s)
    # Return the first number as an integer, or return 0 if there are no numbers
    return int(numbers[0]) if numbers else 0

def plot_something_2locations_relative_to_one_variable_basic():
    loc1=list_csv_files_in_directory_choose('Location 1(dependent)')
    loc2=list_csv_files_in_directory_choose('Location 2(independent)')
    variable=get_common_dictionary_keys_choose(loc1,loc2,'Linking Variable between data Sets')
    # Read the data from the CSV files
    data_df = pd.read_csv(loc1)
    config_df = pd.read_csv(loc2)  # Add subdirectory to file path
    dependent=show_dictionary_keys_from_csv_choose(loc1,'dependent')
    independent=show_dictionary_keys_from_csv_choose(loc2,'independent')
    # Check if the dependent and independent columns exist in the data DataFrame
    independent_columns = [col for col in config_df.columns if col.lower() == independent.lower()]
    dependent_columns = [col for col in data_df.columns if col.lower() == dependent.lower()]

    if not dependent_columns:
        print(f"Column '{dependent}' does not exist in the data CSV file.")
        return
    if not independent_columns:
        print(f"Column '{independent}' does not exist in the data CSV file.")
        return

    dependent_column = dependent_columns[0]
    independent_column = independent_columns[0]

    # Merge the data based on matching names
    merged_df = pd.merge(data_df, config_df, on=variable)

    # Plot the data
    plt.scatter(merged_df[independent_column], merged_df[dependent_column])
    plt.xlabel(independent_column)
    plt.ylabel(dependent_column)
    plot_title = f'{dependent_column} vs {independent_column}'
    plt.title(plot_title)
    
    safe_filename = f'{dependent} vs {independent}.svg'.replace("<", "").replace(">", "").replace(":", "")

    # Save the plot
    plt.savefig(safe_filename)
    plt.close()

def plot_something_2locations_relative_to_one_variable_trend():
    #having 2 locations if I try to read k value and its in both it changes the end to k value_y or _x depending on independent or dependent value 
    #if really care can make exception for it later on to check if its  a linking variable
    loc1=list_csv_files_in_directory_choose('Location 1(dependent)')
    loc2=list_csv_files_in_directory_choose('Location 2(independent)')
    variable=get_common_dictionary_keys_choose(loc1,loc2,'Linking Variable between data Sets')
    # Read the data from the CSV files
    data_df = pd.read_csv(loc1)  # independent side
    config_df = pd.read_csv(loc2)  # dependent side
    dependent=show_dictionary_keys_from_csv_choose(loc1,'dependent')
    independent=show_dictionary_keys_from_csv_choose(loc2,'independent')

    # Merge the data based on matching names
    merged_df = pd.merge(data_df, config_df, on=variable)
    print(merged_df.columns)
    x = merged_df[independent]
    y = merged_df[dependent]

    # Calculate the best-fit line
    slope, intercept = np.polyfit(x, y, 1)
    trendline = slope * x + intercept

    # Calculate the correlation coefficient (R value)
    r, _ = pearsonr(x, y)

    # Plot the data
    plt.scatter(x, y, label='Data')
    plt.plot(x, trendline, color='red', label=f'Trendline: y={slope:.2f}x+{intercept:.2f}, R={r:.2f}')

    plt.xlabel(independent)
    plt.ylabel(dependent)
    plt.title(f'{dependent} vs {independent}')
    plt.legend()

    safe_filename = f'{dependent} vs {independent}.svg'.replace("<", "").replace(">", "").replace(":", "")

    # Save the plot
    plt.savefig(safe_filename)
    plt.close()

def plot_something_location_relative_to_one_variable_basic():
    loc1=list_csv_files_in_directory_choose('Location of csv file')
    data_df = pd.read_csv(loc1)
    dependent=show_dictionary_keys_from_csv_choose(loc1,'dependent')
    independent=show_dictionary_keys_from_csv_choose(loc1,'independent')
    # Check if the dependent and independent columns exist in the DataFrame
    dependent_columns = [col for col in data_df.columns if col.lower() == dependent.lower()]
    independent_columns = [col for col in data_df.columns if col.lower() == independent.lower()]

    if not dependent_columns:
        print(f"Column '{dependent}' does not exist in the CSV file.")
        return
    if not independent_columns:
        print(f"Column '{independent}' does not exist in the CSV file.")
        return

    dependent_column = dependent_columns[0]
    independent_column = independent_columns[0]

    # Plot the data
    plt.scatter(data_df[independent_column], data_df[dependent_column])
    plt.xlabel(independent_column)
    plt.ylabel(dependent_column)
    plot_title = f'{dependent_column} vs {independent_column}'
    plt.title(plot_title)

    safe_filename = f'{dependent} vs {independent}.svg'.replace("<", "").replace(">", "").replace(":", "")

    # Save the plot
    plt.savefig(safe_filename)
    plt.close()

def plot_something_location_relative_to_one_variable_trend():
    loc1=list_csv_files_in_directory_choose('Location of csv file')
    dependent=show_dictionary_keys_from_csv_choose(loc1,'dependent')
    independent=show_dictionary_keys_from_csv_choose(loc1,'independent')
    
    # Read the data from the CSV file
    data_df = pd.read_csv(loc1)

    # Check if the dependent and independent columns exist in the DataFrame
    dependent_columns = [col for col in data_df.columns if col.lower() == dependent.lower()]
    independent_columns = [col for col in data_df.columns if col.lower() == independent.lower()]

    if not dependent_columns:
        print(f"Column '{dependent}' does not exist in the CSV file.")
        return
    if not independent_columns:
        print(f"Column '{independent}' does not exist in the CSV file.")
        return

    dependent_column = dependent_columns[0]
    independent_column = independent_columns[0]

    x = data_df[independent_column]
    y = data_df[dependent_column]

    # Calculate the best-fit line
    slope, intercept = np.polyfit(x, y, 1)
    trendline = slope * x + intercept

    # Calculate the correlation coefficient (R value)
    r, _ = pearsonr(x, y)

    # Plot the data
    plt.scatter(x, y, label='Data')
    plt.plot(x, trendline, color='red', label=f'Trendline: y={slope:.2f}x+{intercept:.2f}, R={r:.2f}')

    plt.xlabel(independent_column)
    plt.ylabel(dependent_column)
    plt.title(f'{dependent_column} vs {independent_column}')
    plt.legend()
    
    safe_filename = f'{dependent} vs {independent}_trend.svg'.replace("<", "").replace(">", "").replace(":", "")

    # Save the plot
    plt.savefig(safe_filename)
    plt.close()

def plot_something_many_labels():
    # Read the CSV file into a pandas DataFrame
    loc1 = list_csv_files_in_directory_choose('Location of csv file')
    dependent = show_dictionary_keys_from_csv_choose(loc1, 'Dependent')
    independent = show_dictionary_keys_from_csv_choose(loc1, 'Independent')
    label = show_dictionary_keys_from_csv_choose(loc1, 'Label')
    df = pd.read_csv(loc1)

    # Evaluate values and convert to floats
    df[dependent] = df[dependent].apply(lambda x: eval(x) if isinstance(x, str) else x).astype(float)
    df[independent] = df[independent].apply(lambda x: eval(x) if isinstance(x, str) else x).astype(float)

    # Group the data by item name
    grouped_data = df.groupby(label)

    # Create a larger figure to plot all items
    fig, ax = plt.subplots(figsize=(12, 8))

    # Generate plots for each item with connected lines
    for item_name, group in grouped_data:
        # Scatter plot
        ax.scatter(group[independent], group[dependent], label=item_name)

        # Connect data points with lines
        ax.plot(group[independent], group[dependent], label='Line for {}'.format(item_name))

    ax.set_xlabel(independent)
    ax.set_ylabel(dependent)
    legend = ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    legend.get_frame().set_linewidth(0.5)  # Set legend frame linewidth
    for text in legend.get_texts():
        text.set_fontsize('small')  # Set legend font size to small

    safe_filename = f'{dependent}_vs_{independent}_labeled_by_{label}.svg'.replace("<", "").replace(">", "").replace(":", "").replace("*", "").replace("?", "").replace("\"", "").replace("\\", "").replace("/", "").replace("|", "")
    plt.savefig(safe_filename, bbox_inches='tight')
    plt.close()

def sort_csv():
    
    csv_name=list_csv_files_in_directory_choose('Location of csv file')
    column_name=show_dictionary_keys_from_csv_choose(csv_name,'Varaibles to choose from')
    # Read the CSV file into a pandas DataFrame
    df = pd.read_csv(csv_name)


    # Validate if all entries in the column are numeric
    if not pd.to_numeric(df[column_name], errors='coerce').notnull().all():
        sorted_df = df.sort_values(by=column_name)

    # Sort the DataFrame based on the specified column
    if column_name == 'name':
        name_values = df[column_name]
        sorted_indices = sorted(range(len(name_values)), key=lambda i: extract_number(name_values[i]) if re.search(r'\d+', name_values[i]) else float('inf'))
        sorted_df = df.iloc[sorted_indices]

        # Move sorted column to the front
        sorted_df = sorted_df[[column_name] + sorted_df.columns.difference([column_name]).tolist()]

    else:
        sorted_df = df.sort_values(by=column_name)

        # Move sorted column to the front
        sorted_df = sorted_df[[column_name] + sorted_df.columns.difference([column_name]).tolist()]

    # Write the sorted DataFrame to a new CSV file
    sorted_csv_name = f"sorted_{csv_name}"
    sorted_df.to_csv(sorted_csv_name, index=False)

def ask_user_actions():
    # Get all definitions in the module
    definitions = get_definitions_useful()
    while True:
        while True:
            # Display the available actions
            print("Please choose an action:")
            for index, (name, _) in enumerate(definitions.items(), start=1):
                print(f"{index}. {name}")
            action_choice = input("Enter the number corresponding to the action: ").strip()
            if action_choice == '0':
                print("Exiting...")
                return
    
            # Validate user's action choice
            if not action_choice.isdigit() or int(action_choice) not in range(1, len(definitions) + 1):
                print("Invalid choice. Please try again.")
                break

            # Get the chosen definition
            chosen_name = list(definitions.keys())[int(action_choice) - 1]
            chosen_definition = definitions[chosen_name]

            # Get the required inputs for the chosen definition
            input_prompts = chosen_definition.__code__.co_varnames[:chosen_definition.__code__.co_argcount]
            inputs = []
            for prompt in input_prompts:
                user_input = input(f"Enter the value for '{prompt}': ").strip()
                if user_input == '0':
                    print("Exiting...")
                    return
                inputs.append(user_input)

            # Execute the chosen definition with the inputs
            exec(f"{chosen_name}(*inputs)")
 
def merge_csv_files():
    file1=list_csv_files_in_directory_choose('File Location 1')
    file2=list_csv_files_in_directory_choose('File Location 2')
    file1_name=file1.replace('.csv','')
    file2_name=file2.replace('.csv','')
    output_file=f'{file1_name}_merged-with_{file2_name}.csv'
    
    # Get the common keys from both CSV files
    common_keys = get_common_dictionary_keys(file1, file2)

    # Read the data from the CSV files
    data1 = pd.read_csv(file1)
    data2 = pd.read_csv(file2)
    
    if not common_keys:
        print("No common dictionary keys found in the CSV files. Unable to merge.")
        return
    else:
        common_keys_set = set(common_keys)
        if len(data1.columns) == len(common_keys_set) and len(data2.columns) == len(common_keys_set):
            print("Both CSV files have only common keys. Unable to merge.")
            return
        elif len(data1.columns) == len(common_keys_set) or len(data2.columns) == len(common_keys_set):
            print("One of the CSV files has only matching keys, So no point in merging")
            if len(data1.columns) == len(common_keys_set):
                print(f"{file2}, contains more keys. You should probably use that one for data.")
            else:
                print(f"{file1}, contains more keys. You should probably use that one for data.")
            return


    # Find matching rows based on the common keys
    merged_data = pd.merge(data1, data2, on=common_keys)

    if merged_data.empty:
        print("No matching rows found based on the common keys. Unable to merge.")
        return

    # Export the combined CSV file
    merged_data.to_csv(output_file, index=False)
    print(f"CSV files merged successfully. Merged data exported to '{output_file}'.")

def append_csv_files():
    file1=list_csv_files_in_directory_choose('File Location 1')
    file2=list_csv_files_in_directory_choose('File Location 2')
    file1_name=file1.replace('.csv','')
    file2_name=file2.replace('.csv','')
    output_file=f'{file1_name}_merged-with_{file2_name}.csv'
    # Read the data from the CSV files
    data1 = pd.read_csv(file1)
    data2 = pd.read_csv(file2)

    # Get the common keys from both CSV files
    common_keys = get_common_dictionary_keys(file1, file2)

    if not common_keys:
        print("No common dictionary keys found in the CSV files. Unable to append.")
        return

    # Check if all common keys exist in both files
    if not all(key in set(data1.columns) for key in common_keys) or not all(key in set(data2.columns) for key in common_keys):
        print("Not all common keys exist in both CSV files. Unable to append.")
        return

    # Check if any row in data1 has all keys matching another row in data2
    match_count = 0
    for _, row1 in data1.iterrows():
        for _, row2 in data2.iterrows():
            if all(row1[key] == row2[key] for key in common_keys):
                match_count += 1
                break

    # Append the data from both CSV files
    appended_data = pd.concat([data1, data2], ignore_index=True)
    appended_data.drop_duplicates(subset=common_keys, inplace=True)
    # Export the appended CSV file
    appended_data.to_csv(output_file, index=False)
    print(f"CSV files appended successfully. Appended data exported to '{output_file}'.")
    print(f"Number of matches found: {match_count}")

def compare_csv_column():
    file1 = list_csv_files_in_directory_choose('True value')
    file2 = list_csv_files_in_directory_choose('Experimental value')
    column_name = get_common_dictionary_keys_choose(file1, file2, 'Linking Variable between data sets')

    df1 = pd.read_csv(file1)
    df2 = pd.read_csv(file2)

    if column_name not in df1.columns or column_name not in df2.columns:
        print(f"Column '{column_name}' not found in both CSV files.")
        return

    comp_name = get_common_dictionary_keys_choose(file1, file2, 'Variable desired to compare between data sets')

    file1_list_of_column_name_entries = [str(item).strip() for item in df1[column_name].copy()]
    file1_list_of_comp_name_entries = [round(float(str(item).strip()),3) for item in df1[comp_name].copy()]
    
    file1_dictionary = dict(zip(file1_list_of_column_name_entries, file1_list_of_comp_name_entries))

    file2_list_of_column_name_entries = [str(item).strip() for item in df2[column_name].copy()]
    file2_list_of_comp_name_entries = [round(float(str(item).strip()),3) for item in df2[comp_name].copy()]
    
    file2_dictionary = dict(zip(file2_list_of_column_name_entries, file2_list_of_comp_name_entries))
    print(file1_dictionary)
    print(file2_dictionary)
    rows = []
    distance_from_true_list = []
    distance_from_true_abs_list = []

    for key in file1_list_of_column_name_entries:
        distance_from_true_list.append(float(float(file1_dictionary[key]) - float(file2_dictionary[key])))
        distance_from_true_abs_list.append(abs(float(float(file1_dictionary[key]) - float(file2_dictionary[key]))))
        print(f'{file1_dictionary[key]}   {file2_dictionary[key]}')
    print(distance_from_true_list)
    print(distance_from_true_abs_list)
    

    distance_from_true_list_size = len(distance_from_true_list)
    distance_from_true_abs_list_size = len(distance_from_true_abs_list)

    avg_distance_from_true = round(np.nanmean(distance_from_true_list),3) if distance_from_true_list_size > 0 else np.nan
    avg_distance_from_true_abs = round(np.nanmean(distance_from_true_abs_list),3) if distance_from_true_abs_list_size > 0 else np.nan
    std_distance_from_true = round(np.nanstd(distance_from_true_list),3) if distance_from_true_list_size > 0 else np.nan
    std_distance_from_true_abs = round(np.nanstd(distance_from_true_abs_list),3) if distance_from_true_abs_list_size > 0 else np.nan

    for key in file1_list_of_column_name_entries:
        if key in file1_dictionary and key in file2_dictionary:
            true_value = file1_dictionary[key]
            exp_value = file2_dictionary[key]
            distance_from_true = round(file1_dictionary[key] - file2_dictionary[key],3)
            distance_from_true_abs = round(abs(distance_from_true),3)
            distance_from_true_rel2_avg = round(distance_from_true / avg_distance_from_true_abs,3)
            distance_from_true_abs_rel2_avg = round(distance_from_true_abs / avg_distance_from_true_abs,3)
            rows.append([key, file1_dictionary[key], file2_dictionary[key], distance_from_true, distance_from_true_abs, distance_from_true_rel2_avg, distance_from_true_abs_rel2_avg])

    columns = [f'{column_name}', f'{comp_name}_file1', f'{comp_name}_file2', 'Difference', 'Absolute Value Difference', 'Difference/(AVG Difference)', '(ABS Difference)/(ABS AVG Difference)']
    df_results = pd.DataFrame(rows, columns=columns)

    safe_filename = f'{comp_name}_comparison_via_{column_name}.csv'.replace("<", "").replace(">", "").replace(":", "")
    df_results.to_csv(safe_filename, index=False)
    
    print(f'STD of Difference: {std_distance_from_true}')
    print(f'STD of ABS Difference: {std_distance_from_true_abs}')





#plot_something_location_relative_to_one_variable_trend('<Rg>','K value')
#sort_csv('config_stat_results.csv','k Value')
#get_def_other_files doesnt work great 
ask_user_actions()