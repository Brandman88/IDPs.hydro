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
            if name=='ask_user_actions':
                definitions['Previous Menu'] = obj
            else:
                definitions[name] = obj
    del definitions['pearsonr']
    del definitions['spec_from_file_location']
    del definitions['module_from_spec']
    del definitions['get_definitions']
    del definitions['extract_number']
    del definitions['get_definitions_useful']
    del definitions['show_dictionary_keys_from_csv_choose']
    del definitions['show_dictionary_keys_from_csv']
    del definitions['list_csv_files_in_directory_choose']
    del definitions['list_py_files_in_directory_choose']
    del definitions['list_dat_files_in_directory_choose']
    del definitions['list_directories_in_directory_choose']
    del definitions['get_common_dictionary_keys_choose']
    del definitions['get_common_dictionary_keys']
    del definitions['search_csv_sequence_return_name']
    del definitions['read_entire_csv_return_dict']
    del definitions['choose_csv_headers_as_int']
    del definitions['pick_from_list_static']
    del definitions['pick_from_list_dynamic']
    del definitions['scale_axes']
    del definitions['choose_category_return_dict']
    return definitions    

def choose_category_return_dict():
    list_of_def=get_definitions_useful()
    list_of_categories=['CSV','Graphing','All']
    chosen_cat=pick_from_list_static(list_of_categories,'Category Type')
    if chosen_cat == 1:
        new_dict = {}
        for name, definition in list_of_def.items():
            if 'csv' in name or 'Previous Menu' in name:  # Check if 'csv' is in the name
                new_dict[name] = definition
        return new_dict
    elif chosen_cat == 2:
        new_dict = {}
        for name, definition in list_of_def.items():
            if 'plot' in name or 'Previous Menu' in name:  # Check if 'plot' is in the name
                new_dict[name] = definition
        return new_dict
    elif chosen_cat == 3:
        return list_of_def    
    
def pick_from_list_static(options_list,option_type,index_start_at=1):
    option_type=option_type.title()
    # Display the available scaling options
    print(f"{option_type}:")
    for index, option in enumerate(options_list, index_start_at):
        print(f"{index}. {option}")

    while True:
        user_input = input(f"Enter the option you want for {option_type} (enter '0' to exit): ").strip()
        if user_input == '0':
            print("Exiting...")
            exit()
        elif not user_input.isdigit() or int(user_input) not in range(index_start_at, len(options_list) + index_start_at):
            print("Invalid choice. Please try again.")
        else:
            chosen_input = options_list[int(user_input) - index_start_at]
            print(f"Chosen option: {chosen_input}")
            user_input = int(user_input)
            break
    return user_input

def pick_from_list_dynamic(options_list, option_type, index_start_at=1):
    option_type = option_type.title()
    remaining_options = options_list.copy()
    chosen_options = []

    print(f"Available {option_type}:")
    for index, option in enumerate(remaining_options, index_start_at):
        print(f"{index}. {option}")

    while True:
        user_input = input(f"{option_type} (enter '0' to exit and '' to stop): ").strip()
        if user_input == '0':
            print("Exiting...")
            exit()
        elif user_input == '':
            print("Exiting selection...")
            break
        elif not user_input.isdigit() or int(user_input) not in range(index_start_at, len(remaining_options) + index_start_at):
            print("Invalid choice. Please try again.")
        else:
            chosen_input = remaining_options[int(user_input) - index_start_at]
            print(f"Chosen option: {chosen_input}")
            chosen_options.append(chosen_input)
            remaining_options.remove(chosen_input)

            if len(remaining_options) == 0:
                print("No more options available.")
                break

            print(f"Remaining {option_type}:")
            for index, option in enumerate(remaining_options, index_start_at):
                print(f"{index}. {option}")

    return chosen_options

def scale_axes(fig):
    options_list = ['Keep scaling linear', 'Change scaling to log (affects all plots)','Change scaling to log (goes through each plot)']
    options_1plot_list = ['Keep scaling linear', 'Change scaling to log']
    # Display the available scaling options
    num_plots=len(fig.get_axes())
    if num_plots ==1:
        scaling_user_selection=pick_from_list_static(options_1plot_list,'Scaling Options')
    else:
        scaling_user_selection=pick_from_list_static(options_list,'Scaling Options')
    
    
    if scaling_user_selection == 1:
        return  # No scaling
    else:
        scaling_options_list = ['Change x-axis scale to log', 'Change y-axis scale to log','Change both axes scale to log']
        if scaling_user_selection == 2:
            scaling_user_selection=pick_from_list_static(scaling_options_list,'Scaling Options')
            if scaling_user_selection == 1:
                for ax in fig.get_axes():
                    ax.set_xscale('log')
            elif scaling_user_selection == 2:
                for ax in fig.get_axes():
                    ax.set_yscale('log')
            elif scaling_user_selection == 3:
                for ax in fig.get_axes():
                    ax.set_xscale('log')
                    ax.set_yscale('log')
        elif scaling_user_selection == 3:
            for ax in fig.get_axes():
                scaling_user_selection=pick_from_list_static(scaling_options_list,f'Scaling Options for {ax}')
                if scaling_user_selection == 1:
                    ax.set_xscale('log')
                elif scaling_user_selection == 2:
                    ax.set_yscale('log')
                elif scaling_user_selection == 3:
                    ax.set_xscale('log')
                    ax.set_yscale('log')

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
    chosen_key_num=pick_from_list_static(common_keys,f'Possbile {name}')
    chosen_key=common_keys[chosen_key_num-1]
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
    chosen_key_num=pick_from_list_static(keys,f'Select {name}')
    chosen_key=keys[chosen_key_num-1]
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

    chosen_file_num=pick_from_list_static(csv_files,'CSV files')
    chosen_file=csv_files[chosen_file_num-1]
    return chosen_file

def list_py_files_in_directory_choose(name):
    # Get the current directory
    current_dir = os.getcwd()

    # Get all files in the directory ending with ".py"
    py_files = [file for file in os.listdir(current_dir) if file.endswith(".py")]

    if not py_files:
        print("No Python files found in the current directory.")
        return None

    chosen_file_num=pick_from_list_static(py_files,'Python files')
    chosen_file=py_files[chosen_file_num-1]
    return chosen_file

def list_dat_files_in_directory_choose():
    # Get the current directory
    current_dir = os.getcwd()

    # Get all files in the directory ending with ".dat"
    dat_files = [file for file in os.listdir(current_dir) if file.endswith(".dat")]

    if not dat_files:
        print("No .dat files found in the current directory.")
        return None
    
    chosen_file_num=pick_from_list_static(dat_files,'Dat files')
    chosen_file=dat_files[chosen_file_num-1]
    return chosen_file

def list_directories_in_directory_choose():
    # Get the current directory
    current_dir = os.getcwd()

    # Get all directories in the current directory
    directories = [dir_name for dir_name in os.listdir(current_dir) if os.path.isdir(os.path.join(current_dir, dir_name))]

    if not directories:
        print("No directories found in the current directory.")
        return None
    chosen_dir_num=pick_from_list_static(directories,'Directories')
    chosen_dir=directories[chosen_dir_num-1]
    return chosen_dir

def extract_number(s):
    # Find the first number in the string
    numbers = re.findall('\d+', s)
    # Return the first number as an integer, or return 0 if there are no numbers
    return int(numbers[0]) if numbers else 0

def filter_dataframe_by_group_selection(grouped_dataframe,label):
    list_of_all_groups = list(grouped_dataframe.groups.keys())
    choosing_type = ['All groups','All but certain groups','Select specific groups']
    choosing_type_num=pick_from_list_static(choosing_type,'What groups do you want in your data (filtering)')
    
    if choosing_type_num == 1:
        # Return all groups
        list_of_selected_groups=list_of_all_groups
    elif choosing_type_num  == 2:
        # Prompt the user for groups to exclude
        list_of_excluded_groups=pick_from_list_dynamic(list_of_all_groups , f'Select a {label} to be excluded from plotted data')
        # Remove the excluded groups from the list of all groups
        list_of_selected_groups = [group for group in list_of_all_groups if group not in list_of_excluded_groups]    
    elif choosing_type_num  == 3:
        # Prompt the user for specific groups
        list_of_selected_groups=pick_from_list_dynamic(list_of_all_groups , f'Select a {label} to be in plotted data')

    # Get the DataFrame for each selected group
    filtered_dataframe = pd.concat([grouped_dataframe.get_group(group) for group in list_of_selected_groups])#for some reason it needs a memory allocation that only comes from grouping doing anything to it removes that
    filtered_grouped_dataframe=filtered_dataframe.groupby(label)
    return filtered_grouped_dataframe

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

def plot_something_many_labels_typical():
    # Read the CSV file into a pandas DataFrame
    loc1 = list_csv_files_in_directory_choose('Location of csv file')
    dependent = show_dictionary_keys_from_csv_choose(loc1, 'Dependent')
    independent = show_dictionary_keys_from_csv_choose(loc1, 'Independent')
    label = show_dictionary_keys_from_csv_choose(loc1, 'Label')
    df = pd.read_csv(loc1)
    print(df)
    # Evaluate values and convert to floats
    df[dependent] = df[dependent].apply(lambda x: eval(x) if isinstance(x, str) else x).astype(float)
    df[independent] = df[independent].apply(lambda x: eval(x) if isinstance(x, str) else x).astype(float)

    # Group the data by item name
    grouped_data = df.groupby(label)
    # Convert grouped_data to a list of groups
    filtered_grouped_data=filter_dataframe_by_group_selection(grouped_data,label)
    
    # Create a larger figure to plot all items
    fig, ax = plt.subplots(figsize=(16, 12))
    scale_axes(fig)  # Use logarithmic scale for y-axis
    # Generate plots for each item with connected lines
    for item_name, group in filtered_grouped_data:
        # Sort the group by the dependent variable
        group.sort_values(by=independent, inplace=True)
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

def plot_something_many_labels_typical_and_relative():
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
    filtered_grouped_data=filter_dataframe_by_group_selection(grouped_data,label)
    # Create a larger figure to plot all items
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(16, 18), sharex=True)
    scale_axes(fig)

    # Generate plots for each item with connected lines
    for item_name, group in filtered_grouped_data:
        # Sort the group by the dependent variable
        group.sort_values(by=independent, inplace=True)

        # Scatter plot
        ax1.scatter(group[independent], group[dependent], label=item_name)

        # Connect data points with lines
        ax1.plot(group[independent], group[dependent], label='Line for {}'.format(item_name))

        # Plot relative to the first variable in the group
        relative_values = group[dependent] / group[dependent].iloc[0]
        ax2.plot(group[independent], relative_values, label=item_name)

    ax1.set_ylabel(dependent)
    ax1.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    ax2.set_xlabel(independent)
    ax2.set_ylabel('Relative Values')
    ax2.legend(loc='center left', bbox_to_anchor=(1, 0.5))

    safe_filename = f'{dependent}_vs_{independent}_labeled_by_{label}.svg'.replace("<", "").replace(">", "").replace(":", "").replace("*", "").replace("?", "").replace("\"", "").replace("\\", "").replace("/", "").replace("|", "")
    plt.savefig(safe_filename, bbox_inches='tight')
    plt.close()

def plot_something_many_labels_relative():
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
    filtered_grouped_data=filter_dataframe_by_group_selection(grouped_data,label)
    # Create a larger figure to plot all items
    fig, ax = plt.subplots(figsize=(16, 12))
    scale_axes(fig)

    # Generate plots for each item with connected lines
    for item_name, group in filtered_grouped_data:
        # Sort the group by the dependent variable
        group.sort_values(by=independent, inplace=True)
        # Plot relative to the first variable in the group
        relative_values = group[dependent] / group[dependent].iloc[0]
        ax.plot(group[independent], relative_values, label=item_name)

    ax.set_xlabel(independent)
    ax.set_ylabel('Relative Values')
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))

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

def read_entire_csv_return_dict(dest_required_file_csv):
    '''Define function to read entire CSV and return a dictionary.'''
    # Using pandas to read the CSV file located at the destination provided.
    df = pd.read_csv(dest_required_file_csv)

    # Converting the dataframe into a list of dictionaries where each dictionary represents a row of data.
    data = df.to_dict('records')

    # Return the list of dictionaries
    return data

def search_csv_sequence_return_name(file_location_protein_sequence,dest_required_file_csv):
    '''Define function to grab specific column information from a specific row in the CSV.'''
    
    file = open(file_location_protein_sequence, 'r')
    read_sequence=file.readline()
    file.close()
    csv_dictionary=read_entire_csv_return_dict(dest_required_file_csv)
    i=0
    matching_rows=[]
    for row in csv_dictionary:
        if row['Sequence']==read_sequence:
            protein_name=row['Name']
            matching_rows.append(i)
        i+=1
    return protein_name,matching_rows

def choose_csv_headers_as_int(csv_data):
    """
    Asks the user if they want any information to be treated as an integer from the CSV headers.

    :param csv_data: DataFrame containing the CSV data
    :return: List of chosen headers to be treated as integers
    """
    headers = csv_data.columns.tolist()
    chosen_headers = []

    while True:
        print("Do you want any information to be treated as an integer?")
        response = input("Enter 'yes' or 'no': ")

        if response.lower() == 'no':
            break

        if response.lower() != 'yes':
            print("Invalid response. Please enter 'yes' or 'no'.")
            continue

        print("Choose from the following CSV headers:")
        for index, header in enumerate(headers, start=1):
            print(f"{index}. {header}")

        while True:
            choice = input("Enter the number corresponding to the header you want to choose (or '0' to finish): ")

            if choice == '0':
                break

            if not choice.isdigit() or int(choice) not in range(1, len(headers) + 1):
                print("Invalid choice. Please try again.")
            else:
                chosen_header = headers[int(choice) - 1]
                chosen_headers.append(chosen_header)
                headers.remove(chosen_header)
                print(f"Chosen header: {chosen_header}")
                break

        if choice == '0':
            break

    return chosen_headers

def read_dat_files_data_merger_create_csv():
    """
    This function reads .dat files from the given directory, merges them,
    and then merges the resulting DataFrame with data from a CSV file.

    :param csv_file: Location of the CSV file to merge data with
    :param dir_comp: The directory where the .dat files are located
    """
    dest_file="1LC.txt"
    csv_file = list_csv_files_in_directory_choose('Parent csv')
    dir_comp=list_directories_in_directory_choose()
    print(dir_comp)
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
        protein_name, matching_rows = search_csv_sequence_return_name(seq_file_loc,csv_file)

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
        df_subset = df[list(non_matching_headers)].copy()  # Select only the non-matching columns from df
        df_subset.index = [matching_rows[matching_array_value]]  # Wrap the value in a list to create a single-item collection

        # Format the values as floats with a maximum of 3 decimal places
        df_subset = df_subset.applymap(lambda x: int(x) if round(float(x), 3).is_integer() else round(float(x), 3))
        dataframes.append(df_subset)  # Append the subset DataFrame to the list of dataframes

    # Concatenate all the dataframes together
    merged_data = pd.concat(dataframes)

    # Merge with the original csv_data
    merged_data = pd.merge(csv_data, merged_data, left_index=True, right_index=True, how='left')
    int_list=choose_csv_headers_as_int(merged_data)
    for item in int_list:
        merged_data[item] = merged_data[item].fillna(0).astype(int)
    # Write it out
    merged_data.to_csv('merged_config_stat_results.csv', index=False)

def ask_user_actions():
    # Get all definitions in the module
    definitions = choose_category_return_dict()
    while True:
        while True:
            # Display the available actions
            print("Please choose an action:")
            for index, (name, _) in enumerate(definitions.items(), start=1):
                print(f"{index}. {name}")
            action_choice = input("Enter the number corresponding to the action (enter '0' to exit): ").strip()
            if action_choice == '0':
                print("Exiting...")
                exit()
            # Validate user's action choice
            elif not action_choice.isdigit() or int(action_choice) not in range(1, len(definitions) + 1):
                if action_choice.lower() == 'back':
                    print("Going back to the previous menu...")
                    definitions = choose_category_return_dict()
                    break
                else:
                    print("Invalid choice. Please try again.")
                    break
            elif list(definitions.keys())[int(action_choice) - 1]=='Previous Menu':
                print("Going back to the previous menu...")
                definitions = choose_category_return_dict()
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


#plot_something_location_relative_to_one_variable_trend('<Rg>','K value')
#sort_csv('config_stat_results.csv','k Value')
#get_def_other_files doesnt work great 
ask_user_actions()