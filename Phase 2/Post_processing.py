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
from bokeh.plotting import figure, show
from bokeh.models import ColumnDataSource, HoverTool
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.patches import Rectangle


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
    del definitions['list_out_files_in_directory_choose']
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
    del definitions['aggregate_group']
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

    chosen_file_num=pick_from_list_static(py_files,f'{name}')
    chosen_file=py_files[chosen_file_num-1]
    return chosen_file

def list_out_files_in_directory_choose(name):
    # Get the current directory
    current_dir = os.getcwd()

    # Get all files in the directory ending with ".py"
    out_files = [file for file in os.listdir(current_dir) if file.endswith(".out")]

    if not out_files:
        print("No Python files found in the current directory.")
        return None

    chosen_file_num=pick_from_list_static(out_files,f'{name}')
    chosen_file=out_files[chosen_file_num-1]
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
    while True:
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
        grouped_dataframe=filtered_grouped_dataframe
        keys=['Yes','No']
        chosen_key_num=pick_from_list_static(keys,f'Would like to filter data further?')
        if chosen_key_num==2:
           return filtered_dataframe

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
    
    safe_filename = f'{dependent} vs {independent}.svg'.replace("<", "").replace(">", "").replace(":", "").replace("*", "").replace("?", "").replace("\"", "").replace("\\", "").replace("/", "").replace("|", "")

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

    safe_filename = f'{dependent} vs {independent}.svg'.replace("<", "").replace(">", "").replace(":", "").replace("*", "").replace("?", "").replace("\"", "").replace("\\", "").replace("/", "").replace("|", "")

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

    safe_filename = f'{dependent} vs {independent}.svg'.replace("<", "").replace(">", "").replace(":", "").replace("*", "").replace("?", "").replace("\"", "").replace("\\", "").replace("/", "").replace("|", "")

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
    
    safe_filename = f'{dependent} vs {independent}_trend.svg'.replace("<", "").replace(">", "").replace(":", "").replace("*", "").replace("?", "").replace("\"", "").replace("\\", "").replace("/", "").replace("|", "")

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

def find_common_csv_name_merge_all():
    """
    Find all CSV files with a specific name and merge them into a single DataFrame.
    The merged data is then exported to a CSV file named 'merged_data.csv'.

    Returns:
        pandas.DataFrame: The merged data.
    """

    csv_filename = 'merged_config_stat_results.csv'
    merged_data = pd.DataFrame()

    for root, dirs, files in os.walk(os.getcwd()):
        for file in files:
            # Check if the file is a CSV file and has the desired name
            if file.endswith('.csv') and file == csv_filename:
                # Get the full path of the CSV file
                file_path = os.path.join(root, file)

                # Read the CSV file into a DataFrame
                data = pd.read_csv(file_path)

                # Merge the data with the existing merged data
                merged_data = pd.concat([merged_data, data], ignore_index=True)

    # Export the merged data to a CSV file
    merged_data.to_csv(f'{os.getcwd()}/merged_data.csv', index=False)

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


def sort_csv_filter():
    
    csv_name=list_csv_files_in_directory_choose('Location of csv file')
    column_name=show_dictionary_keys_from_csv_choose(csv_name,'Varaibles to choose from')
    # Read the CSV file into a pandas DataFrame
    df = pd.read_csv(csv_name)
    
    label = show_dictionary_keys_from_csv_choose(column_name, 'Label')
    
    # Group the data by item name
    grouped_data = df.groupby(label)
    filtered_grouped_data=filter_dataframe_by_group_selection(grouped_data,label)


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

def group_csv_by_range():
    csv_name=list_csv_files_in_directory_choose('Location of csv file')
    while True:
        user_input = input("Please specify the range you would like to group by\n (if you would like to group by 5 then an example would be 0-5,6-10,.. and so on)").strip()
        if user_input == '0':
            print("Exiting...")
            exit()
        elif not user_input.isdigit() or int(user_input) < 0:
            print("Invalid choice. Please try again.")
        else:
            user_input = int(user_input)
            break
    column_name=show_dictionary_keys_from_csv_choose(csv_name,'Varaibles to choose from')
    # Read the CSV file into a pandas DataFrame
    df = pd.read_csv(csv_name)
    max_value = df[column_name].max()
    start_range = 0
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
    
    
    while start_range <= max_value:
        end_range = start_range + user_input
        
        # Create a subset of the DataFrame based on the current range
        subset_df = sorted_df[(sorted_df[column_name] >= start_range) & (sorted_df[column_name] < end_range)]
        
        # Write the subset to a new CSV file
        subset_csv_name = f"subset_{start_range}-{end_range}_{csv_name}"
        subset_df.to_csv(subset_csv_name, index=False)
        
        # Update the start_range and end_range for the next iteration
        start_range += user_input
        end_range += user_input
    
    print("Data has been grouped and saved to separate CSV files.")

def aggregate_and_save_csv_data():
    # 1. Ask the user for the CSV file location
    csv_name = list_csv_files_in_directory_choose('Location of csv file')
    
    # 2. Ask the user for the column to be used as the unique identifier for aggregation
    unique_identifier = show_dictionary_keys_from_csv_choose(csv_name, 'Select the unique identifier for aggregation')
    
    # 3. Read the CSV file into a pandas DataFrame
    df = pd.read_csv(csv_name)

    grouped = df.groupby(unique_identifier).apply(aggregate_group)
    grouped = grouped.reset_index(drop=True)
    num_transactions = df.groupby(unique_identifier).size()
    grouped["number_of_transactions"] = grouped[unique_identifier].map(num_transactions)
    
    # 5. Save the aggregated DataFrame to a new CSV file
    aggregated_csv_name = f"aggregated_{csv_name}"
    grouped.to_csv(aggregated_csv_name, index=False)
    print(f"Aggregated data saved to {aggregated_csv_name}.")


def aggregate_group(group):
    agg_dict = {col: 'mean' for col in group if group[col].dtype in ['int64', 'float64']}
    for col in group:
        if group[col].dtype == 'object' and col != "EOM_TRANS_DATE":
            agg_dict[col] = lambda x: x.iloc[0]
    agg_dict["EOM_TRANS_DATE"] = lambda x: f"{x.min()}-{x.max()}" if x.nunique() > 1 else x.iloc[0]
    return group.agg(agg_dict)


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

def convert_out2csv():
    """
    Convert a .out file into a .csv file using pandas.
    Returns:
    None
    """
    # Reading the .out file
    filestuff = list_out_files_in_directory_choose("Out file you want to convert to csv")
    
    with open(filestuff, "r") as file:
        lines = file.readlines()

    # Parsing the headers
    headers = lines[0].strip().replace('"', '').split('  ')
    
    # Parsing the data rows
    data = []
    for line in lines[1:]:
        row = line.strip().replace('"', '').split('  ')
        data.append(row)

    # Convert data into a pandas DataFrame
    df = pd.DataFrame(data, columns=headers)

    # Letting the user specify the name of the output .csv file
    output_name = input("Please enter a name for the output CSV file (without the .csv extension): ")
    output_path = f"{output_name}.csv"
    
    # Export the DataFrame as a .csv
    df.to_csv(output_path, index=False)
    print(f"File has been converted and saved as: {output_path}")

def merge_out_with_csv():
    # Reading the .out file
    filestuff = list_out_files_in_directory_choose("Out file you want to convert to csv")
    
    with open(filestuff, "r") as file:
        lines = file.readlines()

    # Parsing the headers
    headers = lines[0].strip().replace('"', '').split('  ')
    
    # Parsing the data rows
    data = []
    for line in lines[1:]:
        row = line.strip().replace('"', '').split('  ')
        data.append(row)

    # Convert data into a pandas DataFrame
    out_df = pd.DataFrame(data, columns=headers)

    # Prompt user to choose a CSV file
    csv_file = list_csv_files_in_directory_choose("Choose a CSV file to merge with the .out DataFrame")
    csv_data = pd.read_csv(csv_file)

    # Perform the merge based on common keys
    common_keys = set(out_df.columns).intersection(set(csv_data.columns))
    merged_data = pd.merge(out_df, csv_data, on=list(common_keys))

    if merged_data.empty:
        print("No matching rows found based on the common keys. Unable to merge.")
        return

    # Letting the user specify the name of the output .csv file
    output_name = input("Please enter a name for the merged CSV file (without the .csv extension): ")
    output_path = f"{output_name}.csv"
    
    # Export the merged DataFrame as a .csv
    merged_data.to_csv(output_path, index=False)
    print(f"DataFrames merged successfully. Merged data exported to '{output_path}'.")

def merge_dat_with_out():
    # Read .out file
    out_file_path = list_out_files_in_directory_choose("Out file you want to convert to csv")

    with open(out_file_path, "r") as file:
        out_lines = file.readlines()

    # Parsing the headers from .out file
    out_headers = out_lines[0].strip().split()

    # Read .dat file
    dat_file_path = list_dat_files_in_directory_choose()
    dat_df = pd.read_csv(dat_file_path)

    # Calculate the offset for merging
    offset = len(out_lines) - len(dat_df)

    # Create a DataFrame for combined data
    combined_data = []
    for i, out_line in enumerate(out_lines[offset:], offset):
        out_values = out_line.strip().split()
        dat_values = dat_df.iloc[i - offset].values.tolist()
        combined_values = dat_values + out_values
        combined_data.append(combined_values)

    # Create merged DataFrame with appropriate columns
    merged_df = pd.DataFrame(combined_data, columns=dat_df.columns.tolist() + out_headers)

    # Letting the user specify the name of the output .csv file
    output_name = input("Please enter a name for the merged CSV file (without the .csv extension): ")
    output_path = f"{output_name}.csv"

    # Export the merged DataFrame as a .csv
    merged_df.to_csv(output_path, index=False)
    print(f"Merged data exported to '{output_path}'.")

def format_data_with_commas():
    dat_file_path = list_dat_files_in_directory_choose()
    input_file, output_file=dat_file_path,dat_file_path
    
    with open(input_file, 'r') as f:
        data = f.read()

    lines = data.strip().split('\n')
    header = lines[0].split(', ')
    values = [line.split() for line in lines[1:]]

    formatted_lines = [', '.join(header)]
    for value_set in values:
        formatted_values = [f'{float(value):,.3f}' for value in value_set]
        formatted_lines.append(', '.join(formatted_values))

    formatted_data = '\n'.join(formatted_lines)

    with open(output_file, 'w') as f:
        f.write(formatted_data)

def merge_out_dat_and_convert_to_csv():
    """
    Merge an .out file and a .dat file and save the result as a .csv file using pandas.
    Returns:
    None
    """
    # Read and convert the .out file to a DataFrame
    filestuff = list_out_files_in_directory_choose("Out file you want to merge and convert to csv")
    
    with open(filestuff, "r") as file:
        lines = file.readlines()

    # Parsing the headers
    headers = lines[0].strip().replace('"', '').split('  ')
    
    # Parsing the data rows
    data = []
    for line in lines[1:]:
        row = line.strip().replace('"', '').split('  ')
        data.append(row)

    # Convert data into a pandas DataFrame
    out_df = pd.DataFrame(data, columns=headers)

    # Read the .dat file into a DataFrame
    dat_df = pd.read_csv("your_dat_file.csv")  # Replace with the correct path

    # Calculate the starting index for merging
    start_index = len(out_df) - len(dat_df)

    # Slice the "out" dataframe and merge it with the "dat" dataframe
    merged_df = pd.concat([out_df.iloc[start_index:], dat_df], axis=1)

    # Letting the user specify the name of the output .csv file
    output_name = input("Please enter a name for the merged output CSV file (without the .csv extension): ")
    output_path = f"{output_name}.csv"

    # Export the merged DataFrame as a .csv
    merged_df.to_csv(output_path, index=False)
    print(f"Files have been merged and saved as: {output_path}")

def fix_running_out_file():
    out_file="Running_Config.out"
    with open(out_file,'r') as file:
        lines = []
        for line in file:
            lines.append(line.replace('"',''))
    file.close()
    with open(out_file,'w') as file:
        for filing in lines:
            file.write(f'{filing}')
    file.close()
    
def merge_for_histogram():
    fix_running_out_file()
    out_file="Running_Config.out"
    dat_file="running_stat.dat"
    dat_df = pd.read_csv(dat_file)
    out_df = pd.read_csv(out_file)
    # Calculate the offset for merging
    offset = len(out_df) - len(dat_df)
    interest_out_df=out_df.tail(len(dat_df)).reset_index(drop=True)
    print(interest_out_df)
    
    # Create merged DataFrame with appropriate columns
    merged_df = dat_df.join(interest_out_df)

    # Letting the user specify the name of the output .csv file
    output_name = input("Please enter a name for the merged CSV file (without the .csv extension): ")
    output_path = f"{output_name}.csv"

    # Export the merged DataFrame as a .csv
    merged_df.to_csv(output_path, index=False)
    print(f"Merged data exported to '{output_path}'.")

def plot_histogram():
    loc1 = list_csv_files_in_directory_choose('Location of csv file')
    #dependent = show_dictionary_keys_from_csv_choose(loc1, 'Dependent')
    independent = show_dictionary_keys_from_csv_choose(loc1, 'Independent')
    df = pd.read_csv(loc1)

    # Evaluate values and convert to floats
    #df[dependent] = df[dependent].apply(lambda x: eval(x) if isinstance(x, str) else x).astype(float)
    df[independent] = df[independent].apply(lambda x: eval(x) if isinstance(x, str) else x).astype(float)
    # Create a larger figure to plot all items
    fig, ax = plt.subplots(figsize=(16, 12))
    scale_axes(fig)
    bin_amount= int(input("How many bar sections do you care for?"))
    ax.hist(df[independent], bins=bin_amount, edgecolor='black')
    ax.set_xlabel(independent)  # Update the label accordingly
    ax.set_ylabel('Frequency')  # Update the y-axis label
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))

    safe_filename = f'hist_{independent}.svg'.replace("<", "").replace(">", "").replace(":", "").replace("*", "").replace("?", "").replace("\"", "").replace("\\", "").replace("/", "").replace("|", "")
    plt.savefig(safe_filename, bbox_inches='tight')
    plt.close()

def plot_single_box_and_whisker():
    # Step 1: Get the CSV file location from the user
    loc1 = list_csv_files_in_directory_choose('Location of csv file')

    # Step 2: Choose which column from the CSV file to visualize
    column_to_plot = show_dictionary_keys_from_csv_choose(loc1, 'Column to visualize')

    # Step 3: Load the data into a pandas DataFrame
    df = pd.read_csv(loc1)

    # Step 4: Clean the data
    df[column_to_plot] = df[column_to_plot].apply(lambda x: eval(x) if isinstance(x, str) else x).astype(float)

    # Step 5: Plot the box-and-whisker plot using matplotlib
    fig, ax = plt.subplots(figsize=(12, 8))
    bp = ax.boxplot(df[column_to_plot], showfliers=True, vert=True)
    ax.set_ylabel(column_to_plot)
    ax.set_title(f'Box and Whisker Plot for {column_to_plot}')

    # Extracting statistics from boxplot
    whiskers = [whisker.get_ydata()[1] for whisker in bp["whiskers"]]
    caps = [cap.get_ydata()[1] for cap in bp["caps"]]
    medians = [median.get_ydata()[1] for median in bp["medians"]]
    quartiles = [box.get_ydata()[i] for box in bp["boxes"] for i in [0, 2]]

    # Preparing statistics for the legend
    stats = {
        "Minimum": caps[0],
        "Lower Quartile": quartiles[0],
        "Median": medians[0],
        "Upper Quartile": quartiles[1],
        "Maximum": caps[1]
    }

    legend_labels = [f'{stat}: {value:.2f}' for stat, value in stats.items()]

    # Creating a custom legend
    custom_lines = [Line2D([0], [0], color='blue', lw=4) for _ in stats]
    ax.legend(custom_lines, legend_labels, loc='upper left')

    # Step 6: Allow user to scale the axes
    scale_axes(fig)
    
    safe_filename = f'box_and_whisker_{column_to_plot}.svg'.replace("<", "").replace(">", "").replace(":", "").replace("*", "").replace("?", "").replace("\"", "").replace("\\", "").replace("/", "").replace("|", "")
    plt.savefig(safe_filename, bbox_inches='tight')
    plt.close()

def plot_box_and_whisker_grouped():
    # Step 1: Get the CSV file location from the user
    loc1 = list_csv_files_in_directory_choose('Location of csv file')

    # Step 2: Choose which column from the CSV file to visualize
    column_to_plot = show_dictionary_keys_from_csv_choose(loc1, 'Column to visualize')

    # Step 3: Choose the label column
    label = show_dictionary_keys_from_csv_choose(loc1, 'Label for grouping')

    # Step 4: Load the data into a pandas DataFrame
    df = pd.read_csv(loc1)

    # Step 5: Clean the data
    df[column_to_plot] = df[column_to_plot].apply(lambda x: eval(x) if isinstance(x, str) else x).astype(float)

    unique_labels = df[label].unique()

    # Step 6: Group the data by the label and create separate plots for each group
    fig, axes = plt.subplots(len(unique_labels), 1, figsize=(12, 8 * len(unique_labels)))

    for ax, unique_label in zip(axes, unique_labels):
        data = df[df[label] == unique_label][column_to_plot].dropna().values
        bp = ax.boxplot(data, vert=True)
        ax.set_ylabel(column_to_plot)
        ax.set_title(f'{unique_label}: Box and Whisker Plot for {column_to_plot}')

        # Extract statistics and add to legend
        whiskers = [whisker.get_ydata()[1] for whisker in bp["whiskers"]]
        caps = [cap.get_ydata()[1] for cap in bp["caps"]]
        medians = [median.get_ydata()[1] for median in bp["medians"]]
        quartiles = [box.get_ydata()[i] for box in bp["boxes"] for i in [0, 2]]
        stats = {
            "Minimum": caps[0],
            "Lower Quartile": quartiles[0],
            "Median": medians[0],
            "Upper Quartile": quartiles[1],
            "Maximum": caps[1]
        }
        legend_labels = [f'{stat}: {value:.2f}' for stat, value in stats.items()]
        custom_lines = [Line2D([0], [0], color='blue', lw=4) for _ in stats]
        ax.legend(custom_lines, legend_labels, loc='upper right')

        # Allow user to scale the axes
        scale_axes(fig)

    safe_filename = f'box_and_whisker_{column_to_plot}_grouped_by_{label}.svg'.replace("<", "").replace(">", "").replace(":", "").replace("*", "").replace("?", "").replace("\"", "").replace("\\", "").replace("/", "").replace("|", "")
    plt.tight_layout()
    plt.savefig(safe_filename, bbox_inches='tight')
    plt.close()

def plot_whiskers_with_data():
    base_dir = list_directories_in_directory_choose()

    # Use glob to get all the Histogram.csv files in subdirectories
    histogram_files = glob.glob(f'{base_dir}/**/completed_run/**/Histogram.csv', recursive=True)

    # Ask the user for the column to plot using the first Histogram.csv file
    column_to_plot = show_dictionary_keys_from_csv_choose(histogram_files[0], 'Column to visualize')

    # Ask the user if they want to display fliers or not
    fliers_option_list = ['Display fliers', 'Do not display fliers']
    fliers_user_selection = pick_from_list_static(fliers_option_list, 'Fliers Option')
    show_fliers = True if fliers_user_selection == 1 else False
    
    dataframes = []  # to hold data for box plots
    labels = []  # to hold labels for the box plots
    all_legend_labels = []  # to hold legend labels for each box

    for histogram_file in histogram_files:
        protein_name = os.path.basename(os.path.dirname(histogram_file))
        matched_name = re.match(r"(.*?)(?:_\d+)?$", protein_name).group(1)  # Extract the base protein name

        # Fetch data from the parent directory of the current Histogram.csv (two levels up)
        parent_dir = os.path.abspath(os.path.join(os.path.dirname(histogram_file), '..', '..'))

        if os.path.exists(os.path.join(parent_dir, "data_multi.csv")):
            data_file = os.path.join(parent_dir, "data_multi.csv")
        elif os.path.exists(os.path.join(parent_dir, "data.csv")):
            data_file = os.path.join(parent_dir, "data.csv")

        df_data = pd.read_csv(data_file)
        relevant_row = df_data[df_data["Name"] == matched_name].iloc[0]

        df_histogram = pd.read_csv(histogram_file)
        dataframes.append(df_histogram[column_to_plot].dropna().values)
        labels.append(protein_name)

        # Fetch additional stats from the local data_multi.csv or data.csv
        stats = {
            protein_name: "",  # The protein name on its own line
            f"    Summary Stats": f"Min: {df_histogram[column_to_plot].min():.2f}, Q1: {df_histogram[column_to_plot].quantile(0.25):.2f}, Med: {df_histogram[column_to_plot].median():.2f}, Q3: {df_histogram[column_to_plot].quantile(0.75):.2f}, Max: {df_histogram[column_to_plot].max():.2f}",
            f"    Concentration & Temp": f"Monovalent Concentration: {relevant_row['Monovalent Concentration']}, Absolute Temperature: {relevant_row['Absolute Temperature']}",
            f"    Steps & Forfeiture": f"Number of Steps: {relevant_row['Number of Steps']}, Equilibrium Data Forfeiture: {relevant_row['Equilibrium Data Forfeiture']}%"
        }

        legend_labels = [f"{stat}: {value}" for stat, value in stats.items()]
        all_legend_labels.extend(legend_labels)

    # Plotting
    fig, ax = plt.subplots(figsize=(26, 12))
    bp = ax.boxplot(dataframes, vert=True, labels=labels, showfliers=show_fliers)

    # Create invisible rectangles as handles for the legend
    invisible_rect = Rectangle((0,0), 0, 0, visible=False)
    handles = [invisible_rect] * len(all_legend_labels)
    
    ax.legend(handles=handles, labels=all_legend_labels, loc='upper left', bbox_to_anchor=(1, 1))

    ax.set_xlabel("Protein Names")
    ax.set_ylabel(column_to_plot)

    safe_filename = f'grouped_box_and_whisker_{column_to_plot}.svg'.replace("<", "").replace(">", "").replace(":", "").replace("*", "").replace("?", "").replace("\"", "").replace("\\", "").replace("/", "").replace("|", "")
    plt.tight_layout()
    plt.savefig(safe_filename, bbox_inches='tight')
    plt.close()





def plot_3D_related2_time():
    loc1 = list_csv_files_in_directory_choose('Location of csv file')
    data_df = pd.read_csv(loc1)
    
    dependent = show_dictionary_keys_from_csv_choose(loc1, 'dependent')
    independent = show_dictionary_keys_from_csv_choose(loc1, 'independent')
    
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
    
    fig = plt.figure(figsize=(30, 30))
    ax = fig.add_subplot(111, projection='3d')
    
    y = data_df[independent_column]
    z = data_df[dependent_column]
    x = data_df.index
    
    ax.scatter(x, y, z, c=z, cmap='viridis')  # Color points based on Z values

    ax.set_ylabel(independent_column)
    ax.set_zlabel(dependent_column)
    ax.set_xlabel('Row Index')

    plt.title(f'3D Scatter Plot: {dependent_column} vs {independent_column}')

    # Save the plot as an image
    safe_filename = f'3D_scatter_plot_{dependent_column}_vs_{independent_column}.svg'.replace("<", "").replace(">", "").replace(":", "").replace("*", "").replace("?", "").replace("\"", "").replace("\\", "").replace("/", "").replace("|", "")
    plt.savefig(safe_filename, bbox_inches='tight')

    # Close the plot
    plt.close()



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