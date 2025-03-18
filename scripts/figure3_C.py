import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from sklearn.metrics import classification_report
import os

##########################################################################################################
# load file
scims_file = "./04_chicken/results/chicken_metadata_scims_updated.txt"
rxry_file = "./04_chicken/chicken_rxry_output.txt"
bexy_file = "./04_chicken/bexy_output/chicken_bexy_output.txt"

scims = pd.read_csv(scims_file, sep="\t")
scims.rename(columns={'Run': 'Sample'}, inplace=True)

sub_cols = scims[['Sample', 'SCiMS_reads_mapped', 'host_sex']]

rxry = pd.read_csv(rxry_file, sep="\t")
rxry.rename(columns={'Run': 'Sample'}, inplace=True)
rxry_merged = pd.merge(sub_cols, rxry, on="Sample")

# Modify the 'Sample' column: remove '.sorted' and replace '_' with ''
bexy = pd.read_csv(bexy_file, sep="\t")
bexy.rename(columns={'sample': 'Sample'}, inplace=True)
bexy['Sample'] = bexy['Sample'].str.replace('.sorted', '')
bexy_merged = pd.merge(sub_cols, bexy, on="Sample")


# Add read bins
# Define bins for reads mapped
bins = [1, 500, 1000, 10000, float('inf')]
bin_labels = ['1-500','500-1K','1K-10K', '>10K']

# Classify reads mapped to X chromosome into bins
scims['read_bins'] = pd.cut(scims['SCiMS_reads_mapped'], bins=bins, labels=bin_labels, include_lowest=True)
rxry_merged['read_bins'] = pd.cut(rxry_merged['SCiMS_reads_mapped'], bins=bins, labels=bin_labels, include_lowest=True)
bexy_merged['read_bins'] = pd.cut(bexy_merged['SCiMS_reads_mapped'], bins=bins, labels=bin_labels, include_lowest=True)

##########################################################################################################

# Infer sex based on Rx
inferred_sex_Rx = []
for index, row in rxry_merged.iterrows():
    try:
        # Check if the value is a string before applying strip
        if isinstance(row['Rx 95% CI'], str):
            ci_low_rx, ci_high_rx = map(float, row['Rx 95% CI'].strip('()').split(', '))
        # If it's already a float (or any numeric type), assume it's a single value
        elif isinstance(row['Rx 95% CI'], (int, float)):
            ci_low_rx = ci_high_rx = row['Rx 95% CI']  # Assuming single value CI
        else:
            raise ValueError("Unexpected data type for 'Rx 95% CI'")

        if ci_low_rx > 0.8:
            inferred_sex_Rx.append('female')
        elif ci_high_rx < 0.6:
            inferred_sex_Rx.append('male')
        else:
            inferred_sex_Rx.append('uncertain')
    except ValueError:
        inferred_sex_Rx.append('uncertain')  # Handle missing or malformed CI values
rxry_merged['inferred_sex_Rx'] = inferred_sex_Rx

# Infer sex based on Ry
inferred_sex_Ry = []
for index, row in rxry_merged.iterrows():
    try:
        # Check if the value is a string before applying strip
        if isinstance(row['Ry 95% CI'], str):
            ci_low_ry, ci_high_ry = map(float, row['Ry 95% CI'].strip('()').split(', '))
        # If it's already a float (or any numeric type), assume it's a single value
        elif isinstance(row['Ry 95% CI'], (int, float)):
            ci_low_ry = ci_high_ry = row['Ry 95% CI']  # Assuming single value CI
        else:
            raise ValueError("Unexpected data type for 'Ry 95% CI'")

        if ci_low_ry > 0.077:
            inferred_sex_Ry.append('male')
        elif ci_high_ry < 0.016:
            inferred_sex_Ry.append('female')
        else:
            inferred_sex_Ry.append('uncertain')
    except ValueError:
        inferred_sex_Ry.append('uncertain')  # Handle missing or malformed CI values
rxry_merged['inferred_sex_Ry'] = inferred_sex_Ry

# Infer sex based on bexy
inferred_sex_bexy = []
for index, row in bexy_merged.iterrows():
    if row['sex_karyotype'] == 'XX':
        inferred_sex_bexy.append('female')
    elif row['sex_karyotype'] == 'XY':
        inferred_sex_bexy.append('male')
    else:
        inferred_sex_bexy.append('uncertain')
bexy_merged['inferred_sex_bexy'] = inferred_sex_bexy

##########################################################################################################
# Figure 3C: Bar chart of precision, recall, and F1 score for each method

# Extract the actual and inferred sexes for each method
actual_sex = scims['host_sex']
methods = {
    'SCiMS': scims['SCiMS_sex'],
    'BeXY': bexy_merged['inferred_sex_bexy'],
    'Rx': rxry_merged['inferred_sex_Rx'],
    'Ry': rxry_merged['inferred_sex_Ry']
}

# Initialize an empty DataFrame to hold all metrics
all_reports = pd.DataFrame()

# Extract the actual and inferred sexes for each method
host_sex = scims['host_sex']
methods = {
    'SCiMS': scims['SCiMS_sex'],
    'BeXY': bexy_merged['inferred_sex_bexy'],
    'Rx': rxry_merged['inferred_sex_Rx'],
    'Ry': rxry_merged['inferred_sex_Ry']
}

# Calculate and store the classification report for each method
for method, inferred_sex in methods.items():
    report = classification_report(actual_sex, inferred_sex, labels=['male', 'female'], output_dict=True, zero_division=0)
    report_df = pd.DataFrame(report).transpose()
    report_df['method'] = method  # Add a column for the method name
    all_reports = pd.concat([all_reports, report_df], axis=0)

# Reset the index to properly format the DataFrame
all_reports.reset_index(inplace=True)
all_reports.rename(columns={'index': 'class'}, inplace=True)

# Filter out support and accuracy metrics for plotting
metrics = all_reports[all_reports['class'].isin(['male', 'female'])]

# Plotting
fig, axes = plt.subplots(3, 1, figsize=(5, 5), sharex=True)

# Set colors for methods
colors = {
    'SCiMS': '#1B9E77',  # Green
    'BeXY': '#B31529',   # Red
    'Rx': '#7570B3',     # Purple
    'Ry': '#BD9E39'      # Yellow-Brown
}

color_dict = dict(zip(methods.keys(), colors.values()))

# Precision plot
sns.barplot(x='class', y='precision', hue='method', data=metrics, ax=axes[0], 
            palette=color_dict, edgecolor='black', linewidth=0.8)
axes[0].set_ylabel('Precision', fontsize=12, fontweight='bold')
axes[0].set_xlabel('')
axes[0].legend_.remove()

# Recall plot
sns.barplot(x='class', y='recall', hue='method', data=metrics, ax=axes[1], 
            palette=color_dict, edgecolor='black', linewidth=0.8)
axes[1].set_ylabel('Recall', fontsize=12, fontweight='bold')
axes[1].set_xlabel('')
axes[1].legend_.remove()

# F1 Score plot
sns.barplot(x='class', y='f1-score', hue='method', data=metrics, ax=axes[2], 
            palette=color_dict, edgecolor='black', linewidth=0.8)
axes[2].set_ylabel('F1 Score', fontsize=12, fontweight='bold')
axes[2].set_xlabel('Class', fontsize=12, fontweight='bold')
axes[2].legend_.remove()

# Improve figure readability
for ax in axes:
    ax.grid(axis='y', linestyle='--', linewidth=0.5, alpha=0.5)  # Improve grid visibility
    ax.tick_params(axis='x', labelsize=12, rotation=0)  # Ensure x-axis labels are clear
    ax.tick_params(axis='y', labelsize=12)  # Ensure y-axis labels are readable
    for spine in ax.spines.values():  # Thicken axis borders
        spine.set_linewidth(1)

# Move the legend to a single location
handles, labels = axes[2].get_legend_handles_labels()
fig.legend(handles, labels, loc='upper center', ncol=len(methods), bbox_to_anchor=(0.5, 1.05))

plt.tight_layout(rect=[0, 0, 1, 0.95])
plt.savefig('./figures/fig3C_performance_matrices.eps', dpi=300, bbox_inches='tight', transparent=True, format='eps')
plt.show()

##########################################################################################################
# Figure 3C: Pie chart of accuracy for SCiMS - pie chart is scaled by sample size

# Ensure the save directory exists
save_path = './figures/'
os.makedirs(save_path, exist_ok=True)

# Define the function to categorize predictions
def categorize_predictions(row, method_column):
    if row[method_column] == 'uncertain':
        return 'uncertain'
    elif row[method_column] == row['host_sex']:
        return 'correct'
    else:
        return 'incorrect'

# Group and calculate accuracy metrics
def calculate_accuracy_by_method(df, methods, host_sex_col):
    results = []
    for site, group in df.groupby('Host_site'):
        site_data = {'Host_site': site}
        for method in methods:
            categorized = group.apply(categorize_predictions, axis=1, method_column=method)
            
            # Count the occurrences of each category
            correct = (categorized == 'correct').sum()
            incorrect = (categorized == 'incorrect').sum()
            uncertain = (categorized == 'uncertain').sum()
            total = len(categorized)
            
            # Calculate metrics
            accuracy = correct / (correct + incorrect) if (correct + incorrect) > 0 else 0
            uncertain_rate = uncertain / total if total > 0 else 0
            
            # Store metrics in site_data
            site_data[f'{method}_correct_count'] = correct
            site_data[f'{method}_incorrect_count'] = incorrect
            site_data[f'{method}_uncertain_count'] = uncertain
            site_data[f'{method}_total_count'] = total
            site_data[f'{method}_accuracy'] = accuracy
            site_data[f'{method}_uncertain_rate'] = uncertain_rate
        results.append(site_data)
    return pd.DataFrame(results)


# Define methods and calculate accuracy
methods = ['SCiMS_sex']
accuracy_data = calculate_accuracy_by_method(scims, methods, 'host_sex')

# Add sample size data
sample_counts = scims.groupby('Host_site')['Sample'].count().reset_index()
sample_counts.columns = ['Host_site', 'Sample_size']
accuracy_data = accuracy_data.merge(sample_counts, on='Host_site')

# Save accuracy data to a file
accuracy_data.to_csv(os.path.join(save_path, 'figure3C_accuracy_data.txt'), sep='\t', index=False)

# Display a preview of the accuracy data
accuracy_data.head()

# Normalize sizes
max_size, min_size = 10, 4
sizes = accuracy_data['Sample_size']
normalized_sizes = (sizes - sizes.min()) / (sizes.max() - sizes.min())
scaled_sizes = normalized_sizes * (max_size - min_size) + min_size

# Plot pie charts
for _, row in accuracy_data.iterrows():
    Host_site = row['Host_site']
    for method in methods:
        correct = row[f'{method}_correct_count']
        incorrect = row[f'{method}_incorrect_count']
        uncertain = row[f'{method}_uncertain_count']
        total_definite = correct + incorrect
        total = total_definite + uncertain

        # Calculate rates
        accuracy_rate = correct / total_definite if total_definite > 0 else 0
        inaccuracy_rate = incorrect / total_definite if total_definite > 0 else 0
        uncertainty_rate = uncertain / total if total > 0 else 0
        chart_sizes = [correct, incorrect]
        colors = ['#0099CC', '#FF6633']
        
        # Define labels with percentages for outside annotation
        labels = [f'{accuracy_rate:.1%}', f'{inaccuracy_rate:.1%}']

        # Create the pie chart
        plt.figure(figsize=(6, 6))  # Adjust figure size as needed
        plt.pie(
            chart_sizes, colors=colors, startangle=140,
            labels=labels,  # Add the labels outside the pie chart
            autopct=None,  # Disable the default percentage display on wedges
            wedgeprops={'edgecolor': 'black', 'linewidth': 10},
            textprops={'fontsize': 50}  # Adjust font size for outside labels
        )
        plt.axis('equal')

        # Add the sample size at the top center of the pie chart (specific to SCiMS method)
        if method == 'SCiMS_sex':  # Check if the method is SCiMS
            sample_size = row[f'{method}_total_count']  # Use the SCiMS-specific total count
            plt.text(
                0, 1.3,  # Position the text above the chart
                f'n = {sample_size}', 
                ha='center', fontsize=50
            )

        # Save each chart
        filename = f"{save_path}fig3C_{method}_accuracy_pie_chart.eps"
        plt.savefig(filename, format='eps', dpi=300, bbox_inches='tight', transparent=True)
        plt.close()







