import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from sklearn.metrics import classification_report

scims_file = "./results/scims/simulation_metadata_scims_updated.txt"
rxry_file = "./results/simulation_rxry_output.txt"
bexy_file = "./bexy_output/simulation_bexy_output_0.95.txt"

scims = pd.read_csv(scims_file, sep="\t")
scims.rename(columns={'Run': 'Sample'}, inplace=True)
#scims = scims[scims['Status'] == 'Success']

sub_cols = scims[['Sample', 'SCiMS_reads_mapped', 'actual_sex']]

rxry = pd.read_csv(rxry_file, sep="\t")
rxry.rename(columns={'Run': 'Sample'}, inplace=True)
rxry_merged = pd.merge(sub_cols, rxry, on="Sample")

# Modify the 'Sample' column: remove '.sorted' and replace '_' with ''
bexy = pd.read_csv(bexy_file, sep="\t")
bexy.rename(columns={'sample': 'Sample'}, inplace=True)
bexy['Sample'] = bexy['Sample'].str.replace('.sorted', '')
bexy_merged = pd.merge(sub_cols, bexy, on="Sample")

# Define the categories
categories = [150, 250, 350, 450, 1000, 5000, 10000, 100000, 1000000]
# Function to categorize based on `SCiMS_reads_mapped_y`
def categorize_reads(mapped_y):
    return min(categories, key=lambda x: abs(x - mapped_y))

# Apply the function to create a new `read_bin_category` column
scims['read_bin_category'] = scims['SCiMS_reads_mapped'].apply(categorize_reads)
rxry_merged['read_bin_category'] = rxry_merged['SCiMS_reads_mapped'].apply(categorize_reads)
bexy_merged['read_bin_category'] = bexy_merged['SCiMS_reads_mapped'].apply(categorize_reads)

# Infer sex based on Rx
inferred_sex_Rx = []
for index, row in rxry_merged.iterrows():
    try:
        # Check if the value is a string before applying strip
        if isinstance(row['Rx 95% CI'], str):
            ci_low_rx, ci_high_rx = map(float, row['Rx 95% CI'].strip('()').split(', '))
            if ci_low_rx > 0.8:
                inferred_sex_Rx.append('female')
            elif ci_high_rx < 0.6:
                inferred_sex_Rx.append('male')
            else:
                inferred_sex_Rx.append('uncertain')
        else:
            inferred_sex_Rx.append('uncertain')  # Handle non-string values
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
            if ci_low_ry > 0.077:
                inferred_sex_Ry.append('male')
            elif ci_high_ry < 0.016:
                inferred_sex_Ry.append('female')
            else:
                inferred_sex_Ry.append('uncertain')
        else:
            inferred_sex_Ry.append('uncertain')  # Handle non-string values
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
# Figure 2A: Barchart
# To handle uncertain cases, we'll identify cases where predictions are marked as "uncertain".
# We'll then categorize the predictions into "correct", "incorrect", and "uncertain" for each method.

merged_1 = pd.merge(scims, rxry, on="Sample")
merged_df = pd.merge(merged_1, bexy, on="Sample")

merged_df['read_bin_category'] = merged_df['SCiMS_reads_mapped'].apply(categorize_reads)

def categorize_predictions(row, method_column):
    if row[method_column] == 'uncertain':
        return 'uncertain'
    elif (row['actual_sex_x'] == 'female' and row[method_column] == 'female') or (row['actual_sex_x'] == 'male' and row[method_column] == 'male'):
        return 'correct'
    else:
        return 'incorrect'

inferred_sex_Rx = []
for index, row in merged_df.iterrows():
    try:
        # Check if the value is a string before applying strip
        if isinstance(row['Rx 95% CI'], str):
            ci_low_rx, ci_high_rx = map(float, row['Rx 95% CI'].strip('()').split(', '))
            if ci_low_rx > 0.8:
                inferred_sex_Rx.append('female')
            elif ci_high_rx < 0.6:
                inferred_sex_Rx.append('male')
            else:
                inferred_sex_Rx.append('uncertain')
        else:
            inferred_sex_Rx.append('uncertain')  # Handle non-string values
    except ValueError:
        inferred_sex_Rx.append('uncertain')  # Handle missing or malformed CI values
merged_df['inferred_sex_Rx'] = inferred_sex_Rx

inferred_sex_Ry = []
for index, row in merged_df.iterrows():
    try:
        # Check if the value is a string before applying strip
        if isinstance(row['Ry 95% CI'], str):
            ci_low_ry, ci_high_ry = map(float, row['Ry 95% CI'].strip('()').split(', '))
            if ci_low_ry > 0.077:
                inferred_sex_Ry.append('male')
            elif ci_high_ry < 0.016:
                inferred_sex_Ry.append('female')
            else:
                inferred_sex_Ry.append('uncertain')
        else:
            inferred_sex_Ry.append('uncertain')  # Handle non-string values
    except ValueError:
        inferred_sex_Ry.append('uncertain')  # Handle missing or malformed CI values
merged_df['inferred_sex_Ry'] = inferred_sex_Ry

# Infer sex based on bexy
inferred_sex_bexy = []
for index, row in merged_df.iterrows():
    if row['sex_karyotype'] == 'XX':
        inferred_sex_bexy.append('female')
    elif row['sex_karyotype'] == 'XY':
        inferred_sex_bexy.append('male')
    else:
        inferred_sex_bexy.append('uncertain')
merged_df['inferred_sex_bexy'] = inferred_sex_bexy

# Apply this function for each method and each read depth bin
categorization_by_depth = merged_df.groupby('read_bin_category').apply(lambda x: pd.DataFrame({
    'SCiMS': x.apply(categorize_predictions, axis=1, method_column='SCiMS_sex').value_counts(normalize=True),
    'Rx': x.apply(categorize_predictions, axis=1, method_column='inferred_sex_Rx').value_counts(normalize=True),
    'Ry': x.apply(categorize_predictions, axis=1, method_column='inferred_sex_Ry').value_counts(normalize=True),
    'BeXY': x.apply(categorize_predictions, axis=1, method_column='inferred_sex_bexy').value_counts(normalize=True)
})).reset_index()

sample_counts = merged_df.groupby('read_bin_category')['Sample'].count()

# Prepare the plot for each method, visualizing the fractions of correct, incorrect, and uncertain predictions.
fig, ax = plt.subplots(2, 2, figsize=(10, 10))

methods = ['SCiMS', 'BeXY', 'Rx', 'Ry']
axes = ax.flatten()

colors=['#0099CC', '#FF6633', 'gray']

for i, method in enumerate(methods):
    method_data = categorization_by_depth.pivot(index='read_bin_category', columns='level_1', values=method).fillna(0)
    
    # Plot with heavier edges
    method_data.plot(
        kind='bar',
        stacked=True,
        ax=axes[i],
        color=colors,
        width=0.8,
        edgecolor='black',  # Darker edge color
        linewidth=0.8  # Thicker bar outlines
    )
    
    # Increase spine thickness for the plot
    for spine in axes[i].spines.values():
        spine.set_linewidth(1.5)  # Make axis borders thicker
        spine.set_color("black")  # Darker outline
    
    axes[i].set_title(f'{method}', fontsize=12, fontweight='bold')
    axes[i].set_ylabel('Fraction of Samples', fontsize=12, fontweight='bold')
    axes[i].set_xlabel('Host Reads', fontsize=12, fontweight='bold')
    
    # Customize legend
    axes[i].legend(['Correct', 'Incorrect', 'Uncertain'], loc='lower right', frameon=True, edgecolor='black')

    # Make x-tick labels bold and more readable
    axes[i].set_xticklabels(axes[i].get_xticklabels(), rotation=45, ha='right', fontsize=12)

# Increase figure outline thickness
plt.gca().spines["top"].set_linewidth(1)
plt.gca().spines["bottom"].set_linewidth(1)
plt.gca().spines["left"].set_linewidth(1)
plt.gca().spines["right"].set_linewidth(1)

plt.tight_layout()

# Save with high DPI and heavier outline
plt.savefig('./figures/fig2A_barchart.eps', dpi=300, bbox_inches='tight', transparent=True, format='eps')

plt.show()

##########################################################################################################
# Figure 2B: Accuracy
methods = ['SCiMS', 'BeXY', 'Rx', 'Ry']
markers = ['D', 'h', 'o', 's']


# Initialize lists to store accuracies and SEMs for read bins
read_bins = []
accuracies = {method: [] for method in methods}
sem_accuracies = {method: [] for method in methods}

# Group the data by read bins
read_bin_all = sorted(list(set(scims['read_bin_category'].unique()) | 
                            set(rxry_merged['read_bin_category'].unique()) | 
                            set(bexy_merged['read_bin_category'].unique())))


# Calculate accuracies for each group
## SCiMS
for read_bin in read_bin_all:
    group = scims[scims['read_bin_category'] == read_bin]
    if len(group) == 0:
        accuracies['SCiMS'].append(np.nan)
        sem_accuracies['SCiMS'].append(np.nan)
    else:
        correct_predictions = group[group['SCiMS_sex'] == group['actual_sex']]
        accuracy = len(correct_predictions) / len(group) if len(group) > 0 else 0
        accuracies['SCiMS'].append(accuracy)
        sem = np.sqrt(accuracy * (1 - accuracy) / len(group)) if len(group) > 0 else 0
        sem_accuracies['SCiMS'].append(sem)


## BeXY
for read_bin in read_bin_all:
    group = bexy_merged[bexy_merged['read_bin_category'] == read_bin]
    if len(group) == 0:
        accuracies['BeXY'].append(np.nan)
        sem_accuracies['BeXY'].append(np.nan)
    else:
        correct_predictions = group[group['inferred_sex_bexy'] == group['actual_sex']]
        accuracy = len(correct_predictions) / len(group) if len(group) > 0 else 0
        accuracies['BeXY'].append(accuracy)
        sem = np.sqrt(accuracy * (1 - accuracy) / len(group)) if len(group) > 0 else 0
        sem_accuracies['BeXY'].append(sem)

## Rx
for read_bin in read_bin_all:
    group = rxry_merged[rxry_merged['read_bin_category'] == read_bin]
    if len(group) == 0:
        accuracies['Rx'].append(np.nan)
        sem_accuracies['Rx'].append(np.nan)
    else:
        correct_predictions = group[group['inferred_sex_Rx'] == group['actual_sex_x']]
        accuracy = len(correct_predictions) / len(group) if len(group) > 0 else 0
        accuracies['Rx'].append(accuracy)
        sem = np.sqrt(accuracy * (1 - accuracy) / len(group)) if len(group) > 0 else 0
        sem_accuracies['Rx'].append(sem)

## Ry
for read_bin in read_bin_all:
    group = rxry_merged[rxry_merged['read_bin_category'] == read_bin]
    if len(group) == 0:
        accuracies['Ry'].append(np.nan)
        sem_accuracies['Ry'].append(np.nan)
    else:
        correct_predictions = group[group['inferred_sex_Ry'] == group['actual_sex_x']]
        accuracy = len(correct_predictions) / len(group) if len(group) > 0 else 0
        accuracies['Ry'].append(accuracy)
        sem = np.sqrt(accuracy * (1 - accuracy) / len(group)) if len(group) > 0 else 0
        sem_accuracies['Ry'].append(sem)

# Build the DataFrame with separate columns for accuracy and uncertainty
data = {'read_bin': read_bin_all}
for method in methods:
    data[f'accuracy_{method}'] = accuracies[method]
    data[f'SEM_{method}'] = sem_accuracies[method]

accuracy_df = pd.DataFrame(data)

# Convert read_bin to string (or categorical) for consistent plotting
accuracy_df['read_bin'] = accuracy_df['read_bin'].astype(str)

# Set colors for methods
colors = {
    'SCiMS': '#1B9E77',  # Green
    'BeXY': '#B31529',   # Red
    'Rx': '#7570B3',     # Purple
    'Ry': '#BD9E39'      # Yellow-Brown
}

# Set figure size
plt.figure(figsize=(5, 5))

# Loop over methods in reverse order for plotting order control
for method, marker in zip(methods, markers):
    sns.lineplot(
        x='read_bin',
        y=f'accuracy_{method}',
        data=accuracy_df,
        marker=marker,
        linestyle='-',
        color=colors[method],  # Use pre-defined colors
        label=method,
        linewidth=2,  # Thicker line for better visibility
        markersize=10
    )
    
    # Add error band (SEM)
    plt.fill_between(
        accuracy_df['read_bin'],
        accuracy_df[f'accuracy_{method}'] - accuracy_df[f'SEM_{method}'],
        accuracy_df[f'accuracy_{method}'] + accuracy_df[f'SEM_{method}'],
        color=colors[method],
        alpha=0.2,
        edgecolor=colors[method]  # Add edge color for better visibility
    )

# Customize plot appearance
plt.ylim(0, 1.03)
plt.xlabel('Host Reads', fontsize=12, fontweight='bold')
plt.ylabel('Proportion of Correct Predictions', fontsize=12, fontweight='bold')

# Improve grid visibility
plt.grid(True, which='both', linestyle='--', linewidth=0.5, alpha=0.6)

# Customize tick labels
plt.xticks(rotation=45, fontsize=12)
plt.yticks(fontsize=12)

# Make axis borders darker and thicker
for spine in plt.gca().spines.values():
    spine.set_linewidth(1)  # Thick border
    spine.set_color("black")  # Darker color

# Improve legend appearance
plt.legend(title='Method', fontsize=12, title_fontsize=12, loc='lower right', frameon=True, edgecolor='black')

plt.tight_layout()

# Save with high resolution and heavier outline
plt.savefig('./figures/fig2B_accuracy.eps', dpi=300, bbox_inches='tight', transparent=True, format='eps')

plt.show()


##########################################################################################################
# Figure 2C: Bar chart of precision, recall, and F1 score for each method

# Initialize an empty DataFrame to hold all metrics
all_reports = pd.DataFrame()

# Extract the actual and inferred sexes for each method
actual_sex = merged_df['actual_sex_x']
methods = {
    'SCiMS': merged_df['SCiMS_sex'],
    'BeXY': merged_df['inferred_sex_bexy'],
    'Rx': merged_df['inferred_sex_Rx'],
    'Ry': merged_df['inferred_sex_Ry']
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

# Define a custom color palette
color_dict = dict(zip(methods.keys(), colors.values()))

# Set up the figure with 3 vertically stacked bar charts
fig, axes = plt.subplots(3, 1, figsize=(5, 5), sharex=True)


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

# Move the legend outside the figure
handles, labels = axes[2].get_legend_handles_labels()
fig.legend(handles, labels, loc='upper center', ncol=len(color_dict), bbox_to_anchor=(0.5, 1.1), fontsize=12, title="Method")

plt.tight_layout(rect=[0, 0, 1, 0.95])

# Save with high DPI and heavier outline
plt.savefig('./figures/fig2C_performance_matrices.eps', dpi=300, bbox_inches='tight', transparent=True, format='eps')

plt.show()