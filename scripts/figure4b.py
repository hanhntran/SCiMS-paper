import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from sklearn.metrics import classification_report
import os

##########################################################################################################
# load file
scims_file = "./04_chicken/results/chicken_metadata_scims_updated.txt"
rxry_file = "./04_chicken/results/chicken_rxry_output.txt"
bexy_file = "./04_chicken/bexy_output/chicken_bexy_output.txt"

if not os.path.isdir('./figures'):
    os.mkdir('./figures')

scims = pd.read_csv(scims_file, sep="\t")
scims.rename(columns={'Run': 'Sample'}, inplace=True)

sub_cols = scims[['Sample', 'SCiMS_reads_mapped', 'host_sex', 'Host_site']]


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
        elif isinstance(row['Ry 95% CI'], (int, float)):
            ci_low_ry = ci_high_ry = row['Ry 95% CI']  
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

# Merge data
merged_1 = pd.merge(scims, rxry, on="Sample")
merged_df = pd.merge(merged_1, bexy, on="Sample")

# Define the categories
categories = [150, 250, 350, 450, 1000, 5000, 10000, 100000, 1000000]
# Function to categorize based on `SCiMS_reads_mapped_y`
def categorize_reads(mapped_y):
    return min(categories, key=lambda x: abs(x - mapped_y))

merged_df['read_bin_category'] = merged_df['SCiMS_reads_mapped'].apply(categorize_reads)

def categorize_predictions(row, method_column):
    if row[method_column] == 'uncertain':
        return 'uncertain'
    elif (row['host_sex_x'] == 'female' and row[method_column] == 'female') or (row['host_sex_x'] == 'male' and row[method_column] == 'male'):
        return 'correct'
    else:
        return 'incorrect'

inferred_sex_Rx = []
for index, row in merged_df.iterrows():
    try:
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

def categorize_predictions(row, method_column):
    if row[method_column] == 'uncertain':
        return 'uncertain'
    elif (row['host_sex_x'] == 'female' and row[method_column] == 'female') or (row['host_sex_x'] == 'male' and row[method_column] == 'male'):
        return 'correct'
    else:
        return 'incorrect'

# Apply this function for each method and each read depth bin
categorization_by_host_sites = merged_df.groupby('Host_site_x').apply(lambda x: pd.DataFrame({
    'SCiMS': x.apply(categorize_predictions, axis=1, method_column='SCiMS_sex').value_counts(normalize=True),
    'BeXY': x.apply(categorize_predictions, axis=1, method_column='inferred_sex_bexy').value_counts(normalize=True),
    'Rx': x.apply(categorize_predictions, axis=1, method_column='inferred_sex_Rx').value_counts(normalize=True),
    'Ry': x.apply(categorize_predictions, axis=1, method_column='inferred_sex_Ry').value_counts(normalize=True),
})).reset_index()



sample_counts = merged_df.groupby('Host_site_x')['Sample'].count()

# Configuration
methods = ['SCiMS', 'BeXY', 'Rx', 'Ry']
category_order = ['correct', 'incorrect', 'uncertain']
colors = {
    'correct': '#39c0c8',
    'incorrect': '#f08c2f',
    'uncertain': '#c3c3c3'
}

# Prepare data: ensure category ordering and fill NaNs
categorization_by_host_sites['level_1'] = pd.Categorical(
    categorization_by_host_sites['level_1'],
    categories=category_order,
    ordered=True
)

# Select the columns that contain numerical values before filling NaNs
# Apply fillna(0) only to the columns 'SCiMS', 'BeXY', 'Rx', and 'Ry'
categorization_by_host_sites[methods] = categorization_by_host_sites[methods].fillna(0)

categorization_by_host_sites = (
    categorization_by_host_sites
    .sort_values(by=['Host_site_x', 'level_1'])
)


host_site_label = categorization_by_host_sites['Host_site_x'].unique()[0]

# Build the DataFrame with separate columns for accuracy and uncertainty
classes  = ['male', 'female']
all_rows = []

methods = {
    'SCiMS':  merged_df['SCiMS_sex'],
    'BeXY':   merged_df['inferred_sex_bexy'],
    'Rx':     merged_df['inferred_sex_Rx'],
    'Ry':     merged_df['inferred_sex_Ry'],
}
actual = merged_df['host_sex_x']

for method, preds in methods.items():
    for cls in classes:
        # rows that truly belong to this class
        is_cls = actual == cls

        # define masks
        correct   = (preds == actual) & is_cls & preds.isin(classes)
        incorrect = (preds != actual) & is_cls & preds.isin(classes)
        uncertain = (~preds.isin(classes)) & is_cls

        n_correct   = correct.sum()
        n_incorrect = incorrect.sum()
        n_uncertain = uncertain.sum()

        prec_den  = n_correct + n_incorrect
        rec_den   = n_correct + n_uncertain

        precision = n_correct / prec_den if prec_den > 0 else np.nan
        recall    = n_correct / rec_den  if rec_den  > 0 else np.nan

        if np.isnan(precision) or np.isnan(recall) or (precision+recall) == 0:
            f1 = np.nan
        else:
            f1 = 2 * precision * recall / (precision + recall)

        all_rows.append({
            'method':     method,
            'class':      cls,
            'precision':  precision,
            'recall':     recall,
            'f1-score':   f1,
            'n_correct':   n_correct,
            'n_incorrect': n_incorrect,
            'n_uncertain': n_uncertain,
        })

# tidy DataFrame for plotting
metrics = pd.DataFrame(all_rows)
metrics.to_csv('./figures/metrics_chicken.txt', sep='\t', index=False)

methods_line  = ['SCiMS', 'BeXY', 'Rx', 'Ry']
markers       = ['D', 'h', 'o', 's']
line_colors   = {
    'SCiMS': '#58508f',
    'BeXY':  '#bd5090', 
    'Rx':    '#ac9546', 
    'Ry':    '#eddca5'   
}

# Set up figure
values = metrics[['precision', 'recall', 'f1-score']]
print(values)
# Create a figure with 1 row and 3 columns of subplots
#fig, axes = plt.subplots(1, 3, figsize=(6.69, 3), sharey=True)

long = (
    metrics
    .melt(id_vars=['method', 'class'],
          value_vars=['precision', 'recall', 'f1-score'],
          var_name='metric', value_name='score')
)

# Create a figure with 1 row and 4 columns of subplots
fig, axes = plt.subplots(1, 4, figsize=(6.69, 3)) # Adjust figsize as needed for 4 plots

# --- Plotting the first bar plot (Classification by Library Source) ---
ax1 = axes[0] # Use the first subplot axis

# Redefine the list of methods for plotting purposes
plot_methods = ['SCiMS', 'BeXY', 'Rx', 'Ry']
x = np.arange(len(plot_methods))  # Use the list for x-axis positions

bar_width = 0.6 # Ensure bar_width is defined if not already

bottoms = np.zeros(len(plot_methods)) # bottoms should also match the number of methods

colors = {
    'correct': '#39c0c8',
    'incorrect': '#f08c2f',
    'uncertain': '#c3c3c3'
}
for category in category_order:
    heights = []
    for method in plot_methods: # Iterate through the list of methods
        # Access the data using the method string as the column name
        value = categorization_by_host_sites.loc[
            (categorization_by_host_sites['Host_site_x'] == host_site_label) &
            (categorization_by_host_sites['level_1'] == category),
            method # Use the method string to access the column
        ]
        heights.append(value.values[0])
    ax1.bar(
        x,
        heights,
        bottom=bottoms,
        color=colors[category],
        label=category,
        width=bar_width,
        edgecolor='black'
    )
    bottoms += np.array(heights)

# Formatting for the first plot
ax1.set_xticks(x)
ax1.set_xticklabels(plot_methods, fontsize=10) # Use the list for x-tick labels
ax1.set_ylabel('Fraction of Samples', fontsize=12)
ax1.set_ylim(0, 1.05)
ax1.grid(axis='y', linestyle='--', linewidth=0.5, alpha=0.6)
ax1.set_title('Accuracy', weight='bold', size=10) # Add a title

# Plotting the three metric plots (Precision, Recall, F1-score) 
# The `long` DataFrame and `line_colors` dictionary are correctly defined before this section
# Since axes is 1-dimensional (1 row), access using a single index.
for j, met in enumerate(['precision', 'recall', 'f1-score']):
    ax = axes[j+1]

    sns.barplot(data=long[long['metric'] == met],
                x='class', y='score', hue='method', # 'method' column from long DataFrame
                palette=line_colors, edgecolor='black', linewidth=0.6, ax=ax)

    ax.set_title(met.capitalize(), weight='bold', size=10)
    ax.set_xlabel('')
    ax.set_ylim(0, 1.03)
    ax.grid(axis='y', linestyle='--', linewidth=0.4, alpha=0.5)
    ax.tick_params(axis='x', labelsize=8)
    ax.tick_params(axis='y', labelsize=8)

# Add a single legend for all subplots
handles, labels = axes[1].get_legend_handles_labels()

# Prettify borders on every subplot
for ax in axes:
    for spine in ax.spines.values():
        spine.set_linewidth(0.8)

# Make every axis square
for ax in axes:
    try:  
        ax.set_box_aspect(1)
    except AttributeError:
        ax.set_aspect("equal", adjustable="box")

# Global legend
fig.legend(handles, labels, title='Method',
           ncol=4, loc='upper center', bbox_to_anchor=(0.5, 1.02),
           frameon=False, fontsize=8, title_fontsize=9)

# Adjust layout
fig.tight_layout(w_pad=0.5, h_pad=0.9)
plt.tight_layout(rect=[0, 0.03, 1, 0.95])
plt.savefig('./figures/fig4b.eps', bbox_inches='tight', transparent=True, format='eps', dpi=300)
plt.show()



