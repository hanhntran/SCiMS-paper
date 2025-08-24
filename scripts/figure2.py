import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from sklearn.metrics import classification_report
import os
from numpy import nan

scims_file = "01_simulation/results/simulation_metadata_scims_updated.txt"
rxry_file = "01_simulation/results/simulation_rxry_output.txt"
bexy_file = "01_simulation/bexy_output/simulation_bexy_output_0.95.txt"

if not os.path.isdir('./figures'):
    os.mkdir('./figures')

scims = pd.read_csv(scims_file, sep="\t")
scims.rename(columns={'Run': 'Sample'}, inplace=True)
scims['SCiMS_sex'] = scims['SCiMS_sex'].str.strip().str.capitalize()
#scims = scims[scims['Status'] == 'Success']

sub_cols = scims[['Sample']]

rxry = pd.read_csv(rxry_file, sep="\t")
rxry.rename(columns={'Run': 'Sample'}, inplace=True)
rxry_merged = pd.merge(sub_cols, rxry, on="Sample")

# Modify the 'Sample' column: remove '.sorted' and replace '_' with ''
bexy = pd.read_csv(bexy_file, sep="\t")
bexy.rename(columns={'sample': 'Sample'}, inplace=True)
bexy['Sample'] = bexy['Sample'].str.replace('.sorted', '')
bexy_merged = pd.merge(sub_cols, bexy, on="Sample")

# Infer sex based on Rx
inferred_sex_Rx = []
for index, row in rxry_merged.iterrows():
    try:
        # Check if the value is a string before applying strip
        if isinstance(row['Rx 95% CI'], str):
            ci_low_rx, ci_high_rx = map(float, row['Rx 95% CI'].strip('()').split(', '))
            if ci_low_rx > 0.8:
                inferred_sex_Rx.append('Female')
            elif ci_high_rx < 0.6:
                inferred_sex_Rx.append('Male')
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
                inferred_sex_Ry.append('Male')
            elif ci_high_ry < 0.016:
                inferred_sex_Ry.append('Female')
            else:
                inferred_sex_Ry.append('uncertain')
        else:
            inferred_sex_Ry.append('uncertain')  # Handle non-string values
    except ValueError:
        inferred_sex_Ry.append('uncertain')  # Handle missing or malformed CI values
rxry_merged['inferred_sex_Ry'] = inferred_sex_Ry
rxry_merged_filtered = rxry_merged.filter(['Sample','Rx', 'Rx 95% CI', 'Ry','Ry 95% CI','inferred_sex_Rx', 'inferred_sex_Ry'])

# Infer sex based on bexy
inferred_sex_bexy = []
for index, row in bexy_merged.iterrows():
    if row['sex_karyotype'] == 'XX':
        inferred_sex_bexy.append('Female')
    elif row['sex_karyotype'] == 'XY':
        inferred_sex_bexy.append('Male')
    else:
        inferred_sex_bexy.append('uncertain')
bexy_merged['inferred_sex_bexy'] = inferred_sex_bexy

##########################################################################################################
# Figure 2A: Barchart
# We'll then categorize the predictions into "correct", "incorrect", and "uncertain" for each method.
merged_1 = pd.merge(scims, rxry_merged_filtered, on="Sample")
merged_df = pd.merge(merged_1, bexy_merged, on="Sample")

# Define the categories
categories = [150, 250, 350, 450, 1000, 5000, 10000, 100000, 1000000]
# Function to categorize based on `SCiMS_reads_mapped_y`
def categorize_reads(mapped_y):
    return min(categories, key=lambda x: abs(x - mapped_y))

merged_df['read_bin_category'] = merged_df['SCiMS_reads_mapped'].apply(categorize_reads)

def categorize_predictions(row, method_column):
    if row[method_column] == 'uncertain':
        return 'uncertain'
    elif (row['actual_sex'] == 'Female' and row[method_column] == 'Female') or (row['actual_sex'] == 'Male' and row[method_column] == 'Male'):
        return 'correct'
    else:
        return 'incorrect'

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
plt.savefig('./figures/fig2A_barchart.pdf', dpi=300, bbox_inches='tight', transparent=True, format='pdf')

plt.show()

##########################################################################################################
# Figure 2B: Accuracy
methods = ['SCiMS', 'BeXY', 'Rx', 'Ry']
markers = ['D', 'h', 'o', 's']

colors = {
    'SCiMS': '#58508f',  # blue
    'BeXY': '#bd5090',   # pink
    'Rx': '#ac9546',     # sage green
    'Ry': '#eddca5'      # light yellow
}

correct = categorization_by_depth[categorization_by_depth['level_1'].str.lower() == 'correct'].copy()
correct['read_bin_label'] = correct['read_bin_category'].astype(int).astype(str)

plt.figure(figsize=(5, 5))
# Loop over methods in reverse order for plotting order control
for method, marker in zip(methods, markers):
  sns.lineplot(data=correct, x='read_bin_label', y=method, marker=marker, markersize=10, label=method, color=colors[method],
               linewidth=2)
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

classes  = ['Male', 'Female']
all_rows = []

# Map of method → Series of predictions (same as you had earlier)
methods = {
    'SCiMS':  merged_df['SCiMS_sex'],
    'BeXY':   merged_df['inferred_sex_bexy'],
    'Rx':     merged_df['inferred_sex_Rx'],
    'Ry':     merged_df['inferred_sex_Ry'],
}
actual = merged_df['actual_sex']

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
        rec_den   = n_correct + n_incorrect + n_uncertain

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

# Set up the figure with 3 vertically stacked bar charts
fig, axes = plt.subplots(3, 1, figsize=(5, 5), sharex=True)

colors = {
    'SCiMS': '#58508f',  # blue
    'BeXY': '#bd5090',   # pink
    'Rx': '#ac9546',     # sage green
    'Ry': '#eddca5'      # light yellow
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

# Move the legend outside the figure
handles, labels = axes[2].get_legend_handles_labels()
fig.legend(handles, labels, loc='upper center', ncol=len(color_dict), bbox_to_anchor=(0.5, 1.1), fontsize=12, title="Method")

plt.tight_layout(rect=[0, 0, 1, 0.95])

# Save with high DPI and heavier outline
plt.savefig('./figures/fig2C_performance_matrices.eps', dpi=300, bbox_inches='tight', transparent=True, format='eps')

plt.show()
