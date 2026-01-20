import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from sklearn.metrics import classification_report

scims_file = "01_simulation/results/simulation_metadata_scims_updated.txt"
rxry_file = "01_simulation/results/simulation_rxry_output.txt"
bexy_file = "01_simulation/bexy_output/simulation_bexy_output_0.95.txt"

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

# Initialize lists to store accuracies and SEMs for read bins
read_bins = []
methods = ['SCiMS', 'BeXY', 'Rx', 'Ry']
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

# Convert read_bin to string 
accuracy_df['read_bin'] = accuracy_df['read_bin'].astype(str)

classes  = ['male', 'female']
all_rows = []

# Map of method → Series of predictions 
methods = {
    'SCiMS':  merged_df['SCiMS_sex'],
    'BeXY':   merged_df['inferred_sex_bexy'],
    'Rx':     merged_df['inferred_sex_Rx'],
    'Ry':     merged_df['inferred_sex_Ry'],
}
actual = merged_df['actual_sex_x']

for method, preds in methods.items():
    for cls in classes:
        is_cls = actual == cls

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

metrics = pd.DataFrame(all_rows)

s       = 6.69 / 4          # one-quarter of a 6.69-inch column
fig_w   = s * 6             # 6 squares across
fig_h   = s * 3             # 3 squares down

fig, axes = plt.subplots(nrows=2, ncols=4,
                         figsize=(fig_w, fig_h),
                         dpi=300, constrained_layout=False, sharey=True)

# ---------- 2. top row – stacked bars (correct / incorrect / uncertain) ----
methods_bar  = ['SCiMS', 'BeXY', 'Rx', 'Ry']
stack_cols   = ['#39c0c8', '#f08c2f', '#c3c3c3']   # correct, incorrect, uncertain

for i, meth in enumerate(methods_bar):
    ax = axes[0, i]

    df_m = (
        categorization_by_depth
        .pivot(index='read_bin_category', columns='level_1', values=meth)
        .fillna(0)
        .reindex(columns=['correct', 'incorrect', 'uncertain'])
    )

    df_m.plot(kind='bar', stacked=True, ax=ax,
              color=stack_cols, edgecolor='black', linewidth=0.8, width=0.8)

    ax.set_title(meth, fontsize=10, weight='bold')
    ax.set_xlabel('Host Reads', fontsize=9, weight='bold')
    ax.set_ylabel('Fraction of Samples', fontsize=9, weight='bold')
    ax.set_ylim(0, 1.05)
    ax.tick_params(axis='x', rotation=45, labelsize=8)

    # Thicker border
    for spine in ax.spines.values():
        spine.set_linewidth(0.8)

    # Legend only on first stacked-bar panel
    if i == 0:
        ax.legend(['Correct', 'Incorrect', 'Uncertain'],
                  frameon=False, fontsize=8)
    else:
        ax.get_legend().remove()


methods_line  = ['SCiMS', 'BeXY', 'Rx', 'Ry']
markers       = ['D', 'h', 'o', 's']
line_colors   = {
    'SCiMS': '#58508f',  # blue-violet
    'BeXY':  '#bd5090',  # magenta
    'Rx':    '#ac9546',  # olive
    'Ry':    '#eddca5'   # light khaki
}

# ── Panel: Accuracy vs. read depth
ax_acc = axes[1, 0]
for method, marker in zip(methods_line, markers):
    sns.lineplot(
        data=accuracy_df, x='read_bin', y=f'accuracy_{method}',
        marker=marker, color=line_colors[method], linewidth=1.5,
        markersize=6, ax=ax_acc, label=method
    )
    ax_acc.fill_between(
        accuracy_df['read_bin'],
        accuracy_df[f'accuracy_{method}'] - accuracy_df[f'SEM_{method}'],
        accuracy_df[f'accuracy_{method}'] + accuracy_df[f'SEM_{method}'],
        alpha=0.25, color=line_colors[method]
    )

ax_acc.set_xlabel('Host Reads', weight='bold', size=10)
ax_acc.set_ylabel('Fraction of Samples',    weight='bold', size=10)
ax_acc.set_ylim(0, 1.03)
ax_acc.tick_params(axis='x', rotation=45, labelsize=8)
ax_acc.grid(True, linestyle='--', linewidth=0.4, alpha=0.6)
ax_acc.get_legend().remove()   

# ── Panels: Precision / Recall / F1-score 
long = (
    metrics
    .melt(id_vars=['method', 'class'],
          value_vars=['precision', 'recall', 'f1-score'],
          var_name='metric', value_name='score')
)

for j, met in enumerate(['precision', 'recall', 'f1-score'], start=1):
    ax = axes[1, j]
    sns.barplot(data=long[long['metric'] == met],
                x='class', y='score', hue='method',
                palette=line_colors, edgecolor='black', linewidth=0.6, ax=ax)

    ax.set_title(met.capitalize(), weight='bold', size=10)
    ax.set_xlabel('')
    ax.set_ylim(0, 1.03)
    ax.grid(axis='y', linestyle='--', linewidth=0.4, alpha=0.5)
    ax.tick_params(axis='x', labelsize=8)
    ax.tick_params(axis='y', labelsize=8)
    ax.get_legend().remove()     


# Put one legend centred above the bottom row
handles, labels = ax_acc.get_legend_handles_labels()

# prettify borders on every subplot
for ax_row in axes:
    for ax in ax_row:
        for spine in ax.spines.values():
            spine.set_linewidth(0.8)

for ax in axes.flat:
    try:                              
        ax.set_box_aspect(1)
    except AttributeError:            
        ax.set_aspect("equal", adjustable="box")

# global legend
handles, labels = axes[1,0].get_legend_handles_labels()
fig.legend(handles, labels, title='Method',
           ncol=4, loc='upper center', bbox_to_anchor=(0.5, 1.02),
           frameon=False, fontsize=8, title_fontsize=9)

fig.tight_layout(w_pad=0.5, h_pad=0.9)

# save
fig.savefig('./figures/fig2.eps',
             bbox_inches='tight', transparent=True, format='eps', dpi=300)

plt.show()