import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from sklearn.metrics import classification_report
import matplotlib.gridspec as gridspec
import os

##########################################################################################################
# load file
scims_file = './02_hmp/results/dbGap_metadata_filtered_scims_updated.txt'
bexy_file = './02_hmp/results/bexy/hmp_bexy_output_0.95.txt'
rxry_file = './02_hmp/results/hmp_rxry_output.txt'

if not os.path.isdir('./figures'):
    os.mkdir('./figures')

scims = pd.read_csv(scims_file, sep="\t")
scims.rename(columns={'Run': 'Sample'}, inplace=True)
scims['SCiMS_sex'] = scims['SCiMS_sex'].str.strip().str.capitalize()
scims['host_sex'] = scims['host_sex'].str.strip().str.capitalize()
sub_cols = scims[['Sample']]

rxry = pd.read_csv(rxry_file, sep="\t")
rxry.rename(columns={'Run': 'Sample'}, inplace=True)
rxry_merged = pd.merge(sub_cols, rxry, on="Sample")

# Modify the 'Sample' column: remove '.sorted' and replace '_' with ''
bexy = pd.read_csv(bexy_file, sep="\t")
bexy.rename(columns={'Run': 'Sample'}, inplace=True)
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

merged_1 = pd.merge(scims, rxry_merged_filtered, on="Sample")
merged_df = pd.merge(merged_1, bexy_merged, on="Sample")

classes  = ['Male', 'Female']
all_rows = []

methods = {
    'SCiMS':  merged_df['SCiMS_sex'],
    'BeXY':   merged_df['inferred_sex_bexy'],
    'Rx':     merged_df['inferred_sex_Rx'],
    'Ry':     merged_df['inferred_sex_Ry'],
}
actual = merged_df['host_sex']

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
metrics.to_csv('./figures/fig3_hmp_metrics.txt', sep='\t', index=False)

def calculate_accuracy(
    df: pd.DataFrame,
    methods: list[str],
    group_col: str,
    sex_col:   str
) -> pd.DataFrame:
    """
    For each group in df[group_col], and for each method in methods,
    count correct/incorrect/uncertain and compute rates.
    """
    records = []
    for grp_val, grp in df.groupby(group_col):
        total = len(grp)
        row = {group_col: grp_val, 'total_count': total}
        for m in methods:
            # Count categories
            uncertain = (grp[m] == 'uncertain').sum()
            # only count as correct if it's not 'uncertain' AND matches sex_col
            correct   = ((grp[m] == grp[sex_col]) & (grp[m] != 'uncertain')).sum()
            incorrect = total - correct - uncertain

            # save
            row.update({
                f'{m}_correct_count':    correct,
                f'{m}_incorrect_count':  incorrect,
                f'{m}_uncertain_count':  uncertain,
                f'{m}_accuracy':         correct/total if total else 0,
                f'{m}_uncertain_rate':   uncertain/total if total else 0,
            })
        records.append(row)
    return pd.DataFrame(records)

# 1. SCiMS (group by Host_site; truth in column 'host_sex')
scims_methods = ['SCiMS_sex']
acc_scims = calculate_accuracy(merged_df, scims_methods, group_col='Host_site', sex_col='host_sex')

# merge sample size
sample_size = (
    scims
    .groupby('Host_site')['Sample']
    .count()
    .rename('Sample_size')
    .reset_index()
)
acc_scims = acc_scims.merge(sample_size, on='Host_site')
acc_scims.to_csv('./figures/hmp_scims_accuracy_data.txt', sep='\t', index=False)

# 2. Bexy
bexy_methods = ['inferred_sex_bexy']
acc_bexy = calculate_accuracy(
    merged_df,
    bexy_methods,
    group_col='Host_site',
    sex_col='host_sex'
)
# merge sample size
sample_size = (
    merged_df
    .groupby('Host_site')['Sample']
    .count()
    .rename('Sample_size')
    .reset_index()
)
acc_bexy = acc_bexy.merge(sample_size, on='Host_site')
acc_bexy.to_csv('./figures/hmp_bexy_accuracy_data.txt', sep='\t', index=False)

# 3. Rx
rxry_methods = ['inferred_sex_Rx']
acc_rx = calculate_accuracy(
    merged_df,
    rxry_methods,
    group_col='Host_site',
    sex_col='host_sex'
)
# merge sample size
sample_size = (
    merged_df
    .groupby('Host_site')['Sample']
    .count()
    .rename('Sample_size')
    .reset_index()
)
acc_rx = acc_rx.merge(sample_size, on='Host_site')
acc_rx.to_csv('./figures/hmp_rx_accuracy_data.txt', sep='\t', index=False)

# 4. Ry
ry_methods = ['inferred_sex_Ry']
acc_ry = calculate_accuracy(
    merged_df,
    ry_methods,
    group_col='Host_site',
    sex_col='host_sex'
)
# merge sample size
sample_size = (
    merged_df
    .groupby('Host_site')['Sample']
    .count()
    .rename('Sample_size')
    .reset_index()
)
acc_ry = acc_ry.merge(sample_size, on='Host_site')
acc_ry.to_csv('./figures/hmp_ry_accuracy_data.txt', sep='\t', index=False)

###############################################################################
# 4.  plotting figure 3
###############################################################################
# Colors
colors = {
    'Accuracy': '#39c0c8',
    'Misclassification': '#f08c2f',
    'Uncertainty': '#c3c3c3'
}

labels = acc_scims['Host_site']

y_pos = range(len(labels))

# Set up figure
fig = plt.figure(figsize=(6.69, 5.5))
gs = gridspec.GridSpec(2, 3, height_ratios=[1, 2])  # 2 rows, 3 columns

# === Panel A: Stacked Horizontal Bar Chart ===
ax_top = plt.subplot(gs[0, :])  # Top row spans all 3 columns
bar_height = 0.5

# Plot bar chart
for _, row in acc_scims.iterrows():
    host_site = row['Host_site']
    for method in scims_methods:
        correct = row[f'{method}_correct_count']
        incorrect = row[f'{method}_incorrect_count']
        uncertain = row[f'{method}_uncertain_count']
        total_definite = correct + incorrect
        total = total_definite + uncertain

        # Calculate rates
        accuracy_rate = correct / total if total > 0 else 0
        inaccuracy_rate = incorrect / total if total > 0 else 0
        uncertainty_rate = uncertain / total if total > 0 else 0
        chart_sizes = [correct, incorrect, uncertain]

        # Colors
        colors = {
            'Accuracy': '#39c0c8',
            'Misclassification': '#f08c2f',
            'Uncertainty': '#c3c3c3'
        }        
                
        # Define labels with percentages for outside annotation
        #labels = [f'{accuracy_rate:.1%}', f'{inaccuracy_rate:.1%}']

        ax_top.barh(host_site, accuracy_rate, color=colors['Accuracy'], edgecolor='black', height=bar_height)
        ax_top.barh(host_site, inaccuracy_rate, left=accuracy_rate, color=colors['Misclassification'], edgecolor='black', height=bar_height)
        ax_top.barh(host_site, uncertainty_rate, left=accuracy_rate + inaccuracy_rate, color=colors['Uncertainty'], edgecolor='black', height=bar_height)

    # Annotate percentages
    ax_top.text(accuracy_rate/2, host_site, f"{accuracy_rate*100:.1f}%", ha='center', va='center', fontsize=9, color='black')
    if inaccuracy_rate > 0:
        ax_top.text(accuracy_rate + inaccuracy_rate/2, host_site, f"{inaccuracy_rate*100:.1f}%", ha='center', va='center', fontsize=9, color='black')
    if uncertainty_rate > 0:
        ax_top.text(accuracy_rate + inaccuracy_rate + uncertainty_rate/2, host_site, f"{uncertainty_rate*100:.1f}%", ha='center', va='center', fontsize=9, color='black')

    # Add horizontal gridlines
    ax_top.grid(axis='x', linestyle='--', linewidth=0.5, alpha=0.7)

# Labels and formatting
ax_top.set_yticks(y_pos)
ax_top.set_yticklabels(labels, fontsize=10)
ax_top.set_xlim(0, 1.05)
ax_top.set_xlabel('')
ax_top.invert_yaxis()  # Match figure orientation
ax_top.spines[['top', 'right', 'left', 'bottom']].set_visible(False)
ax_top.tick_params(left=False, bottom=False)

# === Panel B: Metrics ===
long = metrics.melt(id_vars=['method', 'class'],
                    value_vars=['precision', 'recall', 'f1-score'],
                    var_name='metric', value_name='score')
line_colors = {
    'SCiMS': '#58508f',
    'BeXY': '#bd5090',
    'Rx': '#ac9546',
    'Ry': '#eddca5'
}
for i, met in enumerate(['precision', 'recall', 'f1-score']):
    ax = plt.subplot(gs[1, i])
    sns.barplot(data=long[long['metric'] == met],
                x='class', y='score', hue='method',
                palette=line_colors, edgecolor='black', linewidth=0.6, ax=ax)
    ax.set_title(met.capitalize(), fontsize=10, weight='bold')
    ax.set_ylim(0, 1.03)
    ax.grid(axis='y', linestyle='--', linewidth=0.5, alpha=0.5)
    ax.tick_params(axis='x', labelsize=8)
    ax.tick_params(axis='y', labelsize=8)
    ax.set_xlabel('')
    if i != 0:
        ax.get_legend().remove()

    # ── make every axis square ───────────────────────────────────────
    try:                              # Matplotlib ≥3.4
        ax.set_box_aspect(1)
    except AttributeError:            # fallback for older versions
            ax.set_aspect("equal", adjustable="box")

# Shared legend
handles, labels = ax.get_legend_handles_labels()
fig.legend(handles, labels, loc='upper center', title='Method', ncol=4, bbox_to_anchor=(0.5, 0.98))

fig.tight_layout(w_pad=0.5, h_pad=0.9)

plt.tight_layout(rect=[0, 0, 1, 0.95])
plt.savefig('./figures/fig3.eps', dpi=300, bbox_inches='tight', transparent=True, format='eps')
plt.show()