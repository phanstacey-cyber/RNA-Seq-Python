Python 3.13.0 (tags/v3.13.0:60403a5, Oct  7 2024, 09:38:07) [MSC v.1941 64 bit (AMD64)] on win32
Type "help", "copyright", "credits" or "license()" for more information.
>>> # Import required libraries
... import pandas as pd
... import matplotlib.pyplot as plt
... import numpy as np
... 
... # Read the Excel file
... df = pd.read_excel('102212 Col vs rus+gene descriptions.xlsx')
... 
... # Create volcano plot
... plt.figure(figsize=(12, 8))
... 
... # -log10 transform the p-values
... df['log10_padj'] = -np.log10(df['padj'])
... 
... # Define significance thresholds
... padj_threshold = 0.05
... fc_threshold = 2  # log2 fold change threshold of 2 (4-fold change)
... 
... # Create boolean masks for different categories
... significant = (df['padj'] < padj_threshold) & (abs(df['log2 col vs rus1-2 FoldChange']) > fc_threshold)
... 
... # Plot all points first
... plt.scatter(df.loc[~significant, 'log2 col vs rus1-2 FoldChange'], 
...            df.loc[~significant, 'log10_padj'],
...            color='grey',
...            alpha=0.5,
...            label='Not Significant')
... 
... # Plot significant points
... plt.scatter(df.loc[significant, 'log2 col vs rus1-2 FoldChange'],
...            df.loc[significant, 'log10_padj'],
...            color='red',
...            alpha=0.7,
...            label='Significant')
... 
... # Annotate significant genes
for idx in df[significant].index:
    gene_id = df.loc[idx, 'id']
    x = df.loc[idx, 'log2 col vs rus1-2 FoldChange']
    y = df.loc[idx, 'log10_padj']
    plt.annotate(df.loc[idx, 'Unnamed: 0'].split('.')[0],  # Use gene ID without the .1 suffix
                xy=(x, y),
                xytext=(5, 5),
                textcoords='offset points',
                fontsize=8,
                alpha=0.7)

# Add labels and title
plt.xlabel('log2 Fold Change')
plt.ylabel('-log10(Adjusted p-value)')
plt.title('Volcano Plot: Col vs rus1-2 with Differentially Expressed Genes')

# Add threshold lines
plt.axhline(y=-np.log10(padj_threshold), color='r', linestyle='--', alpha=0.3)
plt.axvline(x=-fc_threshold, color='r', linestyle='--', alpha=0.3)
plt.axvline(x=fc_threshold, color='r', linestyle='--', alpha=0.3)

plt.legend()
plt.tight_layout()
plt.show()

# Print summary statistics
sig_count = significant.sum()
print(f"Number of differentially expressed genes: {sig_count}")
