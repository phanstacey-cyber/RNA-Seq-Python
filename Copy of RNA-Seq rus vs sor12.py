Python 3.13.0 (tags/v3.13.0:60403a5, Oct  7 2024, 09:38:07) [MSC v.1941 64 bit (AMD64)] on win32
Type "help", "copyright", "credits" or "license()" for more information.
>>> import pandas as pd
... 
... # Load the Excel file to inspect its contents
... file_path = "102212 rus vs sor12 differential.xlsx" 
... df = pd.read_excel(file_path)
... 
... # Display the first few rows of the dataframe to understand its structure
... print(df.head())
... print(df.columns)
... 
... ########################################
... 
... import matplotlib.pyplot as plt
... import numpy as np
... 
... # Check for non-finite values in the relevant columns
... #nonfinitelogfoldchange-gene has zero counts in one condition and non-zero in another, log transformation can’t be calculated straightforwardly
... non_finite_log2fc = df['log2 rus1-2 vs sor12 FoldChange'].isna().sum()
... #non finite padj-no statistical basis to calculate an adjusted p-value for some genes.
... non_finite_padj = df['padj'].isna().sum()
... 
... print('Non-finite values in log2 Fold Change:', non_finite_log2fc)
... print('Non-finite values in padj:', non_finite_padj)
... 
... # Remove rows with non-finite values in these columns
... df_clean = df.dropna(subset=['log2 rus1-2 vs sor12 FoldChange', 'padj'])
... 
... print('Cleaned dataframe shape:', df_clean.shape)
... 
... ###############################################
... import matplotlib.pyplot as plt
... import numpy as np
... 
... # Define significance thresholds
... padj_threshold = 0.05
log2fc_threshold = 1

# Calculate -log10(padj) for the y-axis
neg_log10_padj = -np.log10(df_clean['padj'])

# Create a new column for significance
# Significantly upregulated
upregulated = (df_clean['padj'] < padj_threshold) & (df_clean['log2 rus1-2 vs sor12 FoldChange'] > log2fc_threshold)
# Significantly downregulated
downregulated = (df_clean['padj'] < padj_threshold) & (df_clean['log2 rus1-2 vs sor12 FoldChange'] < -log2fc_threshold)

# Plot
plt.figure(figsize=(12, 8))
plt.scatter(df_clean['log2 rus1-2 vs sor12 FoldChange'], neg_log10_padj, color='grey', alpha=0.5, label='Not significant')
plt.scatter(df_clean.loc[upregulated, 'log2 rus1-2 vs sor12 FoldChange'], neg_log10_padj[upregulated], color='red', label='Upregulated')
plt.scatter(df_clean.loc[downregulated, 'log2 rus1-2 vs sor12 FoldChange'], neg_log10_padj[downregulated], color='blue', label='Downregulated')

# Annotate top significant genes (to avoid overcrowding)
significant_genes = df_clean[upregulated | downregulated].sort_values('padj').head(20)
for _, gene in significant_genes.iterrows():
    plt.annotate(gene['Unnamed: 0'], 
                xy=(gene['log2 rus1-2 vs sor12 FoldChange'], -np.log10(gene['padj'])),
                xytext=(5, 5), textcoords='offset points', 
                fontsize=8, alpha=0.7)

plt.axhline(y=-np.log10(padj_threshold), color='black', linestyle='--', linewidth=0.8)
plt.axvline(x=log2fc_threshold, color='black', linestyle='--', linewidth=0.8)
plt.axvline(x=-log2fc_threshold, color='black', linestyle='--', linewidth=0.8)

plt.xlabel('Log2 Fold Change (rus1-2 vs sor12)')
plt.ylabel('-Log10 Adjusted P-value')
plt.title('Volcano Plot of Differential Gene Expression')
plt.legend()

# Print summary statistics
print('Number of significantly upregulated genes:', sum(upregulated))
print('Number of significantly downregulated genes:', sum(downregulated))

