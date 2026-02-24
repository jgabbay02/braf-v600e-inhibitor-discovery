# validation_test.py
# Plots top 100 docking scores comparing known ChEMBL binders vs random
# ZINC molecules. Used to validate that the docking pipeline recovers
# known good binders. Output saved to zinc_6p7g.png.

import pandas as pd
import matplotlib.pyplot as plt

if __name__ == '__main__':
    chembls = pd.read_csv('6p7g/docked_molecules/revised_docked_scores.csv')
    zincs = pd.read_csv('6p7g/zinc/docked_molecules/docking_scores.csv')
    zincs.columns = ['name', 'score']
    zincs = zincs.dropna(subset=['score'])
    lowest_scores = zincs.loc[zincs.groupby('name')['score'].idxmin()]

    chembls['Source'] = 'CHEMBL'
    zincs['Source'] = 'ZINC'
    combined = pd.concat([chembls, zincs], ignore_index=True)
    combined = combined.sort_values(by='score', ascending=True)
    top_100 = combined.head(100)

    plt.figure(figsize=(12, 6))
    colors = top_100['Source'].map({'CHEMBL': 'blue', 'ZINC': 'orange'})
    plt.scatter(range(len(top_100)), top_100['score'], c=colors, alpha=0.7)
    plt.xticks(
        ticks=range(len(top_100)),
        labels=top_100['name'],
        rotation=90,
        fontsize=8
    )
    plt.ylabel('Score (Most Negative)')
    plt.title('Top 100 Scores: CHEMBL vs ZINC Molecules')
    plt.legend(handles=[
        plt.Line2D([0], [0], color='blue', marker='o', linestyle='', label='CHEMBL'),
        plt.Line2D([0], [0], color='orange', marker='o', linestyle='', label='ZINC')
    ])
    plt.tight_layout()
    plt.savefig('zinc_6p7g.png')
    print("Saved validation plot to zinc_6p7g.png")
