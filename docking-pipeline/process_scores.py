# process_scores.py
# Step 4 of the docking pipeline: take raw docking scores, keep only the
# best score per molecule, and plot binding affinity vs pChEMBL value.

import pandas as pd
from smiles import plotting_scores_comp

if __name__ == '__main__':
    scores = pd.read_csv('docked_molecules/docking_scores.csv')
    scores.columns = ['name', 'score']
    scores = scores.dropna(subset=['score'])

    lowest_scores = scores.loc[scores.groupby('name')['score'].idxmin()]
    lowest_scores.to_csv('docked_molecules/revised_docked_scores.csv', index=False)

    chembl_info = pd.read_csv('chembl_info.csv')
    docking_scores_path = 'docked_molecules/revised_docked_scores.csv'
    plotting_scores_comp(docking_scores_path, chembl_info)
