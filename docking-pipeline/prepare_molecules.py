# prepare_molecules.py
# Step 1 of the docking pipeline: read ChEMBL SMILES, generate stereoisomers
# and conformers, and save all conformer molecules to an SDF file for
# downstream MMFF optimization and docking.
#
# Note: output paths (conformer_mols.sdf, chembl_info.csv, all_smiles_df.csv)
# were written to a working directory on USC's CARC HPC cluster during execution.
# Update these paths to match your environment before running.

import numpy as np
import pandas as pd
import argparse
from rdkit import Chem

from smiles import get_chembl_info
from smiles import get_all_smiles
from smiles import groupby_stereo_conformers
from smiles import write_mols_sdf

parser = argparse.ArgumentParser(
    description="Generate stereoisomers and conformers from ChEMBL SMILES CSV"
)
parser.add_argument('--sp', dest='smiles_path', type=str, required=True,
                    help="Path to ChEMBL SMILES CSV")
parser.add_argument('--conf', dest='conformer_num', type=int, required=True,
                    help="Number of conformers per stereoisomer")
parser.add_argument('--out_dir', dest='out_dir', type=str, required=True,
                    help="Directory to save output CSV and SDF files")
args = parser.parse_args()

if __name__ == '__main__':
    chembl_info = get_chembl_info(args.smiles_path)
    chembl_info.to_csv(f'{args.out_dir}/chembl_info.csv', index=False)

    conformer_num = args.conformer_num
    all_smiles_df = get_all_smiles(chembl_info, conformer_num)
    all_smiles_df = all_smiles_df.drop_duplicates(
        subset=['CHEMBL_ID', 'CHEMBL_SMILES', 'STEREOISOMERS_INDEX',
                'STEREOISOMERS_SMILES', 'CONFORMER_INDEX']
    )
    all_smiles_df.to_csv(f'{args.out_dir}/all_smiles_df.csv', index=False)

    grouped_conf = groupby_stereo_conformers(all_smiles_df)

    raw_mols = []
    stereo_props = []
    id_props = []
    conf_props = []

    for name, group in grouped_conf:
        chembl_id = name[0]
        stereo_smiles = group['STEREOISOMERS_SMILES'].iloc[0]
        for idx, row in group.iterrows():
            mol = Chem.RemoveHs(row['CONFORMER'], sanitize=False)
            raw_mols.append(mol)
            id_props.append(chembl_id)
            stereo_props.append(stereo_smiles)
            conf_props.append(idx[2])  # CONFORMER_INDEX

    conformer_mols_path = f'{args.out_dir}/conformer_mols.sdf'
    write_mols_sdf(raw_mols, id_props, stereo_props, conf_props, conformer_mols_path)
    print(f"Wrote {len(raw_mols)} conformer molecules to {conformer_mols_path}")
