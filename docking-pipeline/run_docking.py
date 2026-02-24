# run_docking.py
# Step 3 of the docking pipeline: dock MMFF-optimized molecules against a
# receptor pocket using AutoDock QVina and save binding affinity scores.
# Runs in parallel using multiprocessing.
#
# Example usage:
#   python run_docking.py \
#     --clp crystal_ligand.sdf \
#     --p pocket.pdbqt \
#     --conf 20 \
#     --op docked_molecules \
#     --dsf docking_scores.csv

import pandas as pd
import argparse
import sys
import os
import time
from rdkit import Chem

from molecules import get_ref_centroid
from molecules import molecule_docking

if __name__ == '__main__':
    t_0 = time.time()

    parser = argparse.ArgumentParser(
        description="Dock optimized ChEMBL molecules against a BRAF V600E pocket"
    )
    parser.add_argument('--clp', dest='crystal_ligand_path', type=str, required=True,
                        help="Path to crystal ligand SDF (used to center docking box)")
    parser.add_argument('--p', dest='binding_pocket', type=str, required=True,
                        help="Path to binding pocket PDBQT")
    parser.add_argument('--conf', dest='conformer_num', type=int, required=True,
                        help="Number of conformers per molecule")
    parser.add_argument('--op', dest='output_path', type=str, default='docked_molecules',
                        help="Directory to save docked PDBQT files")
    parser.add_argument('--dsf', dest='docking_scores_file', type=str,
                        default='docking_scores.csv',
                        help="Filename for docking scores CSV")
    parser.add_argument('--debug', action='store_true',
                        help="Debug mode: process only first 5 molecules")
    args = parser.parse_args()

    unique_control_mols = Chem.SDMolSupplier('optimized_mols.sdf', sanitize=False)
    if unique_control_mols is None:
        print("Error: could not read optimized_mols.sdf")
        sys.exit(1)

    if args.debug:
        unique_control_mols = [
            unique_control_mols[i] for i in range(5)
            if unique_control_mols[i] is not None
        ]
        print("Running in debug mode (5 molecules)")
    else:
        unique_control_mols = [mol for mol in unique_control_mols if mol is not None]

    # Filter to drug-like size
    unique_control_mols = [
        mol for mol in unique_control_mols if mol.GetNumHeavyAtoms() <= 30
    ]
    print(f"Docking {len(unique_control_mols)} molecules")

    chembl_ids = [mol.GetProp('CHEMBL_ID') for mol in unique_control_mols]
    ref_centroid = get_ref_centroid(args.crystal_ligand_path)

    os.makedirs(args.output_path, exist_ok=True)

    output_path_list = [
        f'{args.output_path}/{mol.GetProp("CHEMBL_ID")}_{mol.GetProp("CONF_ID")}_out.pdbqt'
        for mol in unique_control_mols
    ]
    ref_centroid_list = [ref_centroid] * len(unique_control_mols)
    pocket_list = [args.binding_pocket] * len(unique_control_mols)

    results = molecule_docking(
        unique_control_mols, output_path_list,
        ref_centroid_list, pocket_list, chembl_ids
    )

    docking_scores_path = os.path.join(args.output_path, args.docking_scores_file)
    results_df = pd.DataFrame(results, columns=['CHEMBL_ID', 'Binding Affinity (kcal/mol)'])
    results_df.to_csv(docking_scores_path, index=False)

    print(f"Saved docking scores to {docking_scores_path}")
    print(f"Total time: {time.time() - t_0:.1f}s")

    if args.debug:
        print(results_df)
