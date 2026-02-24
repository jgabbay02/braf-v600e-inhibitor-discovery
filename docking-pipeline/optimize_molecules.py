# optimize_molecules.py
# Step 2 of the docking pipeline: read raw conformer molecules from SDF,
# apply MMFF force field optimization, cluster by RMSD, and save the
# cluster-representative molecules for docking.
#
# Note: input/output SDF paths were on USC's CARC HPC cluster during execution.
# Update --sdf_in and --out_dir to match your environment before running.

import pandas as pd
import argparse
from rdkit import Chem

from smiles import write_mols_sdf_propless
from molecules import optimized_molec

parser = argparse.ArgumentParser(
    description="MMFF optimize conformers and cluster by RMSD"
)
parser.add_argument('--sdf_in', dest='sdf_in', type=str, required=True,
                    help="Path to input conformer_mols.sdf")
parser.add_argument('--conf', dest='conformer_num', type=int, required=True,
                    help="Number of conformers per molecule")
parser.add_argument('--out_dir', dest='out_dir', type=str, required=True,
                    help="Directory to save optimized_mols.sdf and energies.csv")
args = parser.parse_args()

if __name__ == '__main__':
    conformer_num = args.conformer_num

    loaded_raw_mols = Chem.SDMolSupplier(args.sdf_in, sanitize=False)
    print(f"Loaded {len(loaded_raw_mols)} molecules from {args.sdf_in}")

    unique_control_mols = []
    unique_control_mols, energies_df = optimized_molec(
        loaded_raw_mols, conformer_num, unique_control_mols
    )

    unique_control_mols = [Chem.AddHs(mol) for mol in unique_control_mols if mol]

    energies_df.to_csv(f'{args.out_dir}/energies.csv', index=True)
    write_mols_sdf_propless(unique_control_mols, f'{args.out_dir}/optimized_mols.sdf')

    print(f"Wrote {len(unique_control_mols)} optimized molecules to "
          f"{args.out_dir}/optimized_mols.sdf")
