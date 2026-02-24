#Takes in csv file from CHEMBL and synthesizes important info
def get_chembl_info(smiles_path):
    import pandas as pd
    df = pd.read_csv(smiles_path, delimiter=';', quotechar = '"')
    filtered_df = df[(df.iloc[:, 33].str.strip() == 'V600E')]
    result_df = filtered_df.iloc[:, [7, 0, 12]]
    result_df.columns = ["SMILES", "IDS", "pChEMBL_VALS"]
    result_df.reset_index(drop = True, inplace = True)
    return result_df

#Takes in smiles and gives conformer smiles
def get_conformers(smiles, number_of):
    from rdkit import Chem
    from rdkit.Chem import rdDistGeom
    conf_list = []
    for i in range(number_of):
        mol = Chem.MolFromSmiles(smiles)
        mol = Chem.AddHs(mol)
        Chem.rdDistGeom.EmbedMolecule(mol)
        conf_list.append(mol)
    return conf_list

#Takes in smiles and gives stereoisomer smiles
def get_stereoisomers(smiles):
    from rdkit import Chem
    from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers, StereoEnumerationOptions
    opts = StereoEnumerationOptions(tryEmbedding=True, unique=True, onlyUnassigned=True)
    m = Chem.MolFromSmiles(smiles)
    isomers = list(EnumerateStereoisomers(m, options=opts))
    return list(sorted([Chem.MolToSmiles(x, isomericSmiles=True) for x in isomers]))

#Takes in smiles and gives tautomers smiles
def get_tautomers(smiles):
    from rdkit import Chem
    from rdkit.Chem.MolStandardize.rdMolStandardize import TautomerEnumerator
    taute = TautomerEnumerator()
    m = Chem.MolFromSmiles(smiles)
    res = taute.Enumerate(Chem.RemoveHs(m))
    return list(sorted(res.smiles))

#Puts all smiles (CHEMBL, STEREO, TAUT, and confomer mols) into a df
def get_all_smiles(chembl_info, conf_num):
    import pandas as pd
    all_smiles = []
    for chembl_idx, smiles_large in enumerate(chembl_info["SMILES"]):
        for strism_idx, smiles_medium in enumerate(get_stereoisomers(smiles_large)):
            for conf_idx, conf_mol in enumerate(get_conformers(smiles_medium, conf_num)):
                all_smiles.append({"CHEMBL_ID" : chembl_info['IDS'][chembl_idx],
                                   "CHEMBL_SMILES" : smiles_large,
                                   "STEREOISOMERS_INDEX" : strism_idx,
                                   "STEREOISOMERS_SMILES" : smiles_medium,
                                   "CONFORMER_INDEX" : conf_idx,
                                   "CONFORMER" : conf_mol})
    all_smiles_df = pd.DataFrame(all_smiles)
    return(all_smiles_df)

#Creates a groupby object of all conformers per tautomer
def groupby_taut_conformers(all_smiles_df):
    import copy
    import pandas as pd
    
    #Creating a MultiIndex df
    all_smiles_df.set_index(["CHEMBL_ID", "STEREOISOMERS_INDEX", "TAUTOMERS_INDEX", "CONFORMER_INDEX"], inplace=True)
    
    #Dataframe that keeps track of TAUT_SMILES and CONFORMERs
    all_smiles_conf = all_smiles_df[["TAUTOMERS_SMILES", "CONFORMER"]].copy()

    #GroupBy object created for all conformers (molecules) with the same TAUT_INDEX
    grouped_conf = all_smiles_conf.groupby(level=["CHEMBL_ID", "STEREOISOMERS_INDEX", "TAUTOMERS_INDEX"])

    return(grouped_conf)
    

#Creates a groupby object of all conformers per stereoisomer
def groupby_stereo_conformers(all_smiles_df):
    import copy
    import pandas as pd

    #Creating a MultiIndex df
    all_smiles_df_mi = all_smiles_df.set_index(["CHEMBL_ID", "STEREOISOMERS_INDEX", "CONFORMER_INDEX"])
    #Dataframe that keeps track of STEREO_SMILES and CONFORMERs
    all_smiles_conf = all_smiles_df_mi[["STEREOISOMERS_SMILES", "CONFORMER"]]

    #GroupBy object created for all conformers (molecules) with the same STEREOISOMERS_INDEX
    grouped_conf = all_smiles_conf.groupby(level=["CHEMBL_ID", "STEREOISOMERS_INDEX"])

    return(grouped_conf)

#Plots the socres of docking vs standard vals
def plotting_scores_comp(docking_path, chembl_info, code):
    import numpy as np
    import pandas as pd
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    chembl_ba = np.loadtxt(docking_path, delimiter=',', dtype=str)
    names = chembl_ba[:,0]
    scores = chembl_ba[:,1]
    ba_vs_sv = []
    for index, element in enumerate(names):
        for idx, id in enumerate(chembl_info["IDS"]):
            if element.strip() == id.strip():
                ba_vs_sv.append({"Binding affinity" : scores[index].strip(),
                                 "pChEMBL value" : chembl_info.loc[idx, "pChEMBL_VALS"]})
                print(scores[index].strip())
                print(chembl_info.loc[idx, "pChEMBL_VALS"])
    ba_vs_sv = pd.DataFrame(ba_vs_sv)
    ba_vs_sv['Binding affinity'] = pd.to_numeric(ba_vs_sv['Binding affinity'], errors='coerce')
    ba_vs_sv['pChEMBL value'] = pd.to_numeric(ba_vs_sv['pChEMBL value'], errors='coerce')
    ba_vs_sv.dropna(inplace=True)
    plt.scatter(ba_vs_sv['pChEMBL value'], ba_vs_sv['Binding affinity'])
    #plt.xscale('log')
    plt.xlabel('pChEMBL Values (nM)')
    plt.ylabel('Binding Affinity (kcal/mol)')
    plt.tight_layout()
    plt.savefig(f'/home1/gabbayj/braf/{code}/docked_molecules/ba_vs_sm_{code}.png')

#Writes a list of rdkit mols to a .sdf file
def write_mols_sdf(mols, id_props, stereo_props, conf_props, file_path, remove_hs=True, append=False):
    from rdkit import Chem
    import os
    from rdkit.Chem import rdmolfiles
    #isinstance checks?
    os.makedirs(os.path.dirname(file_path), exist_ok=True)
    mode = 'w+' if not append else 'a+'
    with open(file_path, mode=mode) as f, rdmolfiles.SDWriter(f) as writer:
        writer.SetKekulize(False)
        for mol, ids, smile, conf in zip(mols, id_props, stereo_props, conf_props):
            if mol is None:
                mol = Chem.MolFromSmiles('C') #dummy mol
            if remove_hs:
                mol = Chem.RemoveHs(mol, updateExplicitCount=True, sanitize=False)
            try:
                mol.SetProp('CHEMBL_ID', ids)
                mol.SetProp('STEREO_SMILE', smile)
                mol.SetProp('CONF_ID', str(conf))
                writer.write(mol)
            except Exception as e:
                if isinstance(e, RuntimeError):
                    print('Writing mol failed, dummy mol created')
                    mol = Chem.MolFromSmiles('C')
                    mol.SetProp('CHEMBL_ID', ids)
                    mol.SetProp('STEREO_SMILE', smile)
                    mol.SetProp('CONF_ID', str(conf))
                    writer.write(mol)
                else:
                    raise e

def write_mols_sdf_propless(mols, file_path, remove_hs=True, append=False):
    from rdkit import Chem
    import os
    from rdkit.Chem import rdmolfiles
    #isinstance checks?
    os.makedirs(os.path.dirname(file_path), exist_ok=True)
    mode = 'w+' if not append else 'a+'
    with open(file_path, mode=mode) as f, rdmolfiles.SDWriter(f) as writer:
        writer.SetKekulize(False)
        for mol in mols:
            if mol is None:
                mol = Chem.MolFromSmiles('C') #dummy mol
            if remove_hs:
                mol = Chem.RemoveHs(mol, updateExplicitCount=True, sanitize=False)
            try:
                writer.write(mol)
            except Exception as e:
                if isinstance(e, RuntimeError):
                    print('Writing mol failed, dummy mol created')
                    mol = Chem.MolFromSmiles('C')
                    writer.write(mol)
                else:
                    raise e

#Puts all smiles (ZINC, STEREO, TAUT, and confomer mols) into a df
def get_all_zinc_smiles(zinc_info, conf_num):
    import pandas as pd
    all_smiles = []
    for zinc_idx, smiles_large in enumerate(zinc_info["SMILES"]):
        for strism_idx, smiles_medium in enumerate(get_stereoisomers(smiles_large)):
            for conf_idx, conf_mol in enumerate(get_conformers(smiles_medium, conf_num)):
                all_smiles.append({"ZINC_ID" : zinc_info['IDS'][zinc_idx],
                                   "ZINC_SMILES" : smiles_large,
                                   "STEREOISOMERS_INDEX" : strism_idx,
                                   "STEREOISOMERS_SMILES" : smiles_medium,
                                   "CONFORMER_INDEX" : conf_idx,
                                   "CONFORMER" : conf_mol})
    all_smiles_df = pd.DataFrame(all_smiles)
    return(all_smiles_df)

#Creates a groupby object of all conformers per stereoisomer
def groupby_zinc_stereo_conformers(all_smiles_df):
    import copy
    import pandas as pd

    #Creating a MultiIndex df
    all_smiles_df_mi = all_smiles_df.set_index(["ZINC_ID", "STEREOISOMERS_INDEX", "CONFORMER_INDEX"])
    #Dataframe that keeps track of STEREO_SMILES and CONFORMERs
    all_smiles_conf = all_smiles_df_mi[["STEREOISOMERS_SMILES", "CONFORMER"]]

    #GroupBy object created for all conformers (molecules) with the same STEREOISOMERS_INDEX
    grouped_conf = all_smiles_conf.groupby(level=["ZINC_ID", "STEREOISOMERS_INDEX"])

    return(grouped_conf)

