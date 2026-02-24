
#Takes in crystal ligand path and returns the center (for moving molecules to correct place)
def get_ref_centroid(crystal_ligand_path):
    from rdkit import Chem
    from rdkit.Chem import AllChem
    ref_ligand = Chem.MolFromMolFile(crystal_ligand_path)
    ref_conf = ref_ligand.GetConformer() #access 3d points
    ref_centroid = AllChem.ComputeCentroid(ref_conf) #get centroid
    return(ref_centroid)

#Takes in LISTS OF output_path, the reference spot, and the pocket and docks the molecule into a pdbqt file
def molecule_docking(mols, output_paths, ref_centroids, pockets, chembl_ids):
    from rdkit import Chem
    import tempfile
    import os
    from rdkit.Chem import SDMolSupplier
    from subprocess_calls import obabel_call
    from subprocess_calls import qvina_call
    import multiprocessing
    from multiprocessing import Pool
    
#for every molecule in list, create pdbqt file to get a list of temp pdbqt files. then can call pool.starmap(qvina_call, zip(ALL LISSTS IN ORDER OF INPUTS)
    #Create a Named Temporary SDF File for each optimized molecule
    pdbqts = []
    sdfs = []
    for index, (mol, ref_centroid, output_path, pocket, ids) in enumerate(zip(mols, ref_centroids, output_paths, pockets, chembl_ids)):
        with tempfile.NamedTemporaryFile(mode='w+', delete=False, suffix='.sdf') as temp_ligand_file:
            writer = Chem.SDWriter(temp_ligand_file.name)
        
            #Adjusting location of molecule
            control_mol = mol.GetConformer()
            for i in range(mol.GetNumAtoms()):
                pos_1 = control_mol.GetAtomPosition(i)
                new_pos = pos_1 + ref_centroid
                control_mol.SetAtomPosition(i, new_pos)
            writer.write(mol)
            writer.close()
            sdfs.append(temp_ligand_file.name)
            print(f'temp file for {mol.GetProp("CHEMBL_ID")} written') 
        #Convert sdf file to a pdbqt file and dock
        with tempfile.NamedTemporaryFile(mode = 'w+', delete=False, suffix = '.pdbqt') as temp_pdbqt_file:
            obabel_call(temp_ligand_file, temp_pdbqt_file)
            pdbqts.append(temp_pdbqt_file.name)
    #qvina call needs to return two column df)
    num_procs = multiprocessing.cpu_count()
    num_procs = num_procs-2
    with Pool(num_procs) as pool: 
        affinity_score = pool.starmap(qvina_call, zip(mols, pdbqts, ref_centroids, pockets, output_paths, chembl_ids))
        pool.close()
        pool.join()
    for temp_ligand_file_name in sdfs:
        os.remove(temp_ligand_file_name) 
    for temp_pdbqt_file_name in pdbqts:
        os.remove(temp_pdbqt_file_name)
    return(affinity_score)
            
#Takes in a GROUPBY OBJECT (with chembl, stereo, taut, and conformers) and returns the MMFF optimized molecule
def optimized_molec(mols, conformer_num, unique_control_mols):
    from clustering import mmff_opt
    from clustering import full_clustering
    from rdkit import Chem
    from rdkit.Chem import AllChem
    import pandas as pd
    ids = []
    before_energy = []
    after_energy = []

    #By group of twenty molecules from mols (make group = each chunk of 20 mols): before energy with h, remove h, MMFF opt, then save after energy with h, then cluster  
    
    for i in range(0, len(mols), conformer_num):
        group = [mols[j] for j in range(i, min(i + conformer_num, len(mols))) if mols[j].GetNumHeavyAtoms() <= 35]
        if not group:
            continue;
        group = [Chem.AddHs(mol, addCoords = True) for mol in group if mol]
        mmff_mol = []
        for mol in group:
            if mol:
                try:
                    Chem.SanitizeMol(mol)
                except Exception as e:
                    print(f"Sanitization failed for molecule {mol.GetProp('CHEMBL_ID')}: {e}")
                try:
                    Chem.Kekulize(mol, clearAromaticFlags=True)
                except Exception as e:
                    print(f"Kekulization failed for molecule: {e}")
                    ids.append(mol.GetProp('CHEMBL_ID') if mol.HasProp('CHEMBL_ID') else 'Unknown')
                    before_energy.append(None)
                    after_energy.append(None)
                    continue
                #try:
                    #Chem.Kekulize(mol, clearAromaticFlags=True)
                ids.append(mol.GetProp("CHEMBL_ID"))
                before_energy.append(AllChem.MMFFGetMoleculeForceField(mol, AllChem.MMFFGetMoleculeProperties(mol)).CalcEnergy())
                print(AllChem.MMFFGetMoleculeForceField(mol, AllChem.MMFFGetMoleculeProperties(mol)).CalcEnergy())
                new_mol = mmff_opt(mol, conformer_num)
                mmff_mol.append(new_mol)
                if new_mol != None:
                    after_energy.append(AllChem.MMFFGetMoleculeForceField(new_mol, AllChem.MMFFGetMoleculeProperties(new_mol)).CalcEnergy())
                    print(AllChem.MMFFGetMoleculeForceField(new_mol, AllChem.MMFFGetMoleculeProperties(new_mol)).CalcEnergy())
                else:
                    after_energy.append(None)
                    print('Optimization failed and returned a blank mol')
                #except Chem.KekulizeException:
                    #print(f"Kekulization failed for mol: {mol.GetProp('CHEMBL_ID')}. Skipping energy calculation.")
                    #ids.append(mol.GetProp("CHEMBL_ID"))
                    #before_energy.append(None)
                    #print("before_energy added for kekulize")
                    #print(f'Kekulize fail')
                    #mmff_mol.append(None)
                    #after_energy.append(None)
                    #print("after_energy added for kekulize")
            else:
                ids.append(f'No molecule at {index}')
                before_energy.append(None)
                print("before_energy added no mol")
                print(f'No molecule at {index}')
                mmff_mol.append(None)
                after_energy.append(None)
                print("after_energy added for no mol")
            print(f"Current lengths - ids: {len(ids)}, before_energy: {len(before_energy)}, after_energy: {len(after_energy)}")
        mmff_mol = [Chem.RemoveHs(mol, sanitize=False) for mol in mmff_mol if mol]
        unique_mols = []
        for mol_2 in mmff_mol:
            unique_mols.append(mol_2)
            print(mol_2.GetProp('STEREO_SMILE'))
        #If you want to cluster:
        unique_control_mols.extend(full_clustering(unique_mols, conformer_num))
        
    #Turning energies into a DF
    print(len(ids), len(before_energy), len(after_energy))
    energies_df = pd.DataFrame({"MOLEC #" : ids,
                                "Before Energy" : before_energy,
                                "After Energy" : after_energy})
    return(unique_control_mols, energies_df)
