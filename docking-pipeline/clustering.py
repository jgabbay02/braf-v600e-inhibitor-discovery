#Takes in distance matrix and conformer numbers and returns the clsuters (by index)
def butina_clustering(dists, group_size):
    from rdkit.ML.Cluster import Butina
    clusts = Butina.ClusterData(dists, group_size, 1.25, isDistData=True, reordering=True)
    return(clusts)

#MMFF Optimization (takes in centroids indices)
def mmff_opt(mol, num_conf):
    from rdkit import Chem
    from rdkit.Chem import AllChem     
    #if index < num_conf:
    if mol == None:
        print(f'No Mol')
        return None
    if mol.GetNumConformers() == 0:
        print("not 3d")
    molh = mol
         
    Chem.RemoveHs(molh, sanitize=False)
    Chem.AddHs(molh, addCoords=True)
        #removing hydrogens, sanitize, readd
        
        
    try:
        Chem.SanitizeMol(molh)
    except ValueError as e:
        print(f"Sanitization failed: {e}")
        return None

        #Checking MMFF optimization worked
    try:
        n_attempts = 0
        max_attempts = 5
        success = False
        while not success and (n_attempts < max_attempts):
            success = AllChem.MMFFOptimizeMolecule(molh, maxIters = 200) == 0
            n_attempts+=1
    except Exception as e:
        print("MMFF Optimization failed")
        return None
    return molh
    
#Creates distance matrix (group is a group of conformers of the same molecule)
def build_rmsd_matrix(conformer_num, group):
    from rdkit import Chem
    from rdkit.Chem import rdMolAlign
    
    #Create RMSD matrix of aligned conformers
    dists = []
    for i in range(len(group)):
        for j in range(i):
            rmsd = rdMolAlign.GetBestRMS(group[i],group[j])
            dists.append(rmsd)
    return(dists)

