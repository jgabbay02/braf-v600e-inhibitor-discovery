#Obabel call to convert ligand sdf file to a pdbqt file
def obabel_call(temp_ligand_file, temp_pdbqt_file):
    import subprocess
    import tempfile
    obabel_cmd = ["obabel", temp_ligand_file.name, "-O", temp_pdbqt_file.name, "-ph", "7.4"]
    subprocess.run(obabel_cmd, check=True)
    return

#Qvina call to dock pdbqt file and save to output_path
def qvina_call(mol, temp_pdbqt_file, ref_centroid, pocket, output_path, chembl_id):
    import tempfile
    import subprocess
    import time
    qvina_cmd = ["qvina", "--receptor", pocket, "--ligand", temp_pdbqt_file, "--num_modes", "1", "--center_x", str(ref_centroid.x),
                "--center_y", str(ref_centroid.y), "--center_z", str(ref_centroid.z), "--size_x", "20", "--size_y", "20", "--size_z",
                "20", "--out", output_path]
    t_init = time.time()
    try: #to catch errors in smina call
        process_result = subprocess.run(qvina_cmd, capture_output=True, text=True)
        #print(f'{time.time() - t_init} seconds to dock')
        output_lines = process_result.stdout.split('\n') #split up all lines
    except subprocess.CalledProcessError as e:
        print(f"QVina command failed with error: {e}")
        return[]
    #print("QVina output lines:", output_lines)
    score_line_index = 0
    for idx, line in enumerate(output_lines):
        if line.startswith('mode'):  #check if the line starts with 'mode'
            score_line_index = idx+3
            break
    score_breakdown = output_lines[score_line_index].split()
    if len(score_breakdown) > 1:
        print(f'{chembl_id} , {score_breakdown[1]}')
        return (chembl_id, score_breakdown[1])
    else:
        print("Unexpected score breakdown format:", score_breakdown)
        return (chembl_id, None)  

