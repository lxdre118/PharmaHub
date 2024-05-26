import os
from flask import redirect, request, url_for, send_from_directory, jsonify, render_template
from werkzeug.utils import secure_filename
from rdkit import Chem

ALLOWED_EXTENSIONS = {'sdf', 'smi', 'txt'}

def allowed_file(filename):
    """
    Check if the uploaded file has an allowed extension.

    Parameters:
    filename (str): The name of the file to check.

    Returns:
    bool: True if the file has an allowed extension, False otherwise.
    """
    return '.' in filename and filename.rsplit('.', 1)[1].lower() in ALLOWED_EXTENSIONS

def process_structure_file(filepath):
    """
    Process the structure file to extract stepwise rings.

    Parameters:
    filepath (str): The path to the uploaded structure file.

    Returns:
    str: The name of the output file containing the stepwise ring extraction results.
    """
    output_filename = os.path.splitext(os.path.basename(filepath))[0] + "_processed.sdf"
    output_path = os.path.join('static/downloads', output_filename)
    
    # Determine the file type and create the appropriate RDKit supplier
    suppl = None
    file_extension = os.path.splitext(filepath)[1].lower()
    
    if file_extension == ".sdf":
        suppl = Chem.SDMolSupplier(filepath)
    elif file_extension in [".smi", ".txt"]:
        suppl = Chem.SmilesMolSupplier(filepath, delimiter="\t", titleLine=False)
    else:
        raise ValueError("Unsupported file format. Only SDF and SMILES files are supported.")
    
    # Create an SDWriter to write the processed molecules to the output file
    writer = Chem.SDWriter(output_path)
    
    # Process each molecule in the supplier and add the stepwise rings information
    for i, mol in enumerate(suppl):
        if mol is None:
            continue
        smiles = Chem.MolToSmiles(mol)
        steps = extract_stepwise_rings(smiles)
        
        # Add stepwise rings information as a property to the molecule
        for j, step in enumerate(steps):
            mol.SetProp(f"Step {j + 1}", ', '.join(step))
        
        # Write the processed molecule to the output file
        writer.write(mol)
    
    # Close the writer
    writer.close()
    
    return output_filename

def extract_rings(mol):
    """
    Extract rings from a molecule.

    Parameters:
    mol (Mol): An RDKit molecule object.

    Returns:
    list: A list of SMILES strings representing the rings in the molecule.
    """
    ssr = Chem.GetSymmSSSR(mol)
    rings = []
    for ring in ssr:
        ring_smiles = Chem.MolFragmentToSmiles(mol, atomsToUse=list(ring))
        rings.append(ring_smiles)
    return rings

def remove_ring_atoms(mol, ring_atoms):
    """
    Remove atoms that are part of a ring from a molecule.

    Parameters:
    mol (Mol): An RDKit molecule object.
    ring_atoms (list): A list of atom indices that are part of the ring to be removed.

    Returns:
    Mol: The molecule with the ring atoms removed.
    """
    for atom_idx in ring_atoms:
        mol.GetAtomWithIdx(atom_idx).SetAtomicNum(0)
    mol = Chem.DeleteSubstructs(mol, Chem.MolFromSmarts('[#0]'))
    mol = Chem.RemoveHs(mol)
    return mol

def extract_stepwise_rings(smiles):
    """
    Extract rings from a molecule in a stepwise manner until no more rings are left.

    Parameters:
    smiles (str): The SMILES string of the molecule.

    Returns:
    list: A list of lists, where each inner list contains the SMILES strings of rings extracted in one step.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError("Invalid SMILES string")
    
    steps = []
    
    while True:
        rings = extract_rings(mol)
        if not rings:
            break
        steps.append(rings)
        
        ssr = Chem.GetSymmSSSR(mol)
        largest_ring = max(ssr, key=len)
        
        mol = remove_ring_atoms(mol, largest_ring)
    
    return steps

def handle_file_upload(app, request):
    """
    Handle file upload and initiate processing if a valid file is uploaded.

    Parameters:
    app (Flask): The Flask application instance.
    request (Request): The request object.

    Returns:
    Response: The response object, either rendering the upload form or returning JSON with the result.
    """
    if request.method == 'POST':
        if 'file' not in request.files:
            return redirect(request.url)
        file = request.files['file']
        if file.filename == '':
            return redirect(request.url)
        if file and allowed_file(file.filename):
            filename = secure_filename(file.filename)
            filepath = os.path.join(app.config['UPLOAD_FOLDER'], filename)
            file.save(filepath)
            output_filename = process_structure_file(filepath)
            return jsonify({'filename': output_filename})
    return render_template('index.html')

def handle_file_download(app, filename):
    """
    Handle file download by serving the processed file from the downloads directory.

    Parameters:
    app (Flask): The Flask application instance.
    filename (str): The name of the file to be downloaded.

    Returns:
    Response: The response object for sending the file.
    """
    return send_from_directory(app.config['DOWNLOAD_FOLDER'], filename)
