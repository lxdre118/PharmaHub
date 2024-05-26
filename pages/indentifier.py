import streamlit as st
from streamlit_ketcher import st_ketcher
import pubchempy as pcp
import requests
import subprocess
import tempfile
import os
from rdkit import Chem

# Function to fetch identifiers and properties using PubChemPy and CACTUS NCI service
def fetch_identifiers_and_properties(smiles):
    try:
        compound = pcp.get_compounds(smiles, 'smiles')[0]
        inchi = compound.inchi
        inchikey = compound.inchikey
        cas = get_cas_from_cactus(smiles)
        molecular_weight = compound.molecular_weight
        molecular_formula = compound.molecular_formula
        iupac_name = compound.iupac_name
        return inchi, inchikey, cas, molecular_weight, molecular_formula, iupac_name
    except IndexError:
        return ("Not available",) * 6

def get_cas_from_cactus(smiles):
    url = f"http://cactus.nci.nih.gov/chemical/structure/{smiles}/cas"
    response = requests.get(url)
    if response.status_code == 200:
        return response.text.strip()
    else:
        return "Not available"

# Open Babel conversion function
def convert_molecule(input_mol, input_format, output_format, add_hydrogen, generate_3d):
    with tempfile.NamedTemporaryFile(delete=False, suffix=f".{input_format}") as temp_input:
        temp_input.write(input_mol.encode())
        temp_input.flush()

        output_suffix = output_format.lower()
        output_file = temp_input.name.replace(f".{input_format}", f".{output_suffix}")
        
        command = ["obabel", temp_input.name, f"-O{output_file}"]

        if add_hydrogen:
            command.append("-h")
        if generate_3d:
            command.append("--gen3d")

        subprocess.run(command, check=True)

        with open(output_file, "r") as file:
            output_data = file.read()

        os.remove(output_file)
        return output_data

# Function to fetch compound information by CAS number
def fetch_by_cas(cas_number):
    try:
        compound = pcp.get_compounds(cas_number, 'name')[0]
        return compound
    except IndexError:
        return None

# Function to fetch compound information by InChIKey
def fetch_by_inchikey(inchikey):
    try:
        compound = pcp.get_compounds(inchikey, 'inchikey')[0]
        return compound
    except IndexError:
        return None

# Function to fetch compounds by molecular formula
def fetch_by_formula(formula):
    try:
        compounds = pcp.get_compounds(formula, 'formula')
        return compounds
    except IndexError:
        return []

# Function to display compound information
def display_compound_info(compound):
    st.sidebar.markdown(f"**IUPAC Name:** {compound.iupac_name}")
    st.sidebar.markdown(f"**Molecular Formula:** {compound.molecular_formula}")
    st.sidebar.markdown(f"**Molecular Weight:** {compound.molecular_weight}")
    st.sidebar.markdown(f"**Canonical SMILES:** {compound.canonical_smiles}")
    st.sidebar.markdown(f"**InChI:** {compound.inchi}")
    st.sidebar.markdown(f"**InChIKey:** {compound.inchikey}")
    
# Function to convert SMILES to canonical SMILES using RDKit
def convert_to_canonical_smiles(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError("Invalid SMILES string")
    canonical_smiles = Chem.MolToSmiles(mol, canonical=True)
    return canonical_smiles

# Default molecule structure
DEFAULT_MOL = 'C1=CC=C(C(=C1)C(=O)[O-])C(=O)[O-]'

st.title("Molecule Identifier and Converter")

# Initialize session state
if 'drawn_structure' not in st.session_state:
    st.session_state['drawn_structure'] = DEFAULT_MOL
if 'canonical_smiles' not in st.session_state:
    st.session_state['canonical_smiles'] = ''
if 'properties' not in st.session_state:
    st.session_state['properties'] = {}

# Sidebar for input and fetching properties
st.sidebar.header("Input and Fetch Properties")

# Input SMILES and draw structure
smiles_input = st.sidebar.text_input("Input SMILES", DEFAULT_MOL)
if st.sidebar.button("Update Session's Compound"):
    st.session_state['drawn_structure'] = smiles_input

# Canonical SMILES button
if st.sidebar.button("Get Canonical SMILES"):
    try:
        st.session_state['canonical_smiles'] = convert_to_canonical_smiles(st.session_state['drawn_structure'])
    except IndexError:
        st.session_state['canonical_smiles'] = 'Conversion failed.'
        st.error(f"Error converting to canonical SMILES: {e}")

# Display Canonical SMILES if available
if st.session_state['canonical_smiles']:
    st.sidebar.markdown(f"Canonical SMILES: `{st.session_state['canonical_smiles']}`")

# Get Properties button
if st.sidebar.button("Get Properties"):
    try:
        molecule_smiles = pcp.get_compounds(st.session_state['drawn_structure'], 'smiles')[0].canonical_smiles
        properties = fetch_identifiers_and_properties(molecule_smiles)
        st.session_state['properties'] = {
            'Canonical SMILES': molecule_smiles,
            'InChI': properties[0],
            'InChIKey': properties[1],
            'CAS Number': properties[2],
            'Molecular Weight': properties[3],
            'Molecular Formula': properties[4],
            'IUPAC Name': properties[5]
        }
    except IndexError:
        st.session_state['properties'] = {
            'Canonical SMILES': 'Conversion failed',
            'InChI': 'Not available',
            'InChIKey': 'Not available',
            'CAS Number': 'Not available',
            'Molecular Weight': 'Not available',
            'Molecular Formula': 'Not available',
            'IUPAC Name': 'Not available'
        }

# Display properties if available
if st.session_state['properties']:
    st.sidebar.header("Molecule Identifiers and Properties")
    for key, value in st.session_state['properties'].items():
        st.sidebar.markdown(f"**{key}:** `{value}`")

# Ketcher drawing tool
st.header("Molecule Drawing Tool")
drawn_structure = st_ketcher(st.session_state['drawn_structure'], key="ketcher_editor")

# Button to convert drawing to SMILES
if st.button("Convert Drawing to SMILES"):
    try:
        compound = pcp.get_compounds(drawn_structure, 'smiles')
        if compound:
            st.session_state['canonical_smiles'] = compound[0].canonical_smiles
            st.success("Conversion to SMILES successful!")
        else:
            raise IndexError
    except IndexError:
        try:
            st.session_state['canonical_smiles'] = convert_to_canonical_smiles(st.session_state['drawn_structure'])
            if st.session_state['canonical_smiles']:
                st.success("Conversion to SMILES successful!")
            else:
                raise ValueError("Conversion failed.")
        except ValueError as e:
            st.session_state['canonical_smiles'] = 'Conversion failed.'
            st.error(f"Error converting drawing to SMILES: {e}")


# Display the canonical SMILES after conversion
if st.session_state['canonical_smiles']:
    st.markdown(f"Canonical SMILES: `{st.session_state['canonical_smiles']}`")

# Conversion options
st.header("Convert Molecule")
input_format = "smi"  # SMILES format as input
output_format = st.selectbox("Output format", ["mol", "mol2", "sdf"])
add_hydrogen = st.checkbox("Add hydrogen atoms")
generate_3d = st.radio("Generate structure", ["2D", "3D"]) == "3D"

# Convert and display molecule
if st.button("Convert Molecule"):
    converted_molecule = convert_molecule(st.session_state['canonical_smiles'], input_format, output_format, add_hydrogen, generate_3d)
    st.markdown(f"Converted molecule ({output_format}):")
    st.code(converted_molecule)

    # Provide download link for the converted molecule
    st.download_button(
        label="Download converted molecule",
        data=converted_molecule,
        file_name=f"converted_molecule.{output_format}",
        mime="chemical/x-mdl-molfile"
    )

# CAS number input
st.sidebar.header("Fetch by CAS Number")
cas_input = st.sidebar.text_input("Enter CAS Number")
if st.sidebar.button("Fetch by CAS"):
    cas_compound = fetch_by_cas(cas_input)
    if cas_compound:
        display_compound_info(cas_compound)
    else:
        st.sidebar.error("No compound found with this CAS number.")

# InChIKey input
st.sidebar.header("Fetch by InChIKey")
inchikey_input = st.sidebar.text_input("Enter InChIKey")
if st.sidebar.button("Fetch by InChIKey"):
    inchikey_compound = fetch_by_inchikey(inchikey_input)
    if inchikey_compound:
        display_compound_info(inchikey_compound)
    else:
        st.sidebar.error("No compound found with this InChIKey.")

# Molecular formula input
st.sidebar.header("Fetch by Molecular Formula")
formula_input = st.sidebar.text_input("Enter Molecular Formula (e.g., C4H5)")
if st.sidebar.button("Fetch by Formula"):
    formula_compounds = fetch_by_formula(formula_input)
    if formula_compounds:
        st.sidebar.header("Possible Compounds")
        for compound in formula_compounds:
            with st.sidebar.expander(f"Compound ID: {compound.cid}"):
                display_compound_info(compound)
    else:
        st.sidebar.error("No compounds found with this molecular formula.")