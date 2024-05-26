import pandas as pd
import streamlit as st
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit.Chem import MolFromSmiles, MolToSmiles
import plotly.express as px
import mols2grid
import streamlit.components.v1 as components
from rdkit.Chem.rdMolDescriptors import CalcMolFormula

# 页面标题和描述
# 隐藏made with streamlit
hide_streamlit_style = """
            <style>
            #MainMenu {visibility: hidden;}
            footer {visibility: hidden;}
            </style>
            """
st.markdown(hide_streamlit_style, unsafe_allow_html=True)

st.title("ChemInsight: Specialized software for analyzing and visualizing user-upload drug data")
st.markdown("""
<style>
    .sidebar .sidebar-content {
        background-color: #f0f2f6;
    }
    .reportview-container .main .block-container {
        max-width: 80%;
        padding: 1rem;
        background-color: #ffffff;
        border-radius: 10px;
    }
</style>
""", unsafe_allow_html=True)
# Example data
example_data = {
    "smiles": ["CCO", "CC(=O)OC1=CC=CC=C1C(=O)O", "C1CCCCC1"],
    "generic_name": ["Ethanol", "Aspirin", "Cyclohexane"],
    "Pharma_activity(optional)": ["0.79", "2.25", "0.03"]
}
example_df = pd.DataFrame(example_data)

# Convert example data to CSV for download
example_csv = example_df.to_csv(index=False)


# Step 1: File upload
uploaded_file = st.file_uploader("Please upload a dataset (CSV format)", type=["csv"])


if uploaded_file is not None:
    try:
        df = pd.read_csv(uploaded_file)
        # Check if required columns are present
        if "smiles" in df.columns and "generic_name" in df.columns:
            df.rename(columns={"smiles": "SMILES", "generic_name": "Name"}, inplace=True)
            st.success("Dataset uploaded successfully!")
            df = df[df["SMILES"].notna()]    
        
            # Step 2: Calculate molecular descriptors for original molecules
            def calc_descriptors(smiles_string):
                mol = Chem.MolFromSmiles(smiles_string)
                descriptors = {
                    "MW": round(Descriptors.MolWt(mol), 3),
                    "LogP": round(Descriptors.MolLogP(mol), 3),
                    "NumHDonors": round(Descriptors.NumHDonors(mol), 3),
                    "NumHAcceptors": round(Descriptors.NumHAcceptors(mol), 3),
                    "TPSA": round(Descriptors.TPSA(mol), 3),
                    "NumRotatableBonds": round(Descriptors.NumRotatableBonds(mol), 3),
                    "NumAromaticRings": round(Descriptors.NumAromaticRings(mol), 3),
                    "NumAliphaticRings": round(Descriptors.NumAliphaticRings(mol), 3),
                    "NumSaturatedRings": round(Descriptors.NumSaturatedRings(mol), 3),
                    "NumHeteroatoms": round(Descriptors.NumHeteroatoms(mol), 3),
                    "NumValenceElectrons": round(Descriptors.NumValenceElectrons(mol), 3),
                    "NumRadicalElectrons": round(Descriptors.NumRadicalElectrons(mol), 3),
                    "QED": round(Descriptors.qed(mol), 3),
                    "MolMR": round(Descriptors.MolMR(mol), 3),
                    "FractionCSP3": round(Descriptors.FractionCSP3(mol), 3),
                    "RingCount": round(Descriptors.RingCount(mol), 3),
                    "HeavyAtomCount": round(Descriptors.HeavyAtomCount(mol), 3),
                    "NHOHCount": round(Descriptors.NHOHCount(mol), 3),
                    "NOCount": round(Descriptors.NOCount(mol), 3),
                    "NumAliphaticCarbocycles": round(Descriptors.NumAliphaticCarbocycles(mol), 3),
                    "NumAliphaticHeterocycles": round(Descriptors.NumAliphaticHeterocycles(mol), 3),
                    "NumAromaticCarbocycles": round(Descriptors.NumAromaticCarbocycles(mol), 3),
                    "NumAromaticHeterocycles": round(Descriptors.NumAromaticHeterocycles(mol), 3),
                    "NumSaturatedCarbocycles": round(Descriptors.NumSaturatedCarbocycles(mol), 3),
                    "NumSaturatedHeterocycles": round(Descriptors.NumSaturatedHeterocycles(mol), 3),
                    "PEOE_VSA1": round(Descriptors.PEOE_VSA1(mol), 3),
                    "PEOE_VSA2": round(Descriptors.PEOE_VSA2(mol), 3)
                }
                return pd.Series(descriptors)
            
            # Calculate descriptors for original molecules
            descriptors_df = df["SMILES"].apply(calc_descriptors)
            df = pd.concat([df, descriptors_df], axis=1)
            
            # Calculate descriptors for scaffold molecules
            unique_scaffold_smiles_set = set()
            scaffold_smiles_list = []
            scaffold_descriptors = []
            for smiles in df["SMILES"]:
                scaffold_smiles = MurckoScaffold.MurckoScaffoldSmilesFromSmiles(smiles)
                if scaffold_smiles not in unique_scaffold_smiles_set:
                    unique_scaffold_smiles_set.add(scaffold_smiles)
                    scaffold_smiles_list.append(scaffold_smiles)
                    scaffold_descriptors.append(calc_descriptors(scaffold_smiles))
            scaffold_df = pd.DataFrame(scaffold_descriptors)
            scaffold_df["SMILES"] = scaffold_smiles_list
            scaffold_df["Name"] = scaffold_df["SMILES"].apply(lambda x: CalcMolFormula(MolFromSmiles(x)))
            columns = ["Name"] + [col for col in scaffold_df.columns if col != "Name"]
            scaffold_df = scaffold_df[columns]
            scaffold_df = scaffold_df[scaffold_df["MW"] != 0]    
            # Step 3: User selection for dataset analysis
            st.markdown("### Choose Dataset for Analysis")
            option = st.selectbox("Select dataset to analyze:", ("Original Molecules", "Scaffold Molecules"), key="dataset_select")
            
            if option == "Original Molecules":
                selected_df = df
                st.markdown("### Original Molecules")
            else:
                selected_df = scaffold_df
                st.markdown("### Scaffold")
            
            # Step 4: Sidebar settings for filtering
            st.sidebar.header('Parameter Settings')
            st.sidebar.write('*Note: Display compounds with values less than the following thresholds*')
            weight_cutoff = st.sidebar.slider(
                label="Molecular Weight",
                min_value=0,
                max_value=1000,
                value=500,
                step=10,
                key="weight_cutoff"
            )
            logp_cutoff = st.sidebar.slider(
                label="LogP",
                min_value=-10,
                max_value=10,
                value=5,
                step=1,
                key="logp_cutoff"
            )
            NumHDonors_cutoff = st.sidebar.slider(
                label="Number of Hydrogen Donors",
                min_value=0,
                max_value=15,
                value=5,
                step=1,
                key="NumHDonors_cutoff"
            )
            NumHAcceptors_cutoff = st.sidebar.slider(
                label="Number of Hydrogen Acceptors",
                min_value=0,
                max_value=20,
                value=10,
                step=1,
                key="NumHAcceptors_cutoff"
            )

            # Filter dataset based on user input
            filtered_df = selected_df[selected_df["MW"] < weight_cutoff]
            filtered_df = filtered_df[filtered_df["LogP"] < logp_cutoff]
            filtered_df = filtered_df[filtered_df["NumHDonors"] < NumHDonors_cutoff]
            filtered_df = filtered_df[filtered_df["NumHAcceptors"] < NumHAcceptors_cutoff]

            st.write(f"Number of compounds meeting the criteria: {filtered_df.shape[0]}")
            st.write(filtered_df)

            # Display results using mols2grid
            raw_html = mols2grid.display(filtered_df)._repr_html_()
            components.html(raw_html, width=800, height=650, scrolling=False)

            # Visualization options
            st.header("Visualization Options")


            # Automatically detect column names, excluding 'smiles' (case insensitive) and unnamed columns
            column_options = [col for col in filtered_df.columns if col.lower() != 'smiles' and not col.startswith('Unnamed')]

            # Select chart type
            chart_type = st.selectbox(
                "Select Chart Type",
                ["Scatter Plot", "Histogram", "3D Scatter Plot", "Heatmap", "Bar Chart", "Box Plot", "Line Chart", "Pie Chart"],
                key="chart_type"
            )

            # Select factors to display
            x_axis = st.selectbox(
                "Select X-axis Variable",
                column_options,
                key="x_axis"
            )

            if chart_type == "Scatter Plot":
                y_axis = st.selectbox(
                    "Select Y-axis Variable",
                    column_options,
                    key="y_axis"
                )
                fig = px.scatter(filtered_df, x=x_axis, y=y_axis, color=filtered_df[x_axis], hover_data=["Name"], title=f"Scatter Plot of {x_axis} vs {y_axis}")
                st.plotly_chart(fig, use_container_width=True)
            elif chart_type == "Histogram":
                fig = px.histogram(filtered_df, x=x_axis, color=filtered_df[x_axis], title=f"Histogram of {x_axis}")
                st.plotly_chart(fig, use_container_width=True)
            elif chart_type == "3D Scatter Plot":
                y_axis = st.selectbox(
                    "Select Y-axis Variable",
                    column_options,
                    key="y_axis_3d"
                )
                z_axis = st.selectbox(
                    "Select Z-axis Variable",
                    column_options,
                    key="z_axis_3d"
                )
                fig = px.scatter_3d(filtered_df, x=x_axis, y=y_axis, z=z_axis, color=filtered_df[x_axis], hover_data=["Name"], title=f"3D Scatter Plot of {x_axis}, {y_axis}, and {z_axis}")
                st.plotly_chart(fig, use_container_width=True)
            elif chart_type == "Heatmap":
                y_axis = st.selectbox(
                    "Select Y-axis Variable",
                    column_options,
                    key="y_axis_heatmap"
                )
                fig = px.density_heatmap(filtered_df, x=x_axis, y=y_axis, title=f"Heatmap of {x_axis} vs {y_axis}")
                st.plotly_chart(fig, use_container_width=True)
            elif chart_type == "Bar Chart":
                y_axis = st.selectbox(
                    "Select Y-axis Variable",
                    column_options,
                    key="y_axis_bar"
                )
                fig = px.bar(filtered_df, x=x_axis, y=y_axis, color=filtered_df[x_axis], title=f"Bar Chart of {x_axis} vs {y_axis}")
                st.plotly_chart(fig, use_container_width=True)
            elif chart_type == "Box Plot":
                y_axis = st.selectbox(
                    "Select Y-axis Variable",
                    column_options,
                    key="y_axis_box"
                )
                fig = px.box(filtered_df, x=x_axis, y=y_axis, color=filtered_df[x_axis], title=f"Box Plot of {x_axis} vs {y_axis}")
                st.plotly_chart(fig, use_container_width=True)
            elif chart_type == "Line Chart":
                y_axis = st.selectbox(
                    "Select Y-axis Variable",
                    column_options,
                    key="y_axis_line"
                )
                fig = px.line(filtered_df, x=x_axis, y=y_axis, color=filtered_df[x_axis], title=f"Line Chart of {x_axis} vs {y_axis}")
                st.plotly_chart(fig, use_container_width=True)
            elif chart_type == "Pie Chart":
                fig = px.pie(filtered_df, names=x_axis, title=f"Pie Chart of {x_axis}")
                st.plotly_chart(fig, use_container_width=True)

        else:
            st.error("The dataset must contain 'smiles' and 'generic_name' columns.")
    except Exception as e:
        st.error(f"Failed to read the file: {e}")
else:
    st.markdown("### Example Data Format")
    st.write("Please upload a CSV file containing the following columns:")
    st.table(example_df)
    st.image('./source/chemInsight.webp', use_column_width=True)
    
# Footer implementation
footer = """
<style>
.footer {
    position: fixed;
    left: 0;
    bottom: 0;
    width: 100%;
    background-color: #2c3e50;
    color: white;
    text-align: center;
    padding: 10px 0;
    font-family: 'Arial', sans-serif;
    box-shadow: 0 -1px 10px rgba(0, 0, 0, 0.1);
    border-top: 1px solid #bdc3c7;
}
.footer p {
    margin: 0;
    font-size: 14px;
}
.footer a {
    color: #ecf0f1;
    text-decoration: none;
    margin: 0 10px;
}
.footer a:hover {
    color: #3498db;
}
</style>
<div class='footer'>
  <p>&copy; 2024 Wang Shizun. All rights reserved | <a href="https://beian.miit.gov.cn/" target="_blank">鲁ICP备2022040118号-1</a> |  shizunwang@foxmail.com</a></p>
</div>
"""
st.markdown(footer, unsafe_allow_html=True)
