import mols2grid
import pandas as pd
import streamlit as st
import streamlit.components.v1 as components
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem.Lipinski import NumRotatableBonds
import plotly.express as px
from streamlit_echarts import st_echarts

# 页面标题和描述
st.title("FDA ChemInsight: Specialized software for analyzing and visualizing FDA-approved drug data")
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

# 缓存数据集
@st.cache_data
def load_dataset(dataset_path):
    """加载并缓存数据集"""
    df = pd.read_csv(dataset_path, sep="\t")
    return df

# 选择数据集
selected_dataset = st.radio(
    "Select Dataset",
    ("Full Structure", "Molecular Scaffold"),
    help="Choose the dataset you want to analyze."
)

# 根据选择的数据集加载不同的文件
if selected_dataset == "Full Structure":
    dataset_path = "./data/fda_drug.txt"  # 使用全结构数据集文件路径
    st.write("Analyzing the structure and property distribution of FDA-approved drugs.")
elif selected_dataset == "Molecular Scaffold":
    dataset_path = "./data/scaffold.txt"  # 使用分子环状骨架数据集文件路径
    st.write("Analyzing the molecular scaffolds of drugs to inspire scaffold hopping.")


# 加载数据集
df = load_dataset(dataset_path)

# 确认列名后进行重命名（假设 'smiles' 和 'generic_name' 是数据集中的列）
df.rename(columns={"smiles": "SMILES", "generic_name": "Name"}, inplace=True)

from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski

def calc_descriptors(smiles_string):
    if pd.isnull(smiles_string):
        return pd.Series({
            "MW": float('nan'),
            "LogP": float('nan'),
            "NumHDonors": float('nan'),
            "NumHAcceptors": float('nan'),
            "TPSA": float('nan'),
            "NumRotatableBonds": float('nan'),
            "NumAromaticRings": float('nan'),
            "NumAliphaticRings": float('nan'),
            "NumSaturatedRings": float('nan'),
            "NumHeteroatoms": float('nan'),
            "NumValenceElectrons": float('nan'),
            "NumRadicalElectrons": float('nan'),
            "QED": float('nan'),
            "MolMR": float('nan'),  # Molar Refractivity
            "HeavyAtomCount": float('nan'),
            "NHOHCount": float('nan'),
            "NOCount": float('nan'),
            "FractionCSP3": float('nan'),
            "RingCount": float('nan'),
        })
    
    mol = Chem.MolFromSmiles(smiles_string)
    if mol is None:
        return pd.Series({
            "MW": float('nan'),
            "LogP": float('nan'),
            "NumHDonors": float('nan'),
            "NumHAcceptors": float('nan'),
            "TPSA": float('nan'),
            "NumRotatableBonds": float('nan'),
            "NumAromaticRings": float('nan'),
            "NumAliphaticRings": float('nan'),
            "NumSaturatedRings": float('nan'),
            "NumHeteroatoms": float('nan'),
            "NumValenceElectrons": float('nan'),
            "NumRadicalElectrons": float('nan'),
            "QED": float('nan'),
            "MolMR": float('nan'),  # Molar Refractivity
            "HeavyAtomCount": float('nan'),
            "NHOHCount": float('nan'),
            "NOCount": float('nan'),
            "FractionCSP3": float('nan'),
            "RingCount": float('nan'),
        })
    
    descriptors = {
        "MW": round(Descriptors.MolWt(mol), 3),
        "LogP": round(Descriptors.MolLogP(mol), 3),
        "NumHDonors": round(Descriptors.NumHDonors(mol), 3),
        "NumHAcceptors": round(Descriptors.NumHAcceptors(mol), 3),
        "TPSA": round(Descriptors.TPSA(mol), 3),
        "NumRotatableBonds": round(Lipinski.NumRotatableBonds(mol), 3),
        "NumAromaticRings": round(Descriptors.NumAromaticRings(mol), 3),
        "NumAliphaticRings": round(Descriptors.NumAliphaticRings(mol), 3),
        "NumSaturatedRings": round(Descriptors.NumSaturatedRings(mol), 3),
        "NumHeteroatoms": round(Descriptors.NumHeteroatoms(mol), 3),
        "NumValenceElectrons": round(Descriptors.NumValenceElectrons(mol), 3),
        "NumRadicalElectrons": round(Descriptors.NumRadicalElectrons(mol), 3),
        "QED": round(Descriptors.qed(mol), 3),
        "MolMR": round(Descriptors.MolMR(mol), 3),  # Molar Refractivity
        "HeavyAtomCount": round(Descriptors.HeavyAtomCount(mol), 3),
        "NHOHCount": round(Descriptors.NHOHCount(mol), 3),
        "NOCount": round(Descriptors.NOCount(mol), 3),
        "FractionCSP3": round(Descriptors.FractionCSP3(mol), 3),
        "RingCount": round(Descriptors.RingCount(mol), 3),
    }
    return pd.Series(descriptors)


# Calculate molecular descriptors
descriptors_df = df["SMILES"].apply(calc_descriptors)
df = pd.concat([df, descriptors_df], axis=1)

# Sidebar settings
st.sidebar.header('Parameter Settings')
st.sidebar.write('*Note: Display compounds with values less than the following thresholds*')
weight_cutoff = st.sidebar.slider(
    label="Molecular Weight",
    min_value=0,
    max_value=1000,
    value=500,
    step=10,
    key='weight_cutoff'
)
logp_cutoff = st.sidebar.slider(
    label="LogP",
    min_value=-10,
    max_value=10,
    value=5,
    step=1,
    key='logp_cutoff'
)
NumHDonors_cutoff = st.sidebar.slider(
    label="Number of Hydrogen Donors",
    min_value=0,
    max_value=15,
    value=5,
    step=1,
    key='NumHDonors_cutoff'
)
NumHAcceptors_cutoff = st.sidebar.slider(
    label="Number of Hydrogen Acceptors",
    min_value=0,
    max_value=20,
    value=10,
    step=1,
    key='NumHAcceptors_cutoff'
)

# Filter dataset
df_result = df[df["MW"] < weight_cutoff]
df_result2 = df_result[df_result["LogP"] < logp_cutoff]
df_result3 = df_result2[df_result2["NumHDonors"] < NumHDonors_cutoff]
df_result4 = df_result3[df_result3["NumHAcceptors"] < NumHAcceptors_cutoff]

st.write(f"Number of compounds meeting the criteria: {df_result4.shape[0]}")
st.write(df_result4)

# Display results
raw_html = mols2grid.display(df_result4,
                            subset=["img", "Name"],
                            mapping={"SMILES": "smiles", "Name": "generic_name"})._repr_html_()
components.html(raw_html, width=800, height=650, scrolling=False)

# Horizontal line
st.write('<hr style="border: 1px solid #f0f2f6;">', unsafe_allow_html=True)

# Visualization options
st.header("Visualization Options")

# Automatically detect column names, excluding 'smiles' (case insensitive) and unnamed columns
column_options = [col for col in df_result4.columns if col.lower() != 'smiles' and not col.startswith('Unnamed')]

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
    fig = px.scatter(df_result4, x=x_axis, y=y_axis, color=df_result4[x_axis], hover_data=["Name"], title=f"Scatter Plot of {x_axis} vs {y_axis}")
    st.plotly_chart(fig, use_container_width=True)
elif chart_type == "Histogram":
    fig = px.histogram(df_result4, x=x_axis, color=df_result4[x_axis], title=f"Histogram of {x_axis}")
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
    fig = px.scatter_3d(df_result4, x=x_axis, y=y_axis, z=z_axis, color=df_result4[x_axis], hover_data=["Name"], title=f"3D Scatter Plot of {x_axis}, {y_axis}, and {z_axis}")
    st.plotly_chart(fig, use_container_width=True)
elif chart_type == "Heatmap":
    y_axis = st.selectbox(
        "Select Y-axis Variable",
        column_options,
        key="y_axis_heatmap"
    )
    fig = px.density_heatmap(df_result4, x=x_axis, y=y_axis, title=f"Heatmap of {x_axis} vs {y_axis}")
    st.plotly_chart(fig, use_container_width=True)
elif chart_type == "Bar Chart":
    y_axis = st.selectbox(
        "Select Y-axis Variable",
        column_options,
        key="y_axis_bar"
    )
    fig = px.bar(df_result4, x=x_axis, y=y_axis, color=df_result4[x_axis], title=f"Bar Chart of {x_axis} vs {y_axis}")
    st.plotly_chart(fig, use_container_width=True)
elif chart_type == "Box Plot":
    y_axis = st.selectbox(
        "Select Y-axis Variable",
        column_options,
        key="y_axis_box"
    )
    fig = px.box(df_result4, x=x_axis, y=y_axis, color=df_result4[x_axis], title=f"Box Plot of {x_axis} vs {y_axis}")
    st.plotly_chart(fig, use_container_width=True)
elif chart_type == "Line Chart":
    y_axis = st.selectbox(
        "Select Y-axis Variable",
        column_options,
        key="y_axis_line"
    )
    fig = px.line(df_result4, x=x_axis, y=y_axis, color=df_result4[x_axis], title=f"Line Chart of {x_axis} vs {y_axis}")
    st.plotly_chart(fig, use_container_width=True)
elif chart_type == "Pie Chart":
    fig = px.pie(df_result4, names=x_axis, title=f"Pie Chart of {x_axis}")
    st.plotly_chart(fig, use_container_width=True)


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
