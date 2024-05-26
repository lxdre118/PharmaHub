import streamlit as st
import py3Dmol
from stmol import showmol

# Set page configuration
st.set_page_config(layout="wide", page_title="PharmaHub", page_icon="icon.jpg")

# Hide Streamlit style elements
hide_streamlit_style = """
    <style>
    MainMenu {visibility: hidden;}
    footer {visibility: hidden;}
    body {font-size: 100px;}
    </style>
"""
st.markdown(hide_streamlit_style, unsafe_allow_html=True)

def main():
    st.sidebar.success("Choose the tool you want to use")
    st.sidebar.info("Note: The platform is currently experiencing limited network bandwidth, so access may be slower than usual. Please be patient.")
    st.sidebar.info("Note: The platform is in the trial stage. Please avoid uploading sensitive information as security measures may not be fully implemented.")

    col1, col2 = st.columns(2)

    # Load and display custom protein structure
    pdb_cartoon = open("5v3y.pdb", 'r').read()
    view = py3Dmol.view(width=400, height=400)
    view.addModel(pdb_cartoon, "pdb")
    view.setStyle({"cartoon": {"color": "lightblue"}})
    view.addStyle({"elem": 'C', "hetflag": True}, {"stick": {"color": "orange", "radius": 0.15}})
    view.addStyle({"elem": 'O'}, {"stick": {"color": "red", "radius": 0.15}})
    view.addStyle({"elem": 'N'}, {"stick": {"color": "blue", "radius": 0.15}})
    view.zoomTo()
    view.setBackgroundColor("#f0f0f0")

    with col1:
        showmol(view, height=400, width=500)
        st.markdown("**Crystal Structure of Mtb Pks13 Thioesterase domain** (5V3Y) by [RCSBPDB](https://www.rcsb.org/structure/5V3Y)")

    with col2:
        st.header('PharmaHub', divider='rainbow')
        st.markdown("### A Comprehensive, Free, and Fast-Responsive Ligand-Based Drug Design Platform")
        st.markdown("Created by Wang Shizun, Shenyang Pharmaceutical University")

    st.markdown("***")
    st.markdown('<style>body { font-size: 16px; }</style>', unsafe_allow_html=True)
    st.markdown("Welcome to our molecular design and data analysis platform! We provide researchers and data scientists with a comprehensive set of advanced tools designed to streamline your workflow, enhance data accessibility, and achieve more accurate and in-depth research results. By leveraging our platform, you can significantly improve the focus and clarity of your research work, optimize workflows, and attain more insightful research outcomes. We are dedicated to offering indispensable support for drug development, pharmacological research, and molecular simulation. We invite you to experience and utilize these powerful features to elevate your research to the next level.")
    st.markdown("**Key Features:**")
    st.markdown("1. **FDA_ChemInsight:**  The app offers powerful data loading and caching for efficient analysis, flexible dataset selection, and automatic calculation of essential molecular descriptors. Interactive filtering and advanced data visualization options, including scatter plots and 3D scatter plots, enhance the focus and clarity of research. The tool's molecular grid display and customizable visuals provide a user-friendly interface, making it indispensable for drug development, pharmacology research, and molecular simulation. By utilizing FDA ChemInsight, researchers can significantly streamline their workflow, improve data accessibility, and achieve more accurate and insightful results.")
    st.markdown("2. **Pharmafilter:** The app provides specialized tools for analyzing and visualizing user-uploaded drug data. It offers efficient file upload and validation, automatic calculation of molecular descriptors, and scaffold molecule analysis. Interactive filtering and advanced data visualization options, including scatter plots, histograms, and 3D scatter plots, enhance the clarity and focus of research. By using ChemInsight, researchers can streamline their workflow, improve data accessibility, and achieve more accurate and insightful results, making it indispensable for drug development, pharmacology research, and molecular simulation.")
    st.markdown("3. **Drug Machine Learning:** The Dataset Analysis platform provides a zero-code solution for machine learning modeling, enabling users to upload datasets and automatically calculate molecular descriptors or molecular fingerprints. It supports various regression models, including linear, tree-based, and neural network approaches, offering a comprehensive set of tools for data analysis and visualization. With features like interactive filtering, advanced visualizations, and easy-to-use model evaluation metrics, this platform significantly streamlines the workflow, making it essential for researchers and data scientists looking to perform sophisticated data analysis without extensive coding.")
    st.markdown("Ligand-Based Drug Design (LBDD) saves time and costs by utilizing known active ligands to reduce experimental screening, increases drug candidate success rates while lowering clinical trial failure risks, is versatile enough for targets without crystal structures, and optimizes compound selectivity and specificity to minimize side effects.")
    st.markdown("***")

if __name__ == '__main__':
    main()


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