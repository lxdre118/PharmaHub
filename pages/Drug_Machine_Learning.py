import streamlit as st
import pandas as pd
import numpy as np
import altair as alt
from rdkit import Chem
from rdkit.Chem import Descriptors
from sklearn.linear_model import LinearRegression, Ridge, Lasso, ElasticNet
from sklearn.ensemble import RandomForestRegressor, GradientBoostingRegressor
from sklearn.neighbors import KNeighborsRegressor
from sklearn.tree import DecisionTreeRegressor
from sklearn.neural_network import MLPRegressor
from xgboost import XGBRegressor
from lightgbm import LGBMRegressor
from sklearn.svm import SVR
from sklearn.metrics import mean_squared_error

# 页面标题和描述
# 隐藏made with streamlit
hide_streamlit_style = """
            <style>
            #MainMenu {visibility: hidden;}
            footer {visibility: hidden;}
            </style>
            """
st.markdown(hide_streamlit_style, unsafe_allow_html=True)

st.header('Dataset Analysis: Zero-code machine learning modeling platform', divider='rainbow')
st.header('Dataset Analysis :blue[cool] :sunglasses:')

# st.title("Dataset Analysis: Zero-code machine learning modeling platform")

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

# Define function to calculate 1D representation
def calculate_1dmr_repr(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return [np.nan] * 13  # Return NaNs if SMILES cannot be parsed
    else:
        mol_weight = Descriptors.MolWt(mol)
        log_p = Descriptors.MolLogP(mol)
        num_h_donors = Descriptors.NumHDonors(mol)
        num_h_acceptors = Descriptors.NumHAcceptors(mol)
        tpsa = Descriptors.TPSA(mol)
        num_rotatable_bonds = Descriptors.NumRotatableBonds(mol)
        num_aromatic_rings = Descriptors.NumAromaticRings(mol)
        num_aliphatic_rings = Descriptors.NumAliphaticRings(mol)
        num_saturated_rings = Descriptors.NumSaturatedRings(mol)
        num_heteroatoms = Descriptors.NumHeteroatoms(mol)
        num_valence_electrons = Descriptors.NumValenceElectrons(mol)
        num_radical_electrons = Descriptors.NumRadicalElectrons(mol)
        qed = Descriptors.qed(mol)
        return [mol_weight, log_p, num_h_donors, num_h_acceptors, tpsa, num_rotatable_bonds, num_aromatic_rings,
                num_aliphatic_rings, num_saturated_rings, num_heteroatoms, num_valence_electrons, num_radical_electrons, qed]



# Example data
example_data = {
    "smiles": ["CCO", "CC(=O)OC1=CC=CC=C1C(=O)O", "C1CCCCC1"],
    "Pharma_activity(optional)": ["0.79", "2.25", "0.03"]
}
example_df = pd.DataFrame(example_data)

uploaded_file = st.file_uploader("Upload your dataset (CSV file)", type="csv")
if uploaded_file is not None:
    try:
        data = pd.read_csv(uploaded_file)
        st.subheader("Original data")
        st.write(data)

        # 重命名列
        data.columns = ["SMILES", "TARGET"]
        # Sidebar for dataset split ratio and proceeding with calculations
        st.sidebar.header("Machine Learning Settings")
        train_fraction = st.sidebar.slider("Select training set fraction", 0.1, 0.9, 0.8, 0.05)
        
        selected_method = st.sidebar.radio("Select Method", ["Molecular Descriptors", "Molecular Fingerprints"])


        proceed = st.sidebar.checkbox("Proceed with calculations")


        # Calculate 1D-QSAR representations
        st.write("Calculating molecular descriptors...")
        data["1d_mr"] = data["SMILES"].apply(calculate_1dmr_repr)

        # 将数据集分为训练数据集和测试数据集
        train_data = data.sample(frac=train_fraction, random_state=1)
        st.subheader("Train data")
        st.write(train_data)
        train_data.to_csv("hERG_train.csv", index=False)
        test_data = data.drop(train_data.index)
        st.subheader("Test data")
        st.write(test_data)
        test_data.to_csv("hERG_test.csv", index=False)

        # 设定训练/测试目标
        train_y = np.array(train_data["TARGET"].values.tolist())
        test_y = np.array(test_data["TARGET"].values.tolist())

        # 两个数据集分布可视化结果
        st.subheader("Visualization")
        hist_train = alt.Chart(train_data).mark_bar().encode(
            alt.X("TARGET", bin=alt.Bin(maxbins=20), title="pIC50"),
            alt.Y("count()", title="Count"),
            color=alt.value("blue")
        ).properties(
            width=500,
            height=300,
            title="Distribution of pIC50 in Train Data"
        )

        hist_test = alt.Chart(test_data).mark_bar().encode(
            alt.X("TARGET", bin=alt.Bin(maxbins=20), title="pIC50"),
            alt.Y("count()", title="Count"),
            color=alt.value("orange")
        ).properties(
            width=500,
            height=300,
            title="Distribution of pIC50 in Test Data"
        )

        st.altair_chart(hist_train)
        st.altair_chart(hist_test)

        if proceed:
            # Interface for method selection
            # Apply selected method

            #当选用分子描述符
            if selected_method == "Molecular Descriptors":
                st.write("You selected: Molecular Descriptors")

                # Convert to numpy arrays for training
                train_x = np.array(train_data["1d_mr"].tolist())
                train_y = np.array(train_data["TARGET"].tolist())
                test_x = np.array(test_data["1d_mr"].tolist())
                test_y = np.array(test_data["TARGET"].tolist())

                # Define regressors
                regressors = [
                    ("Linear Regression", LinearRegression()),
                    ("Ridge Regression", Ridge(random_state=42)),
                    ("Lasso Regression", Lasso(random_state=42)),
                    ("ElasticNet Regression", ElasticNet(random_state=42)),
                    ("Support Vector", SVR()),
                    ("K-Nearest Neighbors", KNeighborsRegressor()),
                    ("Decision Tree", DecisionTreeRegressor(random_state=42)),
                    ("Random Forest", RandomForestRegressor(random_state=42)),
                    ("Gradient Boosting", GradientBoostingRegressor(random_state=42)),
                    ("XGBoost", XGBRegressor(random_state=42)),
                    ("LightGBM", LGBMRegressor(random_state=42)),
                    ("Multi-layer Perceptron", MLPRegressor(
                        hidden_layer_sizes=(128, 64, 32),
                        learning_rate_init=0.0001,
                        activation='relu', solver='adam',
                        max_iter=10000, random_state=42)),
                ]

                results = {}

                # Train and evaluate models
                for name, regressor in regressors:
                    # 训练回归器
                    regressor.fit(train_x, train_y)
                    # 预测训练和测试数据
                    pred_train_y = regressor.predict(train_x)
                    pred_test_y = regressor.predict(test_x)
                    # 将预测结果添加到数据中
                    train_data[f"1D-QSAR-{name}_pred"] = pred_train_y
                    test_data[f"1D-QSAR-{name}_pred"] = pred_test_y
                    # 计算评估指标
                    mse = mean_squared_error(test_y, pred_test_y)
                    se = abs(test_y - pred_test_y)
                    results[f"{name}"] = {"MSE": mse, "error": se}
                    # st.write(f"[{name}]\tMSE:{mse:.4f}")

                # 创建一个DataFrame来存储结果
                results_df = pd.DataFrame.from_dict(results, orient='index').reset_index()
                results_df.rename(columns={'index': 'Regressor', 'MSE': 'Mean Squared Error'}, inplace=True)

                # 使用Altair创建柱状图
                chart = alt.Chart(results_df).mark_bar().encode(
                    x=alt.X('Regressor', sort=None),
                    y='Mean Squared Error',
                    tooltip=['Regressor', 'Mean Squared Error']
                ).properties(
                    width=600,
                    height=400,
                    title='Mean Squared Error (MSE) of Regressors'
                )
                # 在Streamlit中显示柱状图
                st.altair_chart(chart)
            #当选用分子指纹时
            elif selected_method == "Molecular Fingerprints":
                st.write("You selected: Molecular Fingerprints")
                # Continue with molecular fingerprints method
                st.write("Model training and evaluation will be performed using Molecular Fingerprints.")
                # Train models, etc.
                # For demonstration, we are using a simple placeholder for model training
                # In practice, you would train a model and evaluate its performance
                results = {
                    "Fingerprint Model": {"error": np.random.randn(100)}
                }


##########################################################################################################################################
                # Prepare data for plotting
                residuals_data = []
                for name, result in results.items():
                    if name.startswith("Fingerprint"):
                        model_residuals = pd.DataFrame({"Model": name, "Error": result["error"]})
                        residuals_data.append(model_residuals)

                residuals_df = pd.concat(residuals_data, ignore_index=True)
                residuals_df.sort_values(by="Error", ascending=True, inplace=True)
                model_order = residuals_df.groupby("Model")["Error"].median().sort_values(ascending=True).index

                # Plot using Streamlit
                st.subheader("Residuals Plot")
                residuals_chart = alt.Chart(residuals_df).mark_boxplot().encode(
                    y=alt.Y("Model:N", sort=model_order),
                    x="Error:Q"
                ).properties(
                    width=600,
                    height=400
                )

                st.altair_chart(residuals_chart)
    except Exception as e:
        st.error(f"Failed to read the file: {e}")
else:
    st.markdown("### Example Data Format")
    st.write("Please upload a CSV file containing the following columns:")
    st.dataframe(example_df)
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
