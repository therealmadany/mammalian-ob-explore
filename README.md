# How to Use This Streamlit App

## Getting Started with Streamlit
Streamlit is a Python framework that allows you to build interactive web applications for data science and machine learning projects. To run a Streamlit app, follow these steps:

1. **Install Streamlit** (if not already installed):
   ```bash
   pip install streamlit
   ```
2. **Navigate to the directory** containing the Streamlit app (`app.py`):
   ```bash
   cd /path/to/your/project
   ```
3. **Run the application**:
   ```bash
   streamlit run app.py
   ```
4. A local web server will start, and a URL (typically `http://localhost:8501/`) will be provided. Open this URL in your browser to interact with the app.

---

## Using This Application

This application allows users to manipulate parameters in the left-hand menu to explore the decomposition of three different compounds from a set of 143 compounds. The decomposition is analyzed along four different decoding alphabets with varying alphabet sizes of **k = 2, 4, 6, and 9**.

### Features:
- **Interactive Sidebar Controls**: Adjust parameters to modify the decomposition.
- **Real-time Analysis**: Observe how different settings affect the breakdown of selected compounds.
- **Visualization & Interpretation**: The app provides an intuitive interface to understand how the decoding process varies with different alphabet sizes.

### Steps to Use:
1. **Select the compounds**: Choose up to three compounds from the dataset of 143 compounds.
2. **Modify parameters**: Adjust settings in the sidebar to explore different decompositions.
3. **Analyze the results**: View the changes in decomposition across four different decoding alphabets.
4. **Experiment with alphabet sizes**: Change `k` values (2, 4, 6, or 9) to observe variations in the output.

### Scientific Background
This application illustrates the decomposition of bulbar activity images measured during exposure to single odorant chemicals (odorants) into a linear combination of fundamental activity patterns. This decomposition is achieved using **non-negative matrix factorization (NMF)**, a mathematical approach that ensures all components and their combinations remain positive. This method allows us to reconstruct the original data using optimized, non-negative building blocks, akin to mathematically optimized **Lego bricks**, which collectively resemble the original measured activity patterns. The resulting decomposition forms a generic alphabet representing these neural activities in a representative manner. 

This tool is useful for understanding how different alphabet structures influence the decomposition of chemical compounds and can be a valuable asset for researchers analyzing molecular structures through computational methods.

For any issues or further information, refer to the documentation or reach out to the developer team.
