# Importing necessary libraries
import streamlit as st
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.colors import ListedColormap
import matplotlib as mpl
import pict_with_mask as pic 

# Is only loaded and displayed once for selection purposes
H_sparse = pd.read_csv(f'./CSVs/matrix_H_9_sparse.csv', index_col = 0)
# Variable assignment

# The columns chosen here from H_sparse (= odorants) were manually selected by me
# In an embedded program, the IDs of the odorants would come, for example, from user selection
# Column numbers (Note: the first column has column ID 0)
ind1 = 116 # Corresponds to the first column
ind2 = 136
# Variable assignment
ind3 = 19
# Variable assignment
# The row numbers (from H_sparse) that are nonzero for the respective column, i.e., the module rows that form the codeword (codeword IDs)
# Important to note: Python starts counting from 0, so if the module belonging to the second row is included in the codeword, it corresponds to ID = 1
cw_ids1 = [1,3] # Odorant in the first column has nonzero values in rows 2 and 4
cw_ids2 = [0,1,4]
# Variable assignment
cw_ids3 = [0,1,3,6]
# Variable assignment
# [x+2 for x in [2,6]]
k = 9
# Variable assignment

# lade odorant matrix (Originalmatrix mit nur positiven Einträgen und unique Odorants)
X = pd.read_csv('./CSVs/X_uniplus.csv', index_col = 0)
# Variable assignment

# bestimme die Namen der Odorants zugehörig zur gewählten Spalte aus H_sparse
odor_names = pd.read_csv('./CSVs/odorant_ids.csv', index_col = 0)
# Variable assignment

# die Reihenfolge der Odorants in Matrix X (index) entspricht der Reihenfolge der Spalten in H_sparse
name1 = odor_names.loc[X.index[ind1]]['odorant']
# Variable assignment
name2 = odor_names.loc[X.index[ind2]]['odorant']
# Variable assignment
name3 = odor_names.loc[X.index[ind3]]['odorant']
# Variable assignment

odor_names_list = odor_names.loc[X.index]['odorant']
# Variable assignment

# Create a mapping of index to item
index_to_odorant = {i: item for i, item in enumerate(odor_names_list)}
# Variable assignment
# Reverse mapping (to get index from the selection)
odorant_to_index = {v: k for k, v in index_to_odorant.items()}
# Variable assignment

# preload the NMF matrices for the alphabet sizes with k from  [2, 4, 6, 9] 
W9 = pd.read_csv(f'./CSVs/matrix_W_9.csv', index_col = 0)
H9 = pd.read_csv(f'./CSVs/matrix_H_9.csv', index_col = 0)

W6 = pd.read_csv(f'./CSVs/matrix_W_6.csv', index_col = 0)
H6 = pd.read_csv(f'./CSVs/matrix_H_6.csv', index_col = 0)

W4 = pd.read_csv(f'./CSVs/matrix_W_4.csv', index_col = 0)
H4 = pd.read_csv(f'./CSVs/matrix_H_4.csv', index_col = 0)

W2 = pd.read_csv(f'./CSVs/matrix_W_2.csv', index_col = 0)
H2 = pd.read_csv(f'./CSVs/matrix_H_2.csv', index_col = 0)

def build_dataset(ind, X, W, H):
# Function definition
    # Funktion baut ein numpy array (Dimension k+2 x 2121) auf als Datensatz zur Sammlung der Vektoren, die im Subplot später eine gemeinsame Zeile bilden
    # Reihenfolge der Vektoren:
    # erster Vektor: Originalbild der Substanz (odorant) aus Datenmatrix X
    # zweiter Vektor: nmf Approximation W*h für die Substanz (h sind die Koeffizienten aus Matrix H, die zur Substanz gehören)
    # dritter Vektor bis zu k+2 (letzter) Vektor: nmf Komponenten multipliziert mit passendem Koeffizient aus Matrix H (Komponenten ergeben als lineare Kombination die Approximation der Substanz)
    dataset = np.array([])
    # Variable assignment
    dataset = np.append(dataset, X.iloc[[ind]].values) # füge ersten Vektor hinzu, Dimension (1,2121)
    h = H.iloc[:, ind].values # Koeffizienten aus Matrix H für die Substanz aus Spalte 'ind'
    approx = W.values.dot(h) # Berechne Approximation, Dimension (2121,)
    dataset = np.append([dataset], [approx], axis = 0)  # durch axis = 0 werden die Vektoren als Zeilen angefügt
    nmf = W.values * h # diese Matrix enthält die nmf Komponenten als Spalten, hat also Dimension (2121,k)
    dataset = np.append(dataset, nmf.T, axis = 0) # transponiere nmf, damit die Komponenten als Zeilen an Datensatz angefügt werden können
    h_string = ['{:.2f}'.format(x) for x in h] # speichere die Koeffizienten als String mit 2 Nachkommastellen als Format
    return dataset, h_string

# Data choice initialization: Replace this with actual dataset
codebook_options = {"k=2": ["Component A", "Component B"],
                    "k=4": ["Component I", "Component II", "Component III", "Component IV"],
                    "k=6": ["I.1", "I.2", "III.a", "III.b", "IV.1", "IV.2"],
                    "k=9": ["I.1a", "I.1b", "III.a", "III.b", "IV.1a", "IV.1b", "V", "VI", "VII"]}

odorant_list = ["1-decanol", "quinoline", "2-undecanone"]  # Example compounds

def plot_decomposition(selected_k, selected_odorant1, selected_odorant2, selected_odorant3):
    # choose the current codebook
    allKs = {"k=2": 2, "k=4": 4, "k=6": 6, "k=9": 9}
    allWs = {2: W2, 4: W4, 6: W6, 9: W9}
    allHs = {2: H2, 4: H4, 6: H6, 9: H9}
    
    # set k as the codebook size. 
    if selected_k in allKs:
        k=allKs[selected_k]
    else:
        k=0
        raise ValueError("Invalid value of selected_k. Expected selected_k in [0, 1, 2, 3].") 
    W = allWs[k]
    H = allHs[k] 

    # generate the new figure according to the actual choices
    fig2, axes = plt.subplots(1, k, figsize=(15, 5))

    # generate linear decomposition of all three compounds x1, x2, x3 with W and H from the current alphabet
    dataset1, h1 = build_dataset(selected_odorant1, X, W, H)
    dataset2, h2 = build_dataset(selected_odorant2, X, W, H)
    dataset3, h3 = build_dataset(selected_odorant3, X, W, H)
    # store the encoding vectors W also for displaying them in the first line of the figure
    datasetW = W.values.T
    
    # black boxes position to mark the non-zero words (for k=9 only!) 
    if k==9 :
        cw_ids1 = H_sparse.iloc[:, selected_odorant1].to_numpy().nonzero()[0].tolist()
        cw_ids2 = H_sparse.iloc[:, selected_odorant2].to_numpy().nonzero()[0].tolist()
        cw_ids3 = H_sparse.iloc[:, selected_odorant3].to_numpy().nonzero()[0].tolist()
    else :
        cw_ids1 = [] 
        cw_ids2 = []
        cw_ids3 = []

    # for this specific figure, the first two columns are used differently, thus there is an offset of 2 for the boxes
    cw_ids1 = [x+2 for x in cw_ids1]
    cw_ids2 = [x+2 for x in cw_ids2]
    cw_ids3 = [x+2 for x in cw_ids3]

    sub_titles = np.array(['2DG-uptake','decoded \n approximation','','','','','','','','',''])
    # Subplot mit 4 Zeilen und k+2 Spalten, 'constrained' layout, um automatische Anordnung zu erreichen, figsize etwas Ausprobieren
    fig, axs = plt.subplots(4,k+2, layout = 'constrained', figsize=(k*2,12))
    fig.suptitle(f'NMF encoding with k = {k}') # superior title, also Gesamt-Titel des Plots

    ##############################
    ## plot first line in figure
    ##############################
    for i in range(k+2):
        if i in range(2): # in den ersten 2 Spalten gibt es keine Plots
            axs[0,i].set_yticks([])
            axs[0,i].set_xticks([])
            axs[0,i].spines['top'].set_visible(False)
            axs[0,i].spines['right'].set_visible(False)
            axs[0,i].spines['bottom'].set_visible(False)
            axs[0,i].spines['left'].set_visible(False)
        else:
            img, cmp = pic.pict(datasetW[i-2]) # plotte die NMF Komponenten
            axs[0,i].imshow(img, cmap = cmp)
            axs[0,i].set_yticks([])
            axs[0,i].set_xticks([])
            axs[0,i].spines['top'].set_visible(False)
            axs[0,i].spines['right'].set_visible(False)
            axs[0,i].spines['bottom'].set_visible(False)
            axs[0,i].spines['left'].set_visible(False)
    axs[0,2].set_ylabel('nmf components') # setze im 3. Plot dieser Zeile ein y-label ein
    ##############################
    ## plot second line in figure
    ##############################
    # die Darstellung der Linearkombination über alle NMF Komponenten benötigt eine gemeinsame Colormap für alle Bilder derselben Zeile, 
    # daher wird hier eine gemeinsame Norm über alle Vektoren des Datensatzes berechnet
    norm1 = mpl.colors.Normalize(vmin=np.min(dataset1), vmax=np.max(dataset1))
    xlabel1 = ['',''] + h1 # Koeffizienten unter den NMF Komponenten als x-label, die ersten zwei Plots der Zeile sind uptake-image und decoded approximation
    for i in range(k+2):
        img, cmp = pic.pict(dataset1[i])
        axs[1,i].imshow(img, norm = norm1, cmap = cmp)
        axs[1,i].set_title(sub_titles[i]) # setze hier Untertitel für die Plots ein, um Spalte 2DG-uptake und decoded approximation zu benennen
        axs[1,i].set_xlabel(xlabel1[i])
        axs[1,i].set_yticks([])
        axs[1,i].set_xticks([])
        # Plots, die NMF Komponenten zeigen, die zum Codeword der Substanz gehören, behalten ihren Achsenrahmen (zur Hervorhebung)
        # Plots, die nicht zum Codeword gehören (also deren ids nicht in cw_ids sind), wird der Achsenrahmen entfernt
        if i not in cw_ids1:
            axs[1,i].spines['top'].set_visible(False)
            axs[1,i].spines['right'].set_visible(False)
            axs[1,i].spines['bottom'].set_visible(False)
            axs[1,i].spines['left'].set_visible(False)
    axs[1,0].set_ylabel(index_to_odorant[selected_odorant1]) # Name der Substanz als y-label des ersten Plots in dieser Zeile
    # gemeinsame Colormap wird am Ende der Zeile, also zum letzten Plot zugehörig angezeigt
    fig.colorbar(plt.cm.ScalarMappable(norm=norm1, cmap=cmp), ax=axs[1, k+1], pad = 0.2, orientation='vertical', fraction = 0.1)
    ##############################
    ## plot third line in figure
    ##############################
    norm2 = mpl.colors.Normalize(vmin=np.min(dataset2), vmax=np.max(dataset2))
    xlabel2 = ['',''] + h2
    for j in range(k+2):
        img, cmp = pic.pict(dataset2[j])
        axs[2,j].imshow(img, norm = norm2, cmap = cmp)
        axs[2,j].set_xlabel(xlabel2[j])
        axs[2,j].set_yticks([])
        axs[2,j].set_xticks([])
        if j not in cw_ids2: 
            axs[2,j].spines['top'].set_visible(False)
            axs[2,j].spines['right'].set_visible(False)
            axs[2,j].spines['bottom'].set_visible(False)
            axs[2,j].spines['left'].set_visible(False)
    axs[2,0].set_ylabel(index_to_odorant[selected_odorant2])
    fig2.colorbar(plt.cm.ScalarMappable(norm=norm2, cmap=cmp), ax=axs[2, k+1], pad = 0.2, orientation='vertical', fraction = 0.1)
    ##############################
    ## plot fourth line in figure
    ##############################
    norm3 = mpl.colors.Normalize(vmin=np.min(dataset3), vmax=np.max(dataset3))
    xlabel3 = ['',''] + h3
    for l in range(k+2):
        img, cmp = pic.pict(dataset3[l])
        axs[3,l].imshow(img, norm = norm3, cmap = cmp)
        axs[3,l].set_xlabel(xlabel3[l])
        axs[3,l].set_yticks([])
        axs[3,l].set_xticks([])
        if l not in cw_ids3: # remove axes frames for subplots that are no nmf features from the codeword for this odorant
            axs[3,l].spines['top'].set_visible(False)
            axs[3,l].spines['right'].set_visible(False)
            axs[3,l].spines['bottom'].set_visible(False)
            axs[3,l].spines['left'].set_visible(False)
    axs[3,0].set_ylabel(index_to_odorant[selected_odorant3])
    fig2.colorbar(plt.cm.ScalarMappable(norm=norm3, cmap=cmp), ax=axs[3, k+1], pad = 0.2, orientation='vertical', fraction = 0.1)

    st.pyplot(fig)

# Streamlit UI
st.logo("inb_logo.png", size="large", link="https://www.inb.uni-luebeck.de/en")
st.title("Exploring the Modular Decomposition of Olfactory Bulbar Activity")
st.sidebar.header("Select Options")

# define sidebar selectors
selected_k = st.sidebar.selectbox("Select Codebook Size:", list(codebook_options.keys()), index=3)
selected_odorant1 = st.sidebar.selectbox("Select Odorant Compound 1:", odor_names_list, key="chem1", index=ind1)
selected_odorant2 = st.sidebar.selectbox("Select Odorant Compound 2:", odor_names_list, key="chem2", index=ind2)
selected_odorant3 = st.sidebar.selectbox("Select Odorant Compound 3:", odor_names_list, key="chem3", index=ind3)

st.markdown(
    """This interactive tool allows users to explore how the spatial activity patterns in the 
    olfactory bulb are decomposed into different modular representations. The study behind 
    this tool utilized **Non-Negative Matrix Factorization (NMF)** to identify a compact coding 
    alphabet for odor perception in rodents. 
    Select different codebook sizes to observe how increasing the number of modules refines 
    the spatial encoding of odors.""")

st.write(f"Currently, the alphabet with **{k}** components is used and displayed are the three compounds **{selected_odorant1}** (ID: {odorant_to_index[selected_odorant1]}), **{selected_odorant2}** (ID: {odorant_to_index[selected_odorant2]}), and **{selected_odorant3}** (ID: {odorant_to_index[selected_odorant3]}).")
# textual sanity check

# generate the main figure and display it here
plot_decomposition(selected_k, odorant_to_index[selected_odorant1], odorant_to_index[selected_odorant2], odorant_to_index[selected_odorant3])

# Add some project information
# Explanatory Section on NMF-Codebook Generation
st.markdown("How the NMF Codebook is Generated: The **Non-Negative Matrix Factorization (NMF)** technique decomposes complex olfactory uptake images into a set of fundamental components known as the **codebook**.")
st.markdown("Here's a simplified breakdown:")
st.markdown("1. **Data Matrix Construction:** The 2DG uptake images of the olfactory bulb for different odorants are represented as a large matrix, where each row corresponds to spatial glomerular activity, and each column represents an odorant.")
st.markdown("2. **Matrix Decomposition:** NMF factorizes this matrix into two smaller matrices: **W (Basis Matrix):** Contains the spatial modules (or components) representing distinct activity patterns. **H (Coefficient Matrix):** Indicates how much each module contributes to reconstructing the original uptake image.")
st.markdown("3. **Choosing Codebook Size (k):** The parameter **k** defines the number of components (modules) in the codebook. Optimal values of k (e.g., 2, 4, 6, 9) are identified based on reconstruction accuracy and biological interpretability.")
st.markdown("4. **Reconstruction:** Each odorant's spatial pattern can be reconstructed as a linear combination of these components, forming its unique **odor code**.")
st.markdown("This modular representation helps understand how complex olfactory information is encoded and processed in the brain.")

