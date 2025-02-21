# Imports
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.colors import ListedColormap
import matplotlib as mpl

# file path to X_mask need to be changed if needed
X_mask = pd.read_csv('./CSVs/X_mask.csv', index_col = 0).to_numpy()

def pict_show(vec, name = 'test', weight = 1, mask = X_mask):
    # Funktion zum Plotten des Vectors vec (zB Uptake image oder Module) der Größe 2121x1 gegeben die Maske X_mask zum plotten.
    # optionale Parameter:
    # name: Name, unter dem das Bild abgespeichert wird, default: 'test'
    # weight: Gewichtung zwischen [0,1], Wert gibt an, bis zu welchem Wert die Colormap die Werte anzeigt, default: 1
    # mask: Maske X_mask

    picture = np.zeros(mask.shape)
    np.copyto(picture, mask)

    jet = mpl.colormaps['jet'].resampled(256) # Colormap jet
    if (weight > 1 or weight < 0):
        raise Exception('weight has to be between 0 and 1')
    
    num_colors = round(weight*256)
    myjet = ListedColormap(jet(np.linspace(0, weight, num_colors))) # passe Colormap an
    myjet.set_bad(color='w') # setze maskierte Werte für diese Colormap auf weiß (= Hintergrund)

    if ((picture == 1).sum()!=2121):
        raise Exception('something is wrong with the mask')
    
    np.place(picture, picture == 0, -0.1) # definiere Hintergrund durch eindeutigen Wert (0 könnte durch Modulwerte besetzt werden)
    np.place(picture.T, picture.T == 1, vec) # Picture Array muss transformiert werden, da Python erst zeilenweise vorgeht, die Daten in Vec waren durch Matlab aber erst Spaltenweise gefüllt
    picture_masked = np.ma.masked_where(picture == -0.1, picture) # maskiere das Bild durch den Hintergrund
    plt.imshow(picture_masked, cmap=myjet)
    plt.axis('off')
    plt.savefig(f'Images/{name}.png', bbox_inches='tight')
    return


def pict(vec, weight = 1, mask = X_mask):
    # Funktion bestimmt Bildmatrix (80x44) des Vectors vec (zB Uptake image oder Module) mit Maske X_mask und gibt diese als return-Wert wieder.
    # optionale Parameter:
    # weight: Gewichtung zwischen [0,1], Wert gibt an, bis zu welchem Wert die Colormap die Werte anzeigt, default: 1
    # mask: Maske X_mask

    picture = np.zeros(mask.shape)
    np.copyto(picture, mask)

    jet = mpl.colormaps['jet'].resampled(256) # Colormap jet
    if (weight > 1 or weight < 0):
        raise Exception('weight has to be between 0 and 1')
    
    num_colors = round(weight*256)
    myjet = ListedColormap(jet(np.linspace(0, weight, num_colors))) # passe Colormap an
    myjet.set_bad(color='w') # setze maskierte Werte für diese Colormap auf weiß (= Hintergrund)

    if ((picture == 1).sum()!=2121):
        raise Exception('something is wrong with the mask')
    
    np.place(picture, picture == 0, -0.1) # definiere Hintergrund durch eindeutigen Wert (0 könnte durch Modulwerte besetzt werden)
    np.place(picture.T, picture.T == 1, vec) # Picture Array muss transformiert werden, da Python erst zeilenweise vorgeht, die Daten in Vec waren durch Matlab aber erst Spaltenweise gefüllt
    picture_masked = np.ma.masked_where(picture == -0.1, picture) # maskiere das Bild durch den Hintergrund
    return picture_masked, myjet

