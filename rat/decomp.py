
from sklearn import decomposition
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.markers import MarkerStyle

def simple_pca(mat, label=False, **sets):

    pca = decomposition.PCA(n_components=2)
    pca.fit(mat)
    X = pd.DataFrame(pca.transform(mat))
    X.index = mat.index
    plt.scatter(X[0], X[1], marker='o', c='k', s=50)
    for i, (n, sel) in enumerate(sets.items()):
        if i == 0:
            plt.scatter(X.loc[sel][0], X.loc[sel][1], s=50, c='coral', label=n)
        elif i == 1:
            plt.scatter(X.loc[sel][0], X.loc[sel][1],
                        edgecolors='black',
                        marker=(0, 3, 0), facecolors='none',
                        s=100, label=n)
    if len(sets) > 0:
        plt.legend(loc='best')

    if label:
        for n, r in X.iterrows():
            if isinstance(n, list) or isinstance(n, tuple):
                n = "_".join(n)
            plt.text(r[0], r[1], " " + str(n))
        
