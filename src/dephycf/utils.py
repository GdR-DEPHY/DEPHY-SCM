import numpy as np
from scipy.stats import norm

#### fonctions donné par Anne
def gauss_weight(x,i,sd):
    '''Renvoie un tableau de coefficients de la taille de x dont la somme est 1
    Fenetre gaussienne centrée sur i écart type sd
    '''
    w=norm.pdf(x,loc=i,scale=sd)
    return w/w.mean()

def kernel_mean_gauss(x,y,sd,yout=None):
    '''moyenne mobile de x pondérée par une fenetre gaussienne sur y avec une
    largeur de fenetre sd
    sortie pour les valeurs de y=yout si présent sinon y
    '''
    if yout is None: yout=y
    return [np.average(x,weights=gauss_weight(y,h,sd)) for h in yout]
