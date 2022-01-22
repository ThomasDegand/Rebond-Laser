##################
## Importations ##
##################

import numpy as np
import matplotlib.pyplot as plt
from random import randint
from math import sqrt
from PIL import Image


###############
## Fonctions ##
###############

def m(X):
    if len(X) != 0:
        return sum(X)/len(X)

def power(n):
    def pown(x):
        return x**n
    return pown

def precision(f,X,Y): #corrélation
    fX = v(f)(X)
    EX = m(fX)
    EY = m(Y)
    XE = fX-EX
    YE = Y-EY
    cov = m(XE*YE)
    SX = sqrt(m(XE**2))
    SY = sqrt(m(YE**2))
    return (cov/(SX*SY))**2
    
#################
## Acquisition ##
#################

jet = np.array(Image.open("jet.png"))
def pointer(image, r0, g0, b0):
    p, q, m = np.shape(image)
    X = []
    Y = []
    for j in range(q):
        I = []
        for i in range(p):
            r, g, b = image[i][j][0], image[i][j][1], image[i][j][2]
            if r == r0 and g == g0 and b == b0:
                I.append(i)
        if len(I) != 0:
            X.append(j)
            Y.append(p-sum(I)/len(I))
    return np.array(X, dtype=float), np.array(Y, dtype=float)

Xbot, Ybot = pointer(jet, 0, 0, 255)
Xtop, Ytop = pointer(jet, 255, 0, 0)


##################
## Modélisation ##
##################

def v(f): #Polynome
    n = len(f)
    def vf(x):
        return sum([f[i]*(x**i) for i in range(n)])
    return vf

def model(n,X,Y): #Regression polynomiale recursive
    X = np.array(X, dtype=float)
    Y = np.array(Y, dtype=float)
    if n==1:
        pente = (m(X*Y)-m(X)*m(Y))/(m(X*X)-m(X)*m(X))
        origine = m(Y) - pente*m(X)
        f = [origine,pente]
        return f
    L = len(X)
    Xbis = np.array([(X[i+1]+X[i])/2 for i in range(L-1)])
    Ybis = np.array([(Y[i+1]-Y[i])/(X[i+1]-X[i]) for i in range(L-1)])
    fo = model(n-1,Xbis,Ybis)
    f = [0] + [fo[i]/(i+1) for i in range(n)]
    #print("f;", f, "f[0]:", Y-v(f)(X))
    #print("Y:", Y, "X:", X)
    f[0] = m(Y-v(f)(X))
    return f

def all(n,X,Y): #Affichage de l'approche polynomiale
    p = 200
    Xl = np.linspace(X[0],X[-1],int(p*(X[-1]-X[0])))
    f = model(n,X,Y)
    p = precision(f,X,Y)
    print("f:", f, "precision:", p)
    plt.scatter(X,Y)
    plt.plot(Xl,v(f)(Xl))
    plt.show()
    return p, f




fbot=all(2,Xbot,Ybot)[1]
ftop=all(2,Xtop,Ytop)[1]
Xl = np.linspace(max(Xbot[0],Xtop[0]), max(Xbot[-1],Xtop[-1]), int(200*(Xbot[-1]-Xbot[0])))
plt.plot(Xl,v(fbot)(Xl))
plt.plot(Xl,v(ftop)(Xl))
plt.show()
