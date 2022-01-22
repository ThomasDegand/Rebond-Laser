####################
##  Importations  ##
####################
import numpy as np
import matplotlib.pyplot as plt



##################
##  Constantes  ##
##################
g = 9.81
n0 = 1
n1 = 1.33

##############
##  Outils  ##
##############
def r2d(a):
    return np.rad2deg(a)

def d2r(a):
    return np.deg2rad(a)

def jet(v0, a, z, D):
    u0 = z + D/2
    u1 = np.tan(a)
    u2 = -g/(2*(v0**2)*((np.cos(a))**2))
    jetUp = (u0, u1, u2)
    d0 = u0 - D
    d1 = u1
    d2 = u2
    jetDown = (d0, d1, d2)
    return (jetUp, jetDown)
    
##########

def derive(P):
    n = len(P)
    dP = [0]*(n-1)
    for i in range(1,n):
        dP[i-1] = i*P[i]
    return tuple(dP)
    
def rayon(M, a):
    (x, y) = M
    r1 = np.tan(a) 
    r0 = y - r1*x
    return (r0, r1)
    
def f(P, x):
    n = len(P)
    return sum([P[i]*(x**i) for i in range(n)])

def intersection(P, d):
    a = P[2]
    b = P[1] - d[1]
    c = P[0] - d[0]
    w = (((b**2)-4*a*c)**0.5)/(2*a)
    r = -b/(2*a)
    x1 = r-w
    y1 = d[1]*x1 + d[0]
    x2 = r+w
    y2 = d[1]*x2 + d[0]
    return ((x1,y1), (x2,y2))

def angle(u):
    return np.arctan(u[1]/u[0])

def sens(u):
    if u[0] < 0:
        x = -1
    else:
        x = 1
    if u[1] < 0:
        y = -1
    else:
        y = 1
    return (x, y)

def vecteur(M1, M2):
    x = M2[0] - M1[0]
    y = M2[1] - M1[1]
    return (x, y)

def norme(u):
    return (u[0]**2 + u[1]**2)**0.5

##############
##  Projet  ##
##############

def trajet(M0, u, F):
    a = angle(u)
    d = rayon(M0, a)
    (Mu1, Mu2) = intersection(F[0], d)
    (Md1, Md2) = intersection(F[1], d)
    listM = [(Mu1,0), (Mu2,0), (Md1,1), (Md2,1)]
    #print("M0:",M0,"; listM:",listM)
    I = []
    i = 0
    j = 0
    N = 999999
    for M in listM:
        v = vecteur(M0, M[0])
        #print(u,sens(v),sens(u))
        if sens(v) != sens(u):
            I.append(i)
            j = 1
        if norme(v) < 0.000001 and j != 1:
            I.append(i)
        j = 0
        i += 1
    for i in I:
        listM[i] = 0
    modif = 0
    for i in range(4):
        if listM[i-modif] == 0:
            listM.pop(i-modif)
            modif += 1
    for M in listM:
        v = vecteur(M0, M[0])
        n = norme(v)
        if n < N:
            N = n
            M1 = M
    if not('M1' in locals()):
        u = ((-1)*u[0], (-1)*u[1])
        a = angle(u)
        d = rayon(M0, a)
        (Mu1, Mu2) = intersection(F[0], d)
        (Md1, Md2) = intersection(F[1], d)
        listM = [(Mu1,0), (Mu2,0), (Md1,1), (Md2,1)]
        #print("M0:",M0,"; listM:",listM)
        I = []
        i = 0
        j = 0
        N = 999999
        for M in listM:
            v = vecteur(M0, M[0])
            #print(u,sens(v),sens(u))
            if sens(v) != sens(u):
                I.append(i)
                j = 1
            if norme(v) < 0.000001 and j != 1:
                I.append(i)
            j = 0
            i += 1
        for i in I:
            listM[i] = 0
        modif = 0
        for i in range(4):
            if listM[i-modif] == 0:
                listM.pop(i-modif)
                modif += 1
        for M in listM:
            v = vecteur(M0, M[0])
            n = norme(v)
            if n < N:
                N = n
                M1 = M
    return (M1[0], vecteur(M0, M1[0]), F[M1[1]])

def rebond(M, u, P):
    dP = derive(P)
    t = (1, f(dP, M[0]))
    #tX = np.linspace(0, 5, 300)
    #tY = f(t, tX)
    #axes.plot(tX, tY)
    #print("u:",u)
    #print("t:",t)
    aU = angle(u)
    #print("aU:",r2d(aU))
    aT = angle(t)
    #print("aT:",r2d(aT))
    aRel = np.pi/2 - abs(aU - aT)
    if aRel <= aLim:
        print("OZBIUFBZPFBVUZJFOZUBOFBZUFBZOFBZUFB")
    print("aRel:",r2d(aRel))
    aR = aT - (aU - aT)
    print("aR:",r2d(aR))
    r = (np.cos(aR), np.sin(aR))
    return (M, r)

def rebondOut(M, u, P):
    dP = derive(P)
    t = (1, f(dP, M[0]))
    aU = angle(u)
    aT = angle(t)
    aN = np.pi/2 - abs(aT - aU)
    aR = np.arcsin((n0/n1)*np.sin(aN))
    if aU < aT:
        aRes = (aN - aU) - aR
    else:
        aRes = np.pi/2 - aT + aR
    r = (np.cos(aRes), np.sin(aRes))
    return (M, r)

def parcours(M0, a, v0, z, D):
    ## Le jet et initialisation
    F = (jetUp, jetDown)
    M00 = M0
    u = (1, np.tan(a))
    points = []
    
    ## Position
    dO = (M0[1], 0)
    E = (intersection(F[0], dO), intersection(F[1], dO))
    if M0[0] < E[0][1][0] or M0[0] > E[1][1][0]:
        etape1 = trajet(M0, u, F)
        points.append(etape1[0])
        M = etape1[0]
        u = etape1[1]
        P = etape1[2]
        etape2 = rebondOut(M, u, P)
        M0 = etape2[0]
        u = etape2[1]
        
    ## Début des rebonds
    #F = jet(v0, a, z, D)
    for i in range(57):
        ## Trajet
        etape1 = trajet(M0, u, F)
        points.append(etape1[0])
        #print("trajet:",etape1)
        if points[-1][1] < 0:
            break
        
        ## Rebond
        M = etape1[0]
        u = etape1[1]
        P = etape1[2]
        etape2 = rebond(M, u, P)
        #print("rebond:",etape2)
        #Trajet
        M0 = etape2[0]
        u = etape2[1]
        
    ## Affichage
    X = [M00[0]]
    Y = [M00[1]]
    for point in points:
        X.append(point[0])
        Y.append(point[1])
    #print(X,Y)
    axes.plot(X, Y, color="r")
    axes.plot(np.linspace(-1,1600,160000), f(F[0], np.linspace(-1,1600,160000)), color="b")
    axes.plot(np.linspace(-1,1600,160000), f(F[1], np.linspace(-1,1600,160000)), color="b")
    axes.plot(np.linspace(-1,1600,160000), f(dO, np.linspace(-1,1600,160000)), color="g")
    plt.show()
    return points
    
##########

#Initialisation affichage
figure = plt.figure()
axes = figure.add_subplot(111)

aLim = np.arcsin(n0/n1)
print(r2d(aLim))

#Polynomes issus de la photo
jetUp = (1.9999999999999696, 1.496679292383839, -0.0009806366362053933)
jetDown = (-72.3726962681906, 1.548029126479227, -0.001011504022308292)
parcours((0.1, 0.05), d2r(70), 6.5, 0, 0.01)


