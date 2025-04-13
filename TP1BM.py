import matplotlib.pyplot as plt
import TP2BM as TP2
from math import pi

# Question 2 : Tracé d’une ligne
def ligne(n, xmin, xmax):
    W = []
    dx = (xmax - xmin) / (n - 1)
    for i in range(n):
        x = xmin + i * dx
        W.append([x, 0, 0])
    return TP2.transpose(W)

# Question 3 : Carré vide
def carrevide(n, a):
    W = []
    c = a / 2
    p = int(n // 4)
    for i in range(p):
        x = -c + i * (a / (p - 1))
        W.append([x, -c, 0])
        W.append([x, c, 0])
    for i in range(p):
        y = -c + i * (a / (p - 1))
        W.append([-c, y, 0])
        W.append([c, y, 0])
    return TP2.transpose(W[:n])  # au cas où dépasse

# Question 4 : Carré plein
def carreplein(n, a):
    W = []
    c = a / 2
    k = int(n**0.5)
    for i in range(k):
        for j in range(k):
            x = -c + i * (a / (k - 1))
            y = -c + j * (a / (k - 1))
            W.append([x, y, 0])
    return TP2.transpose(W[:n])

# Question 5 : Pavé plein
def pave_plein(n, a, b, c):
    W = []
    nx = int(n ** (1/3))
    for i in range(nx):
        for j in range(nx):
            for k in range(nx):
                x = -a/2 + i * (a / (nx - 1))
                y = -b/2 + j * (b / (nx - 1))
                z = -c/2 + k * (c / (nx - 1))
                W.append([x, y, z])
    return TP2.transpose(W[:n])


# Question 6c : Cercle plein
def cercle_plein(n, R):
    W = []
    k = int(n**0.5)
    for i in range(k):
        for j in range(k):
            r = (i / (k - 1)) * R
            t = (2 * pi * j) / (k - 1)
            x = r * TP2.cosinus(t)
            y = r * TP2.sinus(t)
            W.append([x, y, 0])
    return TP2.transpose(W[:n])

# Question 7 : Cylindre plein
def cylindre_plein(n, R, h):
    W = []
    k = int(n ** (1/3))
    for i in range(k):
        for j in range(k):
            for l in range(k):
                r = (i / (k - 1)) * R
                t = (2 * pi * j) / (k - 1)
                z = -h/2 + l * (h / (k - 1))
                x = r * TP2.cosinus(t)
                y = r * TP2.sinus(t)
                W.append([x, y, z])
    return TP2.transpose(W[:n])

# ============ EXEMPLES DE TRACÉ ============= #

def plot3D(X, Y, Z, color='ocean'):
    fig = plt.figure()
    ax = plt.axes(projection='3d')
    ax.scatter3D(X, Y, Z, c=Z, cmap=color)
    plt.show()

# Test ligne
X, Y, Z = ligne(20, -5, 5)
plot3D(X, Y, Z)

# Test carré vide
X, Y, Z = carrevide(100, 4)
plot3D(X, Y, Z)

# Test carré plein
X, Y, Z = carreplein(400, 4)
plot3D(X, Y, Z)

# Test pavé plein
X, Y, Z = pave_plein(1000, 4, 3, 2)
plot3D(X, Y, Z)

# Test cercle plein
X, Y, Z = cercle_plein(1000, 4)
plot3D(X, Y, Z)

# Test cylindre plein
X, Y, Z = cylindre_plein(2000, 2, 5)
plot3D(X, Y, Z)

