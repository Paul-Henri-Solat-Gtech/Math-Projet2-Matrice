import matplotlib.pyplot as plt
import math

def prodmat(A, B):
    n = len(A)            # nombre de lignes de A
    p = len(A[0])         # nombre de colonnes de A = nombre de lignes de B
    m = len(B[0])         # nombre de colonnes de B

    C = []
    for i in range(n):
        D = []
        for j in range(m):
            somme = 0
            for k in range(p):
                somme += A[i][k] * B[k][j]
            D.append(somme)
        C.append(D)
    return C

    # Test Q1
# T1 = [[1, 0],
#       [-1, 2],
#       [0, 3]]

# T2 = [[-1, 0, 1],
#       [2, 1, -1]]

# print(prodmat(T1, T2))

def pop(A, i, j):
    B = []
    n = len(A)
    for k in range(n):
        if k != i:
            C = []
            for l in range(n):
                if l != j:
                    C.append(A[k][l])
            B.append(C)
    return B 

def det(A):
    n = len(A)            
    p = len(A[0])         
    if n == 1:
        return A[0][0]
    else:
        somme = 0
        j = 0
        for i in range(n):
            somme += (-1) ** (i + j) * A[i][j] * det(pop(A, i, j))
        return somme

    # Test Q2
# T1 = [[1, 0, -2],
#       [2, 1, 1],
#       [-1, 0, 1]]

# print(det(T1))

def com(A):
    n = len(A)            
    p = len(A[0])         
    B = []
    for i in range(n):
        C = []
        for j in range(p):
            C.append(((-1) ** (i + j)) * det(pop(A, i, j)))
        B.append(C)
    return B

    
    # Test Q3
# T1 = [[1, 0, -2],
#       [2, 1, 1],
#       [-1, 0, 1]]

# print(com(T1))

def prodmatscal(A,k):
    n = len(A)            
    p = len(A[0])         
    B = []
    for i in range(n):
        C = []
        for j in range(p):
            C.append(k * A[i][j])
        B.append(C)
    return B

def transpose(A):
    n = len(A)            
    p = len(A[0])         
    B = []
    for j in range(p):
        C = []
        for i in range(n):
            C.append(A[i][j])
        B.append(C)
    return B

def inverse(A): #nécéssite prodmatscal et transpose
    return prodmatscal(transpose(com(A)), 1/det(A))
   

# Test Q
# T1 = [[1, 0, -2],
#       [2, 1, 1],
#       [-1, 0, 1]]

# print(inverse(T1))

def factoriel(n):
    if n == 0:
        return 1
    else:
        return n * factoriel(n - 1)

def cosinus(x):
    return sum(((-1)**k * x**(2*k)) / factoriel(2*k) for k in range(30))

def sinus(x):
    return sum(((-1)**k * x**(2*k+1)) / factoriel(2*k+1) for k in range(30))

# Question 2 : Tracé d’une ligne
def ligne(n, xmin, xmax):
    W = []
    dx = (xmax - xmin) / (n - 1)
    for i in range(n):
        x = xmin + i * dx
        W.append([x, 0, 0])
    return transpose(W)

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
    return transpose(W[:n])  # au cas où dépasse

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
    return transpose(W[:n])

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
    return transpose(W[:n])


# Question 6c : Cercle plein
def cercle_plein(n, R):
    W = []
    k = int(n**0.5)
    for i in range(k):
        for j in range(k):
            r = (i / (k - 1)) * R
            t = (2 * math.pi * j) / (k - 1)
            x = r * cosinus(t)
            y = r * sinus(t)
            W.append([x, y, 0])
    return transpose(W[:n])

# Question 7 : Cylindre plein
def cylindre_plein(n_r, n_theta, n_z, R, h):
    W = []
    for i in range(n_r):
        r = (i / (n_r - 1)) * R if n_r > 1 else R
        for j in range(n_theta):
            t = (2 * math.pi * j) / n_theta
            for l in range(n_z):
                z = -h/2 + l * (h / (n_z - 1))
                x = r * cosinus(t)
                y = r * sinus(t)
                W.append([x, y, z])
    return transpose(W)

def cylindre_plein(n_r, n_theta, n_z, R, h):
  from math import pi, cos, sin
  pts = []

  for i in range(n_r):
      r = R * i / (n_r - 1) if n_r > 1 else R
      for j in range(n_theta):
          theta = 2 * pi * j / n_theta
          for k in range(n_z):
              z = -h / 2 + (h * k / (n_z - 1)) if n_z > 1 else 0
              x = r * cos(theta)
              y = r * sin(theta)
              pts.append([x, y, z])

  return transpose(pts)

def plot3D(X, Y, Z, color='ocean'):
    fig = plt.figure()
    ax = plt.axes(projection='3d')
    ax.scatter3D(X, Y, Z, c=Z, cmap=color)
    ax.set_xlim(-5, 5)
    ax.set_ylim(-5, 5)
    ax.set_zlim(-5, 5)
    plt.show()
