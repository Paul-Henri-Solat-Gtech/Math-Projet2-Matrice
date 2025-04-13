import matplotlib.pyplot as plt
import numpy as np
from math import pi

def factorial(n):
  fact = 1
  for i in range(n):
    fact *= n - i
  return (fact)

def cos(x):
  signe = -1
  fcos = 0
  for n in range(100):
    if n % 2 != 0:  # Si n n'est pas un multiple de 2
      continue
    signe *= -1
    fcos += (signe * ((x**n) / factorial(n)))

  return fcos

def sin(x):
  return cos(pi / 2 - x)

def transpose(M):
  Mt = [[0 for _ in range(len(M))] for _ in range(len(M[0]))]

  for i in range(len(M)):  # Parcours les lignes de M
    for j in range(len(M[0])):  # Parcours les colonnes de M
      Mt[j][i] = M[i][j]  # Échange les éléments
  return Mt

def pop(A, i, j):
  n = len(A)
  b = []
  for ligne in range(n):
    if ligne != i:
      c = []
      for colonne in range(n):
        if colonne != j:
          c.append(A[ligne][colonne])
      b.append(c)
  return b

def det(A):
  n = len(A)
  d = 0
  j = 0
  i = 0
  #det = (-1)**i + j
  if len(A) == 1:
    return (A[0][0])
  else:
    for i in range(n):
      d += (-1)**(i + j) * A[i][j] * det(pop(A, i, j))

    return (d)

def com(A):
  j = 0
  n = len(A)
  b = []
  for i in range(n):
    c = []
    for j in range(n):
      c.append(((-1)**(i + j)) * det(pop(A, i, j)))
    b.append(c)
  return b

def inverse(A):
  return (1 / det(A) * np.transpose(com(A)))


def cylindre_plein(n, R, h):
  pts = []
  m = int(n**(1 / 3))
  if m < 1:
    m = 1
  for i in range(m):
    # rayon de 0 à R
    r = R * i / (m - 1) if m > 1 else R
    for j in range(m):
      theta = 2 * pi * j / m
      for k in range(m):
        # z de -h/2 à h/2
        z = -h / 2 + (h * k / (m - 1)) if m > 1 else 0
        x = r * cos(theta)
        y = r * sin(theta)
        pts.append([x, y, z])
  return transpose(pts)


def translation(m, F, G, vG, h):
  new_g = 0
  new_vg = 0

  return new_g, new_vg

def rotation(I_inv,F,G,teta,omega,h):
  new_teta = 0
  new_omega = 0

  return

# Tracé cylindre plein
(X, Y, Z) = cylindre_plein(1000, 1, 9)
(Xa, Ya, Za) = cylindre_plein(1000, 2, 2)

x_offset = 0   # décalage en X
y_offset = 0  # décalage en Y
z_offset = 2   # décalage en Z

X, Y, Z = np.array(X), np.array(Y), np.array(Z) #conversion requise

X_translated = X + x_offset
Y_translated = Y + y_offset
Z_translated = Z + z_offset

fig = plt.figure()
ax = plt.axes(projection='3d')
ax.set_title("The star sword !")

ax.scatter3D(X_translated, Y_translated, Z_translated, c=Z_translated, cmap='autumn', s=5)
ax.scatter3D(Xa, Ya, Za)

plt.show()
