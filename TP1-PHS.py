import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D

def transpose(M):
  Mt = [[0 for _ in range(len(M))] for _ in range(len(M[0]))]

  for i in range(len(M)):  # Parcours les lignes de M
    for j in range(len(M[0])):  # Parcours les colonnes de M
      Mt[j][i] = M[i][j]  # Échange les éléments

  return Mt

def ligne(nbPoints, xmin, xmax):
  W = []
  dx = (xmax - xmin) / (nbPoints - 1)
  for x in range(nbPoints):
    W.append([x * dx, 0, 0])
  return (transpose(W))

def carrevide(nbPoint, longueur):
  W = []
  # Utiliser la longueur pour définir le demi-côté
  demi_cote = longueur / 2  
  # On répartit nbPoint de manière égale sur les 4 côtés
  points_par_cote = nbPoint // 4  
  if points_par_cote < 2:
      points_par_cote = 2

  # Côté bas (de -demi_cote à demi_cote, y constant = -demi_cote)
  for i in range(points_par_cote):
      x = -demi_cote + i * (longueur / (points_par_cote - 1))
      y = -demi_cote
      W.append([x, y, 0])

  # Côté droit (x constant = demi_cote, de -demi_cote à demi_cote)
  for i in range(points_par_cote):
      x = demi_cote
      y = -demi_cote + i * (longueur / (points_par_cote - 1))
      W.append([x, y, 0])

  # Côté haut (de demi_cote à -demi_cote, y constant = demi_cote)
  for i in range(points_par_cote):
      x = demi_cote - i * (longueur / (points_par_cote - 1))
      y = demi_cote
      W.append([x, y, 0])

  # Côté gauche (x constant = -demi_cote, de demi_cote à -demi_cote)
  for i in range(points_par_cote):
      x = -demi_cote
      y = demi_cote - i * (longueur / (points_par_cote - 1))
      W.append([x, y, 0])

  return W

# ##### CARRE VIDE #####
n = 40  # Nombre total de points sur la frontière
a = 10  # Longueur du côté du carré
points = carrevide(n, a)

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# Extraire X, Y, Z
X = [p[0] for p in points]
Y = [p[1] for p in points]
Z = [p[2] for p in points]

# Tracer les points
ax.scatter(X, Y, Z, color='red')
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
ax.set_title('Carré vide 3D')

plt.show()