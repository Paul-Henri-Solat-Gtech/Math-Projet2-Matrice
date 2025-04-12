import numpy as np
from math import pi
import matplotlib.pyplot as plt
#from mpl_toolkits.mplot3d import Axes3D

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

def rectangle_plein(nbPoints, L, l, z=0):
    """
    Génère un maillage de points pour un rectangle plein.
    
    nbPoints : nombre approximatif total de points souhaité
    L : dimension en x (longueur) du rectangle
    l : dimension en y (largeur) du rectangle
    z : coordonnée z constante pour le plan du rectangle
    """
    # Détermination du nombre de points sur chaque direction en fonction de nbPoints
    # Ici on choisit une répartition proportionnelle aux dimensions L et l
    n_x = int(np.sqrt(nbPoints * (L / (L + l))))
    n_y = int(np.sqrt(nbPoints * (l / (L + l))))
    # S'assurer d'avoir au moins 2 points pour définir la grille dans chaque direction
    if n_x < 2:
        n_x = 2
    if n_y < 2:
        n_y = 2

    # Générer une répartition régulière de points
    x_coords = np.linspace(-L/2, L/2, n_x)
    y_coords = np.linspace(-l/2, l/2, n_y)
    
    pts = []
    for x in x_coords:
        for y in y_coords:
            pts.append([x, y, z])
    return np.array(pts)

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

# Paramètres pour le cylindre
nbPointsCylindre = 10000
R = 0.5
h = 4
cyl_pts = cylindre_plein(nbPointsCylindre, R, h)

# Paramètres pour le rectangle (la lame)
# Par exemple, on choisit une lame de 10 unités de long et 1 unité de large,
# et on la place en z légèrement au-dessus du centre du cylindre (ici, z = h/2 + 0.5)
nbPointsRectangle = 500  # nombre approximatif de points pour la lame
L = 10  # longueur de la lame
l_rect = 1  # largeur de la lame
z_lame = h/2 + 0.5  # position en z, au-dessus du cylindre

rect_pts = rectangle_plein(nbPointsRectangle, L, l_rect, z=z_lame).T  # 3 x n_points

# Tracé
fig = plt.figure()
ax = plt.axes(projection='3d')

# Tracer le cylindre
ax.scatter(cyl_pts[0], cyl_pts[1], cyl_pts[2], c=cyl_pts[2], cmap='autumn', s=1, label="Cylindre plein")

# Tracer le rectangle
ax.scatter(rect_pts[0], rect_pts[1], rect_pts[2], c='blue', s=20, label="Rectangle plein (lame)")

ax.set_title("Cylindre plein avec rectangle plein placé au-dessus")
ax.legend()
plt.show()

(X, Y, Z) = cylindre_plein(10000, 0.5, 4)
fig = plt.figure()
ax = plt.axes(projection='3d')
ax.scatter3D(X, Y, Z, c=Z, cmap='autumn')
ax.set_title("Cylindre plein")

# plt.show()