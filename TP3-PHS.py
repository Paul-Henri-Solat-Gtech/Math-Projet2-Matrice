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

# variables locales pour le tp
F=[]

# Tracé d'un cylindre plein
(X, Y, Z) = cylindre_plein(10000, 0.5, 4)
fig = plt.figure()
ax = plt.axes(projection='3d')
ax.scatter3D(X, Y, Z, c=Z, cmap='autumn')
ax.set_title("Cylindre plein")

plt.show()
