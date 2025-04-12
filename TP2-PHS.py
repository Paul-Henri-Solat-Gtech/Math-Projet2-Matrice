import matplotlib.pyplot as plt
import numpy as np

tabMatrice = ([5, 0], [-2, 1])


# a = matrice
# i = ligne
# j = colonne
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


print(pop(tabMatrice, 1, 1))


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


print(tabMatrice)
print(inverse(tabMatrice))
