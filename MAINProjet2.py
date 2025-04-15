import LIBProjet2 as LIB
import math
import matplotlib.pyplot as plt


# Somme de deux vecteurs
def somme_vect(A, B):
    return [A[i] + B[i] for i in range(len(A))]

# Produit vecteur * scalaire
def prod_vect_scal(A, k):
    return [k * A[i] for i in range(len(A))]

# Vecteur B - A
def vect(A, B):
    return [B[i] - A[i] for i in range(len(A))]

# Produit vectoriel A ^ B
def prod_vect(A, B):
    return [
        A[1]*B[2] - A[2]*B[1],
        A[2]*B[0] - A[0]*B[2],
        A[0]*B[1] - A[1]*B[0]
    ]

#Calcule la nouvelle position et nouvelle vitesse de G
def translation(m,F,G,vG,h):
    SF=[0,0,0]
    for i in range(len(F)):
        SF=somme_vect(SF,F[i][0])
        aG=prod_vect_scal(SF,1/m)
        newG=somme_vect(G,prod_vect_scal(vG,h))
        newvG=somme_vect(vG,prod_vect_scal(aG,h))
    return(newG,newvG)

#translate un solide W d'un vecteur GnewG
def translate(W,G,newG):
    newW=[]
    GnewG=vect(G,newG)
    for i in range(len(W[0])):
        newW.append(somme_vect([W[0][i],W[1][i],W[2][i]],GnewG))
    return(LIB.transpose(newW))
#Test tracé translation
# R=0.5
# H=4
# W=LIB.cylindre_plein(15,15,15,R,H)
# m=20
# G=[0,0,0]
# vG=[0,0,0]
# F=[[[0,0,20],[0,0,0]]]

# tmax=5

# n=10
# h=tmax/10
# newW=W
# newG=G
# newvG=vG

# for i in range(n):
#     (X,Y,Z)=newW
#     LIB.plot3D(X,Y,Z, 'autumn')
#     (newG,newvG)=translation(m,F,newG,newvG,h)
#     newW=translate(W,G,newG)



def rotation(I_inv, F, G, theta, omega, h):
    # 1. Calcul du moment résultant M_O
    M = [0, 0, 0]
    for force, point in F:
        OG = vect(G, point)               # vecteur OG = point d'application - G
        moment = prod_vect(OG, force)     # moment = OG ^ force
        M = somme_vect(M, moment)         # somme des moments

    # 2. Calcul de l’accélération angulaire alpha = I_inv * M
    # Transformer M en vecteur colonne [[x], [y], [z]]
    M_col = [[x] for x in M]
    alpha_col = LIB.prodmat(I_inv, M_col)  # produit matriciel I_inv * M
    alpha = [a[0] for a in alpha_col]      # retransforme en liste

    # 3. Mise à jour de la vitesse angulaire : omega' = omega + h * alpha
    omega_new = somme_vect(omega, prod_vect_scal(alpha, h))

    # 4. Mise à jour de l’angle : theta' = theta + h * |omega|
    norme_omega = math.sqrt(sum(w**2 for w in omega))
    theta_new = theta #+ h * norme_omega

    return theta_new, omega_new

def rotation_solide(W, M):
    newW = []
    for i in range(len(W[0])):
        P = [[W[0][i]], [W[1][i]], [W[2][i]]]
        P_rot = LIB.prodmat(M, P)  # M * P
        newW.append([P_rot[0][0], P_rot[1][0], P_rot[2][0]])
    return LIB.transpose(newW)

def mat_rot_general(omega, h):
    theta_x = omega[0] * h
    theta_y = omega[1] * h
    theta_z = omega[2] * h

    # Matrices élémentaires
    Rx = [
        [1, 0, 0],
        [0, LIB.cosinus(theta_x), -LIB.sinus(theta_x)],
        [0, LIB.sinus(theta_x),  LIB.cosinus(theta_x)]
    ]
    Ry = [
        [LIB.cosinus(theta_y), 0, LIB.sinus(theta_y)],
        [0, 1, 0],
        [-LIB.sinus(theta_y), 0, LIB.cosinus(theta_y)]
    ]
    Rz = [
        [LIB.cosinus(theta_z), -LIB.sinus(theta_z), 0],
        [LIB.sinus(theta_z),  LIB.cosinus(theta_z), 0],
        [0, 0, 1]
    ]

    # Produit matriciel total : R = Rz * Ry * Rx
    RyRx = LIB.prodmat(Ry, Rx)
    R = LIB.prodmat(Rz, RyRx)
    return R


# I = matrice inertie (3x3)
# G = point de base de la matrice I
# A = nouveau point ou la matrice sera déplacé (arrivée)
# M = masse du solide
# Dans la formule a=x b=y c=z


def deplace_mat(I,m,G,A):
  # G = [g₁, g₂, g₃] | A = [a₁, a₂, a₃]
  # d = [g₁ - a₁, g₂ - a₂, g₃ - a₃]
  d = []
  for i in range(3):
    diff = G[i] - A[i]
    d.append(diff)
  diagonal = d[0]**2 + d[1]**2 + d[2]**2

  I_A = []

  for l in range(3):
    new_ligne = []
    for c in range(3):
      # delta = 1 si on est sur la diagonale, sinon 0
      delta = 1 if l == c else 0
      # Application de la formule : I_A = I + m(d² * delta - d[l]d[c])
      new_val = I[l][c] + m * (diagonal * delta - d[l] * d[c])
      new_ligne.append(new_val)
    I_A.append(new_ligne)

  return I_A


def Somme_mat(a,b):
    return [[a[i][j] + b[i][j]  for j in range(len(a[0]))] for i in range(len(a))]



def solide(n):
    r1, h1 = 1, 9    # cylindre 1 (poignée+lame)
    r2, h2 = 2, 1
    W1 = LIB.cylindre_plein(10, 10, 10, r1, h1)
    W2 = LIB.cylindre_plein(10, 10, 10, r2, h2)
    W3 = LIB.cylindre_plein(10, 10, 10, r2, h2)
    W2_translated = translate(W2, [0, 0, 0], [0, 0, -2])
    W3_translated = translate(W3, [0, 0, 0], [0, 0, +2])

    W1 = somme_vect(W1,W2_translated)

    return somme_vect(W1, W3_translated)

W = solide(1000)
(X, Y, Z) = W

#Test tracé rotation
R1=1
H1=9
m1=2269
m2=942.5
m=m1+m2
G=[0,0,0]
vG=[0,0,0]
FT=[[[0,0,0],[0,0,0]], [[0,0,-9.8*m*0],[0,0,0]]]
FR=[[[0,2000,0],[0,0,3]]]

tmax=5

n=200
h=tmax/10
newW=W
newG=G
newvG=vG

I1 = [[((R1**2)/4) + ((H1**2)/12), 0, 0],
     [0, ((R1**2)/4) + ((H1**2)/12), 0],
     [0, 0, (R1**2)/2]]



I1 = LIB.prodmatscal(I1,m1)
R2=2
H2=1

G2=[0,0,-2]
vG2=[0,0,0]



I2 = [[((R2**2)/4) + ((H2**2)/12), 0, 0],
     [0, ((R2**2)/4) + ((H2**2)/12), 0],
     [0, 0, (R2**2)/2]]

I2 = LIB.prodmatscal(I2,m2)

Idm = deplace_mat(I2,m2,G2,G)
I = Somme_mat(I1,Idm)
Iinv = LIB.inverse(I)

# print (I1)
# print (I2)
# print(Idm)
# print (I)
# print (Iinv)






theta = 0
omega = [0, 0, 0]
theta_new = theta
omega_new = omega


plt.ion()  # Active le mode interactif

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

for i in range(n):
        # Efface le tracé précédent
        ax.clear()
        # Affiche le solide en 3D avec les mêmes limites et palette de couleurs
        ax.scatter3D(X, Y, Z, c=Z, cmap='autumn')
        ax.set_xlim(-10, +10)
        ax.set_ylim(-10, +10)
        ax.set_zlim(-10, +10)

        # Actualise l'affichage
        plt.draw()
        plt.pause(0.1)  # Pause de 0.1 seconde (adapter si besoin)

        # Mise à jour des positions et rotations
        (newG, newvG) = translation(m, FT, newG, newvG, h)
        newW = translate(newW, G, newG)
        (theta_new, omega_new) = rotation(Iinv,  FR, G, theta_new, omega, h)
        R = mat_rot_general(omega_new, h)
        newW = rotation_solide(newW, R)

        # Mise à jour des variables X, Y, Z pour le prochain tracé
        (X, Y, Z) = newW

plt.ioff()
plt.show()

#LIB.plot3D(X,Y,Z, 'autumn')

