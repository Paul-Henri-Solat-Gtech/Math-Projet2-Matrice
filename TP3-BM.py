import TP2BM as TP2
import math
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
    return(TP2.transpose(newW))
#Test tracé translation
R=0.5
H=4
W=TP2.cylindre_plein(15,15,15,R,H)
m=20
G=[0,0,0]
vG=[0,0,0]
F=[[[0,0,20],[0,0,0]]]

tmax=5

n=10
h=tmax/10
newW=W
newG=G
newvG=vG

for i in range(n):
    (X,Y,Z)=newW
    TP2.plot3D(X,Y,Z, 'autumn')
    (newG,newvG)=translation(m,F,newG,newvG,h)
    newW=translate(W,G,newG)



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
    alpha_col = TP2.prodmat(I_inv, M_col)  # produit matriciel I_inv * M
    alpha = [a[0] for a in alpha_col]      # retransforme en liste

    # 3. Mise à jour de la vitesse angulaire : omega' = omega + h * alpha
    omega_new = somme_vect(omega, prod_vect_scal(alpha, h))

    # 4. Mise à jour de l’angle : theta' = theta + h * |omega|
    norme_omega = math.sqrt(sum(w**2 for w in omega))
    theta_new = theta + h * norme_omega

    return theta_new, omega_new

def rotation_solide(W, M):
    newW = []
    for i in range(len(W[0])):
        P = [[W[0][i]], [W[1][i]], [W[2][i]]]
        P_rot = TP2.prodmat(M, P)  # M * P
        newW.append([P_rot[0][0], P_rot[1][0], P_rot[2][0]])
    return TP2.transpose(newW)

def mat_rot_general(omega, h):
    theta_x = omega[0] * h
    theta_y = omega[1] * h
    theta_z = omega[2] * h

    # Matrices élémentaires
    Rx = [
        [1, 0, 0],
        [0, TP2.cosinus(theta_x), -TP2.sinus(theta_x)],
        [0, TP2.sinus(theta_x),  TP2.cosinus(theta_x)]
    ]
    Ry = [
        [TP2.cosinus(theta_y), 0, TP2.sinus(theta_y)],
        [0, 1, 0],
        [-TP2.sinus(theta_y), 0, TP2.cosinus(theta_y)]
    ]
    Rz = [
        [TP2.cosinus(theta_z), -TP2.sinus(theta_z), 0],
        [TP2.sinus(theta_z),  TP2.cosinus(theta_z), 0],
        [0, 0, 1]
    ]

    # Produit matriciel total : R = Rz * Ry * Rx
    RyRx = TP2.prodmat(Ry, Rx)
    R = TP2.prodmat(Rz, RyRx)
    return R

#Test tracé rotation
# R=0.5
# H=4
# W=TP2.cylindre_plein(1000,R,H)
# m=20
# G=[0,0,0]
# vG=[0,0,0]
# F=[[[0,0,10],[1,0,-2]],[[0,0,5],[-1,0,-2]],[[-2,0,0],[1,0,0]]]

# tmax=5

# n=20
# h=tmax/10
# newW=W
# newG=G
# newvG=vG

# I = [[((R**2)/4) + ((h**2)/12), 0, 0],
#     [0, ((R**2)/4) + ((h**2)/12), 0],
#    [0, 0, (R**2)/2]]

# Iinv = TP2.inverse(I)

# theta = 0
# omega = [0, 0.0, 0]
# theta_new = theta
# omega_new = omega

# for i in range(n):
#     (X,Y,Z)=newW
#     TP2.plot3D(X,Y,Z)
#     (theta_new,omega_new) = rotation(Iinv, F, G, theta, omega, h)
#     R=mat_rot_general(omega_new,h)
#     newW=rotation_solide(newW,R)
