import TP2BM as TP2
#Somme 2 vecteurs A et B
def somme_vect(A,B):
    C=[]
    for i in range(len(A)):
        C.append(A[i]+B[i])
    return(C)

#Multiplie un vecteur par un scalaire
def prod_vect_scal(A,k):
    B=[]
    for i in range (len(A)):
        B.append(k*A[i])
    return(B)

#Calcul du vecteur AB
def vect(A,B):
    return(somme_vect(prod_vect_scal(A,-1),B))

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
W=TP2.cylindre_plein(10000,R,H)
m=20
G=[0,0,0]
vG=[0,0,0]
F=[[[0,0,10],[1,0,-2]],[[0,0,5],[-1,0,-2]],[[-2,0,0],[1,0,0]]]

tmax=5

n=10
h=tmax/10
newW=W
newG=G
newvG=vG

for i in range(n):
    (X,Y,Z)=newW
    TP2.plot3D(X,Y,Z)
    (newG,newvG)=translation(m,F,newG,newvG,h)
    newW=translate(W,G,newG)
