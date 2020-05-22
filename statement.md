# Bonjour!

Ce programme présente la factorisation par la méthode de Lenstra qui utilise des courbes elliptiques sur les entiers modulaires.  
Cette version travaille en coordonnées projectives afin d'éviter le recours aux inverses modulaires, lourds à calculer.

```python runnable
from time import time
from random import randint
from math import gcd


# test de Miller pour un témoin a
# Retourne faux si n est composé et vrai si n est probablement premier
# d et r doivent vérifier n = 2^r * d + 1 avec d impair   

def millerTest(a,d,n,r):
         
    x = pow(a, d, n) # Calcule a^d % n   
    if (x == 1  or x == n-1): 
       return True

    for _ in range(r):    
        x = (x * x) % n 
        if (x == 1):
            
            return False 
        if (x == n-1):
            return True    
    
    return False 

# Test de primalité de Miller Rabin  
# Si faux alors n est composé et si vrai alors n est probablement premier 
# k determine le niveau de certitude : P(erreur) < 1/4**k

def isPrime(n, k=25): 
       
    if (n <= 1 or n == 4):
        return False 
    if (n <= 5):
        return True   
    
    # Trouver d et r tels que n = 2^r * d + 1 avec d impair 
    d = n - 1 
    r = 0
    while (d&1 == 0): 
        d  >>= 1 
        r += 1 
    
    # Effectuer k tests de Miller
    for i in range(k):
        a = randint(2,n-2) 
        if (not millerTest(a, d, n, r)):
              return False  
    return True

# Premier suivant n

def nextPrime(n): 
    while not isPrime(n):
        n += 1
    return n

# Crible d'Eratosthenes

def sieve(n):
    b = [True] * n
    ps = []
    for p in range(2, n):
        if b[p]:
            ps.append(p)
            for i in range(p, n, p):
                b[i] = False
    return ps

# Duplication d'un point d'une courbe elliptique (version projective)

def doubleECM(p, a, m):
    px,py,pz=p
    if px==0:
        return p
    if py==0:
        return 0, 1, 0  # Infini
    zz = pz*pz
    dy = (3*px*px + a*zz) %m
    dx = (2*py*pz) %m
    dx2 = dx*dx %m
    dx3 = dx2*dx %m
    
    v = 2*dx*px*py %m
    w = dy*dy - 2*v %m
    rx = dx*w
    ry = dy * (v - w) - 2*dx2*py*py
    rz = dx3
    return rx %m, ry %m, rz %m

# Addition de deux points d'une courbe elliptique (version projective)

def addECM(p, q, a, m): 
    
    px,py,pz=p
    qx,qy,qz=q
    # Si un point est infini, retourne l'autre
    if pz == 0: return q
    if qz == 0: return p
    
    y0 = py*qz %m
    y1 = qy*pz %m
    x0 = px*qz %m
    x1 = qx*pz %m
    
    if x0 == x1:
        if (y0 + y1) %m == 0:
            return 0, 1, 0  # Infini        
        return doubleECM(p, a, m)
    
    dy = (y0 - y1) %m
    dx = (x0 - x1) %m
    
    zz = pz*qz %m
    dx2 =dx*dx %m
    dx3 =dx2*dx %m
    x = (dy*dy*zz - (x0 + x1)*dx2) %m
    y = (dy*(x0*dx2 - x)- y0*dx3) %m
    x *= dx
    z = dx3*zz %m
    
    #if gcd(z,m) >1:
        #print('z found',gcd(z,m))
    return x %m, y, z



# Multiplication par un entier (version projective)

def mulECM(k, p, a, m):
    r = (0, 1, 0)  # Infinity
    while k > 0:                
        if k % 2 == 1:
            r = addECM(p, r, a, m)            
        k >>=1
        p = doubleECM(p, a, m)       
    return r



# Algorithme de Lenstra (version projective)
# Limit specifie le maximum des premiers à tester
# Primes est une liste de nombres premiers < limit
def lenstra(n, limit=1000, primes=None):
  
    if primes is None:
        primes = sieve(limit)
   
    g = n
    while g == n:
        # x et y sont aléatoires
        q = randint(0, n - 1), randint(0, n - 1), 1
        # a est aléatoire, b est calculé pour que le point soit sur la courbe
        a = randint(0, n - 1)
        b = (q[1] * q[1] - q[0] * q[0] * q[0] - a * q[0]) % n
        g = gcd(4 * a * a * a + 27 * b * b, n)  # vérification de singularité
    # Par chance, g peut être un facteur
    if g > 1:
        return g
    # Multiplier le point q par le PPCM des nombres (1, ..., limit)
    for p in primes:
        pp = p
        while pp < limit:
            q = mulECM(p, q, a, n)
            
            x,y,z = q
            # Vérification de l'équation de la courbe
            # Les 5 lignes suivantes peuvent être supprimées pour améliorer la performance
            c = (x*x*x + a * x*z*z + b * z*z*z - y*y*z) %n
            if c != 0:  # Rupture de l'arithmétique elliptique
                print("q n'est pas sur la courbe après multiplication par", p)
                g = gcd(z, n)
                return g

            # Recherche d'un facteur dans z
            # Les 5 lignes suivantes peuvent être supprimées pour améliorer la performance  
            if z > 1:                
                g = gcd(z, n)
                if g >1:
                    #print('g=',g, 'n=',n)
                    return g
                                
            pp = p * pp
    
    g = gcd(z, n)
    #print('lenstra : g=',g)
    if g>1:
        return g
    else:
        return n


def factorECM(n, k=25, limit=10000, div=lenstra, primes=None):
    # Décompose n en facteurs probablement premiers
    if primes is None:
         primes = sieve(limit)    
    if n in primes:
       return [n]
    for p in primes:
        if n%p==0:
            return [p]+factorECM(n//p, primes=primes)
    if isPrime(n,k): # on pense que n est premier
        return [n]
    g = 1
    i = 0  
    while g==1 or g==n: # tant qu'on n'a pas de facteur
        # on essaye avec d'autres valeurs
        i +=1
        if i%10==0:
            print("retry",i)
        g=div(n, limit=limit, primes=primes)                  
    return factorECM(g, primes=primes)+factorECM(n//g, primes=primes)
    
print('Les premiers <20 sont :',sieve(20))
print()

for i in range(2):
    a=nextPrime(randint(1e8,1e12))    
    b=nextPrime(randint(1e8,1e10))
    #c=randint(1e8,1e10) 
    n=a*b
    print('n=',n)
    t=time()
    f=factorECM(n)  
    t=time()-t  
    print('les facteurs de n sont',f,'t=',t,'s')
    print()

```

