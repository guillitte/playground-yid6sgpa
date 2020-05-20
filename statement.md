# Bonjour!

Ce programme présente la factorisation par la méthode de Lenstra qui utilise des courbes elliptiques sur les entiers modulaires

```python runnable
from time import time
from random import randint
from math import gcd

def millerTest(a,d,n,r):
    # test de Miller pour un témoin a
    # Retourne faux si n est composé et vrai si n est probablement premier
    # d et r doivent vérifier n = 2^r * d + 1 avec d impair   
           
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

def isPrime(n, k=25): 
    # Test de primalité de Miller Rabin  
    # Si faux alors n est composé et si vrai alors n est probablement premier 
    # k determine le niveau de certitude : P(erreur) < 1/4**k
    
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

def nextPrime(n):
    # premier suivant n
    while not isPrime(n):
        n += 1
    return n

# Sieve of Eratosthenes
def sieve(n):
    b = [True] * n
    ps = []
    for p in range(2, n):
        if b[p]:
            ps.append(p)
            for i in range(p, n, p):
                b[i] = False
    return ps


def doubleECM(p, a, m):
    px,py,pz=p
    if px==0:
        return p
    if py==0:
        return 0, 1, 0  # Infinity 
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


def addECM(p, q, a, m):
    
    px,py,pz=p
    qx,qy,qz=q
    # If one point is infinity, return other one
    if pz == 0: return q
    if qz == 0: return p
    
    y0 = py*qz %m
    y1 = qy*pz %m
    x0 = px*qz %m
    x1 = qx*pz %m
    
    if x0 == x1:
        if (y0 + y1) %m == 0:
            return 0, 1, 0  # Infinity        
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



# Multiplication (repeated addition and doubling)
def mulECM(k, p, a, m):
    r = (0, 1, 0)  # Infinity
    while k > 0:                
        if k % 2 == 1:
            r = addECM(p, r, a, m)            
        k >>=1
        p = doubleECM(p, a, m)       
    return r



# Lenstra's algorithm for factoring
# Limit specifies the amount of work permitted
def lenstra(n, limit=1000, primes=None):
  
    if primes is None:
        primes = sieve(limit)
   
    g = n
    while g == n:
        # Randomized x and y
        q = randint(0, n - 1), randint(0, n - 1), 1
        # Randomized curve coefficient a, computed b
        a = randint(0, n - 1)
        b = (q[1] * q[1] - q[0] * q[0] * q[0] - a * q[0]) % n
        g = gcd(4 * a * a * a + 27 * b * b, n)  # singularity check
    # If we got lucky, return lucky factor
    if g > 1:
        return g
    # increase k step by step until lcm(1, ..., limit)
    for p in primes:
        pp = p
        while pp < limit:
            q = mulECM(p, q, a, n)
            
            x,y,z = q
            c = (x*x*x + a * x*z*z + b * z*z*z - y*y*z) %n
            if c != 0:  # Elliptic arithmetic breaks
                print('q is not on curve', p)
                g = gcd(z, n)
                return g
                
            if z > 1:                
                g = gcd(z, n)
                if g >1:
                    print('g=',g, 'n=',n)
                    return g
                                
            pp = p * pp
    #print('lenstrap : g=',g)
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
s1,s2=0,0
for i in range(2):
    #a=nextPrime(randint(1e8,1e15))    
    #b=nextPrime(randint(1e8,1e16))
    c=randint(1e8,1e10) 
    n=c
    print('n=',n)
    t=time()
    f=factorECM(n)  
    t=time()-t  
    print('les facteurs de n sont',f,'t=',t,'s')
    print()

```

