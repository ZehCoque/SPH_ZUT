import math 
from numpy import array,seterr
import time
seterr(all='raise')
# Function that returns True if n  
# is prime else returns False  
def isPrime(n): 
      
    # Corner cases  
    if(n <= 1): 
        return False
    if(n <= 3): 
        return True
      
    # This is checked so that we can skip  
    # middle five numbers in below loop  
    if(n % 2 == 0 or n % 3 == 0): 
        return False
      
    for i in range(5,int(math.sqrt(n) + 1), 6):  
        if(n % i == 0 or n % (i + 2) == 0): 
            return False
      
    return True
  
# Function to return the smallest  
# prime number greater than N  
def nextPrime(N): 
  
    # Base case  
    if (N <= 1): 
        return 2
  
    prime = N 
    found = False
  
    # Loop continuously until isPrime returns  
    # True for a number greater than n  
    while(not found): 
        prime = prime + 1
  
        if(isPrime(prime) == True): 
            found = True
  
    return prime 

class Hashing:

    def __init__(self,cell_size,N,d={}):
        self.d = d
        self.cell_size = cell_size
        self.table_size = N
        self.p1 = 73856093
        self.p2 = 19349669
        self.p3 = 83492791
    
    def r_c(self,point):

        return [math.floor(point[0]/self.cell_size),
            math.floor(point[1]/self.cell_size),
            math.floor(point[2]/self.cell_size)]

    def _hash(self,point):
        r = array([point[0]*self.p1,
                   point[1]*self.p2,
                   point[2]*self.p3])
        
        return (r[0]^r[1]^r[2]) % self.table_size

    def _add(self,point,obj):
        r_c = self.r_c(point)
        hh = self._hash(r_c)
        if hh not in self.d.keys():
            self.d[hh] = list()
        if obj not in self.d[hh]:
            self.d[hh].append(obj)
            for i in self.d:
                if obj in self.d[i] and i!=hh:
                    for j in self.d[i]:
                        if j == obj:
                            ind = self.d[i].index(j)
                            del self.d[i][ind]
                            return
            #self._del(point,obj)
            
    def _del(self,point,obj):
        i = self.d.get(self._hash(point))
        for j in i:
            if j == obj:
                del i[j]

    def possible_neighbors(self,point):
        point = array(point)
        L = []
        min_point = point - self.cell_size
        max_point = point + self.cell_size

        BBmin = self.r_c(min_point)
        BBmax = self.r_c(max_point)
        x_count = BBmin[0]
        while x_count <= BBmax[0]:
            y_count = BBmin[1]
            while y_count <= BBmax[1]:
                z_count = BBmin[2]
                while z_count <= BBmax[2]:
                    new_point = [x_count,y_count,z_count]
                    i = self.d.get(self._hash(new_point))
                    if i != None:
                        for j in i:
                            L.append(j)
                    z_count += 1
                y_count += 1
            x_count += 1
        
        return L

# h = 10
# N = 1000

# #Hashing(h,N)._add([10,11,32])

# start = time.time()

# print(Hashing(h,N)._add([10,25,32],88))
# print(Hashing(h,N)._add([10,11,32],68))
# print(Hashing(h,N)._add([10,11,32],67))
# print(Hashing(h,N)._add([10,14,48],65))

# x = Hashing(h,N).possible_neighbors([10,24,32])
# print(x)
# print(len(x))
# print(start-time.time())