#Project 2
#part 1
def BackSub(U,b):
    n = len(U)
    y = []
    #start from row n-1 and loop through rows backwards
    for i in reversed(range(n)):
        s = 0
        for j in range(n-i):
            #store corresponding entry from RHS vector into the solution vector
            if i == n-1:
                s = b[i]
            if j == 0:
                s = b[i]
            else:
                #depending on the row, plug in known values from solution vector
                #and subtract to find the unknown entry
                s = s - U[i][i+j]*y[j-1]
        #divide the RHS quantity by the coefficient of the entry of interest
        s = s/U[i][i]
        y.append(s)
    y.reverse()
    return y

def ForwSub(L,b):
    n = len(L)
    y = []
    #loop through rows
    for i in range(n):
            s = 0
            for j in range(i+1):
                #store corresponding entry from RHS vector into solution vector
                if i == 0:
                    s = b[i]
                if j == 0:
                    s = b[i]
                else:
                    #depending on the row, plug in known entries from solution 
                    #vector and subtract to find the currently unknown entry
                    s = s - L[i][j-1]*y[j-1]
            #divide the current RHS by the coefficient of the entry of interest
            s = s/L[i][i]
            y.append(s)
    return y

upper = [[-1,0,7,3],[0,1,2,5],[0,0,4,2],[0,0,0,-7]]
lower = [[-1,0,0,0],[0,1,0,0],[7,2,4,0],[3,5,2,-7]]
bvec = [9,8,6,-7]
print("upper triangular matrix")
for i in upper:
    print(i)
print("vector")
for i in bvec:
    print(i)
print("result of back substituting with the above matrix and vector")
print(BackSub(upper,bvec))
print("\nlower triangular matrix")
for i in lower:
    print(i)
print("vector")
for i in bvec:
    print(i)
print("result of forward substituting with the above matrix and vector")
print(ForwSub(lower, bvec))

#part 2
def MultiplyNNMatrix(A,B):
    n = len(A)
    prod = []
    for i in range(n):
        prod.append([])
        for j in range(n):
            prod[i].append(0.0)
    
    for i in range(n):
        for j in range(n):
            for k in range(n):
                prod[i][j] = prod[i][j] + A[i][k]*B[k][j]
    prod = [[round(p,3) for p in w] for w in prod]
    return prod

def GE0(A):
    U = A
    n = len(U)
    L = []
    #initialize L matrix to the identity matrix
    for i in range(n):
        L.append([])
        for j in range(n):
            if i == j:
                L[i].append(1)
            else:
                L[i].append(0)
    #loop through rows
    for i in range(1,n):
        fix = i-1
        #fix the column to eliminate and loop through all rows beneath
        for j in range(i,n):
            factor = U[j][fix]/U[fix][fix]
            L[j][fix] = factor
            #eliminate entries to make an upper triangular matrix and update U
            U[j] = [a-u*factor for a,u in zip(U[j],U[fix])] 
    return L,U

testmatrix = [[2,-1,0,0],[-1,2,-1,0],[0,-1,2,-1],[0,0,-1,2]]
testmatrix1 = [[2,-1,0,0],[-1,2,-1,0],[0,-1,2,-1],[0,0,-1,2]]
lo, up = GE0(testmatrix)
print("\nA matrix")
for i in testmatrix1:
    print(i)
print("resulting L matrix")
for i in lo:
    print(i)
print("resulting U matrix")
for i in up:
    print(i)
print("Is A = LU?")
print(MultiplyNNMatrix(lo,up)== testmatrix1)

#part 3
def RowExchange(M,k,l):
    #assumes k and l are between 0 and len(M) and M is nxn
    temp = []
    for u in range(len(M)):
        temp.append(M[l][u])
    M[l] = M[k]
    M[k] = temp
    return M
I = [[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1]]

def Lswitch(M,n, row1, row2):
    l = len(M)
    for p in range(n):
        temp = M[row1][p]
        M[row1][p] = M[row2][p]
        M[row2][p] = temp
    return M


def GE1(A):
    U = [[a for a in w] for w in A]
    n = len(A)
    L = []
    onex = False
    P = []
    #initialize L and P matrices to the identity matrix
    for i in range(n):
        L.append([])
        P.append([])
        for j in range(n):
            if i == j:
                L[i].append(1)
                P[i].append(1)
            else:
                L[i].append(0)
                P[i].append(0)
    #loop through rows
    for i in range(1,n):
        fix = i-1
        max = fix
        factor = 0
        #loop through rows under a fixed diagonal entry
        for j in range(i,n):
            #maximize diagonal entries through row switches and keep track
            # of switches using P matrix, and update L and U matrices accordingly
            if(abs(U[j][fix])>abs(U[max][fix])):
                max = j
                RowExchange(U,fix,max)
                RowExchange(P,fix,max)
                Lswitch(L,j-1,fix,max)
            factor = U[i][fix]/U[fix][fix]
            L[j][fix] = factor
            U[i] = [a-u*factor for a,u in zip(U[i],U[fix])]

    #L = [[round(l,3) for l in w] for w in L]
    #U = [[round(u,3) for u in w] for w in U]
    return P,L,U

A = [[0.2,-0.1,0,0],[1,-2,1,0],[0,1,-2,1],[0,0,-10,20]]
A1 = [[0.2,-0.1,0,0],[1,-2,1,0],[0,1,-2,1],[0,0,-10,20]]
print("\nA matrix")
for i in A1:
    print(i)
P, L, U= GE1(A)
print("P matrix")
for i in P:
    print(i)
print("L matrix")
for i in L:
    print(i)
print("U matrix")
for i in U:
    print(i)
print("Is PA = LU?")
print(MultiplyNNMatrix(P,A)==MultiplyNNMatrix(L,U))

#part 4
def MultiplyMatrixVector(M,v):
    #assumes nxn matrix and n tuples vector
    n = len(M)
    prod = []
    for i in range(n):
        prod.append(0)
    for i in range(n):
        for j in range(n):
            prod[i] = prod[i]+M[i][j]*v[j]
    return prod

def Solve(A, RHS):
    # use P,L,U decomposition to solve the given system
    P,L,U = GE1(A)
    r1 = MultiplyMatrixVector(P,RHS)
    y = ForwSub(L,r1)
    x = BackSub(U,y)
    return x

A = [[0.2,-0.1,0,0],[1,-2,1,0],[0,1,-2,1],[0,0,-10,20]]
RHS = [0.1,0,0,10]
print("\nSolve Ax = b system")
print("A matrix")
for a in A:
    print(a)
print("b vector")
for i in RHS:
    print(i)
print("resulting x vector")
print(Solve(A,RHS))