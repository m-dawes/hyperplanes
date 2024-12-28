import itertools
import time

# L is space of quadratic forms in 4 variables over R

########################################################################################## 
#                                                                                        #
#                                                                                        #
#                                    Main function                                       #
#                                                                                        #
#                                                                                        #
########################################################################################## 


def all_invs():
    '''
    Output:  a list of symplectic involutions fixing interior elements  
             of cones in cells_and_boundary list 
    '''

    cells_and_boundary=[
        [dt6m, dt6mb, "Desargues Transverse (rank 6)"],
        [rt6m, rt6mb, "Reye Transverse (rank 6)"],
        [drt5m, drt5mb, "DR Transverse (rank 5)"],
        [rrt5m, rrt5mb,  "RR Transverse (rank 5)"],
        [hnt4m, hnt4mb,  "Hexagon Non-transverse (rank 4)"],
        [st4m, st4mb,    "Square Transverse (rank 4)"],
        [tt4m, tt4mb,    "Triangle Transverse (rank 4)"],
        [vnt3m, vnt3mb,  "Vertebra Non-transverse (rank 3)"],
        [cnt3m, cnt3mb,  "Crystal Non-transverse (rank 3)"],
        [pt3m, pt3mb,    "Pyramid transverse (rank 3)"],
        [rsnt2m, rsnt2mb, "Red Square Non-Transverse (rank 2)"],
        [ltb2m, ltb2mb,  "Lagrangian Triple Boundary (rank 2)"],
        [lpb1m, lpb1mb, "Lagrangian Pair Boundary (rank 1)"],
        [pb0m, pb0mb, "Point Boundary (rank 0)"]
    ]


    out=[]
    for k in cells_and_boundary:
        print("\n\n\n",k[2],"\n\n\n")
        i=k[0]
        x=find_invs(k)
        for j in x:
            if not j in out:
                out.append(j)
    return out



########################################################################################## 
#                                                                                        #
#                                                                                        #
# Cells, automorphisms and fixed loci calculations                                       #
#                                                                                        #
#                                                                                        #
########################################################################################## 



def vv(a,b,c,d):
    '''
    Input:  - a tuple (a,b,c,d)

    Output: - a vector coefficients of expanded polynomial (aw+bx+cy+dz)^2 on ordered basis
              w^2, x^2, y^2, z^2, wx, wy, wz, xy, xz, yz
    '''

# a^2*w^2 + 2*a*b*w*x + b^2*x^2 + 2*a*c*w*y + 2*b*c*x*y + c^2*y^2
#         + 2*a*d*w*z + 2*b*d*x*z + 2*c*d*y*z + d^2*z^2

    v=[]               # the output vector
    v.append(a*a)      # w^2
    v.append(b*b)      # x^2
    v.append(c*c)      # y^2
    v.append(d*d)      # z^2
    v.append(2*a*b)    # wx
    v.append(2*a*c)    # wy
    v.append(2*a*d)    # wz
    v.append(2*b*c)    # xy
    v.append(2*b*d)    # xz
    v.append(2*c*d)    # yz

    return v           


def matrix_to_cone(m):
    '''
    Input:   - a 4xn matrix in MacPherson and McConnell

    Output   - The associated cone in L, as a Sage polyhedron
    '''
    l=[]
    for i in m.columns():
        l.append(tuple(vv(i[0],i[1],i[2],i[3])))
    return Polyhedron(rays=l)


def maps_to_L(g):
    '''
    Input:   - A symmetric matrix g (in L)

    Output:  - A list, corresponding to a vector, expressing g on ordered basis
                w^2,  x^2, y^2, z^2, (1/2)* { wx, wy, wz, xy, xz, yz }.
    '''
    mww=matrix([[1,0,0,0,0,0,0,0,0,0]])
    mxx=matrix([[0,1,0,0,0,0,0,0,0,0]])
    myy=matrix([[0,0,1,0,0,0,0,0,0,0]])
    mzz=matrix([[0,0,0,1,0,0,0,0,0,0]])
    mwx=matrix([[0,0,0,0,1,0,0,0,0,0]])
    mwy=matrix([[0,0,0,0,0,1,0,0,0,0]])
    mwz=matrix([[0,0,0,0,0,0,1,0,0,0]])
    mxy=matrix([[0,0,0,0,0,0,0,1,0,0]])
    mxz=matrix([[0,0,0,0,0,0,0,0,1,0]])
    myz=matrix([[0,0,0,0,0,0,0,0,0,1]])
        
    v=0*mww
    v=v+g[0,0]*mww
    v=v+g[0,1]*mwx
    v=v+g[0,2]*mwy
    v=v+g[0,3]*mwz
    v=v+g[1,1]*mxx
    v=v+g[1,2]*mxy
    v=v+g[1,3]*mxz
    v=v+g[2,2]*myy
    v=v+g[2,3]*myz
    v=v+g[3,3]*mzz

    return list(v[0])

def g_on_L(g):
    '''
    Input:  - a symplectic transformation g

    Output: - the action of g on L
    '''

    # A list of basis elements for L
    gww=matrix([[1,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0]])  # ww
    gxx=matrix([[0,0,0,0],[0,1,0,0],[0,0,0,0],[0,0,0,0]])  # xx
    gyy=matrix([[0,0,0,0],[0,0,0,0],[0,0,1,0],[0,0,0,0]])  # yy
    gzz=matrix([[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,1]])  # zz
    gwx=matrix([[0,1,0,0],[1,0,0,0],[0,0,0,0],[0,0,0,0]])  # wx
    gwy=matrix([[0,0,1,0],[0,0,0,0],[1,0,0,0],[0,0,0,0]])  # wy
    gwz=matrix([[0,0,0,1],[0,0,0,0],[0,0,0,0],[1,0,0,0]])  # wz
    gxy=matrix([[0,0,0,0],[0,0,1,0],[0,1,0,0],[0,0,0,0]])  # xy
    gxz=matrix([[0,0,0,0],[0,0,0,1],[0,0,0,0],[0,1,0,0]])  # xz
    gyz=matrix([[0,0,0,0],[0,0,0,0],[0,0,0,1],[0,0,1,0]])  # yz
    
    # We calculate the action of g on the above basis, storing result in l
    # the elements of l refer to vectors on ordered basis corresponding to
    # w^2, x^2, y^2, z^2, (1/2)*{ wx, wy, wz, xy, xz, yz }
    
    l=[]                      
    x=g.transpose()*gww*g   # calculate action of g on w^2
    l.append(maps_to_L(x))  # store result in l 
    x=g.transpose()*gxx*g
    l.append(maps_to_L(x))
    x=g.transpose()*gyy*g
    l.append(maps_to_L(x))
    x=g.transpose()*gzz*g
    l.append(maps_to_L(x))
    x=g.transpose()*gwx*g
    l.append(maps_to_L(x))
    x=g.transpose()*gwy*g
    l.append(maps_to_L(x))
    x=g.transpose()*gwz*g
    l.append(maps_to_L(x))
    x=g.transpose()*gxy*g
    l.append(maps_to_L(x))
    x=g.transpose()*gxz*g
    l.append(maps_to_L(x))
    x=g.transpose()*gyz*g
    l.append(maps_to_L(x))
    return matrix(l).transpose()  # convert from rows to columns

def one_eigenbasis(g):
    '''
    Input:  - A linear involution acting on L 

    Output: - A basis for the 1-eigenspace of g
    '''
    E=g.eigenspaces_left()
    for v in E:
        if v[0]==1:
            return v[1].basis()

def fixed_cone(g):
    '''
    Input:   - a symplectic involution g

    Output:  - The fixed locus of g on L
    '''
    G=g_on_L(g)                  # calculate the action of g on L
    B=one_eigenbasis(G)          # calculate the associated 1-eigenbasis
    return Polyhedron(lines=B)   # return the fixed space as a polyhedron

def in_boundary(S, l):
    '''
    Input:   - S is a list of basis vectors for a subspace of L
             - l is a list of functionals (usually taken to define the
               boundary of a cone)
    
    Output:  - True if every element of S evaluates to zero on some
               element of l
             - False otherwise
    '''

    for i in l:                 # iterate over functionals
        v=matrix(i).transpose() # a functional  
        y=True                  # Suppose S is contained in Ker(i)
        for j in S:             # iterate over elements of S
            z=matrix(j)
            c=z*v               # evaluate the functional i on j
            if c[0,0]!=0:       # then some element of S does not evaluate to 0 on i
                y=False         # S is not contained in Ker(i)
                break
        if y==True:             # Then we are contained in some boundary component
            return True
    return False                # We got to the end without being contained in any boundary component


def fixed_in_boundary(g, C, f):
    '''
    Input:   - a symplectic involution g
             - a cone C
             - a list of functionals defining the boundary of C

    Output:  - True if the fixed locus of g lies in the boundary of C
             - False otherwise.
    '''
    H=fixed_cone(g)              # fixed locus H  of g
    Q=H.intersection(C)          # intersection C \cap H
    rays=Q.rays()                # rays of intersection
    return in_boundary(rays, f)  # evaluate rays on boundary functionals


def check_pm_list(P, l):
    '''
    Input:  - A matrix P
            - A list l of vectors in ZZ^4 describing extreme rays of cone C

    Output: - True if P*l == l up to changes of sign
            - False otherwise
    '''
    for i in l:                               # iterate over vectors in list
      u=matrix(i).transpose() 
      v=P*u                                   # transform i (as a row vector)
      x0=[v[0,0],v[1,0],v[2,0],v[3,0]]        # x0 = i   as a list
      x1=[-v[0,0],-v[1,0],-v[2,0],-v[3,0]]    # x1 = -i  as a list
      if (not x0 in l) and (not x1 in l):     # P does not preserve l if neither x0 nor x1 in l
          return False
    return True                               # otherwise, l is preserved

def mm_to_list(x):
    '''
    Input:  - A matrix x
    
    Output: - A list of the columns of x as row vectors
              (i.e. a list of the rows in the transpose of x).
    '''
    out=[]
    x=x.transpose()
    for r in x.rows():
        out.append(list(r))
    return out

def sign_change(s,x):
    '''
    Input:  - a 4-tuple s of +/- 1, corresponding to sign changes 
            - a 4x4 matrix x 
    
    Output: - the list [ [ s[i]*x[i][j] ] ] of i-th row of x * s[i]
    '''
    l=[]   
    for i in range(0,4):
        l.append([s[i]*x[i][0], s[i]*x[i][1], s[i]*x[i][2], s[i]*x[i][3]])
    return l

def find_invs(z):
    '''
    Input:   - A triple z = [cone C in 4xn matrix form, boundary functionals in L, "description of cone" ]
 
    Output:  - Text output to console of a list of symplectic involutions preserving the cone C such that the fixed locus of
    	       g on L is not contained in the boundary of C.
             - A list of all such involutions
    '''
    count=0                                                   # a count of cases checked
    Omega=matrix([[0,0,0,1],[0,0,1,0],[0,-1,0,0],[-1,0,0,0]]) # the symplectic form
    out=[]                                                    # the output list
    start=time.time()                                         # starting time
    l=mm_to_list(z[0])                                        # columns of matrix defining cone, as row vectors
    f=z[1]                                                    # the boundary functionals
    C=matrix_to_cone(z[0])                                    # the cone
    signs=[                                                   # a list of sign changes
        [1, 1, 1, 1],
        [1, 1, 1, -1],
        [1, 1, -1, 1],
        [1, 1, -1, -1],
        [1, -1, 1, 1],
        [1, -1, 1, -1],
        [1, -1, -1, 1],
        [1, -1, -1, -1],
        [-1, 1, 1, 1],
        [-1, 1, 1, -1],
        [-1, 1, -1, 1],
        [-1, 1, -1, -1],
        [-1, -1, 1, 1],
        [-1, -1, 1, -1],
        [-1, -1, -1, 1],
        [-1, -1, -1, -1]
    ]

    
    for i in itertools.combinations(l, 4):                   # iterate over 4-tuples of vectors defining extreme rays of C
        P1=matrix(i)
        if abs(P1.det())==1:                                 # we must start with a primitive 4-tuple
            for j0 in itertools.combinations(l,4):           # the vectors we are mapping to (as a set)
                Q=matrix(j0)
                if abs(Q.det())==1:                          # vectors must be a primitive 4-tuple
                    for k in itertools.permutations(j0):     # the vectors we are mapping to (as a tuple)
                        for s in signs:                      # allow sign changes
                            count=count+1 
                            j=sign_change(s, k)              # change signs of k
                            P2=matrix(j)
                            P=P2*P1.inverse()
                            if (P*P==matrix.identity(4) and
                                P.transpose()*Omega*P==Omega and
                                P!=matrix.identity(4) and
                                P!=-matrix.identity(4) ):     # we have a symplectic involution
                                if check_pm_list(P, l)==True: # P is an automorphism of C
                                    if not ((P in out) or (-P in out)):  # Neither P nor -P is in ouput list
                                                                         # (noting that they act identically on M).

                                        # is fixed locus of P contained in boundary of C?
                                        fb=fixed_in_boundary(P, C, f)
                                        # C contained in fixed locus of P?
                                        subfl=C.intersection(fixed_cone(P))==C
                                        print("Fixed in boundary is: ", fb)
                                        print("Subset of fixed locus is:", subfl)
                                        if  subfl:   # cone lies in Humbert divisor
                                            print("lies in fixed locus of involution")
                                            print(P)
                                            return []
                                        if not fb and not subfl:
                                            print("involution")
                                            print(P)
                                            print(P.charpoly())
                                            out.append(P)
                                            end=time.time()
                                            print("Elapsed time = ", end-start, "seconds")
    end=time.time()
    print("did ", count, " cases in ", end-start, "seconds. Outputting ", len(out), " involutions")
    return out


########################################################################################## 
#                                                                                        #
#                                                                                        #
#                  Cell data                                                             #
#                                                                                        #
#                                                                                        #
########################################################################################## 

dt6m=matrix([
    [1,0,0,0,1,1,1,0,1,1],
    [0,1,0,0,1,0,0,1,1,0],
    [0,0,1,0,0,1,0,0,1,1],
    [0,0,0,1,0,0,-1,1,0,-1]
    ])

rt6m=matrix([
    [1,0,0,0,1,1,0,0,1,1,1,0],
    [0,1,0,0,1,0,1,0,1,1,0,1],
    [0,0,1,0,0,-1,0,1,-1,0,-1,1],
    [0,0,0,1,0,0,-1,-1,0,-1,1,-1]
    ])

drt5m=matrix([
    [1,0,0,0,1,1,1,0,1],
    [0,1,0,0,1,0,0,1,0],
    [0,0,1,0,0,1,0,0,1],
    [0,0,0,1,0,0,-1,1,-1]
    ])

rrt5m=matrix([
    [1,0,0,0,1,0,0,1,1],
    [0,1,0,0,1,1,0,1,0],
    [0,0,1,0,0,0,1,-1,-1],
    [0,0,0,1,0,-1,-1,0,1]
    ])

hnt4m=matrix([
    [1,0,0,0,1,0],
    [0,1,0,0,0,1],
    [0,0,1,0,0,1],
    [0,0,0,1,1,0]
    ])

st4m=matrix([
    [1,0,0,0,1,1,0,1],
    [0,1,0,0,1,0,1,0],
    [0,0,1,0,0,0,0,1],
    [0,0,0,1,0,-1,1,-1]
    ])

tt4m=matrix([
    [1,0,0,0,1,1,0,1],
    [0,1,0,0,1,0,1,0],
    [0,0,1,0,0,1,0,1],
    [0,0,0,1,0,0,1,-1]
    ])

vnt3m = matrix([
    [1,0,0,0,1],
    [0,1,0,0,0],
    [0,0,1,0,0],
    [0,0,0,1,1]
    ])
cnt3m= matrix([
    [1,0,0,0,1,0],
    [0,1,0,0,0,1],
    [0,0,1,0,1,0],
    [0,0,0,1,0,1]
    ])
pt3m=matrix([
    [1,0,0,0,1,0,1],
    [0,1,0,0,1,1,0],
    [0,0,1,0,0,0,1],
    [0,0,0,1,0,1,-1]
    ])

rsnt2m=matrix([
    [1,0,0,0],
    [0,1,0,0],
    [0,0,1,0],
    [0,0,0,1]
    ])
ltb2m=matrix([
    [1,0,1],
    [0,1,1],
    [0,0,0],
    [0,0,0]
    ])

lpb1m=matrix([
    [1,0],
    [0,1],
    [0,0],
    [0,0]
])

pb0m=matrix([
    [1],
    [0],
    [0],
    [0]
    ])

########################################################################################## 
#                                                                                        #
#                                                                                        #
#    Boundary Functionals                                                                #
#                                                                                        #
########################################################################################## 



dt6mb=[
    (0, 0, 0, 0, 0, 0, 0, 1, 0, 0),
    (0, 0, 0, 0, 1, 0, 0, -1, 0, 0),
    (0, 0, 0, 0, 0, 0, 0, 0, 0, -1),
    (0, 0, 0, 0, 0, 1, 0, -1, 0, 1),
    (0, 0, 0, 0, 0, 0, -1, 0, 0, 1),
    (2, 0, 0, 0, -1, -1, 1, 1, 0, -1),
    (0, 0, 0, 0, 0, 0, 0, 0, 1, 0),
    (0, 2, 0, 0, -1, 0, 0, 0, -1, 0),
    (0, 0, 2, 0, 0, -1, 0, 0, 0, 0),
    (0, 0, 0, 2, 0, 0, 1, 0, -1, 0) 
]

rt6mb=[
    (0, 0, 0, 0, 1, 0, 1, 0, 0, 0),
    (0, 2, 0, 0, 0, 0, 0, 0, 1, 0),
    (0, 0, 2, 0, 1, 1, 1, 0, 0, 1),
    (0, 0, 0, 0, 0, -1, -1, 0, 0, 0),
    (0, 0, 0, 0, 0, 0, 0, -1, -1, 0),
    (0, 0, 0, 0, 1, 0, 0, 0, 0, 0),
    (0, 0, 0, 0, 0, 0, -1, -1, 0, -1),
    (0, 0, 0, 0, 0, 0, 0, -1, 0, -1),
    (0, 2, 0, 0, -1, 0, -1, -1, 1, -1),
    (0, 0, 2, 0, 0, 1, 0, -1, 0, 0),
    (0, 0, 2, 0, 0, 0, 0, 0, 0, 1),
    (0, 2, 0, 0, -1, -1, -1, 0, 1, 0),
    (0, 0, 0, 0, 0, -1, 0, 0, 0, 0),
    (0, 0, 2, 0, 0, 1, 1, -1, -1, 1),
    (0, 2, 0, 0, -1, 0, 0, -1, 0, 0),
    (0, 0, 0, 0, 0, 0, 1, -1, -1, 0),
    (2, 0, 0, 0, -1, 1, 1, -1, -1, 0),
    (0, 0, 0, 0, 0, 0, 1, 0, -1, 0),
    (0, 0, 0, 2, 0, 0, 1, 0, 0, 1),
    (0, 0, 0, 2, 1, 0, 1, 1, 1, 1),
    (0, 0, 0, 0, 1, 0, 1, 1, 0, 0),
    (2, 0, 0, 0, 0, 1, 1, 0, 0, 0),
    (0, 2, 0, 2, 0, 0, 0, 1, 2, 1),
    (0, 2, 0, 0, 0, 0, 0, 1, 1, 0),
    (2, 2, 0, 0, -1, 1, 0, 0, 1, 0),
    (0, 0, 2, 2, 1, 1, 1, 1, 1, 2),
    (0, 0, 2, 0, 1, 1, 1, 1, 0, 1),
    (2, 0, 2, 0, 0, 2, 1, 0, 0, 1),
    (0, 0, 0, 2, 0, -1, -1, 1, 1, 1),
    (0, 0, 0, 0, 0, -1, -1, 1, 0, 0),
    (2, 0, 0, 0, -1, 0, -1, 0, 0, 0),
    (0, 0, 0, 2, 0, 0, 0, 0, 0, 1),
    (0, 0, 0, 0, 0, 0, 0, 0, -1, 0),
    (2, 0, 0, 0, -1, 1, 0, -1, -1, 0),
    (0, 0, 0, 2, 1, 0, 0, 1, 1, 1),
    (0, 0, 0, 0, 1, 0, 0, 1, 0, 0),
    (2, 0, 0, 0, 0, 1, 0, 0, 0, 0),
    (0, 0, 0, 2, 0, 0, -1, 0, 1, 0),
    (0, 0, 0, 0, 0, 0, -1, 0, 0, -1),
    (2, 0, 0, 0, -1, 1, -1, -1, 0, -1),
    (0, 0, 0, 2, 0, 0, 0, 0, 1, 0),
    (0, 0, 0, 0, 0, 0, 0, 0, 0, -1),
    (2, 0, 0, 0, -1, 1, 0, -1, 0, -1),
    (0, 2, 0, 2, -1, 0, -1, 0, 2, 0),
    (0, 2, 0, 0, -1, 0, -1, 0, 1, -1),
    (2, 2, 0, 0, -2, 1, -1, -1, 1, -1),
    (0, 0, 2, 2, 0, 1, 0, 0, 1, 1),
    (0, 0, 2, 0, 0, 1, 0, 0, 0, 0),
    (2, 0, 2, 0, -1, 2, 0, -1, 0, 0),
    (0, 0, 2, 2, 0, 0, 0, 1, 1, 2),
    (0, 0, 2, 0, 0, 0, 0, 1, 0, 1),
    (2, 0, 2, 0, -1, 1, 0, 0, 0, 1),
    (0, 2, 0, 2, -1, -1, -1, 1, 2, 1),
    (0, 2, 0, 0, -1, -1, -1, 1, 1, 0),
    (2, 2, 0, 0, -2, 0, -1, 0, 1, 0),
    (0, 0, 0, 2, 0, -1, 0, 1, 1, 1),
    (0, 0, 0, 0, 0, -1, 0, 1, 0, 0),
    (2, 0, 0, 0, -1, 0, 0, 0, 0, 0),
    (0, 0, 2, 2, 0, 1, 1, 0, 0, 2),
    (0, 0, 2, 0, 0, 1, 1, 0, -1, 1),
    (2, 0, 2, 0, -1, 2, 1, -1, -1, 1),
    (0, 2, 0, 2, -1, 0, 0, 0, 1, 1),
    (0, 2, 0, 0, -1, 0, 0, 0, 0, 0),
    (2, 2, 0, 0, -2, 1, 0, -1, 0, 0)
]

drt5mb=[
    (0, 0, 0, 0, 1, 0, 0, 0, 0, 0),
    (0, 0, 0, 0, 0, 0, 0, 0, 0, -1),
    (0, 0, 0, 0, 0, 1, 0, 0, 0, 1),
    (0, 0, 0, 0, 0, 0, -1, 0, 0, 1),
    (2, 0, 0, 0, -1, -1, 1, 0, 0, -1),
    (0, 0, 0, 0, 0, 0, 0, 0, 1, 0),
    (0, 2, 0, 0, -1, 0, 0, 0, -1, 0),
    (0, 0, 2, 0, 0, -1, 0, 0, 0, 0),
    (0, 0, 0, 2, 0, 0, 1, 0, -1, 0)
]

rrt5mb=[
    (0, 0, 0, 0, 0, 0, 0, -1, 0, 0),
    (0, 0, 0, 0, 1, 0, 0, 1, 0, 0),
    (0, 0, 0, 0, 0, 0, 1, 0, 0, 0),
    (2, 0, 0, 0, -1, 0, -1, 0, 0, 0),
    (0, 0, 0, 0, 0, 0, 0, 0, -1, 0),
    (0, 2, 0, 0, -1, 0, 0, 0, 1, 0),
    (0, 0, 0, 0, 0, 0, -1, 0, 0, -1),
    (0, 0, 2, 0, 0, 0, 0, 1, 0, 1),
    (0, 0, 0, 2, 0, 0, 0, 0, 1, 1) 
]

hnt4mb=[
    (0, 0, 0, 0, 0, 0, 1, 0, 0, 0) ,
    (2, 0, 0, 0, 0, 0, -1, 0, 0, 0) ,
    (0, 0, 0, 0, 0, 0, 0, 1, 0, 0) ,
    (0, 2, 0, 0, 0, 0, 0, -1, 0, 0) ,
    (0, 0, 2, 0, 0, 0, 0, -1, 0, 0) ,
    (0, 0, 0, 2, 0, 0, -1, 0, 0, 0) 
]

st4mb=[
    (0, 0, 0, 0, 1, 0, 0, 0, 0, 0),
    (0, 0, 0, 0, 0, 0, 0, 0, 0, -1),
    (0, 0, 0, 0, 0, 0, -1, 0, 0, 1),
    (2, 0, 0, 0, -1, 0, 1, 0, 0, 0),
    (0, 0, 0, 0, 0, 0, 0, 0, 1, 0),
    (0, 2, 0, 0, -1, 0, 0, 0, -1, 0),
    (0, 0, 2, 0, 0, 0, 0, 0, 0, 1),
    (0, 0, 0, 2, 0, 0, 1, 0, -1, 0)
]

tt4mb=[
    (0, 0, 0, 0, 1, 0, 0, 0, 0, 0),
    (0, 0, 0, 0, 0, 0, 0, 0, 0, -1),
    (0, 0, 0, 0, 0, 1, 0, 0, 0, 1),
    (2, 0, 0, 0, -1, -1, 0, 0, 0, 0),
    (0, 0, 0, 0, 0, 0, 0, 0, 1, 0),
    (0, 2, 0, 0, -1, 0, 0, 0, -1, 0),
    (0, 0, 2, 0, 0, -1, 0, 0, 0, 0),
    (0, 0, 0, 2, 0, 0, 0, 0, -1, 1)
]

vnt3mb=[
    (0, 0, 0, 0, 0, 0, 1, 0, 0, 0),
    (2, 0, 0, 0, 0, 0, -1, 0, 0, 0),
    (0, 1, 0, 0, 0, 0, 0, 0, 0, 0),
    (0, 0, 1, 0, 0, 0, 0, 0, 0, 0),
    (0, 0, 0, 2, 0, 0, -1, 0, 0, 0)
]

cnt3mb=[
  (0, 0, 0, 0, 0, 1, 0, 0, 0, 0) ,
  (2, 0, 0, 0, 0, -1, 0, 0, 0, 0),
  (0, 0, 0, 0, 0, 0, 0, 0, 1, 0) ,
  (0, 2, 0, 0, 0, 0, 0, 0, -1, 0),
  (0, 0, 2, 0, 0, -1, 0, 0, 0, 0),
  (0, 0, 0, 2, 0, 0, 0, 0, -1, 0)
]

pt3mb=[
    (0, 0, 0, 0, 1, 0, 0, 0, 0, 0) ,
    (0, 0, 0, 0, 0, 0, 0, 0, 0, -1),
    (2, 0, 0, 0, -1, 0, 0, 0, 0, 1),
    (0, 0, 0, 0, 0, 0, 0, 0, 1, 0) ,
    (0, 2, 0, 0, -1, 0, 0, 0, -1, 0),
    (0, 0, 2, 0, 0, 0, 0, 0, 0, 1),
    (0, 0, 0, 2, 0, 0, 0, 0, -1, 1)
]

rsnt2mb=[
    (1, 0, 0, 0, 0, 0, 0, 0, 0, 0),
    (0, 1, 0, 0, 0, 0, 0, 0, 0, 0),
    (0, 0, 1, 0, 0, 0, 0, 0, 0, 0),
    (0, 0, 0, 1, 0, 0, 0, 0, 0, 0) 
]

ltb2mb=[
    (0, 0, 0, 0, 1, 0, 0, 0, 0, 0),
    (2, 0, 0, 0, -1, 0, 0, 0, 0, 0),
    (0, 2, 0, 0, -1, 0, 0, 0, 0, 0)
]

lpb1mb=[
    (1, 0, 0, 0, 0, 0, 0, 0, 0, 0),
    (0, 1, 0, 0, 0, 0, 0, 0, 0, 0)
]

pb0mb=[
    (1, 0, 0, 0, 0, 0, 0, 0, 0, 0)
    ]
