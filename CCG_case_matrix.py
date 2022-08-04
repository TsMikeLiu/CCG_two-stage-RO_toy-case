from unicodedata import name
import gurobipy as gp
from gurobipy import GRB
import numpy as np
import scipy.sparse as sp

# Constant Setting
f = [400, 414, 326]
a = [18, 25, 20]
C = [[22, 33, 24],
    [33, 23, 30],
    [20, 25, 27]]
# D = [206+40, 274+40, 220+40]
dl = [206, 274, 220]
du = [40, 40, 40]
k = 0 # count iteration
############################# MASTER PROBLEM START ###########################################
MP = gp.Model()
x = MP.addVars(3,3, lb=0, vtype=GRB.CONTINUOUS, name='x_0')
y = MP.addVars(len(f), lb=0, ub=1, obj=f, vtype=GRB.BINARY, name='y')
z = MP.addVars(len(a), lb=0, obj=a, vtype=GRB.CONTINUOUS, name='z')
g = MP.addVars(3, lb=0, ub=1, name='g')
d = MP.addVars(3, lb=0, name='d')
eta = MP.addVar(obj=1, lb=0, name='η')

MP_Cons_1 = MP.addConstrs((z[i] <= 800*y[i] for i in range(3)), name='MP_Cons_1')
MP_Cons_2 = MP.addConstr((gp.quicksum(z[i] for i in range(3)) >= 772), name='MP_Cons_2')

# iteration constraints
MP_Cons_3 = MP.addConstrs((gp.quicksum(x[i,j] for j in range(3)) <= z[i] for i in range(3)), name='MP_Cons_3')
MP_Cons_4 = MP.addConstrs((gp.quicksum(x[i, j] for i in range(3)) >= d[j] for j in range(3)), name='MP_Cons_4')
MP_Cons_eta = MP.addConstr(eta >= gp.quicksum(x[i,j]*C[i][j] for i in range(3) for j in range(3)), name='MP_Cons_eta')


# Master-problem uncertainty constraints
MP_Cons_uncertainty_1 = MP.addConstrs((d[i] == dl[i]+du[i]*g[i] for i in range(3)), name='MP_Uncertainty_Cons1')
MP_Cons_uncertainty_2 = MP.addConstr((gp.quicksum(g[i] for i in range(3)) <= 1.8), name='MP_Uncertainty_Cons2')
MP_Cons_uncertainty_3 = MP.addConstr((gp.quicksum(g[i] for i in range(2)) <= 1.2), name='MP_Uncertainty_Cons3')

MP.optimize()

LB = MP.objval


SP = gp.Model()
x_sub = SP.addVars(1, 9, lb=0, vtype=GRB.CONTINUOUS, name='x_sub')
d_sub = SP.addVars(3, 1, lb=0, name='d_sub')
g_sub = SP.addVars(3, lb=0, ub=1, name='g_sub')
pi = SP.addVars(6, lb=0, vtype=GRB.CONTINUOUS, name='pi')
v = SP.addVars(6, vtype=GRB.BINARY, name='v')
w = SP.addVars(3, 3, vtype=GRB.BINARY, name='w')
M = 10000

C_array = np.array(C)
G = np.array([[-1,-1,-1,0,0,0,0,0,0],
                [0,0,0,-1,-1,-1,0,0,0],
                [0,0,0,0,0,0,-1,-1,-1],
                [1,0,0,1,0,0,1,0,0],
                [0,1,0,0,1,0,0,1,0],
                [0,0,1,0,0,1,0,0,1]])
GC = np.array([])
GC = np.vstack((GC, -z[i].x )for i in range(3))
GC = np.vstack(GC, d_sub)

print(GC.shape)

SP.addConstr(G @ x_sub >= GC)