from unicodedata import name
from gurobipy import Model, quicksum, GRB
import numpy as np

# Constant Setting
f = [400, 414, 326]
a = [18, 25, 20]
C = [[22, 33, 24],
    [33, 23, 30],
    [20, 25, 27]]
D = [206+40, 274+40, 220+40]
dl = [206, 274, 220]
du = [40, 40,40]
k = 0 # count iteration

# Model Creating
MP = Model()    # Master Problem
SP = Model()    # Sub-Problem -> Using KKT to Solve
# SDSP = Model()  # Strong duality of Sub-Problem

# Construction of Master Problem 
# add Variables
y = MP.addVars(len(f), lb=0, ub=1, obj=f, vtype=GRB.INTEGER, name='y')
z = MP.addVars(len(a), lb=0, obj=a, vtype=GRB.CONTINUOUS, name='z')
# x_ij = MP.addVars(3,3, lb=0, vtype=GRB.CONTINUOUS, name='x_ij')
# g = MP.addVars(3, lb=0, ub=1, name='g')
eta = MP.addVar(obj=1, lb=0, name='η')

# add Constraints
MP_Cons_1 = MP.addConstrs((z[i] <= 800 * y[i] for i in range(3)), name='MP_Cons_1')
MP_Cons_2 = MP.addConstr(quicksum(z[i] for i in range(3)) >= 772, name='z') # To make sure SP is solveable 772 = 206_274+220+40*1.8
# MP_Cons_3 = MP.addConstrs((eta >= quicksum(x_ij[i,j] for i in range(3)) for j in range(3)), name='MP_Cons_3')
# MP_Cons_4 = MP.addConstrs((quicksum(x_ij[i,j] for j in range(3)) <= z[i] for i in range(3)), name='MP_Cons_4')

# Cons_3 = MP.addConstrs(quicksum(g[i] for i in range(3) <= 1.8), name='MP_Cons_3')
# Cons_4 = MP.addConstrs(quicksum(g[i] for i in range(2) <= 1.2), name='MP_Cons_4')

# Solve Model
MP.optimize()
# get optimal value
LB = MP.objval
