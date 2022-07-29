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
############################# MASTER PROBLEM END ###########################################


############################# SUB PROBLEM START ###########################################

SP = gp.Model()
x_sub = SP.addVars(3, 3, lb=0, vtype=GRB.CONTINUOUS, name='x_sub')
d_sub = SP.addVars(3, lb=0, name='d_sub')
g_sub = SP.addVars(3, lb=0, ub=1, name='g_sub')
pi = SP.addVars(6, lb=0, vtype=GRB.CONTINUOUS, name='pi')
v = SP.addVars(6, vtype=GRB.BINARY, name='v')
w = SP.addVars(3, 3, vtype=GRB.BINARY, name='w')
M = 10000

# Constraints
SP_Cons_1 = SP.addConstrs((gp.quicksum(x_sub[i,j] for j in range(3)) <= z[i].x for i in range(3)), name='SP_Cons_1') #[0:3]
SP_Cons_2 = SP.addConstrs((gp.quicksum(x_sub[i,j] for i in range(3)) >= d_sub[j] for j in range(3)), name='SP_Cons_2') #[3:6]
SP_Cons_3 = SP.addConstrs((-pi[i]+pi[j+3] <= C[i][j] for i in range(3) for j in range(3)), name='SP_Cons_3') #[6:15]

# slack constraints part 1
SP_SLACK_CONS_1 = SP.addConstrs((z[i].x-gp.quicksum(x_sub[i,j] for j in range(3)) <= M*(1-v[i]) for i in range(3)), name='SP_SLACK_CONS_1') #[15:18]
SP_SLACK_CONS_2 = SP.addConstrs((gp.quicksum(x_sub[i,j] for i in range(3))-d_sub[j] <= M*(1-v[j+3]) for j in range(3)), name='SP_SLACK_CONS_2') #[18:21]
SP_SLACK_CONS_3 = SP.addConstrs((pi[i] <= M*v[i] for i in range(6)), name='SP_SLACK_CONS_3') # [21:27]

# slack constraints part 2
SP_SLACK_CONS_4 = SP.addConstrs((C[i][j]+pi[i]-pi[j+3] <= M*(1-w[i,j]) for i in range(3) for j in range(3)), name='SP_SLACK_CONS_4') # [27:36]
SP_SLACK_CONS_5 = SP.addConstrs((x_sub[i,j] <= M*w[i,j] for i in range(3) for j in range(3)), name='SP_SLACK_CONS_5') # [36:45]

# uncertainty
SP_Cons_uncertainty_1 = SP.addConstrs((d_sub[i] == dl[i]+du[i]*g_sub[i] for i in range(3)), name='SP_Uncertainty_Cons1') #[45:48]
SP_Cons_uncertainty_2 = SP.addConstr((gp.quicksum(g_sub[i] for i in range(3)) <= 1.8), name='MP_Uncertainty_Cons2') #[48]
SP_Cons_uncertainty_3 = SP.addConstr((gp.quicksum(g_sub[i] for i in range(2)) <= 1.2), name='MP_Uncertainty_Cons3') #[49]


sub_obj = gp.quicksum(C[i][j]*x_sub[i,j] for i in range(3) for j in range(3))
SP.setObjective(sub_obj, GRB.MAXIMIZE)
SP.optimize()
SP_objval = SP.objval
UB = LB - eta.x + SP_objval
print(UB)
############################# SUB PROBLEM END ###########################################


############################# CCG START ###########################################
while np.abs(UB-LB)>1e-5 :
    k = k+1
    # Master-problem
    x_new = MP.addVars(3,3, lb=0, vtype=GRB.CONTINUOUS, name='x_{0}'.format(k))
    MP_Cons_3_new = MP.addConstrs((gp.quicksum(x_new[i,j] for j in range(3)) <= z[i] for i in range(3)), name='MP_Cons_3_{}'.format(k))
    MP_Cons_4_new = MP.addConstrs((gp.quicksum(x_new[i, j] for i in range(3)) >= d_sub[j].x for j in range(3)), name='MP_Cons_4_{}'.format(k))
    MP_Cons_eta = MP.addConstr(eta >= gp.quicksum(x_new[i,j]*C[i][j] for i in range(3) for j in range(3)), name='MP_Cons_eta_{}'.format(k))
    MP.optimize()
    LB = max(LB, MP.objval)
    print(LB)

    # Sub-problem update
    # delete old constraints which related to z
    SP.remove(SP.getConstrs()[0:3])
    SP.remove(SP.getConstrs()[15:18])
    # add new constraints which related to z
    SP_Cons_1 = SP.addConstrs((gp.quicksum(x_sub[i,j] for j in range(3)) <= z[i].x for i in range(3)), name='SP_Cons_1') #[0:3]
    SP_SLACK_CONS_1 = SP.addConstrs((z[i].x-gp.quicksum(x_sub[i,j] for j in range(3)) <= M*(1-v[i]) for i in range(3)), name='SP_SLACK_CONS_1') #[15:18]
    
    SP.optimize()
    UB = LB - eta.x + SP.objval
    print(UB)
############################# CCG END ###########################################

# Some informations
print("Iteration finished! We found the optimal solution!")
print("Final Objective:{0}".format(LB))
print(y[0],y[1],y[2])
print(z[0],z[1],z[2])
for i in range(3):
    for j in range(3):
        print(x_new[i,j])