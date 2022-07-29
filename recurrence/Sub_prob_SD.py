from Master_prob import *

# Constant Setting
f = [400, 414, 326]
a = [18, 25, 20]
C = [[22, 33, 24],
    [33, 23, 30],
    [20, 25, 27]]
D = [206+40, 274+40, 220+40]
dl = [206, 274, 220]
du = [40, 40,40]


SDSP = Model()  # Strong duality of Sub-Problem

# add Variables
p = SDSP.addVars(3, ub=0, vtype=GRB.CONTINUOUS, name='pi')
t = SDSP.addVars(3, ub=0, vtype=GRB.CONTINUOUS, name='theta')
# d = SDSP.addVars(3, lb=0, vtype=GRB.CONTINUOUS, name='d')
g_new = SDSP.addVars(3, lb=0, ub=1, vtype=GRB.CONTINUOUS, name='g')


SDSP_Cons1 = SDSP.addConstrs((p[i]-t[j]<=C[i][j] for i in range(3) for j in range(3)), name='SDSP_cons1')
# SDSP_Cons2 = SDSP.addConstrs((d[j]==dl[j]+g[j]*du[j] for j in range(3)), name='SDSP_cons2')
SDSP_Cons3 = SDSP.addConstr(quicksum(g_new[i] for i in range(3)) <= 1.8, name='SDSP_cons3')
SDSP_Cons4 = SDSP.addConstr(quicksum(g_new[i] for i in range(2)) <= 1.2, name='SDSP_cons4')

# obj = quicksum(-p[i]*z[i].x for i in range(3)) + quicksum(t[j]*d[j] for j in range(3))
d = [dl[i] + du[i] * g_new[i] for i in range(3)]
obj = quicksum(-t[i]*d[i]+p[i]*z[i].x for i in range(3))

SDSP.setObjective(obj, GRB.MAXIMIZE)
SDSP.optimize()
UB = LB - eta.x + SDSP.objval
print("Information:")
print(UB)
for i in range(3):
    print(z[i].x)

