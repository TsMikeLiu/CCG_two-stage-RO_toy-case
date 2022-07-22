from Master_prob import *

''' ADD Variables '''
# origin
x = SP.addVars(3, 3, lb=0, obj=np.array(C), vtype=GRB.CONTINUOUS, name='x')
g = SP.addVars(3, lb=0, ub=1, name='g')
d = [206+40*g[0], 274+40*g[1], 220+40*g[2]]

# dual
pi = SP.addVars(3, ub=0, vtype=GRB.CONTINUOUS, name='pi')
theta = SP.addVars(3, ub=0, vtype=GRB.CONTINUOUS, name='theta')

# binary of big M
v = SP.addVars(3, vtype=GRB.BINARY, name='v')
w = SP.addVars(3, vtype=GRB.BINARY, name='w')
h = SP.addVars(3, 3, vtype=GRB.BINARY, name='h')
M = 1e6

''' ADD Constraints '''
# origin
SP_C1 = SP.addConstrs(quicksum(x[i,j] for j in range(3)) <= z[i].x for i in range(3))
SP_C2 = SP.addConstrs(quicksum(x[i,j] for i in range(3)) >= d[j] for j in range(3))

# dual
SP_C3 = SP.addConstrs(pi[i]-theta[j] <= C[i][j] for i in range(3) for j in range(3))

# big M
SP_C4 = SP.addConstrs(-pi[i] <= M*v[i] for i in range(3))
SP_C5 = SP.addConstrs(z[i].x - quicksum(x[i,j] for j in range(3)) <= M*(1-v[i]) for i in range(3))

SP_C6 = SP.addConstrs(-theta[j] <= M*w[j] for j in range(3))
SP_C7 = SP.addConstrs(quicksum(x[i,j] for i in range(3))-d[j] <= M*(1-w[j]) for j in range(3))

SP_C8 = SP.addConstrs(x[i,j] <= M*h[i,j] for i in range(3) for j in range(3)) 
SP_C9 = SP.addConstrs(C[i][j]-pi[i]+theta[j] <= M*(1-h[i,j]) for i in range(3) for j in range(3))

# uncertainty
SP_C10 = SP.addConstr(quicksum(g[i] for i in range(2)) <= 1.2)
SP_C11 = SP.addConstr(quicksum(g[i] for i in range(3)) <= 1.8)

SP.write("SP_KKT.lp")
SP.setObjective(GRB.MAXIMIZE)
SP.optimize()
Q = SP.objval
UB = LB - eta.x - Q

print("****************************************************************")
print("*******************      FINISH     ****************************")
print("****************************************************************")

print("UB is %d" % UB)
print("LB is %d" % LB)