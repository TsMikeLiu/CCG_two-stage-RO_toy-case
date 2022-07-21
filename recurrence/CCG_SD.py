from Sub_prob_SD import *
while UB-LB >10e-4:
    # update varibles and add constraints while iterating
    x_ij = MP.addVars(3,3, lb=0, vtype=GRB.CONTINUOUS, name='x_ij')
    d_op = d.x
    MP_Cons_3 = MP.addConstr((eta >= quicksum(x_ij[i,j] for i in range(3)) for j in range(3)), name='MP_Cons_3')
    MP_Cons_4 = MP.addConstrs((quicksum(x_ij[i,j] for j in range(3)) <= z[i] for i in range(3)), name='MP_Cons_4')
    MP_Cons_5 = MP.addConstrs((quicksum(x_ij[i,j] for i in range(3)) >= d_op[j] for j in range(3)), name='MP_Cons_5')

    MP.optimize()
    LB = MP.obj.values()

    obj = quicksum(p[i]*z[i].x for i in range(3)) - quicksum(t[j]*d[j] for j in range(3))
    SDSP.setObjective(obj, GRB.MAXIMIZE)
    SDSP.optimize()
    UB = min(UB, LB - eta.x + SDSP.objval)  # update UB
    k = k + 1  # Iterative counting

print(LB)
print(UB)
