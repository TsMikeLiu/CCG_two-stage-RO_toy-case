import gurobipy as grb

Prob = grb.Model()

Prob.addVar(lb=-3, ub=0, obj=1)
Prob.optimize()
a = Prob.objval