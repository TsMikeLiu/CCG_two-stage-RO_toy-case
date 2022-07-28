from unicodedata import name
import gurobipy as gp
from gurobipy import GRB
import numpy as np
import scipy.sparse as sp

bbb = 2
testModel = gp.Model()
aaa = testModel.addVar(1,obj=1,vtype=GRB.CONTINUOUS)
testModel.addConstr(aaa>=bbb)
testModel.optimize()
print(aaa.x)
bbb = 1
testModel.remove(testModel.getConstrs()[0:1])
testModel.addConstr(aaa>=bbb)
testModel.optimize()
print(aaa.x)
print(bbb)