from gurobipy import *

# arc travel cost scenarios
c = [[4, 2, 4, 2],
     [2, 4, 2, 4]]
l = len(c)

# arc post-interdiction cost increase
d = [100, 100, 100, 100]

# resource to interdict each arc
r = [1, 1, 1, 1]

# number of nodes and arcs
n = 4
m = 4

# total resource budget
r0 = 1

A = [[0, 1],
     [0, 2],
     [1, 3],
     [2, 3]]

model = Model()

# variable pi[i] shortest post-interdiction path to node i
pi = {}
for i in range(n):
    if (i == 0):
        pi[i] = model.addVar(vtype = GRB.CONTINUOUS, lb = 0, obj = -1) 
    elif (i == n-1):
        pi[i] = model.addVar(vtype = GRB.CONTINUOUS, lb = 0, obj = 1)
    else:
        pi[i] = model.addVar(vtype = GRB.CONTINUOUS, lb = 0, obj = 0)

# variable lmbda[q] - convex combination multiplier for scenario q
lmbda = {}
for q in range(l):
    lmbda[q] = model.addVar(vtype = GRB.CONTINUOUS, lb = 0, obj = 0)

# variable x[a] - attacker interdiction decision on arc a
x = {}
for a in range(m): 
    x[a] = model.addVar(vtype = GRB.BINARY, obj = 0)

# convex combination constraint for scenarios
model.addConstr(quicksum(lmbda[q] for q in range(l)) == 1)

# interdiction resource constraint
model.addConstr(quicksum(r[a]*x[a] for a in range(m)) <= r0)

# shortest path arc cost constraints
for a in range(m):
    i = A[a][0]
    j = A[a][1]
    model.addConstr(pi[j] - pi[i] - quicksum(((c[q][a] + d[a]*x[a])*lmbda[q]) for q in range(l)) <= 0)

# fix pi[0] to 0
model.addConstr(pi[0] == 0)
model.ModelSense = -1
model.update()
model.optimize()

# print output
for a in range(m):
    print("x_{}: {}".format(a, x[a].x))
print("Objective value: {}".format(pi[n-1].x))


        
        