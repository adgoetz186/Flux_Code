import random
import docplex.mp.model as cpx
n = 10
m = 5
set_I = range(1, n+1)
set_J = range(1, m+1)
c = {(i,j): random.normalvariate(0,1) for i in set_I for j in set_J}
a = {(i,j): random.normalvariate(0,5) for i in set_I for j in set_J}
l = {(i,j): random.random()*10 for i in set_I for j in set_J}
u = {(i,j): random.random()*10+10 for i in set_I for j in set_J}
b = {j: random.random()*30 for j in set_J}
opt_model = cpx.Model(name="Model")
x_vars  = {(i,j): opt_model.continuous_var(lb=l[i,j], ub= u[i,j], name="x_{0}_{1}".format(i,j)) for i in set_I for j in set_J}
#constraints = {j : opt_model.add_constraint(ct=opt_model.sum(a[i,j] * x_vars[i,j] for i in set_I) == b[j],ctname="constraint_{0}".format(j))for j in set_J}
objective = opt_model.sum(x_vars[i,j]**2 for i in set_I for j in set_J)
opt_model.minimize(objective)
print(opt_model.solve())