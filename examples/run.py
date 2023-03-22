from scipy.integrate import solve_ivp

from haem_kinetics.models.model1 import Model1

# def integrate(t, init):
#     model._set_initial_conc(init=init)
#     return model.differential_eqs
#
#
# model = Model1()
#
# solution = solve_ivp(integrate, [0, 1000], [0.0, 0.0, 0.0, 0.0])

model = Model1()
print(dir(model))
print(model.initial_values)
model.run(t=[0, 1000], init=[0.0, 0.0, 0.0, 0.0])

print(model.solution)
print(model.time)
print(model.concentrations)
