from haem_kinetics.models.model1 import Model1

t_start = 0
t_end = 1760
t_step = 20

model = Model1()
model.run(t=[t_start, t_end], init=[0.0, 0.0, 0.0, 0.0], t_eval=range(t_start, t_end, t_step), plot='test.png')

print(model.initial_values)
print(model.differential_eqs)
print(model.time)
print(model.concentrations)

