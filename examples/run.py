from haem_kinetics.models.model1 import Model1
from haem_kinetics.models.model2 import Model2
from haem_kinetics.models.model3 import Model3
from haem_kinetics.models.model4 import Model4
from haem_kinetics.models.degradation import Degradation


t_start = 0   # min
t_end = 1700  # min
t_step = 20   # min
init = [0.018,  # Hb-DV: conc in DV, corresponds to 2E-4 M in cell
        0.0,    # Fe2PP
        0.0,    # Fe3PP
        0.0]    # Hz

model = Model3()
model.run(t=[t_start, t_end], init=init[:4], t_eval=range(t_start, t_end, t_step), plot='test.png')

# model = Degradation()
# model.run(t=[t_start, t_end], init=init[:2], t_eval=range(t_start, t_end, t_step), plot='test.png')

# print(model.initial_values)
# print(model.differential_eqs)
# print(model.time)
print(model.concentrations)
# model.concentrations.to_csv('concentrations.tsv', sep='\t')
