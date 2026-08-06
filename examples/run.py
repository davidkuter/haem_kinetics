from haem_kinetics.models.model1 import Model1
from haem_kinetics.models.model2 import Model2
from haem_kinetics.models.model3 import Model3
from haem_kinetics.models.model4 import Model4
from haem_kinetics.models.model5 import Model5
from haem_kinetics.models.model6 import Model6
from haem_kinetics.models.model7 import Model7
from haem_kinetics.models.degradation import Degradation


t_start = 0   # min (parasite age offset 16 h applied inside models)
t_end = 1700  # min
t_step = 20   # min

# ---------------------------------------------------------------------------
# Model 5 — concentrations in M (DV basis); 5 states
# [Hb_DV, Fe2, Fe3_aq, Fe3_lip, Hz]
# ---------------------------------------------------------------------------
init_m5 = [0.018, 0.0, 0.0, 0.0, 0.36]
model = Model5()
model.run(
    t=[t_start, t_end],
    init=init_m5,
    t_eval=range(t_start, t_end, t_step),
    plot='model5.png',
)

# ---------------------------------------------------------------------------
# Model 6 / 7 — amounts in fg Fe/cell; 6 states
# [Hb_DV, Fe2, Fe3_aq, Fe3_lip, Hz, Hb_host]
# Host Fe is adjusted so total ≈ budget (~106 fg) inside run() if needed.
# ---------------------------------------------------------------------------
# init_fg = [1.0, 0.0, 0.0, 0.0, 23.0, 82.0]
# model = Model6()
# model = Model7()
# model.run(
#     t=[t_start, t_end],
#     init=init_fg,
#     t_eval=range(t_start, t_end, t_step),
#     plot='model6.png',
# )

# Legacy:
# model = Model3()
# model.run(t=[t_start, t_end], init=[0.018, 0.0, 0.0, 0.36],
#           t_eval=range(t_start, t_end, t_step), plot='test.png')

print(model.concentrations)
