"""Run active haem kinetics models and save comparison plots."""
from haem_kinetics.models.model1 import Model1
from haem_kinetics.models.model2 import Model2
from haem_kinetics.models.model3 import Model3
from haem_kinetics.models.degradation import Degradation


t_start = 0
t_end = 1700
t_step = 20
t_eval = range(t_start, t_end, t_step)
t_span = [t_start, t_end]

runs = [
    ('examples/model1.png', Model1, [0.018, 0.0, 0.0, 0.36], {}),
    ('examples/model2.png', Model2, [0.018, 0.0, 0.0, 0.36], {}),
    ('examples/model3.png', Model3, [0.018, 0.0, 0.0, 0.36], {}),
    ('examples/degradation.png', Degradation, [0.018, 0.0], {}),
]

for plot_name, cls, init, extra in runs:
    model = cls()
    kwargs = dict(t_eval=t_eval, plot=plot_name, **extra)
    print(f'Running {cls.__name__} -> {plot_name} ...')
    model.run(t=t_span, init=init, **kwargs)
    host = model.concentrations['conc_hb_rbc']
    dv_cols = [c for c in model.concentrations.columns
               if c.startswith('conc_') and c != 'conc_hb_rbc'
               and c not in ('conc_hb_dv_obs',)]
    tot = float(model.concentrations[dv_cols].iloc[-1].sum() + host.iloc[-1])
    hz = float(model.concentrations['conc_hz'].iloc[-1]) if 'conc_hz' in model.concentrations else 0.0
    print(f'  end total Fe={tot:.2f} fg  Hz={hz:.2f} fg  host={float(host.iloc[-1]):.2f} fg')

print('Done.')
