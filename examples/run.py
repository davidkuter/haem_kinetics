"""Run active haem kinetics models and save comparison plots."""
from haem_kinetics.models.model1 import Model1
from haem_kinetics.models.model2 import Model2a
from haem_kinetics.models.model2b import Model2b
from haem_kinetics.models.model3 import Model3
from haem_kinetics.models.model4a import Model4a
from haem_kinetics.models.model4b import Model4b
from haem_kinetics.models.model5 import Model5
from haem_kinetics.models.model6 import Model6
from haem_kinetics.models.model7 import Model7
from haem_kinetics.models.model8 import Model8
from haem_kinetics.models.model9a import Model9a
from haem_kinetics.models.model9b import Model9b
from haem_kinetics.models.model9c import Model9c
from haem_kinetics.models.model10 import Model10
from haem_kinetics.models.model11 import Model11
from haem_kinetics.models.model12a import Model12a
from haem_kinetics.models.model12b import Model12b
from haem_kinetics.models.model12c import Model12c
from haem_kinetics.models.model13 import Model13
from haem_kinetics.models.model14a import Model14a
from haem_kinetics.models.model14b import Model14b
from haem_kinetics.models.model99 import Model99
from haem_kinetics.models.degradation import Degradation


t_start = 0
t_end = 1700
t_step = 20
t_eval = range(t_start, t_end, t_step)
t_span = [t_start, t_end]

runs = [
    ('examples/model1.png', Model1, [0.018, 0.0, 0.0, 0.36], {}),
    ('examples/model2a.png', Model2a, [0.018, 0.0, 0.0, 0.36], {}),
    ('examples/model2b.png', Model2b, [0.018, 0.0, 0.0, 0.36], {}),
    ('examples/model3.png', Model3, [0.018, 0.0, 0.0, 0.36], {}),
    ('examples/model4a.png', Model4a, [0.018, 0.0, 0.0, 0.36], {}),
    ('examples/model4b.png', Model4b, [0.018, 0.0, 0.0, 0.36], {}),
    ('examples/model5.png', Model5, [0.018, 0.0, 0.0, 0.36], {}),
    ('examples/model6.png', Model6, [0.018, 0.0, 0.0, 0.36], {}),
    ('examples/model7.png', Model7, [0.018, 0.0, 0.0, 0.36], {}),
    ('examples/model8.png', Model8, [0.018, 0.0, 0.0, 0.36], {}),
    ('examples/model9a.png', Model9a, [0.018, 0.0, 0.0, 0.36], {}),
    ('examples/model9b.png', Model9b, [0.018, 0.0, 0.0, 0.36], {}),
    ('examples/model9c.png', Model9c, [0.018, 0.0, 0.0, 0.36], {}),
    ('examples/model10.png', Model10, [0.018, 0.0, 0.0, 0.36], {}),
    ('examples/model11.png', Model11, [0.018, 0.0, 0.0, 0.36], {}),
    ('examples/model12a.png', Model12a, [0.018, 0.0, 0.0, 0.36], {}),
    ('examples/model12b.png', Model12b, [0.018, 0.0, 0.0, 0.36], {}),
    ('examples/model12c.png', Model12c, [0.018, 0.0, 0.0, 0.36], {}),
    ('examples/model13.png', Model13, [0.018, 0.0, 0.0, 0.36], {}),
    ('examples/model14a.png', Model14a, [0.018, 0.0, 0.0, 0.36], {}),
    ('examples/model14b.png', Model14b, [0.018, 0.0, 0.0, 0.36], {}),
    ('examples/model99.png', Model99, [0.018, 0.0, 0.0, 0.36], {}),
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
               and c not in ('conc_hb_dv_obs', 'conc_hb_assay', 'conc_fe3pp_free')]
    tot = float(model.concentrations[dv_cols].iloc[-1].sum() + host.iloc[-1])
    hz = float(model.concentrations['conc_hz'].iloc[-1]) if 'conc_hz' in model.concentrations.columns else 0.0
    print(f'  end total Fe={tot:.2f} fg  Hz={hz:.2f} fg  host={float(host.iloc[-1]):.2f} fg')
    has_hm = (
        'conc_fe3pp' in model.concentrations.columns
        or 'conc_fe3pp_aq' in model.concentrations.columns
    )
    if has_hm and 'conc_hz' in model.concentrations.columns:
        model.score_vs_experiment()
        print(model.format_fit_metrics())
    else:
        print('  fit metrics skipped (sandbox / incomplete speciation)')

print('Done.')
