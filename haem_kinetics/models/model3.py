from haem_kinetics.models.model2 import Model2
from haem_kinetics.models.helpers import garnie_pm_amount_scale


class Model3(Model2):
    """
    Model 2 + Garnie Fig. 3 plasmepsin amount schedule s_PM(t).

    Uptake stays f_exp. Shared variable_dv_volume bookkeeping is unchanged.
    [E]_i,eff = s_PM(t) · n_E,i / V_DV(t).
    """

    def __init__(self, model_name: str = 'Model 3'):
        super().__init__(model_name=model_name)

    def _enzyme_conc(self, enzyme, t):
        return garnie_pm_amount_scale(t) * super()._enzyme_conc(enzyme, t)
