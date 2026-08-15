from haem_kinetics.models.model2 import Model2a
from haem_kinetics.models.helpers import (
    F_EXP_BREAK_T_MIN,
    fraction_exp_growth_two_phase,
    two_phase_a_late,
)


class Model2b(Model2a):
    """
    Model 2a with a two-phase f_exp, still × remaining host.

    Break at Garnie Fig. 5B 29 h (t = 780 min from 16 h), not fitted.
    f_exp is continuous at the join (a_late from two_phase_a_late).
    a_early, b_early, b_late are empirical (Dd2 DV Fe defaults).
    """

    def __init__(
        self,
        model_name: str = 'Model 2b',
        a_early: float = 0.3224,
        b_early: float = 0.0007607,
        b_late: float = 0.003342,
        t_break_min: float = F_EXP_BREAK_T_MIN,
    ):
        super().__init__(model_name=model_name)
        self.a_early = a_early
        self.b_early = b_early
        self.b_late = b_late
        self.t_break_min = t_break_min
        self.a_late = two_phase_a_late(
            a_early, b_early, b_late, t_break_min=t_break_min,
        )

    def _uptake_dv(self, t):
        host = self._nonneg(self.initial_values[self.HOST_KEY])
        if host <= 0.0:
            return 0.0
        tot_hb_conc = host * self.const.vol_rbc / self._vol_dv(t)
        return fraction_exp_growth_two_phase(
            t,
            self.a_early,
            self.b_early,
            self.a_late,
            self.b_late,
            t_break_min=self.t_break_min,
        ) * tot_hb_conc
