import numpy as np

from scipy.integrate import solve_ivp


class KineticsModel:
    def __int__(self):

        # To be solved
        self.initial_values = {}    # Stores initial concentrations of haem species using in integration
        self.differential_eqs = []  # Stores the differential eqs to be integrated

        # Solved results
        self.concentrations = []    # Stores the solved time-course concentrations of haem species
        self.time = []              # Stores the time-series for the solution
        self.solution = None        # Stores the entire integrated solution (gives access to additional info if needed)

    def _set_initial_conc(self, init):
        """
        Function to set the initial concentration of haem species. Since each model could have different haem species
        or ordering of haem species, this class function must be overwritten for each model
        :param init: Initial values for haem concentrations (Order matters!)
        :return:
        """
        raise NotImplementedError('_set_initial_conc must be overwritten by the model class')

    def _set_diff_eqs(self):
        """
        Function to add differential equations to the class list to be integrated. Since each model has different
        differential equations, this class function must be overwritten for each model.
        :return:
        """
        raise NotImplementedError('_set_diff_eqs must be overwritten by the model class')

    def _integrate(self, t, init) -> np.array[np.array]:
        """
        Function that will integrate differential equations for the model. The differential equations must be set
        by within the model but adding them to the self.differential_eqs list.

        :param t: Time range that will be integrated over
        :param init: Initial values for haem concentrations (Order matters!)
        :return: returns concentrations for each haem species in the model at each time point. E.g.:
                   t0      t1      t2     t3     t4
                 [
                  [0.    0.001    0.01   0.02   0.03]     # Haem species 1
                  [0.    0.       0.001  0.002  0.003]    # Haem species 2
                 ]
        """
        # Set initial concentration values
        self._set_initial_conc(init=init)

        return self.differential_eqs

    def _solve(self, t, init, kwargs):
        """
        Performs integration of differential equations for a model. Reference used for setting this up:
        https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.solve_ivp.html

        :param t: Time range that will be integrated over
        :param init: Initial values for haem concentrations (Order matters!)
        :param kwargs: [Optional] Additional parameters for the integrator
        :return:
        """
        # Ensure there are differential equations to integrate
        self._set_diff_eqs()
        if len(self.differential_eqs) == 0:
            raise ValueError(f'No differential equations have been added to "differential_eqs" variable')

        # Integrate the differential equations
        self.solution = solve_ivp(self._integrate, t, init, **kwargs)
        self.time = self.solution.t
        self.concentrations = self.solution.y

    # ToDo: Implement later
    def _plot(self):
        pass

    def run(self, t, init, kwargs):
        """
        API that solves differential equations and saves output.

        :param t:
        :param init:
        :param kwargs:
        :return:
        """

        raise NotImplementedError('"run" must be overwritten by the model class')
