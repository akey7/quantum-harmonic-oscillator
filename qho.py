import numpy as np
from scipy.integrate import quad
from math import pi, sqrt, exp, factorial


class QuantumHarmonicOscillator:
    """
    This models a harmonic oscillator. Parameters used throughout the methods
    are:

    v: The quantum of the harmonic oscillator
    k: The force constant
    mass_r: The reduced mass of the system. If you need to calculate this for a 
            diatomic molecule, see the diatomic_reduced_mass() method below.

    The nomenclature of variable names follows Atkins and de Paula Physical
    Chemistry 8th ed, pp. 290-297.
    """

    def __init__(self, mass_r, k):
        """
        Setup the values used by all methods in this model.

        Parameters
        ----------
        mass: float
            The reduced mass of the system, in kg

        k: float
            The force constant, in N/m
        """
        self.hbar = 1.054571817e-34
        self.k = k
        self.mass_r = mass_r
        self.omega = sqrt(k / mass_r)

    @staticmethod
    def diatomic_reduced_mass(mass_1, mass_2):
        """
        Computes the reduced mass of a diatomic system. In order to provide
        meaningful results, please use the exact isotopic mass for your molecule.

        Parameters
        ----------
        mass_1: float
            The mass of the first atom in kg.

        mass_2: float
            The mass of the second atom in kg.
        """
        return mass_1*mass_2/(mass_1+mass_2)

    def hermite(self, v, gamma):
        """
        Returns the value of the v-th (nth) Hermite polynomial evaluated on gamma

        The v and gamma notation follows Atkins' Physical Chemistry 8th ed.

        Parameters
        ----------
        v: int
            The v-th (nth) Hermite polynomial

        gamma: float
            The value to calculate with the Hermite polynomial

        Returns
        -------
        float
            The value of the v-th Hermite polynomial evaluated with gamma.

        Raises
        ------
        Exception
            Raises an exception if the nth Hermite polynomial is not
            supported.
        """
        if v == 0:
            return 1
        elif v == 1:
            return 2 * gamma
        elif v == 2:
            return 4 * gamma**2 - 2
        elif v == 3:
            return 8 * gamma**3 - 12 * gamma
        elif v == 4:
            return 16 * gamma**4 - 48 * gamma**2 + 12
        elif v == 5:
            return 32 * gamma**5 - 160 * gamma**3 + 120 * gamma
        elif v == 6:
            return 64 * gamma**6 - 480 * gamma**4 + 720 * gamma**2 - 120
        else:
            raise Exception(f"Hermite polynomial {v} is not supported")

    def max_v(self):
        """
        Returns
        -------
        int
            The maximum v value for the instance.
        """
        return 6

    def energy(self, v):
        """
        Calculate the energy at the given level v of the system

        Parameters
        ----------
        v: int
            The quantum vmber v for the energy level of this system

        Returns
        -------
        float
            Energy of the system in Joules.
        """
        return (v + 0.5) * self.hbar * self.omega

    def energy_separation(self):
        """
        Returns
        -------
        float
            The energy difference between adjacent energy levels in Joules.
        """
        return self.hbar * self.omega

    def wavefunction(self, v, x):
        """
        Returns the value of the wavefunction at energy level v
        at coordinate x.

        Parameters
        ----------
        v: float
            Energy level of the system.

        x: float
            x coordinate of the particle in m.

        Returns
        -------
        float
            Value of the wavefunction v at x.
        """
        alpha = (self.hbar**2 / self.mass_r / self.k) ** 0.25
        gamma = x / alpha
        normalization = sqrt(1 / (alpha * sqrt(pi) * 2**v * factorial(v)))
        gaussian = exp((-(gamma**2)) / 2)
        hermite = self.hermite(v, gamma)
        return normalization * hermite * gaussian

    def wavefunction_across_range(self, v, x_min, x_max, points=100):
        """
        Calculates the wavefunction across a range.

        Parameters
        ----------
        v: int
            The quantum vmber of the system.

        x_min: float
            The minimum x value to calculate.

        x_max: float
            The maximum x value to calculate.

        points: int
            The vmber of points across the range

        Returns
        -------
        np.array, list
            The first array are the x coordinates, the second list are the
            float values of the wavefunction.
        """
        xs = np.linspace(x_min, x_max, points)
        ys = [self.wavefunction(v, x) for x in xs]
        return xs, ys

    def prob_density(self, v, x_min, x_max, points=100):
        """
        Returns the probability density between x_min and x_max for a given
        vmber of points at energy level v.

        Parameters
        ----------
        v: int
            Quantum vmber of the system.

        x_min: float
            Minimum of length being calculated. Probably negative. Units are
            meters.

        x_max: float
            Maximum of length being calculated. Probably positive. Units are
            meters.

        points: int
            The vmber of points to compute the probability density for

        Returns
        -------
        np.array, list
            The first array is the list of x coordinates. The list are the
            corresponding values of the probability density.
        """
        xs = np.linspace(x_min, x_max, points)
        ys = [self.wavefunction(v, x) ** 2 for x in xs]
        return xs, ys

    def integrate_prob_density_between_limits(self, v, x_min, x_max):
        """
        As a way of testing the methods in this class, provide a way to
        integrate across the wavefunction squared between limits.

        Parameters
        ----------
        v: int
            Quantum vmber of the system.

        x_min: float
            lower bound of integration

        x_max: float
            upper bound of integration
        """

        def integrand(x):
            return self.wavefunction(v, x) ** 2

        result, _ = quad(integrand, x_min, x_max)
        return result
