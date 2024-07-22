import numpy as np
import pandas as pd
#from scipy.optimize import curve_fit, differential_evolution, brute, direct, shgo

class Model:
    def __init__(self, beta, k_J_per_K):
        self.beta = beta
        self.k_J_per_K = k_J_per_K
    def n0_N_constraint(self, x):
        for i in range(0, x.size, 4):
            if(x[i+2] > x[i+3]):
                return np.inf
        return 0
    def antiderivative(self, T_prime, E_J):
        return self.k_J_per_K / E_J * (T_prime**2) * np.exp(-E_J / (self.k_J_per_K * T_prime))
    def evaluate(self, T, *args):
        res = np.zeros(T.size)
        for i in range(0, len(args), 4):
            E_ev = args[i]
            S = args[i + 1]
            #nf = args[i + 2]
            n0 = args[i+2]
            N = args[i + 3]
            #N = 1e5
            #n0 = N*nf
            E_J = E_ev * (1.602e-19)
            T0_antiderivative = self.antiderivative(T[0], E_J)
            for j in range(0, T.size):
                integral_result = self.antiderivative(T[j], E_J) - T0_antiderivative
                numerator = n0**2 * S * np.exp(-E_J / (self.k_J_per_K * T[j]))
                denominator = N * (1 + (n0 * S) / (self.beta * N) * integral_result)**2 
                res[j] += numerator / denominator
        return res
    def residual(self, x, T, intensity_values):
        #print("Res ->", x)
        for i in range(0, x.size, 4):
            if(x[i+2] > x[i+3]):
                return np.inf
        res = self.evaluate(T, *x)
        return np.sum((intensity_values - res)**2)
    def std_dev(self, x, T, intensity_values):
        #print("Res ->", x)
        res = self.evaluate(T, *x)
        return np.sum(abs(intensity_values - res)) / T.size
    def vectorized_residual(self, xs, T, intensity_values):
        #print("VR -> ", xs.shape)
        if(xs.ndim == 1):
            return self.residual(xs, T, intensity_values)
        rows, columns = xs.shape
        ans = np.empty(columns)
        for i in range(0, columns):
            x = xs[:, i]
            #print("x -> ", x)
            res = self.evaluate(T, *x)
            ans[i] = np.sum((intensity_values - res)**2)
            for j in range(0, x.size, 4):
                if(x[j+2] > x[j+3]):
                    ans[i] = np.inf
                    break
        return ans