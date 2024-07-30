from scipy.optimize import curve_fit, differential_evolution, brute, direct, shgo, dual_annealing
import numpy as np
import pandas as pd


class Solver:
    def __init__(self, model, x_data, y_data):
        self.model = model
        self.x_data = x_data
        self.y_data = y_data
    def n0_N_constraint(self, x):
        for i in range(0, x.size, 4):
            if(x[i + 2] > x[i + 3]):
                return np.inf
        return 0.5
    def differential_evolution(self, bounds, strategy='best1bin', maxiter=1000, popsize=15, tol=0.01, mutation=(0.5, 1), recombination=0.7, callback=None, disp=False, polish=True, init='latinhypercube', atol=0, updating='immediate', workers=1, vectorized = False, x0 = None, constraints = None):
        if(vectorized):
            if(constraints != None):
                res = differential_evolution(
                    func = self.model.vectorized_residual,
                    args = (self.x_data, self.y_data),
                    bounds = bounds,
                    strategy=strategy,
                    maxiter = maxiter,
                    popsize = popsize,
                    tol = tol,
                    mutation=mutation,
                    recombination=recombination,
                    disp = disp,
                    polish = polish,
                    init = init,
                    updating=updating,
                    workers = workers,
                    vectorized=True,
                    callback = callback,
                    atol = atol,
                    seed = 3,
                    x0 = x0,
                    constraints = constraints
                )
            else:
                res = differential_evolution(
                    func = self.model.vectorized_residual,
                    args = (self.x_data, self.y_data),
                    bounds = bounds,
                    strategy=strategy,
                    maxiter = maxiter,
                    popsize = popsize,
                    tol = tol,
                    mutation=mutation,
                    recombination=recombination,
                    disp = disp,
                    polish = polish,
                    init = init,
                    updating=updating,
                    workers = workers,
                    vectorized=True,
                    callback = callback,
                    atol = atol,
                    seed = 3,
                    x0 = x0
                )
        else:
            if(constraints != None):
                res = differential_evolution(
                    func = self.model.residual,
                    args = (self.x_data, self.y_data),
                    bounds = bounds,
                    strategy=strategy,
                    maxiter = maxiter,
                    popsize = popsize,
                    tol = tol,
                    mutation=mutation,
                    recombination=recombination,
                    disp = disp,
                    polish = polish,
                    init = init,
                    updating=updating,
                    workers = workers,
                    vectorized=False,
                    callback = callback,
                    atol = atol,
                    seed = 3,
                    x0 = x0,
                    constraints=constraints
                )
            else:
                res = differential_evolution(
                    func = self.model.residual,
                    args = (self.x_data, self.y_data),
                    bounds = bounds,
                    strategy=strategy,
                    maxiter = maxiter,
                    popsize = popsize,
                    tol = tol,
                    mutation=mutation,
                    recombination=recombination,
                    disp = disp,
                    polish = polish,
                    init = init,
                    updating=updating,
                    workers = workers,
                    vectorized=False,
                    callback = callback,
                    atol = atol,
                    seed = 3,
                    x0 = x0
                )
        return res.x
    def direct(self, bounds, eps=0.0001, maxfun=None, maxiter=1000, locally_biased=True, f_min=0, f_min_rtol=0.0001, vol_tol=1e-16, len_tol=1e-06, callback=None):
        res = direct(
            func = self.model.residual,
            bounds = bounds,
            args = (self.x_data, self.y_data),
            eps = eps,
            maxfun = maxfun,
            maxiter = maxiter,
            locally_biased = locally_biased,
            f_min = f_min,
            f_min_rtol = f_min_rtol,
            vol_tol = vol_tol,
            len_tol = len_tol,
            callback = callback
        )
        return res.x
    def dual_annealing(self, bounds, maxiter=1000, initial_temp=5230.0, restart_temp_ratio=2e-05, visit=2.62, accept=-5.0, maxfun=10000000.0, no_local_search=False, callback=None, x0 = None):
        res = dual_annealing(
            func = self.model.residual,
            bounds = bounds,
            args = (self.x_data, self.y_data),
            maxiter = maxiter,
            initial_temp = initial_temp,
            restart_temp_ratio = restart_temp_ratio,
            visit = visit,
            accept = accept,
            maxfun = maxfun,
            seed = 3,
            no_local_search = no_local_search,
            callback = callback,
            x0 = x0
            )
        return res.x
    def shgo(self, bounds, n=100, iters=1, method='simplicial', options = {}, workers=1, callback = None):
        res = shgo(
            func = self.model.residual,
            bounds = bounds,
            args = (self.x_data, self.y_data),
            n=n,
            iters = iters,
            sampling_method=method,
            workers = workers,
            options = options,
            callback = callback
        )
        return res.x
    def brute_force(self, bounds, Ns, workers, disp, x0 = None):
        print(bounds, Ns, workers, disp)
        x0 = brute(
            func = self.model.residual, 
            ranges = bounds,
            args = (self.x_data, self.y_data),
            Ns = Ns, 
            workers = workers, 
            disp = disp,
            x0 = x0
        )
        return x0
    def curve_fit(self, initial_guess, bounds, ftol = 0.001, sigma = None, absolute_sigma = False):
        popt, pcov, a, b, c = curve_fit(
            f = self.model.evaluate,
            xdata = self.x_data,
            ydata = self.y_data,
            p0 = initial_guess,
            maxfev = 1000000,
            ftol = ftol,
            sigma = sigma,
            absolute_sigma = absolute_sigma,
            bounds = bounds,
            full_output = True
        )
        #print(a,b,c)
        return popt
    def curve_fit2(self, initial_guess, bounds, ftol = 0.001, sigma = None, absolute_sigma = False):
        popt, pcov, a, b, c = curve_fit(
            f = self.model.evaluate,
            xdata = self.x_data,
            ydata = self.y_data,
            p0 = initial_guess,
            maxfev = 1000000,
            ftol = ftol,
            sigma = sigma,
            absolute_sigma = absolute_sigma,
            bounds = bounds,
            full_output = True
        )
        print(a,b,c)
        return popt