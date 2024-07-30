import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.integrate import quad
from scipy.optimize import curve_fit, differential_evolution, brute, direct, shgo, NonlinearConstraint
import argparse
from plotter import Plotter
from solver import Solver
from model import Model
from cpp_differential_evolution import lib, ffi
from progress_bar import ProgressBar
import csv
from rich.console import Console
from rich.markdown import Markdown

def minimums_maxes(initial_parameters):
    Es = []
    Ss = []
    n0s = []
    Ns = []
    for i in range(0, initial_parameters.size, 3):
        Es.append(initial_parameters[i])
    for i in range(1, initial_parameters.size, 3):
        Ss.append(initial_parameters[i])
    for i in range(2, initial_parameters.size, 3):
        n0s.append(initial_parameters[i])
    for i in range(3, initial_parameters.size, 3):
        Ns.append(initial_parameters[i])
    return ([min(Es), min(Ss), min(n0s), min(Ns)], [max(Es), max(Ss), max(n0s), max(Ns)])

def write_to_file(filename, peaks, x_data, input_y_data, peak_y_datas, fitted_y_data, answer, initial_guess, residual, fom):
    #Format: E, S, n0, N, x_data, input_y_data, fitted_y_data, residual, fom, residual
    with open(filename, 'w') as f:
        csv_writer = csv.DictWriter(f, fieldnames=["E_fitted", "s_fitted", "N_fitted", "n0_fitted", "E_initial", "s_initial", "N_initial", "n0_initial", "temp_K", "int_peak1", "int_peak2", "int_peak3", "int_peak4", "int_peak5", "int_peak6", "int_peak7", "int_peak8", "int_peak9", "int_sum", "temp_data_K", "int_data_cps"])
        csv_writer.writeheader()
        for i in range(0, peaks):
            row = {
                "E_fitted": answer[4*i],
                "s_fitted": answer[4*i + 1],
                "n0_fitted": answer[4*i + 2],
                "N_fitted": answer[4*i + 3],
                "E_initial": initial_guess[4*i],
                "s_initial": initial_guess[4*i + 1],
                "n0_initial": initial_guess[4*i + 2],
                "N_initial": initial_guess[4*i + 3],
                "temp_K": x_data[i],
                "int_sum": fitted_y_data[i],
                "temp_data_K": x_data[i],
                "int_data_cps": input_y_data[i]
            }
            for j in range(0, peaks):
                row[f"int_peak{j+1}"] = peak_y_datas[j][i]
            csv_writer.writerow(row)
        for i in range(peaks, x_data.size):
            row = {
                "temp_K": x_data[i],
                "int_sum": fitted_y_data[i],
                "temp_data_K": x_data[i],
                "int_data_cps": input_y_data[i]
            }
            for j in range(0, peaks):
                row[f"int_peak{j+1}"] = peak_y_datas[j][i]
            csv_writer.writerow(row)

def generate_bounds(N, Eb, Sb, n0b, Nb):
    res = []
    for i in range(0, N):
        res.append(Eb)
        res.append(Sb)
        res.append(n0b)
        res.append(Nb)
    return res


def generate_cpp_bounds(N, Eb, Sb, n0b, Nb):
    res = [Eb[0], Eb[1], Sb[0], Sb[1], n0b[0], n0b[1], Nb[0], Nb[1]]
    return res


def generate_curve_bounds(N, Eb, Sb, n0b, Nb):
    res = []
    lowest = []
    highest = []
    for i in range(0, N):
        lowest.append(Eb[0])
        lowest.append(Sb[0])
        lowest.append(n0b[0])
        lowest.append(Nb[0])
    for i in range(0, N):
        highest.append(Eb[1])
        highest.append(Sb[1])
        highest.append(n0b[1])
        highest.append(Nb[1])
    res.append(lowest)
    res.append(highest)
    return res


if __name__ == "__main__":
    #python main.py three_peak.csv --shgo --bounds 0.07 3 1e8 1e11 1 10 --peaks 3 --maxfev 3000 --f_min 0 --verbose --workers 2
    parser = argparse.ArgumentParser(
        prog='TL_curve_fitting',
        add_help=False
    )
    parser.add_argument('filename')           # positional argument

    parser.add_argument('--initial_guess', nargs = "*", type = float, default = None)

    parser.add_argument('--basinhopping', action='store_true')
    
    parser.add_argument('--brute', action='store_true')
    parser.add_argument('--ns', type = float, default = 100)

    parser.add_argument('--cppdifferential_evolution', action='store_true')

    parser.add_argument('--differential_evolution', action='store_true')
    parser.add_argument('--strategy', type = str, default = "best1bin")#ignored by cppdiff
    parser.add_argument('--popsize', type = int, default = 15)
    parser.add_argument('--tol', type = float, default = 0.01)
    parser.add_argument('--mutation', nargs = 2, type = float, default = (0.5, 1))#ignored by cppdiff
    parser.add_argument('--recombination', type = float, default = 0.7)#ignored by cppdiff
    parser.add_argument('--polish', action='store_true', default = False)#ignored by cppdiff
    parser.add_argument('--init', type = str, default = "latinhypercube")#ignored by cppdiff
    parser.add_argument('--atol', type = float, default = 0)
    parser.add_argument('--updating', type = str, default = "immediate")#ignored by cppdiff
    parser.add_argument('--vectorized', action='store_true', default = False)#ignored by cppdiff

    parser.add_argument('--shgo', action='store_true')
    parser.add_argument('--n', type = int, default = 100)
    parser.add_argument('--iters', type = int, default = 1)
    parser.add_argument('--sampling_method', type = str, default = "simplicial")
    parser.add_argument('--maxfev', type = int) #options {}
    parser.add_argument('--f_min', type = float, default = 0)
    parser.add_argument('--f_tol', type = float)
    parser.add_argument('--ftol', type = float)
    parser.add_argument('--maxiter', type = int)
    parser.add_argument('--maxev', type = int)
    parser.add_argument('--maxtime', type = float)
    parser.add_argument('--minhgrd', type = int)
    parser.add_argument('--no_minimize_every_iter', action='store_true', default = False)
    parser.add_argument('--local_iter', type = int)

    parser.add_argument('--dual_annealing', action='store_true')
    #parser.add_argument('--maxiter', type = int, default = 1000)
    parser.add_argument('--initial_temp', type = float, default = 5230)
    parser.add_argument('--restart_temp_ratio', type = float, default = 2e-5)
    parser.add_argument('--visit', type = float, default = 2.62)
    parser.add_argument('--accept', type = float, default = -5)
    #parser.add_argument('--maxfun', type = int, default = 1e7)
    parser.add_argument('--no_local_search', action='store_true')

    parser.add_argument('--direct', action='store_true')
    parser.add_argument('--eps', type = float, default = 1e-4)
    parser.add_argument('--maxfun', type = int, default = None)
    #parser.add_argument('--maxiter', type = int, default = 1000)
    parser.add_argument('--locally_biased', action='store_true', default = False)
    #parser.add_argument('--f_min', type = float, default = 0)
    parser.add_argument('--f_min_rtol', type = float, default = 1e-4)
    parser.add_argument('--vol_tol', type = float, default = 1e-16)
    parser.add_argument('--len_tol', type = float, default = 1e-6)

    parser.add_argument('--peaks', type = int)
    parser.add_argument('--bounds', nargs = '*', type = str)
    parser.add_argument('--workers', type = int, default = 1)
    parser.add_argument('--sigma', type = float, default = None)
    parser.add_argument('--curve_ftol', type = float, default = 0.001)

    parser.add_argument('--absolute_sigma', action='store_true')
    parser.add_argument('--verbose', action='store_true') #enables disp = True

    parser.add_argument('--output', type = str, default = None)
    parser.add_argument('--beta', type = float, default = None)
    parser.add_argument('-h', '--help', action='store_true')
    args = parser.parse_args()
    #print(args)

    if(args.bounds == None or args.peaks == None):
        if(args.help):
            pass
            with open("USAGE.md", 'r') as f:
                console = Console()
                console.print(Markdown(f.read()))
                exit(1)
        else:
            if(args.peaks == None and args.bounds == None):
                raise RuntimeError("the following arguments are required: --peaks, --bounds")
            elif(args.peaks == None):
                raise RuntimeError("the following arguments are required: --peaks")
            elif(args.bounds == None):
                raise RuntimeError("the following arguments are required: --bounds")

    if(args.initial_guess != None):
        if(len(args.initial_guess) != args.peaks * 4):
            print(args.initial_guess)
            raise RuntimeError(f"{args.peaks*4} parameters for initial_guess expected, {len(args.initial_guess)} provided")
        
    #Process bounds
    if(len(args.bounds) == 1):
        #Read bounds from .csv
        arr = np.loadtxt(args.bounds[0], delimiter=",", skiprows=1, max_rows=args.peaks).flatten()

        #store in bounds
        bounds = []
        cpp_bounds = []
        for i in range(0, arr.size, 2):
            bounds.append([arr[i], arr[i+1]])
            cpp_bounds.append(arr[i])
            cpp_bounds.append(arr[i+1])

        #get lower bounds and higher bounds for curve_bounds
        lowest = []
        highest = []
        for i in range(0, arr.size, 2):
            lowest.append(arr[i])
            highest.append(arr[i+1])
        
        #store in curve_bounds
        curve_bounds = [lowest, highest]

    elif(len(args.bounds) == 8):

        #store bounds based on command line arguments
        bounds = generate_bounds(
            N = args.peaks,
            Eb = [args.bounds[0], args.bounds[1]],
            Sb = [args.bounds[2], args.bounds[3]],
            n0b = [args.bounds[4], args.bounds[5]],
            Nb = [args.bounds[6], args.bounds[7]]
        )

        curve_bounds = generate_curve_bounds(
            N = args.peaks,
            Eb = (args.bounds[0], args.bounds[1]),
            Sb = (args.bounds[2], args.bounds[3]),
            n0b = (args.bounds[4], args.bounds[5]),
            Nb = (args.bounds[6], args.bounds[7])
        )

        cpp_bounds = generate_cpp_bounds(
            N = args.peaks,
            Eb = (args.bounds[0], args.bounds[1]),
            Sb = (args.bounds[2], args.bounds[3]),
            n0b = (args.bounds[4], args.bounds[5]),
            Nb = (args.bounds[6], args.bounds[7])
        )
    else:
        raise RuntimeError("Not enough bounds/no filename provided")

    #Read data from input .csv file
    x_data = []
    y_data = []
    with open(args.filename) as file:
        csv_reader = csv.DictReader(file)
        for row in csv_reader:
            #print(row)
            row = dict(row)
            if "temp_c" in row:
                #print("Celsius")
                x_data.append(float(row["temp_c"]) + 273.15)
            elif "temp_C" in row:
                #print("Celsius")
                x_data.append(float(row["temp_C"]) + 273.15)
            elif "temp_k" in row:
                #print("Kelvin")
                x_data.append(float(row["temp_k"]))
            elif "temp_K" in row:
                #print("Kelvin")
                x_data.append(float(row["temp_K"]))
            else:
                raise RuntimeError(f"Incorrect .csv file format, expected temp_c/temp_k, intensity_cps, got {row.keys()}")
            y_data.append(float(row["intensity_cps"]))
    x_data = np.array(x_data)
    y_data = np.array(y_data)
    #print(x_data, y_data)
    #exit(1)

    #Create Plotter class to plot data
    print("Creating plotter...")
    plotter = Plotter()

    #Create Model class to model curves
    print("Creating model...")
    model = Model(
        #N = 1E5,    # Total concentration of traps
        beta = 2.0 if args.beta == None else args.beta,    # Heating rate
        k_J_per_K = 1.381e-23  # Boltzmann constant in J/K
    )

    #Create Solver class to fit data
    print("Creating solver...")
    solver = Solver(
        model = model,
        x_data = x_data,
        y_data = y_data
    )

    #Constraint ensures that n0 is less than N
    constraints = NonlinearConstraint(fun = solver.n0_N_constraint, lb = 0, ub = 1)
    
    #Clip initial_guess to fit within bounds
    initial_guess = args.initial_guess
    if(initial_guess!=None):
        initial_guess = np.clip(initial_guess, curve_bounds[0], curve_bounds[1])

    #Switch between different types of Global Solver
    if(args.basinhopping):
        pass
        #Basin Hopping Differential Evolution
    elif(args.brute):
        #Brute Global Solver
        if(args.vectorized):
            raise RuntimeError("Vectorization cannot be used with brute...")
        print("The program does not support displaying a progress bar for this solver...")
        initial_parameters = solver.brute_force(
            bounds = bounds,
            workers = args.workers,
            Ns = args.ns,
            disp = True
        )
    elif(args.dual_annealing):
        #Dual Annealing Global Solver
        print("Generating initial parameters with dual_annealing...")
        if(args.workers != 1 or args.vectorized):
            raise RuntimeError("Vectorization and/or Parallelization cannot be used with dual_annealing...")
        if(not args.verbose):
            print("The program does not support displaying a progress bar for this solver...")
        pbar = ProgressBar(total = 1e7 if(args.maxfun == None) else args.maxfun, desc = "Progress", bar = False)
        initial_parameters = solver.dual_annealing(
            bounds = bounds,
            maxiter = 1000 if(args.maxiter == None) else args.maxiter,
            initial_temp = args.initial_temp,
            restart_temp_ratio = args.restart_temp_ratio,
            visit = args.visit,
            accept = args.accept,
            maxfun = 1e7 if(args.maxfun == None) else args.maxfun,
            no_local_search = args.no_local_search,
            callback = pbar.update_da,
            x0 = initial_guess
        )
    elif(args.cppdifferential_evolution):
        #C++ Differential Evolution Global Solver
        print("Generating initial parameters with cpp_differential_evolution...")
        if(args.workers != 1 or args.vectorized):
            raise RuntimeError("Vectorization and/or Parallelization cannot be used with cpp_differential_evolution...")
        x_buf = ffi.new("double[]", x_data.tolist())
        y_buf = ffi.new("double[]", y_data.tolist())
        cpp_buf = ffi.new("double[]", cpp_bounds)
        res = lib.solve(args.verbose, args.peaks, 1000 if args.maxiter == None else args.maxiter, 2.0 if args.beta == None else args.beta, args.atol, args.tol, args.popsize, cpp_buf, x_buf, y_buf, len(x_data))
        #print("HERE\n")
        initial_parameters = np.zeros(args.peaks * 4)
        for i in range(0, args.peaks * 4):
            initial_parameters[i] = ffi.cast("double", res[i + 1])
    elif(args.differential_evolution):
        #Differential Evolution Global Solver
        print("Generating initial parameters with differential_evolution...")
        if(args.workers != 1 and args.vectorized):
            raise RuntimeError("Vectorization and Parallelization cannot be used at the same time...")
        pbar = ProgressBar(total = 100, desc = "Progress", bar = not args.verbose)
        #print(initial_guess)
        initial_parameters = solver.differential_evolution(
            bounds = bounds,
            strategy = args.strategy,
            maxiter = args.maxiter,
            tol = args.tol,
            mutation = args.mutation,
            recombination = args.recombination,
            disp = args.verbose,
            init = args.init,
            atol = args.atol,
            updating = "deferred" if(args.workers!=1 or args.vectorized) else args.updating,
            workers = args.workers,
            vectorized = args.vectorized,
            polish = args.polish,
            callback = pbar.update,
            x0 = initial_guess,
            constraints = constraints
        )
    elif(args.direct):
        #DIRECT Global Solver
        print("Generating initial parameters with DIRECT...")
        if(args.workers != 1):
            raise RuntimeError("DIRECT doesn't support parallelism")
        if(not args.verbose):
            print("The program does not support displaying a progress bar for this solver...")
        if(args.initial_guess != None):
            print("This solver doesn't support initial guesses...")
        pbar = ProgressBar(total = 3000*args.peaks if args.maxfun == None else args.maxfun, desc = "Progress", bar = False)
        initial_parameters = solver.direct(
            bounds = bounds,
            eps = args.eps,
            maxfun = 3000*args.peaks if args.maxfun == None else args.maxfun,
            maxiter = 1000 if args.maxiter == None else args.maxiter,
            locally_biased = args.locally_biased,
            f_min = args.f_min,
            f_min_rtol=args.f_min_rtol,
            vol_tol = args.vol_tol,
            len_tol = args.len_tol,
            callback=pbar.update_xk
        )
    elif(args.shgo):
        #SHGO Global Solver
        if(args.initial_guess != None):
            print("This solver doesn't support initial guesses...")
        print("Generating initial parameters with SHGO...")
        options = {}
        if(args.maxfev != None):
            options["maxfev"] = args.maxfev
        if(args.f_min != None):
            options["f_min"] = args.f_min
        if(args.f_tol != None):
            options["f_tol"] = args.f_tol
        if(args.ftol != None):
            options["ftol"] = args.ftol
        if(args.maxiter != None):
            options["maxiter"] = args.maxiter
        if(args.maxev != None):
            options["maxev"] = args.maxev
        if(args.maxtime != None):
            options["maxtime"] = args.maxtime
        if(args.minhgrd != None):
            options["minhgrd"] = args.minhgrd
        if(args.local_iter != None):
            options["local_iter"] = args.local_iter
        options["minimize_every_iter"] = not args.no_minimize_every_iter
        options["disp"] = False
        if(not args.verbose):
            print("The program does not support displaying a progress bar for this solver...")
        pbar = ProgressBar(total = args.maxfev, desc = "Progress", bar = False)
        initial_parameters = solver.shgo(
            bounds = bounds,
            n = args.n,
            iters = args.iters,
            method = args.sampling_method,
            options = options,
            workers = args.workers,
            callback = pbar.update_xk
        )
        pbar.close()
    else:
        #No Global Solver
        if(args.initial_guess!=None):
            print("Moving directly to curve_fit...")
            initial_parameters = args.initial_guess
        else:
            raise RuntimeError("Specify a curve-fitting method")

    #Local Curve Fit Solver
    print("Fitting curve with initial parameters...")
    initial_parameters = np.clip(initial_parameters, curve_bounds[0], curve_bounds[1])
    answer = solver.curve_fit(
        initial_guess = initial_parameters,
        bounds = curve_bounds,
        ftol = args.curve_ftol,
        sigma = args.sigma,
        absolute_sigma = args.absolute_sigma
    )
    print("Initial parameters (generated by Global Solver) -> ", initial_parameters)
    print("Answer (generated by (locally) curve fitting from Initial parameters) -> ", answer)


    #calculate residuals and figure of merit
    residual = model.residual(answer, x_data, y_data)
    fom = model.figure_of_merit(answer, x_data, y_data)
    initial_parameters_y_data = model.evaluate(x_data, *initial_parameters)
    fitted_y_data = model.evaluate(x_data, *answer)

    #plot data
    print("Plotting...")
    plotter.plot(
        x_data = x_data,
        y_data = y_data,
        alpha = 1,
        label = 'Real data curve'
    )
    plotter.plot(
        x_data = x_data, 
        y_data = initial_parameters_y_data,
        alpha = 0.2,
        label = 'Initial data curve'
    )
    plotter.plot(
        x_data = x_data,
        y_data = fitted_y_data,
        alpha = 1,
        label = 'Fitted data curve'
    )
    peak_y_datas = []
    for i in range(0, int(answer.size / 4)):
        print(answer[i * 4: (i + 1) * 4])
        peak_y_datas.append(model.evaluate(x_data, *(answer[i * 4: (i + 1) * 4])))
        plotter.fill_between(
            x = x_data,
            y1 = peak_y_datas[len(peak_y_datas)-1],
            alpha = 0.1,
            label = f"Peak #{i+1}"
        )
    
    #print residual & figure of merit
    print(f"Residual -> {residual}")
    print(f"FOM -> {fom}%")
    
    #write to file
    if(args.output!=None):
        write_to_file(args.output, args.peaks, x_data, y_data, peak_y_datas, fitted_y_data, answer, initial_parameters, residual, fom)
    
    #display plot
    plt.xlabel('Temperature, $T$ [K]')
    plt.ylabel('Intensity, $I$ [cps]')
    plt.legend(fontsize=7, loc='best')
    plt.show()