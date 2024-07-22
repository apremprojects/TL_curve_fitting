from inspect import getmembers, isfunction
from _pi_cffi import ffi, lib
import numpy as np
#functions_list = getmembers(_pi_cffi)

bounds = [0.5,2,1e6,1e11,1e5,10e5]

arr = np.transpose(np.loadtxt("TL_5peak_a.csv", delimiter=",", skiprows=1))
x_data = ffi.new("double[]", arr[0].tolist())# + (np.ones(arr[0].size) * 273.15)
#above code is for files in kelvin, uncomment if in celcius
y_data = ffi.new("double[]", arr[1].tolist())#

print(x_data)
#double solve(const double _peaks, const double atol, const double tol, const int popsize, double* _bounds, double* _x_data, double* _y_data, size_t _xy_size) {
res = lib.solve(5, 0.0, 10000, 15, bounds, x_data, y_data, len(x_data))
print(res)
print(res[0])
for i in range(1, 1 + 5 * 3):
    print(res[i], end = " ")