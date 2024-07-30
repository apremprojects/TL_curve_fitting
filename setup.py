from cffi import FFI
ffibuilder = FFI()

# cdef() expects a single string declaring the C types, functions and
# globals needed to use the shared object. It must be in valid C syntax.
ffibuilder.cdef(
"""
     double* solve(bool verbose, const int _peaks, const int maxiter, const double beta, const double atol, const double tol, const int popsize, double* _bounds, double* _x_data, double* _y_data, size_t _xy_size);
"""
)

# set_source() gives the name of the python extension module to
# produce, and some C source code as a string.  This C code needs
# to make the declarated functions, types and globals available,
# so it is often just the "#include".
ffibuilder.set_source("cpp_differential_evolution",
     """#include "C++/include/CDifferentialEvolution.h"   // the C header of the library""",
     sources=['C++/src/CDifferentialEvolution.cpp', 'C++/src/CPPDifferentialEvolution.cpp'],
     source_extension='.cpp'
)

if __name__ == "__main__":
    ffibuilder.compile(verbose=True)