#include <R.h>
#include <Rdefines.h>

SEXP affy_residuals (SEXP sets, SEXP vals) {
    unsigned int n_ps   = length (sets), // num. of probesets
                 n_vals = length (vals); // num. of total probes
    double *p_vals = NUMERIC_POINTER (vals);
    SEXP elem;

    // Compute all residuals
    SEXP res;
    PROTECT (res = NEW_NUMERIC (n_vals));
    double *p_res = NUMERIC_POINTER (res);
    for (unsigned int i = 0; i < n_vals; ++i)
        p_res[i] = NA_REAL;
    for (unsigned int ps = 0; ps < n_ps; ++ps) {
        PROTECT (elem = VECTOR_ELT (sets, ps));
        int *p_elem = INTEGER_POINTER (elem);
        unsigned int n_p = length (elem);

        // Compute this probeset's mean value once first
        double mean = 0.0;
        for (unsigned int p = 0; p < n_p; ++p)
            mean += p_vals[p_elem[p] - 1];
        mean /= n_p;

        // Subtract mean from each value in probeset to get residual
        for (unsigned int p = 0; p < n_p; ++p)
            p_res[p_elem[p] - 1] = p_vals[p_elem[p] - 1] - mean;

        UNPROTECT (1);
    }

    UNPROTECT (1);
    return res;
}
