// =====================================================================================
//                              -*- Mode: C -*- 
// pyn.cpp
// Copyright 2010 Laboratoire de Bioinformatique Fonctionnelle et Structurale,
//                Institut de recherche en immunologie et en cancerologie (IRIC),
//                Universite de Montreal.
// Authors:       Philippe Serhal (philippe.serhal@umontreal.ca)
//                Sebastien Lemieux (s.lemieux@umontreal.ca)
// =====================================================================================

#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <algorithm>
#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include "pyn.h"

void grid_neighbours (int *n, int *x, int *y, int *k, int *neighbours) {
    Grid grid = generate_grid (*n, x, y);

    int dx, dy, nx, ny, p;
    double dist;
    int offset = 0;
    for (int i = 1; i <= *n; ++i, offset += *k) {
        int thisX = x[i - 1], // from R, 1-offset indices to C, 0-offset... and back to R, 1-offset
            thisY = y[i - 1];
        int ni = 0;

        if (thisX == NA_INTEGER)
            continue;

        for (int d = 1; ni < *k; ++d) {
            for (int dk = 0; dk < 2 * d; ++dk) {
                dx = -d + dk;
                dy = -d;
                nx = thisX + dx;
                ny = thisY + dy;
                dist = (dx * dx) + (dy * dy);
                if (ni < *k && nx > 0 && ny > 0 && nx <= grid.w && ny <= grid.h && (p = grid_get (&grid, nx, ny)) != -1)
                    neighbours[offset + ni++] = p;
     
                dx = d;
                dy = -d + dk;
                nx = thisX + dx;
                ny = thisY + dy;
                dist = (dx * dx) + (dy * dy);
                if (ni < *k && nx > 0 && ny > 0 && nx <= grid.w && ny <= grid.h && (p = grid_get (&grid, nx, ny)) != -1)
                    neighbours[offset + ni++] = p;

                dx = d - dk;
                dy = d;
                nx = thisX + dx;
                ny = thisY + dy;
                dist = (dx * dx) + (dy * dy);
                if (ni < *k && nx > 0 && ny > 0 && nx <= grid.w && ny <= grid.h && (p = grid_get (&grid, nx, ny)) != -1)
                    neighbours[offset + ni++] = p;

                dx = -d;
                dy = d - dk;
                nx = thisX + dx;
                ny = thisY + dy;
                dist = (dx * dx) + (dy * dy);
                if (ni < *k && nx > 0 && ny > 0 && nx <= grid.w && ny <= grid.h && (p = grid_get (&grid, nx, ny)) != -1)
                    neighbours[offset + ni++] = p;
            }
        }
    }
}

Grid generate_grid (int n, int *x, int *y) {
    // Find grid width and height
    int xMax = 0, yMax = 0;
    for (int i = 0; i < n; ++i) {
        if (x[i] > xMax)
            xMax = x[i];
        if (y[i] > yMax)
            yMax = y[i];
    }

    Grid grid;
    grid.w = xMax;
    grid.h = yMax;
    grid.n = (xMax + 1) * (yMax + 1);
    grid.g = (int *) R_alloc (grid.n, sizeof (int));

    // Initialize all cells to -1
    for (int i = 0; i < grid.n; ++i)
        grid.g[i] = -1;

    // Set each probe's cell to its (R, 1-offset) index
    for (int i = 1; i <= n; ++i)
        if (x[i - 1] != NA_INTEGER)
            grid_set (&grid, x[i - 1], y[i - 1], i);

    return grid;
}

/* This is pretty tricky: you need to figure out in advance the correct length
 * for the output matrix (actually, vector) 'g'; that is max(x) * max(y).
 * TODO: Maybe eventually rewrite this with the .Call interface to avoid this. */
void map_to_grid (int *n, double *vals, int *x, int *y, double *g) {
    int w = 0; /* max(x), i.e. width of output grid */
    for (int i = 0; i < *n; ++i)
        if (x[i] > w)
            w = x[i];
    ++w;

    for (int i = 0; i < *n; ++i)
        g[x[i] + y[i] * w] = vals[i];
}

inline int grid_get (Grid *grid, int i, int j) {
    return grid->g[i + j * grid->w];
}

inline void grid_set (Grid *grid, int i, int j, int val) {
    grid->g[i + j * grid->w] = val;
}

void array_neighbours (int *n, int *pos, int *k, int *neighbours) {
    Array array = generate_array (*n, pos);

    srand (time (NULL));
    int offset = 0;
    for (int i = 1; i <= *n; ++i, offset += *k) {
        int thisX = pos[i - 1]; // from R, 1-offset indices to C, 0-offset... and back to R, 1-offset
        int ni = 0;
        int l = 0, r = 0, rnd = rand () % 2; // Randomly select left or right to explore first
        int p, pi;
        while (ni < *k) {
            if ((l + r + rnd) % 2 == 0) {
                if ((p = thisX - ++l) >= 0 && (pi = array.a[p]) != -1)
                    neighbours[offset + ni++] = pi;
            } else if ((p = thisX + ++r) < array.n && (pi = array.a[p]) != -1)
                neighbours[offset + ni++] = pi;
        }
    }
}

Array generate_array (int n, int *pos) {
    // Find virtual array dimension
    int xMax = 0;
    for (int i = 0; i < n; ++i)
        if (pos[i] > xMax)
            xMax = pos[i];

    Array array;
    array.n = xMax + 1;
    array.a = (int *) R_alloc (array.n, sizeof (int));

    // Initialize all cells to -1
    for (int x = 0; x < array.n; ++x)
        array.a[x] = -1;

    // Set each probe's cell to its (R, 1-offset) index
    for (int i = 1; i <= n; ++i)
        array.a[pos[i - 1]] = i;

    return array;
}

void map_values (int *n, int *indices, double *values, double *mapped) {
    for (int i = 0; i < *n; ++i)
        if (indices[i] != NA_INTEGER)
            mapped[i] = values[indices[i] - 1]; // 1-offset to 0-offset
}

SEXP affy_residuals (SEXP sets, SEXP vals, SEXP _summary_stat) {
    int summary_stat = *INTEGER_POINTER (_summary_stat);
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

        double res_estimate = 0.0;
        switch (summary_stat) {
            case MEAN:
                for (unsigned int p = 0; p < n_p; ++p)
                    res_estimate += p_vals[p_elem[p] - 1];
                res_estimate /= n_p;
                break;
            case MEDIAN:
                {
                    double *p_ps_res = (double *) R_alloc (n_p, sizeof (double));
                    for (unsigned int p = 0; p < n_p; ++p)
                        p_ps_res[p] = p_vals[p_elem[p] - 1];
                    std::nth_element (p_ps_res, p_ps_res + (int) floor ((n_p - 1) / 2), p_ps_res + n_p);
                    double lower_median = p_ps_res[(int) floor ((n_p - 1) / 2)];
                    std::nth_element (p_ps_res, p_ps_res + (int) ceil ((n_p - 1) / 2), p_ps_res + n_p);
                    double upper_median = p_ps_res[(int) ceil ((n_p - 1) / 2)];
                    res_estimate = (lower_median + upper_median) / 2;
                }
                break;
            case MIN:
                {
                    double *p_ps_res = (double *) R_alloc (n_p, sizeof (double));
                    for (unsigned int p = 0; p < n_p; ++p)
                        p_ps_res[p] = p_vals[p_elem[p] - 1];
                    res_estimate = *std::min_element (p_ps_res, p_ps_res + n_p);
                }
                break;
            case MAX:
                {
                    double *p_ps_res = (double *) R_alloc (n_p, sizeof (double));
                    for (unsigned int p = 0; p < n_p; ++p)
                        p_ps_res[p] = p_vals[p_elem[p] - 1];
                    res_estimate = *std::max_element (p_ps_res, p_ps_res + n_p);
                }
                break;
        }

        // Subtract res. estimate from each value in probeset to get residual
        for (unsigned int p = 0; p < n_p; ++p)
            p_res[p_elem[p] - 1] = p_vals[p_elem[p] - 1] - res_estimate;

        UNPROTECT (1);
    }

    UNPROTECT (1);
    return res;
}

