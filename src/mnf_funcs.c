// =====================================================================================
//                              -*- Mode: C -*- 
// mnf_funcs.c
// Copyright 2010 Laboratoire de Bioinformatique Fonctionnelle et Structurale,
//                Institut de recherche en immunologie et en cancerologie (IRIC),
//                Universite de Montreal.
// Authors:       Philippe Serhal (philippe.serhal@umontreal.ca)
//                Sebastien Lemieux (s.lemieux@umontreal.ca)
// =====================================================================================

#include <stdlib.h>
#include <time.h>
#include <R.h>
#include <Rinternals.h>
#include "mnf_funcs.h"

void grid_neighbours (int *n, int *x, int *y, int *k, int *neighbours) {
    Grid grid = generate_grid (*n, x, y);

    int dx, dy, nx, ny, p;
    double dist;
    int offset = 0;
    for (int i = 1; i <= *n; ++i, offset += *k) {
        int thisX = x[i - 1], // from R, 1-offset indices to C, 0-offset... and back to R, 1-offset
            thisY = y[i - 1];
        int ni = 0;

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
    grid.n = xMax * yMax;
    grid.g = (int *) R_alloc (grid.n, sizeof (int));

    // Initialize all cells to -1
    for (int i = 0; i < grid.n; ++i)
        grid.g[i] = -1;

    // Set each probe's cell to its (R, 1-offset) index
    for (int i = 1; i <= n; ++i)
        grid_set (&grid, x[i - 1], y[i - 1], i);

    return grid;
}

inline int grid_get (Grid *grid, int i, int j) {
    return grid->g[(i - 1) + (j - 1) * grid->w];
}

inline void grid_set (Grid *grid, int i, int j, int val) {
    grid->g[(i - 1) + (j - 1) * grid->w] = val;
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

void map_values (int *n, int *indices, int *values, int *mapped) {
    for (int i = 0; i < *n; ++i)
        mapped[i] = values[indices[i] - 1]; // 1-offset to 0-offset
}

