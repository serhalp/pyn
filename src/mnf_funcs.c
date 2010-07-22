// =====================================================================================
//                              -*- Mode: C -*- 
// mnf_funcs.c
// Copyright 2010 Laboratoire de Bioinformatique Fonctionnelle et Structurale,
//                Institut de recherche en immunologie et en cancerologie (IRIC),
//                Universite de Montreal.
// Authors:       Philippe Serhal (philippe.serhal@umontreal.ca)
//                Sebastien Lemieux (s.lemieux@umontreal.ca)
// =====================================================================================

#include <R.h>
#include <Rinternals.h>
#include "mnf_funcs.h"

void grid_neighbours (int *n, int *x, int *y, int *k, int *neighbours) {
    Grid grid = generate_grid (*n, x, y);

    int dx, dy, nx, ny, p;
    double dist;
    int offset = 0;
    for (int i = 1; i <= *n; ++i, offset += *k) {
        int thisX = x[i - 1], // from R, 1-offset indices to C, 0-offset
            thisY = y[i - 1];
        int ni = 0;

        for (int d = 1; ni < *k; ++d) {
            for (int dk = 0; dk < 2 * d; ++dk) {
                dx = -d + dk;
                dy = -d;
                nx = thisX + dx;
                ny = thisY + dy;
                dist = (dx * dx) + (dy * dy);
                if (nx > 0 && ny > 0 && nx <= grid.w && ny <= grid.h && (p = grid_get (&grid, nx, ny)) != -1)
                    neighbours[offset + ni++] = p;
     
                dx = d;
                dy = -d + dk;
                nx = thisX + dx;
                ny = thisY + dy;
                dist = (dx * dx) + (dy * dy);
                if (nx > 0 && ny > 0 && nx <= grid.w && ny <= grid.h && (p = grid_get (&grid, nx, ny)) != -1 && ni < *k)
                    neighbours[offset + ni++] = p;

                dx = d - dk;
                dy = d;
                nx = thisX + dx;
                ny = thisY + dy;
                dist = (dx * dx) + (dy * dy);
                if (nx > 0 && ny > 0 && nx <= grid.w && ny <= grid.h && (p = grid_get (&grid, nx, ny)) != -1 && ni < *k)
                    neighbours[offset + ni++] = p;

                dx = -d;
                dy = d - dk;
                nx = thisX + dx;
                ny = thisY + dy;
                dist = (dx * dx) + (dy * dy);
                if (nx > 0 && ny > 0 && nx <= grid.w && ny <= grid.h && (p = grid_get (&grid, nx, ny)) != -1 && ni < *k)
                    neighbours[offset + ni++] = p;
            }
        }
    }
}

Grid generate_grid (int n, int *x, int *y) {
    Grid grid;
    int i;
    int xMax = 0, yMax = 0;

    // Find grid width and height
    for (i = 0; i < n; ++i) {
        if (x[i] > xMax)
            xMax = x[i];
        if (y[i] > yMax)
            yMax = y[i];
    }

    grid.w = xMax;
    grid.h = yMax;
    grid.n = xMax * yMax;
    grid.g = (int *) R_alloc (grid.n, sizeof (int));

    // Initialize all cells to -1
    for (i = 0; i < grid.n; ++i)
        grid.g[i] = -1;

    // Set each probe's cell to its (R, 1-offset) index
    for (i = 0; i < n; ++i)
        grid_set (&grid, x[i], y[i], i + 1);

    return grid;
}

inline int grid_get (Grid *grid, int i, int j) {
    return grid->g[(i - 1) + (j - 1) * grid->w];
}

inline void grid_set (Grid *grid, int i, int j, int val) {
    grid->g[(i - 1) + (j - 1) * grid->w] = val;
}

