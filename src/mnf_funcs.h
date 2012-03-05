// =====================================================================================
//                              -*- Mode: C -*- 
// mnf_funcs.h
// Copyright 2010 Laboratoire de Bioinformatique Fonctionnelle et Structurale,
//                Institut de recherche en immunologie et en cancerologie (IRIC),
//                Universite de Montreal.
// Authors:       Philippe Serhal (philippe.serhal@umontreal.ca)
//                Sebastien Lemieux (s.lemieux@umontreal.ca)
// =====================================================================================

#ifndef _MNF_FUNCS_H
#define _MNF_FUNCS_H

#include <Rdefines.h>

typedef struct {
    int *g;
    int w, h, n;
} Grid;

typedef struct {
    int *a;
    int n;
} Array;

enum SummaryStat {MEAN, MEDIAN, MIN, MAX};

extern "C" {
    void grid_neighbours (int *, int *, int *, int *, int *);
    Grid generate_grid (int, int *, int *);
    void map_to_grid (int *, double *, int *, int *, double *);
    inline int grid_get (Grid *, int, int);
    inline void grid_set (Grid *, int, int, int);

    void array_neighbours (int *, int *, int *, int *);
    Array generate_array (int, int *);

    void map_values (int *, int *, double *, double *);

    SEXP affy_residuals (SEXP, SEXP, SEXP);
}

#endif

