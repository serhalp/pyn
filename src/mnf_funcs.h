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

typedef struct {
    int *g;
    int w, h, n;
} Grid;

void grid_neighbours (int *, int *, int *, int *, int *);
Grid generate_grid (int, int *, int *);
inline int grid_get (Grid *, int, int);
inline void grid_set (Grid *, int, int, int);

#endif

