/*************************************************************************************************
 * sausage.h
 * sausage-like confusion network is an alternative representation of word lattice, and it merges 
 * those edges into equivalence classes with orders consistent with the original lattice. 
 * Author: Jiada 
 *
 *************************************************************************************************/
#ifndef __SAUSAGE_H__
#define __SAUSAGE_H__

#include "pocketsphinx.h"

/**
 * node_t
 */
typedef struct node_s node_t;
/***
 * node_set_t
 */
typedef struct node_set_s node_set_t;
/**
 * edge_t
 */
typedef struct edge_s edge_t;
/***
 * edge_set_t
 */
typedef struct edge_set_s edge_set_t;
/***
 * sausage_t
 */
typedef struct sausage_s sausage_t;

/**
 * function: convert_lattice_to_sausage()
 * convert incoming lattice to a sausage.
 */
sausage_t* convert_lattice_to_sausage(ps_lattice_t* dag);
/**
 * function: sausage_wirte()
 * sausage_wirte
 */
void sausage_write(sausage_t* s, ps_lattice_t* dag, const char* filename);
void sausage_free(sausage_t* s);
void sausage_last_node_set(sausage_t* s,  ps_lattice_t* dag);

typedef struct lite_node_s lite_node_t;
typedef struct lite_edge_s lite_edge_t;
typedef struct lite_edge_set_s lite_edge_set_t;
typedef struct lite_sausage_s lite_sausage_t;
/**
 * function: sausage_simplify()
 * simplify sausage to a lite one
 */
lite_sausage_t* sausage_simplify(sausage_t* s, ps_lattice_t* dag);
void lite_sausage_write(lite_sausage_t* lite_s, const char* filename);
void lite_sausage_free(lite_sausage_t* lite_s);

typedef struct s_hit_s s_hit_t;
typedef struct s_hits_pos_s s_hits_pos_t;
typedef struct s_hits_word_s s_hits_word_t;
typedef struct dualclue_index_s dualclue_index_t;
//typedef struct dualclue_index_cache_s dualclue_index_cache_t;


dualclue_index_t* dualclue_index_init(const char* filename);
dualclue_index_t* dualclue_index_read(const char* filename);
void dualclue_index_write(dualclue_index_t* index, const char* filename);
void dualclue_index_addhit(dualclue_index_t* index, const char* uttid, lite_sausage_t* lite_s);
void dualclue_index_free(dualclue_index_t* index);
/*
dualclue_index_cache_t* dualclue_index_get_cache(dualclue_index_t* index);*/
void dualclue_index_search(dualclue_index_t* index, char** terms, int n_term);
#endif
