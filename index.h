
#ifndef __INDEX_H__
#define __INDEX_H__

#include "pocketsphinx.h"

/** 
 * hit_t
 */
typedef struct hit_s hit_t;

/** 
 * inverted_index_t 
 */
typedef struct inverted_index_s inverted_index_t;

/**
 * result_list_t
 */
typedef struct result_list_s result_list_t;

/**
 * function: inverted_index_init();
 * Create and Initialize a primitive inverted_index from a file.
 */
inverted_index_t* inverted_index_init(const char* filename);

/**
 * function: inverted_index_free()
 * free the memory of a inverted_index
 */
void inverted_index_free(inverted_index_t* index);

/**
 * function: inverted_index_read()
 * Construct a inverted_index from file
 */
inverted_index_t* inverted_index_read(const char* filename);
 
/**
 * function: inverted_index_write()
 * Save a inverted_index to a file
 */
int inverted_index_write(inverted_index_t* index, const char* filename);

/**
 * function: inverted_index_get_wid()
 * return the index number of **word** in the word_list
 */
int inverted_index_get_wid(inverted_index_t* index, const char* word);

/**
 * function: inverted_index_addhits()
 * Add new hits from a lattice
 */
void inverted_index_addhits(inverted_index_t* index, const char* uttid, ps_lattice_t* lat, float32 ascale);


/**
 * function: inverted_index_search()
 * 
 */ 
void inverted_index_search(inverted_index_t* index, ngram_model_t* lm, float32 ascale, char** terms, int n_term, result_list_t** rl);

#endif
