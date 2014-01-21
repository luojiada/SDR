#include <stdio.h>
#include <string.h>
#include "index.h"

#define SENSCR_SHIFT 10

#define WORD_MAX_LENGTH 15
#define MAX_LINE_LENGTH 256
#define INTERVAL 0.3 /*  */
/** 
 * hit_t
 */
struct hit_s {
    char *uttid;   /** uttarence id */
    int32 norm;  /** utterance normalizer */
   
    int wid;
    char* word;
    char* subseq_word;
    
    int from_id, to_id;
    
    double start_time;   /** start time in seconds */ 
    double end_time;    /** end time in seconds */
    int32 alpha;  /** forward likelihood */
    int32 beta;   /** backward likelihood */
    int32 ascr;   /** acoustic model score */

    struct hit_s *next;   /** pointer to next hit */
};

/** 
 * inverted_index_t 
 */
struct inverted_index_s {
    int n_word; /** total number of real words in dictionary, here word refers to syllable */
    char** word_list;  /** word list */
    hit_t** first_hits; /** each element inside points to the first hit of that WORD*/
    hit_t** last_hits;  /** each element inside points to the last hit of that WORD*/
};


/**
 * partial_path_t
 */
typedef struct partial_path_s {
    hit_t* first_term;
    hit_t* last_term;
    int n_term;
    struct partial_path_s *next;
} partial_path_t;

/**
 * queue_t
 *
 */
typedef struct path_queue_s {
    partial_path_t* head;
    partial_path_t* tail;
    int n_path;
} path_queue_t;


/**
 * result_t
 */
typedef struct result_s {
    char* uttid;
    int32 similarity;
    struct result_s* next;
} result_t;

struct result_list_s {
    int n_result;
    result_t *first;
};

void result_list_free(result_list_t* rl)
{
    if (!rl)
        return;
    result_t* r;
    for (r = rl->first; r; r = r->next) {
        free(r->uttid);
        free(r);
        rl->n_result--;
    }
    
    free(rl);
}

/* =====================================================================
 * partial_path_t's function definitions 
 * ===================================================================== */

/** Initialize an empty path */ 
partial_path_t* partial_path_init() {
    partial_path_t* p;
    p = (partial_path_t*) malloc( sizeof(partial_path_t) );
    p->first_term = NULL;
    p->last_term = NULL;
    p->next = NULL;
    p->n_term = 0;  
    return p;
}

/** Finitialize a path */ 
void partial_path_free(partial_path_t* p) 
{
    if (p) {
        hit_t* h;
        while (p->first_term) {
            h = p->first_term;
            p->first_term = p->first_term->next;
            /** Don't forget to free the **uttid**,**word** and **subseq_word** member in hit_t */
            free(h->uttid); 
            free(h->word);
            free(h->subseq_word);
            free(h);
            p->n_term--;
        }
        free(p);
    }
}

/** Extend the path with a new term */
int partial_path_extend(partial_path_t* p, hit_t* h) 
{
    if (!p) 
        return -1;
    if (!h)
        return -1;

    hit_t* hit;
    hit = (hit_t*) malloc(sizeof(hit_t));
    
    hit->uttid = (char*) malloc(MAX_LINE_LENGTH * sizeof(char));
    memset(hit->uttid, 0, MAX_LINE_LENGTH * sizeof(char));
    strncpy(hit->uttid, h->uttid, strlen(h->uttid));
    
    hit->norm = h->norm;
    
    hit->wid = h->wid;
    hit->word = (char*) malloc( (WORD_MAX_LENGTH+1) * sizeof(char));
    memset(hit->word, 0, (WORD_MAX_LENGTH+1) * sizeof(char));
    strncpy(hit->word, h->word, strlen(h->word));
    
    hit->subseq_word = (char*) malloc( (WORD_MAX_LENGTH+1) * sizeof(char));
    memset(hit->subseq_word, 0, (WORD_MAX_LENGTH+1) * sizeof(char));
    strncpy(hit->subseq_word, h->subseq_word, strlen(h->subseq_word));
    
    hit->start_time = h->start_time;
    hit->end_time = h->end_time;
    hit->alpha = h->alpha;
    hit->beta = h->beta;
    hit->ascr = h->ascr;
    hit->next = NULL;
    
    if (p->n_term == 0) {
        p->first_term = hit;
    } else {
        p->last_term->next = hit;
    }
    p->last_term = hit;
    p->n_term++;
    
    return 0;
    
}

/** Genarete a new path from an existing one */
partial_path_t* partial_path_copy(partial_path_t* p) {
       if (!p) {
           return NULL;
       }
       
       hit_t* h;
       partial_path_t* out = partial_path_init();
       for (h = p->first_term; h; h = h->next) {
            if (partial_path_extend(out, h) != 0) {
                perror("Add hits to path failure");
                partial_path_free(out);
                return NULL;
            }
       }
       out->next = NULL;
       return out;
}


/** Get posterior log-likelihood of the partial path: P(path|O) */
int32 partial_path_get_posterior(partial_path_t* p, ngram_model_t* lm, float32 ascale)
{
    hit_t *hit, *subseq_hit;
    int32 n_used;
    int32 result;
    if (p) {
        // p->n_term >= 0
        switch(p->n_term) {
            case 0:
                return 0;
            case 1:
                return p->first_term->alpha + p->first_term->beta - p->first_term->norm;
            default:
                result = p->first_term->alpha + p->last_term->beta - p->first_term->norm;
                for (hit = p->first_term, subseq_hit = p->first_term->next;
                        hit && subseq_hit; 
                        hit = hit->next, subseq_hit = subseq_hit->next) {
                    result = result 
                            + ngram_score_to_prob(lm, ngram_bg_score(lm, ngram_wid(lm, subseq_hit->word), ngram_wid(lm, hit->word), &n_used))
                            + (subseq_hit->ascr << SENSCR_SHIFT) * ascale;         
                }
                return result;
        }   
    }
    return 0;
}


void partial_path_print(partial_path_t* p, ngram_model_t* lm, float32 ascale)
{
    if (!p)
        return;
    hit_t* h;
    printf("[%.2f-%.2f] ", p->first_term->start_time, p->last_term->end_time);
    for (h = p->first_term; h; h = h->next) {
        printf("%s ", h->word);
    }
    printf("%d\n", partial_path_get_posterior(p, lm, ascale));
}

/* =====================================================================
 * path_queue_t's function definitions 
 * ===================================================================== */ 
 
path_queue_t* path_queue_init() 
{
    path_queue_t* q = (path_queue_t*) malloc( sizeof(path_queue_t) );
    q->head = q->tail = NULL;
    q->n_path = 0;
    return q;
}

void path_queue_free(path_queue_t* q)
{
    if (!q) return;
    partial_path_t* p;
    while (q->head) {
        p = q->head;
        q->head = q->head->next;
        partial_path_free(p);
        q->n_path--;
    }
    free(q);    
}

void path_queue_print(path_queue_t* q, ngram_model_t* lm, float32 ascale) {
    if(!q)
        return;
    partial_path_t* p = q->head;
    while (p) {
        partial_path_print(p, lm, ascale);
        p = p->next;
    }
}

void path_queue_add(path_queue_t* q, partial_path_t* p) {
    if (!q) {
        perror("No queue object");
        return;
    }
    if (!p) {
        perror("No path object");
        return;
    }
    /* Next Path in p Should be NULL, becoz p is newest. */
    if (p->next != NULL) {
        perror("Next Path Should be NULL");
        return;
    }
    if (NULL == q-> head && NULL == q->tail) {  // empty queue
        q->head = q->tail = p;
    } else { 
        // here we should make sure head and tail point to the right memory
        // head -> [path]-> [path] -> NULL
        //                     ^
        //                     |
        //                    tail
        q->tail->next = p;
        q->tail = p;
    }
    q->n_path++;
}
/* =====================================================================
 * inverted_index's functons
 * ===================================================================== */
 
inverted_index_t* inverted_index_init(const char* filename)
{
    int i;
    FILE* fh;
    char s[WORD_MAX_LENGTH + 2] = {'\0',};
    
    fh = fopen(filename, "r");    
    if (fh == NULL) {
        perror("Failed to open Word List file.");
        return NULL;
    }
    
    inverted_index_t* index = (inverted_index_t*) malloc(sizeof(inverted_index_t));
    index->n_word = 0;
    
    while ((fgets(s, WORD_MAX_LENGTH + 2, fh) != '\0')) {
        index->n_word++;
	    //printf("%d: %s", index->n_word, s); 
    }
    index->word_list = (char**) malloc(index->n_word * sizeof(char*));
    for (i = 0; i < index->n_word; i++) {
        index->word_list[i] = (char*) malloc( (WORD_MAX_LENGTH + 1) * sizeof(char));
        memset(index->word_list[i], 0, (WORD_MAX_LENGTH + 1) * sizeof(char));
    }
    fseek(fh, 0, SEEK_SET);
    i = 0;
    while ((fgets(s, WORD_MAX_LENGTH + 2, fh) != '\0')) {
        strncpy(index->word_list[i], strtok(s, "\n"), WORD_MAX_LENGTH);
	    //printf("%s", index->word_list[i]);	    
	    i++;
    }
    fclose(fh);
    
    /** allocate memory for index */
    index->first_hits = (hit_t**) malloc(index->n_word * sizeof(hit_t*));
    for (i = 0; i < index->n_word; i++) {
        index->first_hits[i] = NULL;
    }
    index->last_hits = (hit_t**) malloc(index->n_word * sizeof(hit_t*));
    for (i = 0; i < index->n_word; i++) {
        index->last_hits[i] = NULL;
    }
    return index;
}

void inverted_index_free(inverted_index_t* index)
{
    int i;
    hit_t* p;
    
    for(i = 0; i < index->n_word; i++) {
        free(index->word_list[i]);
    }
        
    for(i = 0; i < index->n_word; i++) {
        if (index->first_hits[i] == NULL && index->last_hits[i] == NULL)
            continue;   // no hits in this WORD, skip to next WORD
        // delete each hit from first hit to last
        while(index->first_hits[i]) {
            p = index->first_hits[i];
            index->first_hits[i] = index->first_hits[i]->next;
            free(p->uttid);
            free(p->word);
            free(p->subseq_word);
            free(p);    
        }
        index->last_hits[i] = NULL;
    }
    
    free(index->word_list);
    free(index->first_hits);
    free(index->last_hits);
    free(index);
    printf("Finialize index Successfully\n");
}

int inverted_index_get_wid(inverted_index_t* index, const char* word)
{
    int i;
    int len = strlen(word);
    
    for (i = 0; i < index->n_word; i++) {
        if ( len == strlen(index->word_list[i]) 
            && strncmp(word, index->word_list[i], len ) == 0 ) 
            break;
    }
    
    if ( i < index->n_word ) 
        return i;
    else
        return -1; /**  **word** not found in the word_list */

}

int inverted_index_write(inverted_index_t* index, const char* filename)
{
    FILE* fp;
    int i;
    hit_t* hit;
    
    if ( (fp = fopen(filename, "w")) == NULL) {
        perror("Failed to open file");
        return -1;
    }
    
    //fprintf(fp, "# Index by Jiada\n");
    fprintf(fp, "# Words: %d\n\n", index->n_word);
    
    for (i = 0; i < index->n_word; i++) {
        fprintf(fp, "%d:%s\n", i, index->word_list[i]);
        for (hit = index->first_hits[i]; hit; hit = hit->next) {
            fprintf(fp, "(%s, %.2f, %.2f, %d, %d, %d, %d, %d, %d, %s)\n",
                hit->uttid, hit->start_time, hit->end_time,
                hit->ascr, hit->alpha, hit->beta, hit->norm, hit->from_id, hit->to_id, hit->subseq_word);
        }
    }
    
    fclose(fp);
    return 0;
}

inverted_index_t* inverted_index_read(const char* filename)
{ 
    FILE* fp;
    int i, k;
    int n_word;
    char line[MAX_LINE_LENGTH] = {'\0',}; 
    inverted_index_t* index;
    
    int wid, from_id, to_id;
    char word[MAX_LINE_LENGTH] = {'\0',};
    char subseq_word[MAX_LINE_LENGTH] = {'\0',};
    char uttid[MAX_LINE_LENGTH] = {'\0',};
    float st, et;
    int32 norm, ascr, alpha, beta;
    
    hit_t* hit;
    
    if ( (fp = fopen(filename, "r")) == NULL) {
        perror("Failed to open file.");
        return NULL;
    }
    /** Get total number of words */
    fscanf(fp, "# Words: %d\n\n", &n_word);
    //printf("words: %d\n", n_word);
    /** allocate space to store index*/
    index = (inverted_index_t*) malloc(sizeof(inverted_index_t));
    index->n_word = n_word;
    /** word list */
    index->word_list = (char**) malloc(index->n_word * sizeof(char*));
    for (i = 0; i < index->n_word; i++) {
        index->word_list[i] = (char*) malloc( (WORD_MAX_LENGTH + 1) * sizeof(char));
        memset(index->word_list[i], 0, (WORD_MAX_LENGTH + 1) * sizeof(char));
    }
    /** allocate memory for hits */
    index->first_hits = (hit_t**) malloc(index->n_word * sizeof(hit_t*));
    for (i = 0; i < index->n_word; i++) {
        index->first_hits[i] = NULL;
    }
    index->last_hits = (hit_t**) malloc(index->n_word * sizeof(hit_t*));
    for (i = 0; i < index->n_word; i++) {
        index->last_hits[i] = NULL;
    }
    
    while ( NULL != fgets(line, MAX_LINE_LENGTH, fp)) {
        if ( ( (k = sscanf(line, "%d:%s\n", &wid, word)) != 2) 
            && ( (k = sscanf(line, "(%[a-z], %f, %f, %d, %d, %d, %d, %d, %d, %s)\n",
                        uttid, &st, &et, &ascr, &alpha, &beta, &norm, &from_id, &to_id, subseq_word)) != 10) ) 
        {
            //printf("k=%d %s", k, line);
            perror("Format Error");
            inverted_index_free(index);
            return NULL;    
        }
        
        if ( k == 2) {
            strncpy(index->word_list[wid], word, WORD_MAX_LENGTH);
            //printf("%s\n", index->word_list[wid]);           
        }
        
        if ( k == 10) {
            //printf("(%d:%s %s, %.2f, %.2f, %d, %d, %d)\n", wid, word, uttid, st, et, ascr, alpha, beta); 
            
            // add a new hit 
            
            hit = (hit_t*) malloc(sizeof(hit_t));
            hit->uttid = (char*) malloc(MAX_LINE_LENGTH * sizeof(char));
            memset(hit->uttid, 0, MAX_LINE_LENGTH * sizeof(char));
            strncpy(hit->uttid, uttid, strlen(uttid));
            
            hit->norm = norm;
            
            hit->wid = wid;
            hit->word = (char*) malloc( (WORD_MAX_LENGTH+1) * sizeof(char));
            memset(hit->word, 0, (WORD_MAX_LENGTH+1) * sizeof(char));
            strncpy(hit->word, word, strlen(word));
            
            hit->from_id = from_id;
            hit->to_id = to_id;
            
            hit->subseq_word = (char*) malloc( (WORD_MAX_LENGTH+1) * sizeof(char));
            memset(hit->subseq_word, 0, (WORD_MAX_LENGTH+1) * sizeof(char));
            strncpy(hit->subseq_word, subseq_word, strlen(subseq_word));
            
            hit->start_time = (double) st;
            hit->end_time = (double) et;
            hit->alpha = alpha;
            hit->beta = beta;
            hit->ascr = ascr;
            hit->next = NULL;
            if (index->first_hits[wid] == NULL) {
                index->first_hits[wid] = hit;
            } else {
                index->last_hits[wid]->next = hit;
            }
        
            index->last_hits[wid] = hit;            
        }
        
    }

    fclose(fp);
    return index;    
}

/** Add new hits from a lattice */
void inverted_index_addhits(inverted_index_t* index, const char* uttid, ps_lattice_t* lat, float32 ascale)
{
    int wid;
    int32 norm;
    int16 sf, ef;
    int32 nfrate;
    int32 ascr, alpha, beta;
    const char* word;
    const char* subseq_word;
    ps_latnode_iter_t* node_iter;
    ps_latlink_iter_t* link_iter;
    ps_latnode_t *d, *to;
    ps_latlink_t* link;
    
    hit_t* hit;
    
    norm = ps_lattice_get_norm(lat);
    
    // Traverse all edges in the lattice to add new hits
    for (node_iter = ps_latnode_iter(lat); node_iter; node_iter = ps_latnode_iter_next(node_iter)) {
        d = ps_latnode_iter_node(node_iter); 
        if (!ps_latnode_reachable(d))
            continue;
        for (link_iter = ps_latnode_exits(d); link_iter; link_iter = ps_latlink_iter_next(link_iter)) {
        
            link = ps_latlink_iter_link(link_iter);
            to = ps_latlink_nodes(link, NULL);
            if ( to == NULL || !ps_latnode_reachable(to) )
                continue;
            if ( (ps_latlink_get_ascr(link) < (int)0xE0000000 ) || (ps_latlink_get_ascr(link) > 0))
                continue; 
        
            /** Extract infomation from each link */
    	    word = ps_latlink_word(lat, link);
    	    wid = inverted_index_get_wid(index, word);
    	    if (wid < 0) {
    	        continue;
    	    }
    	    
    	    subseq_word = ps_latnode_word(lat, ps_latlink_nodes(link, NULL));
    	    ef = ps_latlink_times(link, &sf); 
    	    nfrate = ps_lattice_get_frate(lat);
    	    ascr = (ps_latlink_get_ascr(link) << SENSCR_SHIFT) * ascale;
    	    alpha = ps_latlink_get_alpha(link);
    	    beta = ps_latlink_get_beta(link);
    	    
    	    //printf("%d: %s st:%.2f et:%.2f ascr:%d alpha:%d beta:%d\n", wid, word, (double) sf/nfrate, (double) ef/nfrate, ascr, alpha, beta);
    	    // add a new hit 
    	    hit = (hit_t*) malloc(sizeof(hit_t));
    	    hit->uttid = (char*) malloc(MAX_LINE_LENGTH * sizeof(char));
    	    memset(hit->uttid, 0, MAX_LINE_LENGTH * sizeof(char));
    	    strncpy(hit->uttid, uttid, strlen(uttid));
    	    
    	    hit->norm = norm;
    	    
    	    hit->wid = wid;
            hit->word = (char*) malloc( (WORD_MAX_LENGTH+1) * sizeof(char));
            memset(hit->word, 0, (WORD_MAX_LENGTH+1) * sizeof(char));
            strncpy(hit->word, word, strlen(word));
    	    
    	    hit->subseq_word = (char*) malloc( (WORD_MAX_LENGTH+1) * sizeof(char));
            memset(hit->subseq_word, 0, (WORD_MAX_LENGTH+1) * sizeof(char));
            strncpy(hit->subseq_word, subseq_word, strlen(subseq_word));
    	    
    	    hit->from_id = ps_latnode_get_id(d);
    	    
    	    hit->to_id = ps_latnode_get_id(to);
    	    
    	    hit->start_time = (double) sf / nfrate;
    	    hit->end_time = (double) ef / nfrate;
    	    hit->alpha = alpha;
    	    hit->beta = beta;
    	    hit->ascr = ascr;
    	    hit->next = NULL;
    	    
    	    if (index->first_hits[wid] == NULL) {
    	        index->first_hits[wid] = hit;
    	    } else {
    	        index->last_hits[wid]->next = hit;
    	    }
    	    
    	    index->last_hits[wid] = hit;
    	    
    	}
    }
}  


/**
 * function: inverted_index_search()
 * Return id of utterances in which all query terms are matched 
 */ 
void inverted_index_search(inverted_index_t* index, ngram_model_t* lm, char** terms, int n_term, result_list_t** rl)
{
    if (!index) {
        perror("Index not found");
        return;
    }
    
    if (!lm) {
        perror("lm not found");
        return;
    }

    if (!terms) {
        perror("no query terms");
        return;
    }
    int i, k;
    int rv;
    int wid;
    hit_t* hit;
    partial_path_t *p, *q;
    path_queue_t** queues;
    
    *rl = NULL;
    queues =  (path_queue_t**) malloc( n_term * sizeof(path_queue_t*) );
    for (i = 0; i < n_term; i++) {
        queues[i] = path_queue_init();
    }
    /** Seach candidate partial pathes which match all query terms */
    for (k = 0; k < n_term; k++) {
        if ( (wid = inverted_index_get_wid(index, terms[k])) != -1) {
            hit = index->first_hits[wid];
            if (!hit) {
                perror("No hits on current query term");
                break;
            }
            while (hit) {
                if (k == 0) { /** first query term */
                    q = partial_path_init();  
                    rv = partial_path_extend(q, hit);
                    if ( 0 != rv ) {
                        perror("Error when adding hit to path, skip it");
                        partial_path_free(q);
                        continue;
                    }
                    path_queue_add(queues[k], q);              
                                                          
                } else {
                    if (!queues[k-1]->head) {
                        fprintf(stderr, "No hits on previous term k:%d\n", k);
                        goto exit;
                    }
                    for (p = queues[k-1]->head; p; p = p->next) {
                        if ( (strncmp(p->first_term->uttid, hit->uttid, strlen(hit->uttid)) == 0)
                             && ( ( (p->last_term->end_time) <= hit->start_time) && hit->start_time <= (p->last_term->end_time + INTERVAL)) ) {
                               q = partial_path_copy(p);
                               if ( 0 != partial_path_extend(q, hit) ) {
                                   perror("Error when adding hit to path, skip it");
                                   partial_path_free(q);
                                   continue;
                               }
                               path_queue_add(queues[k], q);  
                        }
                    }                
                }               
                hit = hit->next;
            }
        path_queue_print(queues[k], lm, 1.0/20.0);
        } else {
            perror("Query term not found inside the index");
            break;
        }
    }
    /** */
    if (k == n_term ) {
        k = n_term - 1;
        if (queues[k]->n_path > 0) { /** candidate path exists*/
            /** Compare similiarity between query terms and utterances */
            int32 norm;
            printf("#result: %d\n", queues[k]->n_path);
            /*
            for (i = 0, p = queues[k]->head; (i < queues[k]->n_path) && p; i++, p= p->next) {
                norm = partial_path_get_posterior(p, lm, 1.0/20.0);
                printf("(%s, %.2f, %.2f, %d)\n", p->first_term->uttid, p->first_term->start_time, p->last_term->end_time, norm);
            }
            */
            
        }
    }
    
  
exit:      
    for (i = 0; i < n_term; i++) {
        path_queue_free(queues[i]);
    }
    free(queues);   
}

