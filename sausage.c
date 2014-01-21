#include "sausage.h"

#define MATCH 0
#define MISMATCH 2
#define INSERT 1
#define DELETE 1

#define WORD_MAX_LENGTH 15
#define MAX_LINE_LENGTH 256

struct node_s {
    ps_latnode_t *node;
    node_t *next;
};

struct node_set_s {
    int ns_id;  /** node_set id */

    node_t* nodes; /** nodes in original lattice contained in this node_set */
    node_t* last_node;
    int n_node; /** total number of nodes inside a node_set */
    
    int initial; /** initial node set (0--false )*/
    int final;  /** final node set 0:false */
    
    int t_max, t_min;
    
    struct edge_set_s *entry, *exit;    /** edges connected to this node_set */
    
    struct node_set_s *prev;
    struct node_set_s *next;    
};

struct edge_s {
    ps_latlink_t *edge;
    edge_t *next;   
};

struct edge_set_s {
    int from_ns_id, to_ns_id;   /** edge_set is between <from_ns_id> and <to_ns_id> */
    node_set_t *from, *to;
    int n_edge; /**total number of edges inside a edge_set */
    edge_t* edges;    /** edges in original lattice contained in this edge_set */
    edge_t* last_edge;
};

struct sausage_s {
    node_set_t* nodesets;
    int n_nodeset;
};




int node_set_add(node_set_t* ns, ps_latnode_t* node)
{
    if (!ns || !node)
        return -1;
        
    node_t *n;
    int t = ps_latnode_times(node, NULL, NULL);
    
    n = (node_t*) malloc(sizeof(node_t));
    n->node = node;
    n->next = NULL;
    
    if (!ns->nodes) { // n-o any nodes in node_set
        ns->nodes = n;
        ns->t_min = ns->t_max = t;
    } else {
        ns->last_node->next = n;
        // Upadte t_max t_min
    
        if ( t < ns->t_min) {
            ns->t_min = t;    
        }
        if ( t > ns->t_max ) {
            ns->t_max = t;
        }
    }
    ns->last_node = n;
    ns->n_node++;
     
    return 0;
}

int edge_set_add(edge_set_t* es, ps_latlink_t* link)
{
    if ( !es || !link)
        return -1;
        
    edge_t *e;
    e = (edge_t*) malloc(sizeof(edge_t));
    e->edge = link;
    e->next = NULL;
    if (!es->edges) {
        es->edges = e;
    } else {
        es->last_edge->next = e;
    }
    es->last_edge = e;
    es->n_edge++;
    return 0;
}

/** To tell whether there is a link between node n and any nodes inside node_set ns
 *  Return 1 means YES, 0 means NO
 */
int node_set_has_connection(node_set_t* ns, ps_latnode_t* node)
{
    /** Go through all entries of latnode n to see if there is a link whose from node'id inside node_set */
    ps_latlink_iter_t *iter;
    ps_latlink_t *l;
    ps_latnode_t *d;
    node_t* n;
    for (iter = ps_latnode_entries(node); iter; iter = ps_latlink_iter_next(iter)) {
        l = ps_latlink_iter_link(iter);
        ps_latlink_nodes(l, &d);
        for (n = ns->nodes; n; n = n->next) {
            if ( d == n->node ) {
                return 1;
            }
        }
    }
    return 0;
}

int node_set_contains(node_set_t* ns, ps_latnode_t* node)
{
    node_t* n;
    for (n = ns->nodes; n; n = n->next) {
            if (node == n->node) {
                return 1;
            }
    }
    return 0;
}

/***
 * get similarity score between two words from two latlink
 * compare minimum edit distance of two words, here word also refers to syllable which represents the pronunciation.
 * supposed  we define the cost of insertion and deletion is  1, the cost of subsitution is 2 
 * 
 */

double similarity_of_words(const char* s, const char* t)
{
    int m, n; 
    int i, j, k;
    int cost_m, cost_i, cost_d;
    int cost;
    m = strlen(s);
    n = strlen(t);
    int** vcost = (int**) calloc(m+1, sizeof(int*));
    for (k = 0; k < m+1; k++) {
        vcost[k] = (int*) calloc(n+1, sizeof(int));
    } 
    
    vcost[0][0] = 0;
    for (i = 1; i < m+1; i++) {
        vcost[i][0] = vcost[i-1][0] + DELETE;
    }
    for (j = 1; j < n+1; j++) {
        vcost[0][j] = vcost[0][j-1] + INSERT;
    }
    for (i = 1; i < m+1; i++) {
        for (j = 1; j < n+1; j++) {
           cost_m = vcost[i-1][j-1]  + ((s[i-1] == t[j-1]) ? MATCH : MISMATCH);
           cost_i = vcost[i][j-1] + INSERT;
           cost_d = vcost[i-1][j] + DELETE;
           int smaller = (cost_m < cost_i) ? cost_m : cost_i;
           vcost[i][j] = (smaller < cost_d) ? smaller : cost_d;     
        }
    }
    
    
   cost = vcost[m][n];
   for (k = 0; k < m+1; k++) {
        free(vcost[k]);
    } 
   free(vcost);
   return  1.0/(cost + 1);
}

double overlap(edge_set_t* es, ps_latlink_t* l) 
{
    int16 sf;
    int ef;
    ef = ps_latlink_times(l, &sf);
    
    int es_sf, es_ef;
    es_sf = es->from->t_max;
    es_ef = es->to->t_min;
     
    int overlap_lft, overlap_rt;
    int lft, rt;
    
    if (ef <= es_sf || sf >= es_ef) {
        return 0;
    }
    
    if (sf <= es_sf ) {
        lft = sf;
        overlap_lft = es_sf;
    } else {
        lft = es_sf;
        overlap_lft = sf;
    }
    
    if (ef <= es_ef ) {
        rt = es_ef;
        overlap_rt = ef;
    } else {
        rt = ef;
        overlap_rt = es_ef;
    }
    if (rt - lft == 0 ) {
        perror("Impossible");
    }
    return (double)(overlap_rt - overlap_lft) / (rt - lft);
}


/**
 * PAPER: IMPROVED CONFUSION NETWORK ALGORITHM AND SHORTEST PATH SEARCH
 * 
 */
sausage_t* convert_lattice_to_sausage(ps_lattice_t* dag)
{
    int i;
    int n_node;
    ps_latnode_iter_t* itor;
    ps_latnode_t* node;
    
    n_node = 0;
    for (itor = ps_latnode_iter(dag); itor; itor = ps_latnode_iter_next(itor)){
        n_node++;
    }
    ps_latnode_t** node_stack;
    node_stack = (ps_latnode_t**) malloc( n_node * sizeof(ps_latnode_t*) );
    for (i = 0, itor = ps_latnode_iter(dag); i < n_node && itor; i++, itor = ps_latnode_iter_next(itor)) {
        node_stack[i] = ps_latnode_iter_node(itor);
    }
    
    sausage_t* sausage;
    sausage = (sausage_t*) malloc( sizeof(sausage) );
    sausage->nodesets = NULL;
    sausage->n_nodeset = 0;
    // Assign initial node n_0 to NS_0
    sausage->nodesets = (node_set_t*) malloc( sizeof(node_set_t) );
    sausage->nodesets->ns_id = 0;
    
    sausage->nodesets->nodes = NULL;
    sausage->nodesets->last_node = NULL;
    sausage->nodesets->n_node = 0;
    
    sausage->nodesets->initial = 1;
    sausage->nodesets->final = 0;
    sausage->nodesets->entry = NULL;
    sausage->nodesets->exit = NULL;
    sausage->nodesets->prev = NULL;
    sausage->nodesets->next = NULL;
    sausage->nodesets->t_min = sausage->nodesets->t_max = 0;
    
    sausage->n_nodeset++;
    
    if ( 0 != node_set_add(sausage->nodesets, node_stack[n_node-1]) ) {
        perror("Assign initial node to First Node Set error");
        // need to do some cleaning
        return NULL;
    } 
    node_set_t *ns = sausage->nodesets;
    for ( i = n_node - 2; i >= 0; i--) {
        node = node_stack[i];
        if ( 1 == node_set_has_connection(ns, node)) {
            // Construct a new node set
            node_set_t *new_ns;
            new_ns = (node_set_t*) malloc( sizeof(node_set_t) );
            sausage->n_nodeset++;
            new_ns->ns_id = ns->ns_id + 1;
            
            new_ns->nodes = NULL;
            new_ns->last_node = NULL;
            new_ns->n_node = 0;
            
            new_ns->initial = 0;
            new_ns->final = ((i == 0) ? 1 : 0);
            new_ns->entry = NULL;
            new_ns->exit = NULL;
            new_ns->prev = ns;
            new_ns->next = NULL;
            new_ns->t_min = new_ns->t_max = 0;
            // construct a new edge set between ns and new_ns
            edge_set_t* es;
            es = (edge_set_t*) malloc( sizeof(edge_set_t) );
            es->from_ns_id = ns->ns_id;
            es->to_ns_id = new_ns->ns_id;
            es->from = ns;
            es->to = new_ns;
            es->n_edge = 0;
            es->edges = NULL;
            es->last_edge = NULL;
            
            new_ns->entry = es;
            ns->exit = es;
                        
            // add new node set to the previous one
            ns->next = new_ns;
            // update current node set pointer
            ns = new_ns;
        }
        // add latnode to current node set
        if ( node_set_add(ns, node) != 0) {
            perror("Failed to add latnode to node set");
        }
              
        // process all edges from this node's entries
        ps_latnode_t* from;
        ps_latlink_t* link;
        ps_latlink_iter_t* itor;
        for ( itor = ps_latnode_entries(node); itor; itor = ps_latlink_iter_next(itor)) {
            link = ps_latlink_iter_link(itor);
            ps_latlink_nodes(link, &from);
            node_set_t* locator;
            for (locator = ns->prev; locator; locator = locator->prev) {
            
                if (1 == node_set_contains(locator, from) ) {
                   break;
                }
                
            }
            if (!locator) {
                perror("Impossible: ");
                // ...
            }
            // 
            if ( (ns->ns_id - locator->ns_id) == 1) {
                // assign edge to to the edge set which is between node set <ns> and node set <locator> directly
                if ( 0 != edge_set_add(ns->entry, link) ) {
                    perror("Failed to add latlink to edge set");
                    // need to do some cleaning
                }
                
            } else { // Assign edge to one of edge sets between <ns> and <locator>
                double sim = 0, tsim = 0, bstsim = 0;
                edge_t* e; 
                node_set_t *cur_ns, *bst_ns;
                for (cur_ns = locator; cur_ns != ns; cur_ns = cur_ns->next) {
                    // caculate SIM(E, e)
                    
                    for (e = cur_ns->exit->edges; e; e = e->next) {
                        sim += similarity_of_words(ps_latlink_word(dag, e->edge), ps_latlink_word(dag,link));
                    }
                    tsim =  sim * overlap(cur_ns->exit, link) / cur_ns->exit->n_edge;
                    
                    if (tsim >= bstsim) {
                        bstsim = tsim;
                        bst_ns = cur_ns;
                    }
                }                
                if ( 0 != edge_set_add(bst_ns->exit, link) ) {
                    perror("Failed to add latlink to edge set");
                    // need to do some cleaning
                }
            }
            
        }
        
    }
    free(node_stack);
    return sausage;
}

void sausage_last_node_set(sausage_t* s, ps_lattice_t* dag)
{
    node_set_t* ns;
	for (ns = s->nodesets; ns; ns = ns->next) {
	    if (ns->final == 1) {
	        printf("#n_node:%d\n", ns->n_node);
	        node_t* n;
	        for (n = ns->nodes; n; n = n->next) {
	            printf("%s\n", ps_latnode_word(dag, n->node));
	        }
	    }
	}
}


void sausage_write(sausage_t* s, ps_lattice_t* dag, const char* filename)
{
    FILE* fp;
    if ((fp = fopen(filename, "w")) == NULL) {
        perror("Failed to open file to dump sausage");
        return;
    }
    node_set_t* ns;
    edge_set_t* es;
    edge_t* e;
    for (ns = s->nodesets; ns; ns = ns->next) {
        fprintf(fp,"Node Set #%d %d %d %d:\n", ns->ns_id, ns->n_node, ns->t_min, ns->t_max);
        es = ns->exit;
        if (es) {
            for (e = es->edges; e; e = e->next) {
                fprintf(fp, "(%s, %d)\n", ps_latlink_word(dag, e->edge), ps_latlink_prob(dag, e->edge, NULL) );
            }
        }
    }    
    fclose(fp);
}



void sausage_free(sausage_t* s)
{
    if (!s)
        return;
    node_t* n;
    edge_t* e;
    node_set_t* ns;
    edge_set_t* es;
    while (s->nodesets) {
        ns = s->nodesets;
        s->nodesets = s->nodesets->next;
        
        es = ns->exit;
        if (es) {
            while (es->edges) {
                e = es->edges;
                es->edges = es->edges->next;
                free(e);
                es->n_edge--;
            }
            free(es);
        }
        
        while(ns->nodes) {
            n = ns->nodes;
            ns->nodes = ns->nodes->next;
            free(n);
            ns->n_node--;
        }
        free(ns);
        s->n_nodeset--;
    }
    free(s);
}

struct lite_node_s {
  int id;  
  struct lite_edge_set_s *edge_set; 
};

struct lite_edge_set_s {
    int n_edge;
    lite_edge_t* edges;
    lite_edge_t* last_edge;
};

struct lite_edge_s {
    char* word;  /** WORD */
    int32 post; /**posterior likelihood*/
    struct lite_edge_s* next;
};


struct lite_sausage_s {
    int n_node;
    lite_node_t* nodes;
};


void lite_edge_set_add(lite_edge_set_t* lite_es, edge_t* e, ps_lattice_t* dag)
{
    if (!lite_es || !e) {
        perror("Err: Bad lite_es or e in function lite_edge_set_add()");
        return;
    }
    ps_latlink_t* link = e->edge;
    const char* word = ps_latlink_word(dag, link);
    lite_edge_t* lite_e;
    for (lite_e = lite_es->edges; lite_e; lite_e = lite_e->next) {
        if ( strncmp(lite_e->word, word, strlen(word)) == 0) {
            break;
        }
    }
    if (!lite_e) {
        lite_edge_t* lite_edge;
        lite_edge = (lite_edge_t*) malloc(sizeof(lite_edge_t));
        lite_edge->word = (char*) calloc(WORD_MAX_LENGTH + 1, sizeof(char));
        strncpy(lite_edge->word, word, strlen(word));
        lite_edge->post = ps_latlink_prob(dag, link, NULL);
        lite_edge->next = NULL;
        
        if (!lite_es->edges) {
            lite_es->edges = lite_edge;
        } else {
            lite_es->last_edge->next = lite_edge;
        }
        lite_es->last_edge = lite_edge;
        lite_es->n_edge++;
            
    } else {
        lite_e->post = logmath_add(ps_lattice_get_logmath(dag), lite_e->post, ps_latlink_prob(dag, link, NULL) ); 
		if (lite_e->post > 0) 
			lite_e->post = 0;
    }
    
}



lite_sausage_t* sausage_simplify(sausage_t* s, ps_lattice_t* dag)
{
    if (!s || !(s->n_nodeset > 0)) { 
        return NULL;
    }
    
    int n_node = s->n_nodeset;
    lite_sausage_t* lite_s = (lite_sausage_t*) malloc( sizeof(lite_sausage_t) );
    lite_s->n_node = n_node;    
    lite_s->nodes = (lite_node_t*) calloc( n_node, sizeof(lite_node_t) );
    
    int i;
    node_set_t* ns;
    edge_set_t* es;
    edge_t* e;
    for (i = 0, ns = s->nodesets; i < n_node && ns; i++, ns = ns->next) {
        lite_s->nodes[i].id = ns->ns_id;
        lite_s->nodes[i].edge_set = NULL;
        if (ns->exit) {
            es = ns->exit;
            lite_s->nodes[i].edge_set = (lite_edge_set_t*) malloc( sizeof(lite_edge_set_t) );
            lite_s->nodes[i].edge_set->n_edge = 0;
            lite_s->nodes[i].edge_set->edges = lite_s->nodes[i].edge_set->last_edge = NULL;
            for (e = es->edges; e; e = e->next) {
                lite_edge_set_add(lite_s->nodes[i].edge_set, e, dag);
            }
        }        
    }
    return lite_s;      
}

void lite_sausage_write(lite_sausage_t* lite_s, const char* filename)
{
    FILE* fp;
    if ((fp = fopen(filename, "w")) == NULL) {
        perror("Failed to open file to dump sausage");
        return;
    }
    
    if(!lite_s || !(lite_s->n_node >0)) {
        perror("Bad lite_s");
        return;
    }
    int i;
    lite_node_t* node;
    lite_edge_t* edge;
    for (i = 0; i < lite_s->n_node; i++) {
        node = &(lite_s->nodes[i]);
        fprintf(fp, "Node %d\n", node->id);
        if (node->edge_set) {
            for (edge = node->edge_set->edges; edge; edge = edge->next) {
                fprintf(fp, "(%s, %d)\n", edge->word, edge->post);
            }
        }
    }
    
    fclose(fp);
}

void lite_sausage_free(lite_sausage_t* lite_s)
{
    if (!lite_s)
        return;
        
    int i;
    lite_node_t* node;
    lite_edge_t* edge;
    for (i = 0; i < lite_s->n_node; i++) {
        node = &(lite_s->nodes[i]);  
        if (node->edge_set) {
            while (node->edge_set->edges) {
                edge = node->edge_set->edges;
                node->edge_set->edges = node->edge_set->edges->next;
                free(edge->word);
                free(edge);
                node->edge_set->n_edge--;
            }
        }
    }
    free(lite_s);
}


struct s_hit_s {
    char* uttid;    /** utterance id*/
    int32 post;     /** posterior likelihood */
    struct s_hit_s* next;
};


struct s_hits_pos_s {
	int pos;
    int n_hit;
    s_hit_t* first;
    s_hit_t* last;
    struct s_hits_pos_s* next;
};

struct s_hits_word_s {
    int n_pos;
    s_hits_pos_t* first;
    s_hits_pos_t* last;
};

struct dualclue_index_s {
    int n_word;
    char** word_list;
    s_hits_word_t* s_hits;
};

s_hits_pos_t* s_hits_word_get_pos(s_hits_word_t* hits_word, int pos)
{
	s_hits_pos_t* hits_pos;
	for (hits_pos = hits_word->first; hits_pos; hits_pos = hits_pos->next) {
		if (pos == hits_pos->pos) {
			return hits_pos;
		}
	}
	return NULL;
}

dualclue_index_t* dualclue_index_init(const char* filename)
{
    int i;
    FILE* fp;
    char s[WORD_MAX_LENGTH + 2] = {'\0',};
    if ( (fp = fopen(filename, "r")) == NULL ) {
        perror("Failed to open file to initialize dualclue index");
        return NULL;
    }
    dualclue_index_t* index = (dualclue_index_t*) malloc( sizeof(dualclue_index_t) );
    index->n_word = 0;
    index->word_list = NULL;
    index->s_hits = NULL;
    
    while ( fgets(s, WORD_MAX_LENGTH + 2, fp) != '\0') {
        index->n_word++;
    }
    index->word_list = (char**) malloc(index->n_word * sizeof(char*));
    for (i = 0; i < index->n_word; i++) {
        index->word_list[i] = (char*) malloc( (WORD_MAX_LENGTH + 1) * sizeof(char));
        memset(index->word_list[i], 0, (WORD_MAX_LENGTH + 1) * sizeof(char));
    }
    fseek(fp, 0, SEEK_SET);
    i = 0;
    while ((fgets(s, WORD_MAX_LENGTH + 2, fp) != '\0')) {
        strncpy(index->word_list[i], strtok(s, "\n"), WORD_MAX_LENGTH);
	    //printf("%s", index->word_list[i]);	    
	    i++;
    }
    fclose(fp);
    
    index->s_hits = (s_hits_word_t*) calloc(index->n_word, sizeof(s_hits_word_t));
    return index;
}

int dualclue_index_get_wid(dualclue_index_t* index, const char* word)
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

void dualclue_index_addhit(dualclue_index_t* index, const char* uttid, lite_sausage_t* lite_s)
{
    if (!index) {
        perror("dualclue_index_addhit: Bad index");
        return;
    }
    if(!lite_s || !(lite_s->n_node >0)) {
        perror("dualclue_index_addhit: Bad lite_s");
        return;
    }
    int i;
    int wid, pos;
    lite_node_t* node;
    lite_edge_t* edge;
	s_hits_pos_t* hits_pos;
    for (i = 0; i < lite_s->n_node; i++) {
        node = &(lite_s->nodes[i]);
        pos = node->id;
        if (node->edge_set) {
            for (edge = node->edge_set->edges; edge; edge = edge->next) {
                wid = dualclue_index_get_wid(index, edge->word);
                if (wid == -1) { /** skip incorrect word */
                    continue;
                }   
                /**looking for pointer which points to the position of the word , */
                if (index->s_hits[wid].n_pos == 0) { /** no hits in this word right now, 
					also means there's not any position information inside the word,
					so create the first position information */
					hits_pos = (s_hits_pos_t*) calloc(1, sizeof(s_hits_pos_t));
					hits_pos->pos = pos;
					index->s_hits[wid].n_pos++;
					index->s_hits[wid].first = index->s_hits[wid].last = hits_pos;
				}else { /** look for the wanted position pointer among existing position information*/
					hits_pos = s_hits_word_get_pos(&(index->s_hits[wid]), pos);
					if (NULL == hits_pos){ /** current position is not found among existing position information*/
						hits_pos = (s_hits_pos_t*) calloc(1, sizeof(s_hits_pos_t));
						hits_pos->pos = pos;
						index->s_hits[wid].n_pos++;
						index->s_hits[wid].last->next = hits_pos;
						index->s_hits[wid].last = hits_pos;
					} 
				}
				/**Found pointer **hits_pos** which points to the position of the word , then add hit to that position */
				s_hit_t* hit;
				hit = (s_hit_t*) calloc(1, sizeof(s_hit_t));
				hit->uttid = (char*) calloc(strlen(uttid)+1, sizeof(char));
				strncpy(hit->uttid, uttid, strlen(uttid));
				hit->post = edge->post;
				if (hits_pos->n_hit == 0) {
					hits_pos->first = hit;
				} else {
					hits_pos->last->next = hit;
				}
				hits_pos->last = hit;
				hits_pos->n_hit++;
            }
        }
    }
}

void dualclue_index_free(dualclue_index_t* index)
{
	if(!index) {
		return;
	}
	int i;
	s_hits_word_t* hits_word;
	s_hits_pos_t* hits_pos;
	s_hit_t* hit;
	for (i = 0; i < index->n_word; i++) {
		hits_word = &(index->s_hits[i]);
		/* == free hits inside this word ==*/
		while (hits_word->first) {
			hits_pos = hits_word->first;
			hits_word->first = hits_word->first->next;
			/* == free hits inside this position ==*/
			while (hits_pos->first) {
				hit = hits_pos->first;
				hits_pos->first = hits_pos->first->next;
				free(hit->uttid);
				free(hit);
				hits_pos->n_hit--;			
			}			
			free(hits_pos);
			hits_word->n_pos--;			
		}
	}
	for (i = 0; i < index->n_word; i++) {
		free(index->word_list[i]);
	}
	free(index->word_list);
	free(index->s_hits);
	free(index);
}


void dualclue_index_write(dualclue_index_t* index, const char* filename)
{
	FILE* fp;
	if ( (fp = fopen(filename, "w")) == NULL) {
		perror("dualclue_index_write: BAD filename");
		return;
	}
	
	fprintf(fp, "# Words: %d\n", index->n_word);
    s_hits_word_t* hits_word;
	s_hits_pos_t* hits_pos;
	s_hit_t* hit;
	int i;
    for (i = 0; i < index->n_word; i++) {
		hits_word = &(index->s_hits[i]);
		fprintf(fp, "WORD#%d %s (%d)\n", i, index->word_list[i], hits_word->n_pos);
		for (hits_pos = hits_word->first; hits_pos; hits_pos = hits_pos->next) {
			fprintf(fp, "POS #%d\n", hits_pos->pos);
			for (hit = hits_pos->first; hit; hit = hit->next) {
				fprintf(fp, "(%d, %s)\n", hit->post, hit->uttid);
			}
		}
		
    }
	
	fclose(fp);
}

dualclue_index_t* dualclue_index_read(const char* filename)
{
	FILE* fp;
	if ((fp = fopen(filename, "r")) == NULL) {
		perror("dualclue_index_read: BAD filename");
		return NULL;
	}
	
	int n_word;
	if ( 1 != fscanf(fp,"# Words: %d\n", &n_word)) {
		perror("dualclue_index_read: format error");
		return NULL;
	};
	dualclue_index_t* index = (dualclue_index_t*) malloc( sizeof(dualclue_index_t) );
    index->n_word = n_word;
    index->word_list = NULL;
    index->s_hits = NULL;
	/** word list */
    index->word_list = (char**) calloc(index->n_word, sizeof(char*));
	/** hits */
	index->s_hits = (s_hits_word_t*) calloc(index->n_word, sizeof(s_hits_word_t));
	
	int i, k;
	char line[MAX_LINE_LENGTH] = {'\0',}; 
	char word[WORD_MAX_LENGTH + 1] = {'\0',};
	int n_pos, pos;
	char uttid[MAX_LINE_LENGTH] = {'\0',};
	int32 post;
	s_hits_pos_t* hits_pos;
	while ( NULL != fgets(line, MAX_LINE_LENGTH, fp) ) {
		if ( (( k = sscanf(line, "WORD#%d %s (%d)\n", &i, word, &n_pos) ) != 3) &&
				(( k = sscanf(line, "POS #%d\n", &pos) ) != 1) &&
				(( k = sscanf(line, "(%d, %[^)])\n", &post, uttid) ) != 2) ) {
			dualclue_index_free(index);
			return NULL;
		}
		if ( k == 3) { // new word
			//printf("WORD#%d %s (%d)\n", i, word, n_pos);
			index->word_list[i] = (char*) calloc(WORD_MAX_LENGTH+1, sizeof(char));
			strncpy(index->word_list[i], word, strlen(word));
			index->s_hits[i].n_pos = n_pos;			
		} 
		if ( k == 1) { // new position
			//printf("POS #%d\n", pos);
			hits_pos = (s_hits_pos_t*) calloc(1, sizeof(s_hits_pos_t));
			hits_pos->n_hit = 0;
			hits_pos->pos = pos;
			hits_pos->first = hits_pos->last = NULL;
			if (!index->s_hits[i].first) {
				index->s_hits[i].first = index->s_hits[i].last = hits_pos;
			} else {
				index->s_hits[i].last->next = hits_pos;
				index->s_hits[i].last = hits_pos;
			}
		}
		if ( k == 2 ) { // new hit
			//printf("(%d, %s)\n", post, uttid);
			s_hit_t* hit = (s_hit_t*) calloc(1, sizeof(s_hit_t));
			hit->uttid = (char*) calloc(strlen(uttid)+1, sizeof(char));
			strncpy(hit->uttid, uttid, strlen(uttid));
			hit->post = post;
			if (!hits_pos->first) {
				hits_pos->first = hit;
			}else {
				hits_pos->last->next = hit;
			}
			hits_pos->last = hit;
			hits_pos->n_hit++;
		} 
	}
	
	fclose(fp);
	return index;
}
/*
struct dualclue_index_cache_s {
	int n_word;
	int* n_pos;
	char** word_list;
	s_hit_t*** s_hits;
};

s_hit_t** dualclue_index_cache_get_element(dualclue_index_cache_t* cache, int i, int j)
{
	// i [0, n_word) ; j [0, n_pos[i])
	int k;
	int c = 0;
	for ( k = 0; k < i; k++) {
		c = c + cache->n_pos[k];
	} 
	c += j;
	return ((s_hit_t**)(cache->s_hits) + c);	
}

dualclue_index_cache_t* dualclue_index_get_cache(dualclue_index_t* index)
{
	if (!index) {
		return NULL;
	}
	dualclue_index_cache_t* cache = (dualclue_index_cache_t*) calloc(1, sizeof(dualclue_index_cache_t));
	cache->n_word = index->n_word;
	cache->word_list = index->word_list;
	cache->n_pos = (int*) calloc(cache->n_word, sizeof(int));
	cache->s_hits = (s_hit_t***) calloc(cache->n_word, sizeof(s_hit_t**));
	int i, j;
	int c;
	int n_pos;
	for (i = 0; i < index->n_word; i++) {
		cache->n_pos[i] = index->s_hits[i].n_pos;
		cache->s_hits[i] = (s_hit_t**) calloc(cache->n_pos[i], sizeof(s_hit_t*));
	}
	s_hit_t** tmp;
	s_hits_pos_t* hits_pos;
	for (i = 0; i < index->n_word; i++) {
		for (j = 0; j < cache->n_pos[i]; j++) {
			tmp = dualclue_index_cache_get_element(cache, i, j);
			hits_pos = s_hits_word_get_pos( &(index->s_hits[i]), j);
			if (hits_pos) 			
				*tmp = hits_pos->first;
			else
				*tmp = NULL;
		}
	}
	return cache;
}

void dualclue_index_cache_write(dualclue_index_cache_t* cache, const char* filename)
{
	FILE* fp;
	if ( (fp = fopen(filename, "w")) == NULL) {
		perror("dualclue_index_cache_write: BAD filename");
		return;
	}
	
	fprintf(fp, "# Words: %d\n", cache->n_word);
	s_hit_t** hits_pos;
	s_hit_t* hit;
	int i, j;
    for (i = 0; i < cache->n_word; i++) {
		fprintf(fp, "WORD#%d %s (%d)\n", i, cache->word_list[i], cache->n_pos[i]);
		for (j = 0; j < cache->n_pos[i]; j++) {
			hits_pos = dualclue_index_cache_get_element(cache, i, j);
			if (*hits_pos) {
				fprintf(fp, "POS #%d\n", j);
				hit = *hits_pos;
				while(hit) {
					fprintf(fp, "(%d, %s)\n", hit->post, hit->uttid);
					hit = hit->next;
				}
			}
		}
		
    }
	
	fclose(fp);
}

*/

/**
 * s_partial_path_t
 */
typedef struct s_partial_path_s {
    s_hit_t* first_term;
    s_hit_t* last_term;
    int n_term;
	int pos;
    struct s_partial_path_s *next;
} s_partial_path_t;

/**
 * s_path_queue_t
 *
 */
typedef struct s_path_queue_s {
    s_partial_path_t* head;
    s_partial_path_t* tail;
    int n_path;
} s_path_queue_t; 

/** Initialize an empty path */ 
s_partial_path_t* s_partial_path_init() {
    s_partial_path_t* p;
    p = (s_partial_path_t*) calloc(1, sizeof(s_partial_path_t) );
    return p;
}

/** Finitialize a path */ 
void s_partial_path_free(s_partial_path_t* p) 
{
    if (p) {
        s_hit_t* h;
        while (p->first_term) {
            h = p->first_term;
            p->first_term = p->first_term->next;
            /** Don't forget to free the **uttid** member in hit_t */
            free(h->uttid); 
            free(h);
            p->n_term--;
        }
        free(p);
    }
}

/** Extend the path with a new term */
int s_partial_path_extend(s_partial_path_t* p, s_hit_t* h) 
{
    if (!p) 
        return -1;
    if (!h)
        return -1;

    s_hit_t* hit;
    hit = (s_hit_t*) malloc(sizeof(s_hit_t));
    
    hit->uttid = (char*) malloc(MAX_LINE_LENGTH * sizeof(char));
    memset(hit->uttid, 0, MAX_LINE_LENGTH * sizeof(char));
    strncpy(hit->uttid, h->uttid, strlen(h->uttid));
    hit->post = h->post;
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
s_partial_path_t* s_partial_path_copy(s_partial_path_t* p) {
       if (!p) {
           return NULL;
       }
       s_hit_t* h;
       s_partial_path_t* out = s_partial_path_init();
       for (h = p->first_term; h; h = h->next) {
            if (s_partial_path_extend(out, h) != 0) {
                perror("Add hits to path failure");
                s_partial_path_free(out);
                return NULL;
            }
       }
       out->next = NULL;
       return out;
}


/** Get posterior log-likelihood of the partial path: P(path|O) */
int32 s_partial_path_get_posterior(s_partial_path_t* p)
{
    s_hit_t* hit;
    int32 result = 0;	
	for (hit = p->first_term; hit; hit = hit->next) {
		result += hit->post;
	}
    return result;
}

/* =====================================================================
 * path_queue_t's function definitions 
 * ===================================================================== */ 
 
s_path_queue_t* s_path_queue_init() 
{
    s_path_queue_t* q = (s_path_queue_t*) malloc( sizeof(s_path_queue_t) );
    q->head = q->tail = NULL;
    q->n_path = 0;
    return q;
}

void s_path_queue_free(s_path_queue_t* q)
{
    if (!q) return;
    s_partial_path_t* p;
    while (q->head) {
        p = q->head;
        q->head = q->head->next;
        s_partial_path_free(p);
        q->n_path--;
    }
    free(q);    
}

void s_path_queue_add(s_path_queue_t* q, s_partial_path_t* p) {
    if (!q) {
        perror("s_path_queue_add: No queue object");
        return;
    }
    if (!p) {
        perror("s_path_queue_add: No path object");
        return;
    }
    /* Next Path in p Should be NULL, becoz p is newest. */
    if (p->next != NULL) {
        perror("s_path_queue_add: Next Path Should be NULL");
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

void dualclue_index_search(dualclue_index_t* index, char** terms, int n_term)
{
	int i;
	s_path_queue_t** queues = (s_path_queue_t**) calloc(n_term, sizeof(s_path_queue_t*));
	for (i = 0; i < n_term; i++) {
		queues[i] = s_path_queue_init();
	}
	// Process the 1st query term
	int wid = dualclue_index_get_wid(index, terms[0]);
	if (wid == -1)
		goto exit;
	s_hits_word_t* hits_word = &(index->s_hits[wid]);
	s_hits_pos_t* hits_pos;
	s_hit_t* hit; 
	for (hits_pos = hits_word->first; hits_pos; hits_pos = hits_pos->next) {
		hit = hits_pos->first;
		while (hit) {
			s_partial_path_t* p = s_partial_path_init();
			s_partial_path_extend(p, hit);
			p->pos = hits_pos->pos;
			s_path_queue_add(queues[0], p);
			hit = hit->next;
		}
	}
	for (i = 1; i < n_term; i++ ) {
		if (queues[i-1]->n_path == 0) {
			goto exit;
		}
		wid = dualclue_index_get_wid(index, terms[i]);
		hits_word = &(index->s_hits[wid]);
		s_partial_path_t* p;
		for (p = queues[i-1]->head; p; p = p->next) {
			hits_pos = s_hits_word_get_pos(hits_word, p->pos + 1); // adjust position range
			if (!hits_pos) {
				continue;
			}
			hit = hits_pos->first;
			while (hit) {
				if ( strncmp(hit->uttid, p->first_term->uttid, strlen(hit->uttid)) == 0) {
					s_partial_path_t* q = s_partial_path_copy(p);
					s_partial_path_extend(q, hit);
					q->pos = hits_pos->pos;
					s_path_queue_add(queues[i], q);
				}
				hit = hit->next;
			}
		}	
	}
	
	if ( i == n_term) {
		s_partial_path_t* p;
		for (p = queues[n_term-1]->head; p; p = p->next) {
			printf("%s %d\n", p->first_term->uttid, s_partial_path_get_posterior(p));
		}
	}
	
exit:
	for (i = 0; i < n_term; i++) {
		s_path_queue_free(queues[i]);
	}
	free(queues);
}
