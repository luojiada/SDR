#include <stdio.h>
#include <stdlib.h>

#include "index.h"


int main(int argc, char** argv)
{
    int i;
    int rv;
    FILE* fh;
    char const *hyp, *uttid;
    int32 score;
    ps_decoder_t *ps;
	cmd_ln_t *config;
    ps_lattice_t* dag;

    config = cmd_ln_init(NULL, ps_args(), TRUE,
			     "-hmm", "./hmm/zh_broadcastnews_ptm256_8000",
			     "-lm", "./lm/syllables.lm.DMP",
			     "-dict", "./lm/syllables_sorted.dic",
			     NULL);
	if (config == NULL)
		return 1;
	ps = ps_init(config);
	if (ps == NULL)
		return 1;
   
	fh = fopen(argv[1], "rb");
	if (fh == NULL) {
		perror("Failed to open audio file.");
		return 1;
	}

	rv = ps_decode_raw(ps, fh, "test", -1);
	if (rv < 0)
		return 1;

	hyp = ps_get_hyp(ps, &score, &uttid);
	if (hyp == NULL)
		return 1;
	printf("Recognized: %s\n", hyp);
	inverted_index_t* index = inverted_index_init("./syllable.lst");
    if (index == NULL) {
       exit(1);
    }

	dag = ps_get_lattice(ps);
	if (dag == NULL) {
	    perror("No lattice");
	    return 1;
	}

    /*
    printf("# Total number of words: %d\n", index->n_word);
    for(i = 0; i < index->n_word; i++) {
        printf("%3d: %s\n", i+1, index->word_list[i]);
    }*/
    float32 ascale = cmd_ln_float32_r(config, "-ascale");
    printf("ascale: %f\n", ascale);
    //printf("%d: %s\n", inverted_index_get_wid(index, "ba"), "ba");
    //printf("%d: %s\n", inverted_index_get_wid(index, "bia"), "bia");
    inverted_index_addhits(index, "test", dag, 1.0/ascale);
    inverted_index_write(index, "./index");
    inverted_index_free(index);
    index = inverted_index_read("./index");
    inverted_index_write(index, "./index2");
    
    char* query[] = {"jin", "tian", "jie", "mu"};
    result_list_t* rl;
    inverted_index_search(index, ps_get_lmset(ps), 1.0/ascale, query, 4, &rl);

	inverted_index_free(index);
	return 0;
}
