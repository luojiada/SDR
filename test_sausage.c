#include <stdio.h>
#include <stdlib.h>

#include "sausage.h"


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

	dag = ps_get_lattice(ps);
	if (dag == NULL) {
	    perror("No lattice");
	    return 1;
	}

	sausage_t* s = convert_lattice_to_sausage(dag);
	sausage_write(s, dag, "sausage.txt");
	
    lite_sausage_t* lite_s = sausage_simplify(s,dag);
    lite_sausage_write(lite_s, "simplified_sausage.txt");
	
    dualclue_index_t* index = dualclue_index_init("./syllable.lst");
    dualclue_index_addhit(index, argv[1], lite_s);
	
	dualclue_index_write(index, "dualclue_index.txt");
    dualclue_index_free(index);
	
	index = dualclue_index_read("dualclue_index.txt");
	dualclue_index_write(index, "dualclue_index2.txt");
	
	
	char* query[] = {"jin", "tian", "jie", "mu", "de", "zhu", "yao", "nei", "rong", "you"};
	dualclue_index_search(index, query, 10);
	
	
	dualclue_index_free(index);
    lite_sausage_free(lite_s);
    sausage_free(s);
	return 0;
}
