/*
Regulatory Element Locus Intersection (RELI) analysis
Copyright (C) <2017>  <Xiaoting.Chen@cchmc.org>

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <thread>
#include <mutex>
#include <condition_variable>
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <numeric>
#include <limits.h>
#include <stdlib.h>
#include <stdio.h>
#include <string>
#include <string.h>
#include <time.h>
#include <map>
#include <math.h>
#include <algorithm>
#include <ctime>
#include <unordered_map>
#include <chrono>
#include <cstring>
#include <queue>
#include <random>
#include <unistd.h>
#include "RELI_impl.h"

using namespace std;
using namespace RELI;


void display_RELI(){
	cout << "+---------------------------------------------------------------+" << endl;
	cout << "|                                                               |" << endl;
	cout << "|     Regulatory Element Locus Intersection (RELI) Analysis     |" << endl;
	cout << "|                    Current version: 0.91                      |" << endl;
	cout << "|                                                               |" << endl;
	cout << "+---------------------------------------------------------------+" << endl;

}
void display_help(){

	cout << "Usage:	./RELI [options]" << endl;
	cout << "OPTIONS are:" << endl;
	cout << endl;
	cout << "-snp FILE: Phenotype snp file in 4 column bed format. [required]" << endl;
	cout << "-ld FILE: Phenotype linkage disequilibrium structure for snps, default: no ld file. [optional]" << endl;
	cout << "-index FILE: ChIP-seq index file. [required if -target points to the label]" << endl;
	cout << "-data DIR: Specify directory where ChIP-seq data are stored. [required if -target points to the label]" << endl;
	cout << "-target STRING: Target label of ChIP-seq experiment to be tested from index file. Or path to the file [required]" << endl;
	cout << "-build FILE: Genome build file. [required]" << endl;
	cout << "-null FILE: Null model file. [required]" << endl;
	cout << "-dbsnp FILE: dbSNP table file. [required]" << endl;
	cout << "-out DIR: Specify output directory name under current working folder. [required]" << endl;
    cout << "-prefix STRING: Specify output files prefix. [optional]" << endl;
	cout << "-match: Boolean switch to turn on minor allele frequency based matching, default: off. [optional]" << endl;
	cout << "-rep NUMBER: Number of permutation/simulation to be performed, default: 2000. [optional]" << endl;
	cout << "-corr NUMBER: Bonferroni correction multiplier for multiple test, default: 1 [optional]" << endl;
	cout << "-phenotype STRING: User provided phenotype name, default: \".\". [optional]" << endl;
	cout << "-ancestry STRING: User provided ancestry name, default: \".\". [optional]" << endl;
	cout << endl;
	cout << "EXAMPLE:" << endl;
	cout << "../script/RELI \\" << endl;
	cout << "-snp SLE_EU.snp  \\" << endl;
	cout << "-ld SLE_EU.ld \\" << endl;
	cout << "-index ../data/ChIPseq.index \\" << endl;
	cout << "-data ../data/ChIP-seq \\" << endl;
	cout << "-target hg19_0302 \\" << endl;
	cout << "-build ../data/GenomeBuild/hg19.txt \\" << endl;
	cout << "-null ../data/Null/CommonSNP_MAFmatch \\" << endl;
	cout << "-dbsnp ../data/SNPtable/SNPtable \\" << endl;
	cout << "-out Output   \\" << endl;
	cout << "-match \\" << endl;
	cout << "-rep 2000 \\" << endl;
	cout << "-corr 1544 \\" << endl;
	cout << "-phenotype Systemic_Lupus_Erythematosus \\" << endl;
	cout << "-ancestry EU " << endl;

	exit(-1);
}


int main(int argc, char* argv[]){
	/*
		initilize RELI instance and handle options
	*/
    RELI::RELIobj *RELIinstance = new RELI::RELIobj;
	for (auto i = 1; i < argc; ++i){
		if (strcmp(argv[i], "-null") == 0){  // null model
			RELIinstance->public_ver_null_fname = argv[i + 1];
			RELIinstance->flag_null_file = true;
		}
		if (strcmp(argv[i], "-dbsnp") == 0){  // dbSNP table that contains MAF allele frequency info
			RELIinstance->public_ver_snp_table_fname = argv[i + 1];
			RELIinstance->flag_dbsnp_table = true;
		}
		if (strcmp(argv[i], "-match") == 0){  // snp matching indicator
			RELI::snp_matching = true;
		}
		if (strcmp(argv[i], "-snp") == 0){  // phenotype snp file
			RELIinstance->public_ver_snp_fname = argv[i + 1];
			RELIinstance->flag_input_snp = true;
		}
		if (strcmp(argv[i], "-ld") == 0){  // phenotype ld file
			RELI::ldfile = argv[i + 1];
			RELI::ldfile_flag = true;
			RELIinstance->flag_ld_file = true;
		}
		if (strcmp(argv[i], "-data") == 0){  // directory of target ChIP-seq files
			RELIinstance->public_ver_data_dir= argv[i + 1];
			RELIinstance->flag_chipseq_data_dir = true;
		}
		if (strcmp(argv[i], "-out") == 0){  // output directory
			RELIinstance->public_ver_output_dir = argv[i + 1];
			RELIinstance->flag_output_dir = true;
		}
        if (strcmp(argv[i], "-prefix") == 0){  // output file prefix
            RELIinstance->public_ver_output_prefix = argv[i + 1];
            RELIinstance->flag_output_prefix = true;
        }
		if (strcmp(argv[i], "-rep") == 0){  // replication number
			RELI::repmax = atoi(argv[i + 1]);
		}
		if (strcmp(argv[i], "-corr") == 0){  // bonferroni correction multiplier
			RELI::corr_muliplier = atof(argv[i + 1]);
		}
		if (strcmp(argv[i], "-target") == 0){  // string corresponding to target ChIP-seq data or path to file
			RELIinstance->public_ver_target_label = argv[i + 1];
            if (access(RELIinstance->public_ver_target_label.c_str(), F_OK ) != -1 ) {
				// if -target points to the file, -data and -index are not required anymore and should be set to True
				RELIinstance->flag_chipseq_data_dir = true;
				RELIinstance->flag_chipseq_data_index = true;
                RELIinstance->public_ver_data_dir = "Ignored";
                RELIinstance->public_ver_data_index_fname = "Ignored";
			}
			RELIinstance->flag_target_label = true;
		}
		if (strcmp(argv[i], "-index") == 0){  // ChIP-seq index file
			RELIinstance->public_ver_data_index_fname = argv[i + 1];
			RELIinstance->flag_chipseq_data_index = true;
		}
		if (strcmp(argv[i], "-build") == 0){  // genome build file
			RELI::species_chr_mapping_file = argv[i + 1];
			RELI::using_default_species = false;
			RELIinstance->flag_genome_build = true;
		}
		if (strcmp(argv[i], "-phenotype") == 0){  // phenotype name
			RELIinstance->public_ver_phenotype_name = argv[i + 1];
		}
		if (strcmp(argv[i], "-ancestry") == 0){  // phenotype ancestry
			RELIinstance->public_ver_ancestry_name= argv[i + 1];
		}
	}
	display_RELI();
	if (argc < 2 || !RELIinstance->minimum_check()) {  // check minimum arguments for run
		display_help();
	}

	/*
		load data and pre-processing
	*/
    RELIinstance->create_output_dir();
	RELI::createSpeciesMap(RELI::using_default_species);  // load genome structure
	RELIinstance->public_ver_read_data_index();  // read ChIP-seq index file
	RELIinstance->public_ver_set_target_data();  // set target ChIP-seq file
	RELI::target_bed_file TBF;  // initialize target ChIP-seeq file object
	TBF.readingData(RELIinstance->public_ver_target_data_fname, false);  // load target ChIP-seeq file
	TBF.makeIndex();  // create target ChIP-seeq file index
	RELI::targetbedinfilevec = TBF.myData;  //
	RELI::targetbedfileindex_start = TBF.index;  //
	RELI::binned_null_model_data.loading_null_data(RELIinstance->public_ver_null_fname.c_str());  // load null model
	RELI::loadSnpFile(RELIinstance->public_ver_snp_fname.c_str());  // load phenotype snp file
	RELIinstance->load_snp_table();  // load dbsnp table
	RELIinstance->extract_snp_info(RELIinstance->ATGCmap); // extract snp bin info from dbsnp table
	RELI::SNP_vec_temp = RELI::SNP_vec;  // copy over snp data for loading ld structure
	RELIinstance->load_ld_snps(RELI::ldfile_flag,RELI::ldfile);	 // load phenotype ld structure file

	/*
		permutation/simulation
	*/
	RELIinstance->sim();

	/*
		handle statistics and output
	*/
	RELIinstance->output();


	return 0;
}

