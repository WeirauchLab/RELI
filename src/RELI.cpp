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

#include <iostream>
#include <cstring>      // for strcmp, strncmp
#include <cstdlib>      // for EXIT_SUCCESS, EXIT_FAILURE
#include "RELI_impl.h"
#include "config.h"
#include "ansicolor.hpp"

using namespace std;
using namespace RELI;


void display_banner() {
	cout << "+-------------------------------------------------------------------------+" << endl
	     << "|          Regulatory Element Locus Intersection (RELI) Analysis          |" << endl
	     << "|                           Current version: 0.90                         |" << endl
	     << "+-------------------------------------------------------------------------+" << endl
	     << endl;
}

void display_help(int ret=EXIT_SUCCESS, int level=HELP_BRIEF) {
	display_banner();

	cout << "  " << color::ul << "usage" << color::reset << ":" << endl
	     << "    RELI [-h | -help | --help | -help-all] [-V | -version]" << endl
	     << "    RELI -snp FILE -ld FILE -index FILE -data DIR -target STRING"<< endl
	     << "         -build FILE -null FILE -dnsnp FILE -out DIR [-rep NUMBER]" << endl
	     << "         [-corr NUMBER] [-phenotype STRING] [-ancestry STRING]" << endl
	     << endl;

	if (level < HELP_OPTIONS) {
		cout << "  " << color::bold_green << "Hint" << color::reset
		     << ": Try 'RELI -h' for explanations of the options." << endl
		     << "  Please report issues at " RELI_ISSUE_TRACKER << endl
		     << endl;
		exit(ret);
	}

	cout << "  " << color::ul << "where" << color::reset << ":" << endl
	     << "    -h, -help, --help  displays this help" << endl
	     << "    -help-all          displays full help, including examples" << endl
	     << "    -V, -version       displays the RELI version number and quits" << endl
	     << "    -snp FILE          (required) phenotype SNP file in 4-column BED format" << endl
	     << "    -ld FILE           phenotype linkage disequilibrium structure for SNPs" << endl
	     << "                       [default: no LD file]" << endl
	     << "    -index FILE        (required) ChIP-seq index file." << endl
	     << "    -data DIR          (required) directory where ChIP-seq data are stored." << endl
	     << "    -target STRING     (required) target label of ChIP-seq experiment to be" << endl
	     << "                       tested from index file" << endl
	     << "    -build FILE        (required) specifies the genome build file" << endl
	     << "    -null FILE         (required) specifies the null model file" << endl
	     << "    -dbsnp FILE        (required) specifies the dbSNP table file" << endl
	     << "    -out DIR           (required) output directory name to be created inside" << endl
	     << "                       inside the current working directory" << endl
	     << "    -match             Boolean switch to enable minor allele frequency based" << endl
	     << "                       matching [default: off]" << endl
	     << "    -rep NUMBER        specifies number of permutations/simulations to be" << endl
	     << "                       performed [default: 2000]" << endl
	     << "    -corr NUMBER       Bonferroni correction multiplier for multiple tests" << endl
	     << "                       [default: 1]" << endl
	     << "    -phenotype STRING  user-specified phenotype name [default: \".\"]" << endl
	     << "    -ancestry STRING   user-specified ancestry name [default: \".\"]" << endl
	     << endl;

	if (level < HELP_EXAMPLES) {
		cout << "  " << color::bold_green << "Hint" << color::reset
		     << ": Try 'RELI -help-all' for example usage." << endl
		     << "  Please report issues at " RELI_ISSUE_TRACKER << endl
		     << endl;
		exit(ret);
	}

	cout << "  " << color::ul << "example" << color::reset << ":" << endl
	     << "    RELI \\" << endl
	     << "      -snp SLE_EU.snp \\" << endl
	     << "      -ld SLE_EU.ld \\" << endl
	     << "      -index data/ChIPseq.index \\" << endl
	     << "      -data data/ChIP-seq \\" << endl
	     << "      -target hg19_0302 \\"<< endl
	     << "      -build data/GenomeBuild/hg19.txt \\" << endl
	     << "      -null data/Null/CommonSNP_MAFmatch \\" << endl
	     << "      -dbsnp data/SNPtable/SNPtable \\" << endl
	     << "      -out Output \\" << endl
	     << "      -match \\" << endl
	     << "      -rep 2000 \\" << endl
	     << "      -corr 1544 \\" << endl
	     << "      -phenotype Systemic_Lupus_Erythematosus \\ " << endl
	     << "      -ancestry EU" << endl
	     << endl
	     << "  Please report issues at " RELI_ISSUE_TRACKER << endl
	     << endl;

	exit(ret);
}


int main(int argc, char* argv[]){
	/*
		initilize RELI instance and handle options
	*/
    RELI::RELIobj *RELIinstance = new RELI::RELIobj;

	for (auto i = 1; i < argc; ++i){
		if (strcmp(argv[i], "-help-all") == 0 ||
			strcmp(argv[i], "--help-all") == 0) {
			display_help(EXIT_SUCCESS, HELP_EXAMPLES);
		}
		if (strncmp(argv[i], "-h", 2) == 0 ||
			strcmp(argv[i], "--help") == 0) {
			display_help(EXIT_SUCCESS, HELP_OPTIONS);
		}
        if (strncmp(argv[i], "-V", 2) == 0 ||
            strncmp(argv[i], "-vers", 5) == 0 ||
		    strncmp(argv[i], "--vers", 6) == 0) {
            cout << "RELI v" RELI_VERSION << endl;
            exit(EXIT_SUCCESS);
        }
		if (strcmp(argv[i], "-null") == 0){		// null model
			RELIinstance->public_ver_null_fname = argv[i + 1];
			RELIinstance->flag_null_file = true;
		}
		if (strcmp(argv[i], "-dbsnp") == 0){		// dbSNP table that contains MAF allele frequency info
			RELIinstance->public_ver_snp_table_fname = argv[i + 1];
			RELIinstance->flag_dbsnp_table = true;
		}
		if (strcmp(argv[i], "-match") == 0){		// snp matching indicator
			RELI::snp_matching = true;
		}
		if (strcmp(argv[i], "-snp") == 0){		// phenotype snp file
			RELIinstance->public_ver_snp_fname = argv[i + 1];
			RELIinstance->flag_input_snp = true;
		}
		if (strcmp(argv[i], "-ld") == 0){		// phenotype ld file
			RELI::ldfile = argv[i + 1];
			RELI::ldfile_flag = true;
			RELIinstance->flag_ld_file = true;
		}
		if (strcmp(argv[i], "-data") == 0){		// directory of target ChIP-seq files
			RELIinstance->public_ver_data_dir= argv[i + 1];
			RELIinstance->flag_chipseq_data_dir = true;
		}
		if (strcmp(argv[i], "-out") == 0){		// output directory
			RELIinstance->public_ver_output_dir = argv[i + 1];
			RELIinstance->flag_output_dir = true;
		}
		if (strcmp(argv[i], "-rep") == 0){		// replication number
			RELI::repmax = atoi(argv[i + 1]);
		}
		if (strcmp(argv[i], "-corr") == 0){		// bonferroni correction multiplier
			RELI::corr_muliplier = atof(argv[i + 1]);
		}
		if (strcmp(argv[i], "-target") == 0){		//	string corresponding to target ChIP-seq data
			RELIinstance->public_ver_target_label = argv[i + 1];
			RELIinstance->flag_target_label = true;
		}
		if (strcmp(argv[i], "-index") == 0){		// ChIP-seq index file
			RELIinstance->public_ver_data_index_fname = argv[i + 1];
			RELIinstance->flag_chipseq_data_index = true;
		}
		if (strcmp(argv[i], "-build") == 0){		// genome build file
			RELI::species_chr_mapping_file = argv[i + 1];
			RELI::using_default_species = false;
			RELIinstance->flag_genome_build = true;
		}
		if (strcmp(argv[i], "-phenotype") == 0){		// phenotype name
			RELIinstance->public_ver_phenotype_name = argv[i + 1];
		}
		if (strcmp(argv[i], "-ancestry") == 0){		// phenotype ancestry
			RELIinstance->public_ver_ancestry_name= argv[i + 1];
		}
	}

	if (argc < 2 || !RELIinstance->minimum_check()) {	// check minimum arguments for run
		display_help(EXIT_FAILURE);
	}

	display_banner();
	
	/*
		load data and pre-processing
	*/
	RELI::createSpeciesMap(RELI::using_default_species); //	load genome structure
	RELIinstance->public_ver_read_data_index();	//	read ChIP-seq index file
	RELIinstance->public_ver_set_target_data();	//	set target ChIP-seq file
	RELI::target_bed_file TBF;	// initialize target ChIP-seq file object
	TBF.readingData(RELIinstance->public_ver_target_data_fname, false);	// load target ChIP-seeq file
	TBF.makeIndex();	// create target ChIP-seq file index
	RELI::targetbedinfilevec = TBF.myData;	//
	RELI::targetbedfileindex_start = TBF.index;	//
	RELI::binned_null_model_data.loading_null_data(RELIinstance->public_ver_null_fname.c_str());  //	load null model
	RELI::loadSnpFile(RELIinstance->public_ver_snp_fname.c_str()); //	load phenotype snp file
	RELIinstance->load_snp_table();	//	load dbsnp table
	RELIinstance->extract_snp_info(RELIinstance->ATGCmap);	//	extract snp bin info from dbsnp table
	RELI::SNP_vec_temp = RELI::SNP_vec;   // copy over snp data for loading ld structure
	RELIinstance->load_ld_snps(RELI::ldfile_flag,RELI::ldfile);	//	load phenotype ld structure file

	/*
		permutation/simulation
	*/
	RELIinstance->sim();

	/*
		handle statistics and output
	*/
	RELIinstance->output();

	exit(EXIT_SUCCESS);
}
