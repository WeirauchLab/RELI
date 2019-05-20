#ifndef RELI_impl_h
#define RELI_impl_h
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
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_randist.h>



using namespace std;
extern string buffer;
extern const int bufferSize;
extern char bufferChar[];
extern vector<string> tfbs2SNP_DirAnnoSet;
extern vector<string> tfbs2SNP_RefAnnoSet;

namespace RELI{
	#define BedSigNullModelDir ""
	#define BedSigBedFileDir ""
	#define BedSigFastaDir ""

	extern vector<string> linehandler(string);
	//classes 
	class phastConData{
	public:
		int bin;
		string chr;
		int st;
		int end;
		double sum;
		double sumsqare;
		double avg_score;
		inline bool operator<(const phastConData& rhs)const{
			return this->chr < rhs.chr || (this->chr == rhs.chr && this->st < rhs.st);
		}
		inline bool operator==(const phastConData& rhs)const{
			return  (this->chr == rhs.chr && this->st == rhs.st && this->end == rhs.end);
		}
	};
	class bed3col{  
	public:
		string bed_chr;
		unsigned int bed_start;
		unsigned int bed_end;
		unsigned int length;
		vector<phastConData> overlapped_phastCon_data_vec;
		double avg_peak_length_adjusted_phastCon_score;
		void cal_avg_peak_length_adjusted_phastCon_score();
		void cal_avg_peak_length_adjusted_phastCon_score_ez();
		void cal_avg_peak_length_adjusted_phastCon_score_ez_50bp();
		inline bool operator<(const bed3col& rhs) const{
			return (this->bed_chr < rhs.bed_chr || (this->bed_chr == rhs.bed_chr && this->bed_start < rhs.bed_start));
		}
		vector<double> ez_phastCon_score_vec;
		vector<double> ez_phastCon_score_vec_padding_50bp; 
	};
	class target_bed_file{
	public:
		vector<unsigned int> lengthvec;
		vector<bed3col> myData;
		vector<bed3col> myData_bgnull;
		map<string, int> index;
		map<pair<string, unsigned int>, int> index2;	
		void makeIndex();
		void makeIndex2();
		void readingData(string, bool);
		unsigned int median_data_length;
		vector<double> phastCon_score_vec_from_peaks;
		double avg_phastCons_score;
		double sd_phastCons_score;
		void cal_avg_phastCons_score();
	};
	class SNP{  
	public:
		unsigned int inherited_unique_key_from_LD;
		string snp_chr;
		unsigned int snp_start;
		unsigned int snp_end;
		unsigned int length;		
		string snp_name;			
		string snp_type;    
		string obs_strand;
		string _ref_allele;
		vector<string> _alt_allele;
		double _MAF;
		int _MAF_Bin;
		int _TSS_Bin;
		int cal_MAF_Bin(const double &inVal){
			if (inVal < 0.05){ return 0; }
			if (0.05 <= inVal &&inVal < 0.1){ return 1; }
			if (0.1 <= inVal &&inVal < 0.15){ return 2; }
			if (0.15 <= inVal &&inVal < 0.2){ return 3; }
			if (0.2 <= inVal &&inVal < 0.25){ return 4; }
			if (0.25 <= inVal &&inVal < 0.3){ return 5; }
			if (0.3 <= inVal &&inVal < 0.35){ return 6; }
			if (0.35 <= inVal &&inVal < 0.4){ return 7; }
			if (0.4 <= inVal &&inVal < 0.45){ return 8; }
			if (0.45 <= inVal &&inVal <= 0.5){ return 9; }
            return 0;
		}
		string _upstreamseq_50bp;
		string _downstreamseq_50bp;
		string ref_seq_padded_50bp;
		vector<string> alt_seq_padded_50bp;
		string _MAF_alt_allele_string;
		inline bool operator<(const SNP& rhs) const{
			return (this->snp_chr < rhs.snp_chr ||
				(this->snp_chr == rhs.snp_chr && this->snp_end < rhs.snp_end));
		}
		inline bool operator==(const string& rhs) const{
			return (this->snp_name == rhs);
		}
	};
	class LD_template{  
	public:
		vector<string> mySNP;
		string keySNP;
	};
	class LD{   
	public:
		vector<SNP> mySNP;
		SNP keySNP;
		vector<int> dis2keySNP;
		int max_dis;
		int min_dis;
		// for goshifter
		string LD_chr;
		int LD_left_edge;
		int LD_right_edge;
		LD() :max_dis(0), min_dis(0){};
		target_bed_file features_within_LDblock;
		void get_features_within_LDblock(const target_bed_file& rhs);
		vector<bed3col> goShifting_feature_data();
	};
	class LD_sim :public LD{   
	public:
		bool overlap_sim;
		unsigned int unique_key;
	};
	class LD_shomo_method :public LD{
	public:
		unsigned int extended_LD_start;
		unsigned int extended_LD_end;
		vector<bed3col> intersected_peaks;

		bool overlap_sim;
	};
	class myless{
	public:
		bool operator()(SNP snpA, SNP snpB) const {
			return ((snpA.snp_name<snpB.snp_name) || ((snpA.snp_name == snpB.snp_name) && (snpA.snp_chr<snpB.snp_chr))
				|| ((snpA.snp_name == snpB.snp_name) && (snpA.snp_chr == snpB.snp_chr) && (snpA.snp_start < snpB.snp_start)));
		}
	};
	class mymap :public map < string, string > {
	public:
		static iterator find2(iterator _it1, iterator _it2, string _val){
			for (auto it = _it1; it != _it2; ++it){
				if (it->first == _val){
					return it;
				}
			}
			return _it2;
		}
	};
	class MAF_binned_null_model{
	public:
		vector<unsigned int> bin0;	//0 - 5%
		vector<unsigned int> bin1;	//5 - 10%
		vector<unsigned int> bin2;	//10 - 15%
		vector<unsigned int> bin3;	//15 - 20%
		vector<unsigned int> bin4;	//20 - 25%
		vector<unsigned int> bin5;	//25 - 30%
		vector<unsigned int> bin6;	//30 - 35%
		vector<unsigned int> bin7;	//35 - 40%
		vector<unsigned int> bin8;	//40 - 45%
		vector<unsigned int> bin9;	//45 - 50%
		map<int, vector<unsigned int>*> bin_map;
		MAF_binned_null_model(){
			bin_map[0] = &bin0;
			bin_map[1] = &bin1;
			bin_map[2] = &bin2;
			bin_map[3] = &bin3;
			bin_map[4] = &bin4;
			bin_map[5] = &bin5;
			bin_map[6] = &bin6;
			bin_map[7] = &bin7;
			bin_map[8] = &bin8;
			bin_map[9] = &bin9;
		}
		void loading_null_data(string);
	};
	class MAF_TSS_binned_null_model{
	public:
		map<pair<int, int>, vector<unsigned int>> bin_map;
	};
	class resultClass{	
	public:
		string resultLine;
		string track;
		string cell;
		string tf;
		int overlap;
		int total;
		double ratio;
		double mean;
		double sd;
		double zscore;
		double enrichment;
		double pval;
		double correctedPval;
		string nullModel;
		string species;
		bool randomizedSet;
		bool significant;
		string mappedTf;
		resultClass() :randomizedSet(false), significant(false){};
		inline bool operator<(const resultClass & rhs) const{
			return this->correctedPval < rhs.correctedPval;
		}
	};

	class data_index{public:
		string datalabel;
		string source;
		string cell;
		string tf;
		string cell_label;
		string pmid;
		string group;
		string ebv_status;
		string species;

		inline bool operator==(const string& rhs)const{
			return this->datalabel == rhs;
		}
	};
	class snp_table_data{public:
		string chr;
		string start;
		string end;
		string rsid;
		string obs_strand;
		string ref_allele;
		string alt_alleles;
		string type;
		string alt_allele_info;
		string alt_allele_freq;
		string bin;
	};
	class RELIobj{public:
		static const int lcsize = 500000000;
		char linechar[lcsize];
		string line; 
		string dnaSeqReverse(string inSeq, map<char, char> thismap){
			string tempseq;
			for (string::reverse_iterator sit = inSeq.rbegin(); sit != inSeq.rend(); ++sit){
				tempseq += thismap[*sit];
			}
			return tempseq;
		}
		map<char, char> ATGCmap;

		/*
			variables and data containers
		*/ 
		string public_ver_snp_fname;		// phenotype snp after LD expansion
		string public_ver_ld_fname;			// ld structure
		string public_ver_genome_build_fname;	// build
		string public_ver_null_fname;	// null model	
		string public_ver_null_index_fname;	// null model index
		string public_ver_data_dir;		// chipseq data
		string public_ver_target_label;			// target label
		string public_ver_target_data_fname;	// target file name: dir + target label
		string public_ver_data_index_fname;	// chipseq index
		vector<data_index> dataindexvec;		// index obj vec
		data_index public_ver_selected_data_index;	// selected index obj 
		string public_ver_snp_table_fname;	// snp table
		string public_ver_output_dir;
        string public_ver_output_prefix;    // output files prefix
		string public_ver_output_fname;		// initialized as target label + "RELI.stats"
		string public_ver_output_fname_overlaps;		// initialized as target label + "RELI.overlaps"
		unordered_map<string, RELI::snp_table_data> snptablemap;
		string public_ver_phenotype_name;
		string public_ver_ancestry_name;
		bool public_ver_debug;
		bool flag_input_snp;
		bool flag_ld_file;
		bool flag_dbsnp_table;
		bool flag_null_file;
		bool flag_chipseq_data_index;
		bool flag_target_label;
		bool flag_genome_build;
		bool flag_chipseq_data_dir;
		bool flag_output_dir;
        bool flag_output_prefix;
		/*
			functions
		*/
		void public_ver_read_null(); 
		void public_ver_output();
		void public_ver_set_target_data();
		void public_ver_read_data_index();
		void public_ver_read_snp_table(); 
		void create_output_dir();
		void load_snp_table();
		void extract_snp_info(map<char,char>);
		void load_ld_snps(bool, string);
		void output();
		void sim();

		bool minimum_check();

		// constructor
		RELIobj(){
			this->ATGCmap['A'] = 'T';
			this->ATGCmap['T'] = 'A';
			this->ATGCmap['C'] = 'G';
			this->ATGCmap['G'] = 'C';
			this->ATGCmap['N'] = 'N';
			this->public_ver_phenotype_name = ".";
			this->public_ver_ancestry_name = ".";
			this->public_ver_debug = false;	
			this->flag_input_snp=false;
			this->flag_ld_file = false;
			this->flag_dbsnp_table = false;
			this->flag_null_file = false;
			this->flag_target_label = false;
			this->flag_chipseq_data_index = false;
			this->flag_genome_build = false;
			this->flag_chipseq_data_dir = false;
			this->flag_output_dir = false;
            this->flag_output_prefix = false;
		}
	};
	//variables  
	enum stats_model{ normal, empirical, phasetype, binomial, hypergeometric, fishers_exact }; 
	// external variables 
	extern vector<SNP> SNP_vec; 	
	extern vector<SNP> SNP_vec_temp;   
	extern vector<LD> LD_vec;	
	extern vector<LD_template> LD_template_vec; 	
	extern vector<LD_sim> LD_sim_vec, ldsimvec_after_intersection; 	
	extern vector<bed3col> targetbedinfilevec;
	extern vector<double> statsvec;
	extern vector<pair<string, unsigned int>> chromosome_strucuture;
	extern vector<unsigned int> chromosome_strucuture_val;
	extern vector<int> simulated_number_vec;
	extern vector<string> species_name;
	extern string TFtype;
	extern string targetbedinfilename;
	extern string nullmodel;
	extern string nullmodelinfiledir;
	extern string nullmodelinfilename;
	extern string bgnullmodel;
	extern string bgnullmodelinfilename;
	extern string defaultBedFileDir;
	extern string dnase_coverage_filename;
	extern string snpinfilename;
	extern string temp_nullmodelname, temp_nullmodelcount;
	extern string cmdline;
	extern string outputpath;
	extern string outputpath_simulated_intersection;
	extern string ldfile;
	extern string species;
	extern string default_used_species;
	extern string species_chr_mapping_file;
	extern unsigned int dnase_coverage;
	extern int repmax;
	extern bool ldfile_flag;
	extern bool using_default_species;
	extern bool snp_matching;
	extern bool snp_matching_local_shift;
	extern bool bgnull;
	extern double mu, sd, zscore, pval, corr_pval;
	extern double corr_muliplier;
	extern double sig_pct;
	extern map<string, unsigned int> dnase_coverage_map;
	extern map<string, int> targetbedfileindex_start;
	extern map<SNP, LD_sim, myless> snp2ldsim;
	extern mymap speciesMap;
	extern stats_model used_stats_model;
	extern MAF_binned_null_model binned_null_model_data;
	extern MAF_binned_null_model bg_null_model_data;
	extern MAF_TSS_binned_null_model binned_null_model_data2;
	extern unordered_map<string, vector<unsigned int>> indexing_mapping;
	bool SNPfit(LD, SNP &, unsigned int, vector<pair<string, unsigned int>>, vector<unsigned int>);
	bool SNPfit(LD, SNP &, unsigned int, vector<pair<string, unsigned int>>, vector<unsigned int>, bool);
	bool SNPfit_local(LD &, SNP &, unsigned int, vector<pair<string, unsigned int>>, vector<unsigned int>, bool);
	bool SNPfit_goshift(LD &, SNP &, unsigned int, vector<pair<string, unsigned int>>, vector<unsigned int>, bool);
	bool SNPfit_local_index_only(LD &, SNP &, unsigned int, vector<pair<string, unsigned int>>, vector<unsigned int>, bool);
	bool mybedsort(bed3col, bed3col);
	bool mysort(LD_sim, LD_sim);
	bool myunique(LD_sim, LD_sim);
	bool snpsort(SNP, SNP);
	unsigned int mymax(vector<unsigned int>);
	unsigned int mymin(vector<unsigned int>);
	const int lookback_step = 50; // 
	int lookback_with_zerocheck(int);
	void callSpeciesMap();
	void snpmodifier(SNP &, SNP, int);
	void overlapping2(vector<SNP>, vector<bed3col>);
	void overlapping3(vector<SNP>, vector<bed3col>, vector<unsigned int>&);
	void overlapping_w_index(vector<SNP>, vector<bed3col>, vector<unsigned int>&, map<string, int>);
	void createSpeciesMap(bool);
	void cal_stats(stats_model);
	void read_ld_file(string);
	void ReadBedSigResult(string, vector<RELI::resultClass>&, bool &);
	void ReadTFMapping(string, map<string, string>&);
	void initiate_BedSig();
	int get_index_to_be_used(string, int, map<pair<string, unsigned int>, int>);
	extern int secondary_index_size;
	int binomialCoeff(int, int);
	double binomial_pvalue(int, int, double);
	double binomial_pvalue_appr(int, int, double);
	void loadSnpFile(string);
}
// end of RELI namespace  




#endif
