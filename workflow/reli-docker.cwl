class: CommandLineTool
cwlVersion: v1.0
$namespaces:
  sbg: 'https://www.sevenbridges.com/'
id: reli
baseCommand:
  - RELI
inputs:
  - id: snp_file
    type: File
    inputBinding:
      position: 0
      prefix: '-snp'
    label: SNP file
    doc: Phenotype SNP file in 4 column BED format
  - id: ld_file
    type: File?
    inputBinding:
      position: 1
      prefix: '-ld'
    label: LD file
    doc: 'Phenotype linkage disequilibrium structure for SNPs (default: none)'
  - id: index_file
    type: File
    inputBinding:
      position: 2
      prefix: '-index'
    label: Index file
    doc: ChIP-seq index file
  - id: chipseq_data_dir
    type: Directory
    inputBinding:
      position: 3
      prefix: '-data'
    label: Data directory
    doc: Directory where ChIP-seq data are stored
  - id: target
    type: string
    inputBinding:
      position: 4
      prefix: '-target'
    label: Target
    doc: Target label of ChIP-seq experiment to be tested from index file
  - id: build_file
    type: File
    inputBinding:
      position: 5
      prefix: '-build'
    label: Build file
    doc: Genome build file
  - id: 'null'
    type: File
    inputBinding:
      position: 6
      prefix: '-null'
    label: Null file
    doc: Null model file
  - id: dbsnp_table
    type: File
    inputBinding:
      position: 7
      prefix: '-dbsnp'
    label: dbSNP table
    doc: dbSNP table file
  - 'sbg:toolDefaultValue': 'off'
    id: match
    type: boolean?
    inputBinding:
      position: 9
      prefix: '-match'
    label: Minor AF matching
    doc: 'Enable/disable minor allele frequency based matching (default: off)'
  - id: rep
    type: int?
    inputBinding:
      position: 9
      prefix: '-rep'
    label: Number of simulations
    doc: 'Number of permutations/simulations to be performed (default: 2000)'
  - id: corr
    type: int?
    inputBinding:
      position: 11
      prefix: '-corr'
    label: Bonferroni correction multiplier
    doc: 'Bonferroni correction multiplier for multiple tests (default: 1)'
  - id: phenotype
    type: string?
    inputBinding:
      position: 12
      prefix: '-phenotype'
    label: Phenotype
    doc: 'User provided phenotype name (default: ".")'
  - id: ancestry
    type: string?
    inputBinding:
      position: 13
      prefix: '-ancestry'
    label: Ancestry
    doc: 'User provided ancestry name (default: ".")'
outputs:
  - id: overlapped_snps_file
    doc: Output files from completed RELI analysis
    label: Overlaps output file
    type: File
    outputBinding:
      glob: '*.rsids'
  - id: stats_file
    label: Statistics output file
    type: File
    outputBinding:
      glob: '*.stats'
  - id: stdout_log
    doc: The stderr of the RELI analysis run
    label: stdout log file
    type: File?
    outputBinding:
      glob: '*stdout.txt'
  - id: stderr_log
    doc: The stderr of the RELI analysis run
    label: stderr log file
    type: File?
    outputBinding:
      glob: '*stderr.txt'
label: RELI_public
arguments:
  - position: 8
    prefix: '-out'
    valueFrom: $(runtime.outdir)
requirements:
  - class: DockerRequirement
    dockerPull: 'weirauchlab/reli:latest'
  - class: InlineJavascriptRequirement
stdout: reli-stdout.txt
stderr: reli-stderr.txt
