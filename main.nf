#!/usr/bin/env nextflow

/*
 * gavin.kelly@crick.ac.uk
 * harshil.patel@crick.ac
 * nourdine.bah@crick.ac.uk
 * philip.east@crick.ac.uk
 */

/* ------------------------------------------------------------------------- */
/* ------------------------------------------------------------------------- */
/* --				 ______         ______    _                  				  -- */
/* --				(____  \   /\  (____  \  | |                 				  -- */
/* --				 ____)  ) /  \  ____)  )  \ \                				  -- */
/* --				|  __  ( / /\ \|  __  (    \ \               				  -- */
/* --				| |__)  ) |__| | |__)  )____) )              				  -- */
/* --				|______/|______|______(______/               				  -- */
/* --				                                             				  -- */
/* --				 ______  ______            _                 				  -- */
/* --				(_____ \|  ___ \   /\     | |                				  -- */
/* --				 _____) ) |   | | /  \ ___ \ \   ____ ____   				  -- */
/* --				(_____ (| |   | |/ /\ (___) \ \ / _  ) _  |  				  -- */
/* --				      | | |   | | |__| |_____) | (/ / | | |  				  -- */
/* --				      |_|_|   |_|______(______/ \____)_|| |  				  -- */
/* --				                                        |_|  				  -- */
/* --				 ______                    ___ _             				  -- */
/* --				|  ___ \             _    / __) |            				  -- */
/* --				| |   | | ____ _   _| |_ | |__| | ___  _ _ _ 				  -- */
/* --				| |   | |/ _  | \ / )  _)|  __) |/ _ \| | | |				  -- */
/* --				| |   | ( (/ / ) X (| |__| |  | | |_| | | | |				  -- */
/* --				|_|   |_|\____|_/ \_)\___)_|  |_|\___/ \____|				  -- */
/* --				                                             				  -- */
/* --				 ______ _             _ _                    				  -- */
/* --				(_____ (_)           | (_)                   				  -- */
/* --				 _____) ) ____   ____| |_ ____   ____        				  -- */
/* --				|  ____/ |  _ \ / _  ) | |  _ \ / _  )       				  -- */
/* --				| |    | | | | ( (/ /| | | | | ( (/ /        				  -- */
/* --				|_|    |_| ||_/ \____)_|_|_| |_|\____)       				  -- */
/* --				         |_|                                 				  -- */
/* ------------------------------------------------------------------------- */
/* ------------------------------------------------------------------------- */


// the java modules
import java.nio.file.Paths
import java.nio.file.Files

// exotic source modules
import org.yaml.snakeyaml.Yaml


///////////////////////////////////////////////////////////////////////////////
/* --                                                                     -- */
/* --                             USAGE                                   -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////

def usage() {
	log.info"""
	Usage:

	The workflow can be parametrised with a json or a yaml file or with the
	command line. Here are the parameters with the command line:


		--design "<path>"
		--output_directory "<path>"
		--binomial_nomenclature "<Genus> <Species>"
		--strandedness "<none|forward|reverse>"
		--read_length "<integer>"
		--single_end (run single end mode if present)
		--genome_type (for example for human "<ensembl>")
		--genome_version (for example for human "<GRCh38>")
		--genome_release (for example for human "86")
		--sequencing_run_directory "<path>"
		--sequencing_project_directory "<path>"
		--modules "<path>"
		--genome_sequence_extension "toplevel.fa"
		--fastq_screen_conf "<path>"
		--publish_directory_mode "<symlink|link|copy|move>"
		--publish_directory_overwrite (overwrite if present)
		--r_script "<path>"
		--multiqc_config "<path>"
	"""
}

params.help = false
if (params.help){
	usage()
	exit 0
}


///////////////////////////////////////////////////////////////////////////////
/* --                                                                     -- */
/* --                             METHODS                                 -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////

/*
 * Takes a String, creates the corresponding path after checking if it exists
 * and returns the directory path as a File object.
 */
static File createDirPathFromString(String dirpath) {
	File dir = new File(dirpath)
	if (!dir.exists()) {
		boolean dir_creation = dir.mkdirs()
		String msg =
			"Problem: the " + dirpath +
					" directory couldn't have been created..."
		assert dir_creation : msg
	}
	return dir
}


/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/* are a list elements contained in another list ? */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
def contains_elements(list1, list2) {

	// are the elements of list1 present in list2 ?
	is_present = []
	list1.each{ list2.contains(it) }

	// if at least one is absent
	return ! is_present.contains(false)
}


/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/* gets the columns from the header of a csv file */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
def get_columns(csv_path, sep) {

	// the header
	def header = ""
	new File(csv_path).withReader{ header = it.readLine() }

	// the columns
	def columns = header.split(sep)

	return columns
}


/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/* is the design file conform ? */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
def check_design_file(file_path, single_end_columns, paired_end_columns,
								single_end) {

	// get the columns from the header of the csv design file
	columns = get_columns(file_path, ",")
	
	// the design file must be in accordance with the SINGLE_END parameter
	def is_single_end_conform =
		contains_elements(single_end_columns, columns) && single_end
	def is_paired_end_conform =
		contains_elements(paired_end_columns, columns) && (!single_end)

	// exit if the situation is neither a SINGLE_END or a !SINGLE_END situation
	if ( (!is_single_end_conform) && (!is_paired_end_conform) ) {

		// error message
		msg = 'Error: the design file must contain either the columns ' +
				'["sample,"file"] if it is SINGLE_END, or the columns ' +
				'["sample,"file1","file2"] if it is not SINGLE_END'

		println msg
		exit -1

	} else {

		return true
	}
}


/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/* the conditions are all the columns that are not information columns */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
def get_conditions(file_path, info_columns) {

	// all the columns
	def columns = get_columns(file_path, ",")

	// the non info columns
	def conditions = columns.findAll{ ! info_columns.contains(it) }

	return conditions
}


/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/* maps the colum name to its position in the design file */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
def get_column_map(file_path, info_columns) {

	// associate name to position
	def column_map = [:]

	// the columns of the csv file
	def columns = get_columns(file_path, ",")

	// the conditions of the experiment
	def conditions = get_conditions(file_path, info_columns)

	// the sample name
	column_map["sample"] = columns.findIndexOf{ it == "sample" }

	// the conditions
	conditions.each{ cd -> column_map[cd] = columns.findIndexOf{ it == cd } }

	// the files location
	files = columns.findAll{ ! (conditions+["sample"]).contains(it) }
	files.each{ f -> column_map[f] = columns.findIndexOf{ it == f } }

	return column_map
}


/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/* information in design file is converting into a list of map */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
def get_samples_info(file_path, sep, info_columns) {

	// the list of sample information
	def samples = []

	// the columns and their position
	def column_map = get_column_map(file_path, info_columns)

	new File(file_path).splitEachLine(sep) { fields ->
		def sample = [:]
		def keys = column_map.keySet() as List
		def values = keys.collect{ fields[column_map[it]] }

		// check if it's not the header itself
		if ( ! keys.equals(values) ) {
			keys.eachWithIndex{ k, i -> sample[k] = values[i] }
			samples.add(sample)
		}
	}

	return samples
}


///////////////////////////////////////////////////////////////////////////////
/* --                                                                     -- */
/* --                         INPUT PARAMETERS                            -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////

if (workflow.commandLine.indexOf("-params-file")>0) {

	// design
	DESIGN_FILEPATH = new File(params.design).getAbsolutePath()

	// the output directory
	OUTPUT_DIRPATH = params.output_directory

	// the reads
	SINGLE_END = params.single_end
	STRANDEDNESS = params.strandedness
	READ_LENGTH = params.read_length

	// the binomial nomenclature
	GENUS = params.binomial.genus.toLowerCase()
	SPECIES = params.binomial.species.toLowerCase()

	// annotation
	GENOME_TYPE = params.genome_type
	GENOME_VERSION = params.genome_version
	GENOME_RELEASE = params.genome_release
	GENOME_SEQ_EXTENSION = params.genome_sequence_extension

	// fastq screen configuration file
	FSCREEN_CONF_FILEPATH = new File(params.fastq_screen_conf).getAbsolutePath()

	// modules file
	MODULES_FILEPATH = new File(params.modules).getAbsolutePath()

	// session parameters
	PUBLISHDIR_MODE = params.publish_directory_mode
	PUBLISHDIR_OVERWRITE = params.publish_directory_overwrite

	// r script path
	R_SCRIPT_FILEPATH = new File(params.r_script).getAbsolutePath()

	// multiqc config file
	if (params.multiqc_config) {
		MULTIQC_CONFIG_FILEPATH =
			new File(params.multiqc_config).getAbsolutePath()
	} else {
		MULTIQC_CONFIG_FILEPATH = ""
	}

} else {

	// design
	DESIGN_FILEPATH = new File(params.design).getAbsolutePath()

	// the output directory
	OUTPUT_DIRPATH = params.output_directory

	// the binomial nomenclature
	BINOMIAL = params.binomial_nomenclature.split(' ')
	assert BINOMIAL.size()==2: "Should be --binomial_nomenclature " +
										"<Genus> <Species>"
	GENUS = BINOMIAL[0].toLowerCase()
	SPECIES = BINOMIAL[1].toLowerCase()
	
	// is single end ?
	SINGLE_END = params.single_end ? true : false

	// strandedness
	STRANDEDNESS = params.strandedness.toLowerCase()
	def strandedness_values = ["none", "forward", "reverse"]
	assert STRANDEDNESS in strandedness_values : "Should be --strandedness " +
																"<none|forward|reverse>"

	// read length
	READ_LENGTH = params.read_length

	// annotation
	GENOME_TYPE = params.genome_type
	GENOME_VERSION = params.genome_version
	GENOME_RELEASE = params.genome_release
	GENOME_SEQ_EXTENSION = params.genome_sequence_extension

	// fastq screen configuration file
	FSCREEN_CONF_FILEPATH = new File(params.fastq_screen_conf).getAbsolutePath()

	// modules file
	MODULES_FILEPATH = new File(params.modules).getAbsolutePath()

	// session parameters
	PUBLISHDIR_MODE = params.publish_directory_mode.toLowerCase()
	PUBLISHDIR_OVERWRITE = params.publish_directory_overwrite ? true : false

	// r script path
	R_SCRIPT_FILEPATH = new File(params.r_script).getAbsolutePath()

	// multiqc config file
	if (params.multiqc_config) {
		MULTIQC_CONFIG_FILEPATH =
			new File(params.multiqc_config).getAbsolutePath()
	} else {
		MULTIQC_CONFIG_FILEPATH = ""
	}
}


///////////////////////////////////////////////////////////////////////////////
/* --                                                                     -- */
/* --                              MODULES                                -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////

// the default modules
def MODULE_BAMTOOLS_DEFAULT = "BamTools/2.4.0-foss-2016b"
def MODULE_CUTADAPT_DEFAULT = "cutadapt/1.9.1-foss-2016b-Python-2.7.12"
def MODULE_FSCREEN_DEFAULT = "fastq_screen/0.9.3-2016a-Perl-5.22.1"
def MODULE_MULTIQC_DEFAULT = "multiqc/1.2-2016b-Python-2.7.12"
def MODULE_PICARD_DEFAULT = "picard/2.1.1-Java-1.8.0_92"
def MODULE_R_DEFAULT = "R/3.4.0-intel-2017a-X11-20170314"
def MODULE_RNASEQC_DEFAULT = "RNA-SeQC/1.1.8-Java-1.7.0_80"
def MODULE_RSEM_DEFAULT = "RSEM/1.3.0-foss-2016b"
def MODULE_RSEQC_DEFAULT = "RSeQC/2.6.4-foss-2016b-Python-2.7.12-R-3.3.1"
def MODULE_SAMTOOLS_DEFAULT = "SAMtools/1.3.1-foss-2016b"
def MODULE_STAR_DEFAULT = "STAR/2.5.2a-foss-2016b"

// the modules file
Yaml parser = new Yaml()
def yml = parser.load((MODULES_FILEPATH as File).text)

// get the modules if they are specified
MODULE_BAMTOOLS = yml["bamtools"] ? yml["bamtools"] : MODULE_BAMTOOLS_DEFAULT
MODULE_CUTADAPT = yml["cutadapt"] ? yml["cutadapt"] : MODULE_CUTADAPT_DEFAULT
MODULE_FSCREEN = yml["fscreen"] ? yml["fscreen"] : MODULE_FSCREEN_DEFAULT
MODULE_MULTIQC = yml["multiqc"] ? yml["multiqc"] : MODULE_MULTIQC_DEFAULT
MODULE_PICARD = yml["picard"] ? yml["picard"] : MODULE_PICARD_DEFAULT
MODULE_R = yml["r"] ? yml["r"] : MODULE_R_DEFAULT
MODULE_RNASEQC = yml["rnaseqc"] ? yml["rnaseqc"] : MODULE_RNASEQC_DEFAULT
MODULE_RSEM = yml["rsem"] ? yml["rsem"] : MODULE_RSEM_DEFAULT
MODULE_RSEQC = yml["rseqc"] ? yml["rseqc"] : MODULE_RSEQC_DEFAULT
MODULE_SAMTOOLS = yml["samtools"] ? yml["samtools"] : MODULE_SAMTOOLS_DEFAULT
MODULE_STAR = yml["star"] ? yml["star"] : MODULE_STAR_DEFAULT


///////////////////////////////////////////////////////////////////////////////
/* --                                                                     -- */
/* --                            READ LENGTH                              -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////

// three babs star indices have been built with these read length parameters
def starIndexReadLengths = [50, 75, 100]

// take the index with the closest read length to the experiment's
def diffs = []
starIndexReadLengths.each() { length ->
	diff = (length - READ_LENGTH).abs()
	diffs.add(diff)
}
def index = diffs.findIndexValues() { i -> i == diffs.min() }[0]
def ROUGH_READ_LENGTH = starIndexReadLengths[index.toInteger()]


///////////////////////////////////////////////////////////////////////////////
/* --                                                                     -- */
/* --                     EXPERIMENTAL VARIABLES                          -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////

// the species
def BINOMIAL = GENUS.capitalize() + " " + SPECIES
def BINOMIAL_DIRNAME = BINOMIAL.replace(" ", "_").toLowerCase()
def BINOMIAL_FILENAME = BINOMIAL.replace(" ", "_")


///////////////////////////////////////////////////////////////////////////////
/* --                                                                     -- */
/* --                       OUTPUT DIRECTORIES                            -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////
 
// the annotation
def GENOME_DIRNAME = "release-" + GENOME_RELEASE
def ROUGH_READ_LENGTH_DIRNAME = ROUGH_READ_LENGTH + "bp"

// the babs directories
def BABS_WORKING_DIRPATH = "/camp/stp/babs/working"

// the subdirectories names
def ALIGNMENT_DIRNAME = "alignment"
def ANALYSIS_DIRNAME = "analysis"
def OUTPUT_DIRNAME = "output"
def LOG_DIRNAME = "log"

// the subsubdirectories names (alignment)
def CUTADAPT_DIRNAME = "cutadapt"
def FSCREEN_DIRNAME = "fastq_screen"
def RSEM_DIRNAME = "rsem"
def STAR_DIRNAME = "star"
def STAR_MBSCREEN_DIRNAME = "mbscreen"
def PICARD_DIRNAME = "picard"
def RSEQC_DIRNAME = "rseqc"
def RNASEQC_DIRNAME = "rnaseqc"

// -- ######################## -- //
// -- the subdirectories paths -- //
// -- ######################## -- //

def CUTADAPT_DIRPATH = Paths.get(OUTPUT_DIRPATH, CUTADAPT_DIRNAME).toString()
def FSCREEN_DIRPATH = Paths.get(OUTPUT_DIRPATH, FSCREEN_DIRNAME).toString()
def STAR_DIRPATH = Paths.get(OUTPUT_DIRPATH, STAR_DIRNAME).toString()
def STAR_MBSCREEN_DIRPATH = Paths.get(OUTPUT_DIRPATH,
										STAR_MBSCREEN_DIRNAME).toString()
def PICARD_DIRPATH = Paths.get(OUTPUT_DIRPATH, PICARD_DIRNAME).toString()
def RSEQC_DIRPATH = Paths.get(OUTPUT_DIRPATH, RSEQC_DIRNAME).toString()
def RNASEQC_DIRPATH = Paths.get(OUTPUT_DIRPATH, RNASEQC_DIRNAME).toString()
def ANALYSIS_DIRPATH = Paths.get(OUTPUT_DIRPATH, ANALYSIS_DIRNAME).toString()
def LOG_DIRPATH = Paths.get(OUTPUT_DIRPATH, LOG_DIRNAME).toString()


///////////////////////////////////////////////////////////////////////////////
/* --                                                                     -- */
/* --                     BABS COMMON FILES PATHS                         -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////

// annotation
def ANNOT_EXTENSION = "gtf"
def ANNOT_BASENAME =
	BINOMIAL_FILENAME + "." + GENOME_VERSION + "." + GENOME_RELEASE

// sequence
def SEQ_BASENAME = BINOMIAL_FILENAME + "." + GENOME_VERSION

// directory names
def BABS_DATA_DIRNAME = "data"
def BABS_GENOME_DIRNAME = "genomes"
def INDICE_DIRNAME = "genome_idx"
def ANNOT_DIRNAME = ANNOT_EXTENSION
def RSEM_STAR_MBIOL_DIRNAME = "molecular_biology"
def RSEM_STAR_MBIOL_INDICE_DIRNAME = "all.idx" 
def SEQ_DIRNAME = "genome"

// indice names
def RSEM_STAR_INDICE_NAME = "genome"
def RSEM_STAR_MBSCREEN_INDICE_NAME = "molecular_biology.all"

/////////////////////
// directory paths //
/////////////////////

// the babs shared data directory path
def BABS_DATA_DIRPATH =
	Paths.get(BABS_WORKING_DIRPATH, BABS_DATA_DIRNAME).toString()

// the babs genomes directory path
def BABS_GENOME_DIRPATH =
	Paths.get(BABS_DATA_DIRPATH, BABS_GENOME_DIRNAME).toString()

// the directory path containing the annotation files for this experiment
def GENOME_DIRPATH =
	Paths.get(
		BABS_GENOME_DIRPATH, BINOMIAL_DIRNAME,
		GENOME_TYPE, GENOME_VERSION,
		GENOME_DIRNAME
		).toString()

def ANNOT_DIRPATH = Paths.get(GENOME_DIRPATH, ANNOT_DIRNAME).toString()

def INDICE_DIRPATH = Paths.get(GENOME_DIRPATH, INDICE_DIRNAME).toString()

// the indice file for the alignment with star
def RSEM_STAR_INDICE_DIRPATH =
	Paths.get(
		INDICE_DIRPATH, RSEM_DIRNAME,
		STAR_DIRNAME,
		ROUGH_READ_LENGTH_DIRNAME
		).toString()

// the indice file for the molecular biology alignment with star
def RSEM_STAR_MBSCREEN_INDICE_DIRPATH =
	Paths.get(
		BABS_GENOME_DIRPATH,
		RSEM_STAR_MBIOL_DIRNAME,
		RSEM_STAR_MBIOL_INDICE_DIRNAME,
		RSEM_DIRNAME, STAR_DIRNAME,
		ROUGH_READ_LENGTH_DIRNAME
		).toString()

////////////////
// file names //
////////////////

def ANNOT_FILENAME = ANNOT_BASENAME + "." + ANNOT_EXTENSION
def ANNOT_RNASEQC_FILENAME = ANNOT_BASENAME + ".rnaseqc." + ANNOT_EXTENSION
def ANNOT_BED_FILENAME = ANNOT_BASENAME + ".bed"
def ANNOT_REFFLAT_FILENAME = ANNOT_BASENAME + ".refflat"
def ANNOT_RRNA_FILENAME = ANNOT_BASENAME + ".rRNA.list"
def ANNOT_RRNA_INTERVAL_FILENAME = ANNOT_BASENAME + ".rRNA.interval_list"

// i don't know if this nomenclature can be unified...
// mouse, human --> "primary_assembly.fa"
// yeast, fly --> "toplevel.fa"
def SEQ_FILENAME = SEQ_BASENAME + ".dna_sm." + GENOME_SEQ_EXTENSION

////////////////
// file paths //
////////////////

def ANNOT_FILEPATH = Paths.get(ANNOT_DIRPATH, ANNOT_FILENAME).toString()

def ANNOT_RNASEQC_FILEPATH =
	Paths.get(ANNOT_DIRPATH, ANNOT_RNASEQC_FILENAME).toString()

def ANNOT_REFFLAT_FILEPATH =
	Paths.get(ANNOT_DIRPATH, ANNOT_REFFLAT_FILENAME).toString()

def ANNOT_BED_FILEPATH =
	Paths.get(ANNOT_DIRPATH, ANNOT_BED_FILENAME).toString()

def ANNOT_RRNA_FILEPATH =
	Paths.get(ANNOT_DIRPATH, ANNOT_RRNA_FILENAME).toString()

def ANNOT_RRNA_INTERVAL_FILEPATH =
	Paths.get(ANNOT_DIRPATH, ANNOT_RRNA_INTERVAL_FILENAME).toString()

def SEQ_FILEPATH =
	Paths.get(GENOME_DIRPATH, SEQ_DIRNAME, SEQ_FILENAME).toString()


///////////////////////////////////////////////////////////////////////////////
/* --                                                                     -- */
/* --                   OUTPUT DIRECTORIES CREATION                       -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////

// the folder that contains the cutadapt output
File cutadapt_dir = createDirPathFromString(CUTADAPT_DIRPATH)

// the folder that contains the fastq screen output
File fscreen_dir = createDirPathFromString(FSCREEN_DIRPATH)

// the folder that contains the star output
File star_dir = createDirPathFromString(STAR_DIRPATH)

// the folder that contains the star output
File star_mbscreen_dir = createDirPathFromString(STAR_MBSCREEN_DIRPATH)

// the directory that contains picard's output
File picard_dir = createDirPathFromString(PICARD_DIRPATH)

// the directory that contains picard's output
File rseqc_dir = createDirPathFromString(RSEQC_DIRPATH)

// the directory that contains picard's output
File rnaseqc_dir = createDirPathFromString(RNASEQC_DIRPATH)

// the directory that contains logs
File analysis_dir = createDirPathFromString(ANALYSIS_DIRPATH)

// the directory that contains logs
File log_dir = createDirPathFromString(LOG_DIRPATH)


///////////////////////////////////////////////////////////////////////////////
/* --                                                                     -- */
/* --                        THE INPUT CHANNEL                            -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////

// the information columns of the design file
def single_end_columns = ["file","sample"]
def paired_end_columns = ["file1","file2","sample"]
def info_columns = SINGLE_END ?  single_end_columns : paired_end_columns

// is the design file OK ?
def conform_design_file =
	check_design_file(
		DESIGN_FILEPATH,
		single_end_columns,
		paired_end_columns,
		SINGLE_END
		)

// the information about each sample
def samples_info = get_samples_info(DESIGN_FILEPATH, ",", info_columns)

/* ============================================================= */
/* === is the path of the fastq files absolute or relative ? === */

// get the first file path as a representative
def file_columns = SINGLE_END ? ["file"] : ["file1", "file2"]
def first_filepath = samples_info[0][ file_columns[0] ]

// build the absolute path if it's necessary
if ( ! first_filepath.startsWith("/") ) {
	samples_info.each { sample ->
		file_columns.each { file_col ->
			sample[file_col] = new File(sample[file_col]).getAbsolutePath()
		}
	}
}

/* ============================================================= */

// create a channel from the samples
samples = Channel.from(samples_info)

// files for dgea
r_script = Channel.fromPath(R_SCRIPT_FILEPATH)
design_file = Channel.fromPath(DESIGN_FILEPATH)


///////////////////////////////////////////////////////////////////////////////
/* --                                                                     -- */
/* --                             CUTADAPT                                -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////

process cutadapt {

	// custom label
	tag { name }

	// problem with slurm otherwise
	beforeScript "module purge"

	// modules
	module MODULE_CUTADAPT

	// HPC
	cpus 1
	executor "slurm"
	memory "6000"

	// output directory
	publishDir cutadapt_dir.toString(),
		mode: PUBLISHDIR_MODE,
		overwrite: PUBLISHDIR_OVERWRITE

	input:
		val sample from samples
	
	output:
		set val(name), file('*.cutadapt.fastq.gz') into trimmed_fastq
		set val(name), file('*.log') into cutadapt_log

	////////////////////////////////////////////////////////////////////////////
	shell:

		if (SINGLE_END) {

			// name and location
			name = sample["sample"]
			fastq = sample["file"]

			"""
			cutadapt \
				-a AGATCGGAAGAGC \
				-o ${name}.cutadapt.fastq.gz \
				-e 0.1 \
				-q 10 \
				-m 25 \
				-O 1 \
				${fastq} > ${name}.log
			"""

		} else {

			// name and location
			name = sample["sample"]
			name_1 = sample["sample"] + "_1"
			name_2 = sample["sample"] + "_2"
			fastq_1 = sample["file1"]
			fastq_2 = sample["file2"]

			"""
			cutadapt \
				-a AGATCGGAAGAGC -A AGATCGGAAGAGC \
				-o ${name_1}.cutadapt.fastq.gz \
				-p ${name_2}.cutadapt.fastq.gz \
				-e 0.1 \
				-q 10 \
				-m 25 \
				-O 1 \
				${fastq_1} ${fastq_2} > ${name}.log
			"""
		}
}

// the fastq files
trimmed_fastq_star = Channel.create()
trimmed_fastq_star_mbscreen = Channel.create()
trimmed_fastq_screen = Channel.create()
trimmed_fastq.into(
	trimmed_fastq_star,
	trimmed_fastq_star_mbscreen,
	trimmed_fastq_screen
	)


///////////////////////////////////////////////////////////////////////////////
/* --                                                                     -- */
/* --			                  FASTQ SCREEN                                -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////

process fastq_screen {

	// custom label
	tag { sample }

	// problem with slurm otherwise
	beforeScript "module purge"

	// modules
	module MODULE_FSCREEN

	// HPC
	cpus 32
	executor "slurm"
	memory "6000"

	// output directory
	publishDir fscreen_dir.toString(),
		mode: PUBLISHDIR_MODE,
		overwrite: PUBLISHDIR_OVERWRITE

	input:
		set val(sample), file(fastq) from trimmed_fastq_screen

	output:
		set val(sample), file("*.html") into fastq_screen_html
		set val(sample), file("*.txt") into fastq_screen_txt
	
	shell:
		"""
		fastq_screen \
			--force \
			--outdir ./ \
			--subset 200000 \
			--conf ${FSCREEN_CONF_FILEPATH} \
			--threads ${task.cpus} \
			--aligner bowtie2 \
			${fastq}
		"""
}


///////////////////////////////////////////////////////////////////////////////
/* --                                                                     -- */
/* --                               STAR                                  -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////

process rsem_star {

	// custom label
	tag { sample }

	// problem with slurm otherwise
	beforeScript "module purge"

	// modules
	module MODULE_RSEM
	module MODULE_STAR

	// HPC
	cpus 32
	executor "slurm"
	memory "6000"

	// output directory
	publishDir star_dir.toString(),
		mode: PUBLISHDIR_MODE,
		overwrite: PUBLISHDIR_OVERWRITE

	input:
		set val(sample), file(fastq) from trimmed_fastq_star

	output:
		set val(sample), file("*.STAR.genome.bam") into bam_star
		set val(sample), file("*.transcript.bam") into transcript_bam_star
		set val(sample), file("*.results") into star_results
		set val(sample), file("*.stat") into star_stat

	////////////////////////////////////////////////////////////////////////////
				//--temporary-folder ${star_tmp_dir.toString()} \
				//--calc-ci \
				//--ci-memory 10240 \
	shell:

		if (SINGLE_END) {

			"""
			rsem-calculate-expression \
				--temporary-folder "tmp" \
				--star \
				--num-threads ${task.cpus} \
				--strandedness ${STRANDEDNESS} \
				--estimate-rspd \
				--seed 1 \
				--star-output-genome-bam \
				--star-gzipped-read-file \
				${fastq} \
				${RSEM_STAR_INDICE_DIRPATH}/${RSEM_STAR_INDICE_NAME} \
				${sample}
			"""

		} else {

			"""
			rsem-calculate-expression \
				--temporary-folder "tmp" \
				--star \
				--num-threads ${task.cpus} \
				--strandedness ${STRANDEDNESS} \
				--estimate-rspd \
				--seed 1 \
				--star-output-genome-bam \
				--star-gzipped-read-file \
				--paired-end ${fastq[0]} ${fastq[1]} \
				${RSEM_STAR_INDICE_DIRPATH}/${RSEM_STAR_INDICE_NAME} \
				${sample}
			"""
		}
}

// fork for multiqc
bam_star_sort = Channel.create()
bam_star_multiqc = Channel.create()
bam_star.into(bam_star_sort, bam_star_multiqc)

// fork for multiqc
star_results_dgea = Channel.create()
star_results_multiqc = Channel.create()
star_results.into(star_results_dgea, star_results_multiqc)


///////////////////////////////////////////////////////////////////////////////
/* --                                                                     -- */
/* --                DIFFERENTIAL GENE EXPRESSION ANALYSIS                -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////

process dgea {

	// custom label
	tag { script }

	// just to be sure as it's still a test version
	errorStrategy 'ignore'

	// problem with slurm otherwise
	beforeScript "module purge"

	// modules
	module MODULE_R

	// HPC
	cpus 1
	executor "slurm"
	memory "6000"

	// output directory
	publishDir analysis_dir.toString(),
		mode: PUBLISHDIR_MODE,
		overwrite: PUBLISHDIR_OVERWRITE

	input:
		file script from r_script
		file design from design_file
		file results from star_results_dgea.map{x->x[1]}.collect()

	output:
		file "*" into dgea
	
	shell:
		"""
		Rscript ${script} \
			-r . \
			-d "${design}" \
			-b "${BINOMIAL}" \
			-a "${ANNOT_FILEPATH}" \
			-t 0.05 \
			-o .
		"""
}


///////////////////////////////////////////////////////////////////////////////
/* --                                                                     -- */
/* --                           STAR MBSCREEN                             -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////

def star_mbscreen_indice_path = RSEM_STAR_MBSCREEN_INDICE_DIRPATH +
							"/" + RSEM_STAR_MBSCREEN_INDICE_NAME

process rsem_star_mbscreen {

	/* WARNING!
	 * This process tries to align the reads on a molecular biology contaminants
	 * database, so it can fail sometimes.
	 */
	errorStrategy 'ignore'

	// custom label
	tag { sample }

	// problem with slurm otherwise
	beforeScript "module purge"

	// modules
	module MODULE_RSEM
	module MODULE_STAR

	// HPC
	cpus 32
	executor "slurm"
	memory "6000"

	// output directory
	publishDir star_mbscreen_dir.toString(),
		mode: PUBLISHDIR_MODE,
		overwrite: PUBLISHDIR_OVERWRITE

	input:
		set val(sample), file(fastq) from trimmed_fastq_star_mbscreen

	output:
		set val(sample), file("*.STAR.genome.bam") into bam_star_mbscreen
		set val(sample),
			file("*.transcript.bam") into transcript_bam_star_mbscreen
		set val(sample), file("*.results") into results_star_mbscreen
		set val(sample), file("*.stat") into stat_star_mbscreen

	////////////////////////////////////////////////////////////////////////////
	shell:

		if (SINGLE_END) {

			"""
			rsem-calculate-expression \
				--temporary-folder "tmp" \
				--star \
				--num-threads ${task.cpus} \
				--strandedness ${STRANDEDNESS} \
				--seed 1 \
				--star-gzipped-read-file \
				--no-bam-output \
				${fastq} \
				${star_mbscreen_indice_path} \
				${sample}.mbscreen
			"""

		} else {

			"""
			rsem-calculate-expression \
				--temporary-folder "tmp" \
				--star \
				--num-threads ${task.cpus} \
				--strandedness ${STRANDEDNESS} \
				--seed 1 \
				--star-gzipped-read-file \
				--no-bam-output \
				--paired-end ${fastq[0]} ${fastq[1]} \
				${star_mbscreen_indice_path} \
				${sample}.mbscreen
			"""
		}
}


///////////////////////////////////////////////////////////////////////////////
/* --                                                                     -- */
/* --                     SORTING AND INDEXING STAR                       -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////

process sort_index_star {

	// custom label
	tag { sample }

	// problem with slurm otherwise
	beforeScript "module purge"

	// modules
	module MODULE_SAMTOOLS

	// HPC
	cpus 32
	executor "slurm"
	memory "6000"

	// output directory
	publishDir star_dir.toString(),
		mode: PUBLISHDIR_MODE,
		overwrite: PUBLISHDIR_OVERWRITE

	input:
		set val(sample), file(bam) from bam_star_sort

	output:
		set val(sample), file("*.bam"), file("*.bai") into bam_star_sorted

	shell:
		filename = sample + ".sorted.bam"
		"""
		samtools sort \
			--threads ${task.cpus} \
			-o ${filename} \
			${bam}
		samtools index ${filename}
		"""
}

// fork for multiqc
bam_star_sorted_picard = Channel.create()
bam_star_sorted_multiqc = Channel.create()
bam_star_sorted.into(bam_star_sorted_picard, bam_star_sorted_multiqc)


///////////////////////////////////////////////////////////////////////////////
/* --                                                                     -- */
/* --                            PICARD                                   -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////

process picard_group {

	// custom label
	tag { sample }

	// problem with slurm otherwise
	beforeScript "module purge"

	// modules
	module MODULE_PICARD

	// HPC
	cpus 32
	executor "slurm"
	memory "6000"

	// output directory
	publishDir picard_dir.toString(),
		mode: PUBLISHDIR_MODE,
		overwrite: PUBLISHDIR_OVERWRITE

	input:
		set val(sample), file(bam), file(bai) from bam_star_sorted_picard

	output:
		set val(sample), file("*.bam") into picard_rg

	shell:

		tmp_dirname = "tmp"
		filename = sample + ".rg.bam"

		"""
		java -Xmx10g -Djava.io.tmpdir=${tmp_dirname} \
			-jar \$EBROOTPICARD/picard.jar AddOrReplaceReadGroups \
			VALIDATION_STRINGENCY=SILENT \
			INPUT=${bam} \
			OUTPUT=${filename} \
			RGID=${sample} \
			RGLB=${sample} \
			RGPU=${sample} \
			RGSM=${sample} \
			RGCN=TheFrancisCrickInsitute \
			RGPL=Illumina
		"""
}

// forking for multiqc
picard_rg_duplicate = Channel.create()
picard_rg_multiqc = Channel.create()
picard_rg.into(picard_rg_duplicate, picard_rg_multiqc)

process picard_duplicates {

	// custom label
	tag { sample }

	// problem with slurm otherwise
	beforeScript "module purge"

	// modules
	module MODULE_PICARD

	// HPC
	cpus 32
	executor "slurm"
	memory "6000"

	// output directory
	publishDir picard_dir.toString(),
		mode: PUBLISHDIR_MODE,
		overwrite: PUBLISHDIR_OVERWRITE

	input:
		set val(sample), file(bam) from picard_rg_duplicate

	output:
		set val(sample), file("*.bam") into picard_dupmarked
		set val(sample), file("*.marked_duplicates") into picard_duplicate

	shell:

		tmp_dirname = "tmp"
		filename = sample + ".dupmarked.bam"

		// the metrics file
		metrics_filename = sample + ".marked_duplicates"

		"""
		java -Xmx10g -Djava.io.tmpdir=${tmp_dirname} \
			-jar \$EBROOTPICARD/picard.jar MarkDuplicates \
			VALIDATION_STRINGENCY=SILENT \
			INPUT=${bam} \
			OUTPUT=${filename} \
			METRICS_FILE=${metrics_filename} \
			ASSUME_SORTED=true \
			REMOVE_DUPLICATES=false \
			TMP_DIR=${tmp_dirname}
		"""
}


///////////////////////////////////////////////////////////////////////////////
/* --                                                                     -- */
/* --                            FORKING                                  -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////

picard_dupmarked_indexing = Channel.create()
picard_dupmarked_complexity = Channel.create()
picard_dupmarked_rnaseqmetrics = Channel.create()
picard_dupmarked_multimetrics = Channel.create()
picard_dupmarked_infer_experiment = Channel.create()
picard_dupmarked_junction_annotation = Channel.create()
picard_dupmarked_junction_saturation = Channel.create()
picard_dupmarked_mismatch_profile = Channel.create()
picard_dupmarked_read_distribution = Channel.create()
picard_dupmarked_rnaseqc = Channel.create()

picard_dupmarked.into(
	picard_dupmarked_indexing,
	picard_dupmarked_complexity,
	picard_dupmarked_rnaseqmetrics,
	picard_dupmarked_multimetrics,
	picard_dupmarked_infer_experiment,
	picard_dupmarked_junction_annotation,
	picard_dupmarked_junction_saturation,
	picard_dupmarked_mismatch_profile,
	picard_dupmarked_read_distribution,
	picard_dupmarked_rnaseqc)


///////////////////////////////////////////////////////////////////////////////
/* --                                                                     -- */
/* --                          INDEXING PICARD                            -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////

process index_picard {

	// custom label
	tag { sample }

	// problem with slurm otherwise
	beforeScript "module purge"

	// modules
	module MODULE_SAMTOOLS

	// HPC
	cpus 32
	executor "slurm"
	memory "6000"

	// output directory
	publishDir picard_dir.toString(),
		mode: PUBLISHDIR_MODE,
		overwrite: PUBLISHDIR_OVERWRITE

	input:
		set val(sample), file(bam) from picard_dupmarked_indexing

	output:
		set val(sample), file(bam), file("*.bai") into bai_picard

	shell:
		"""
		samtools index ${bam}
		"""
}

// forking
bai_picard_transcript_integrity_number = Channel.create()
bai_picard_rnaseqc = Channel.create()
bai_picard_multiqc = Channel.create()
bai_picard.into(
	bai_picard_transcript_integrity_number,
	bai_picard_rnaseqc,
	bai_picard_multiqc
	)


///////////////////////////////////////////////////////////////////////////////
/* --                                                                     -- */
/* --                               PICARD                                -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////

process picard_complexity {

	// custom label
	tag { sample }

	// problem with slurm otherwise
	beforeScript "module purge"

	// modules
	module MODULE_PICARD

	// HPC
	cpus 32
	executor "slurm"
	memory "6000"

	// output directory
	publishDir picard_dir.toString(),
		mode: PUBLISHDIR_MODE,
		overwrite: PUBLISHDIR_OVERWRITE

	input:
		set val(sample), file(bam) from picard_dupmarked_complexity

	output:
		set val(sample), file("*.complexity") into picard_complexity

	shell:

		// the temporary dictory and the metrics file
		tmp_dirname = "tmp"
		metrics_filename = sample + ".complexity"

		"""
		java -Xmx10g -Djava.io.tmpdir=${tmp_dirname} \
			-jar \$EBROOTPICARD/picard.jar EstimateLibraryComplexity \
			VALIDATION_STRINGENCY=SILENT \
			INPUT=${bam} \
			OUTPUT=${metrics_filename} \
			TMP_DIR=${tmp_dirname}
		"""
}

process picard_rnaseqmetrics {

	// custom label
	tag { sample }

	// problem with slurm otherwise
	beforeScript "module purge"

	// modules
	module MODULE_PICARD

	// HPC
	cpus 32
	executor "slurm"
	memory "6000"

	// output directory
	publishDir picard_dir.toString(),
		mode: PUBLISHDIR_MODE,
		overwrite: PUBLISHDIR_OVERWRITE

	input:
		set val(sample), file(bam) from picard_dupmarked_rnaseqmetrics

	output:
		set val(sample), file("*.rnaseqmetrics") into picard_rnaseqmetrics

	shell:

		// the temporary dictory and the metrics file
		tmp_dirname = "tmp"
		metrics_filename = sample + ".rnaseqmetrics"

		"""
		java -Xmx10g -Djava.io.tmpdir=${tmp_dirname} \
			-jar \$EBROOTPICARD/picard.jar CollectRnaSeqMetrics \
			VALIDATION_STRINGENCY=SILENT \
			INPUT=${bam} \
			OUTPUT=${metrics_filename} \
			REF_FLAT=${ANNOT_REFFLAT_FILEPATH} \
			STRAND_SPECIFICITY=NONE \
			RIBOSOMAL_INTERVALS=${ANNOT_RRNA_INTERVAL_FILEPATH} \
			TMP_DIR=${tmp_dirname}
		"""
}

process picard_multimetrics {

	// custom label
	tag { sample }

	// problem with slurm otherwise
	beforeScript "module purge"

	// modules
	module MODULE_PICARD
	module MODULE_R

	// HPC
	cpus 32
	executor "slurm"
	memory "6000"

	// output directory
	publishDir picard_dir.toString(),
		mode: PUBLISHDIR_MODE,
		overwrite: PUBLISHDIR_OVERWRITE

	input:
		set val(sample), file(bam) from picard_dupmarked_multimetrics

	output:
		set val(sample), file("*.pdf") into picard_multimetrics_pdf
		set val(sample), file("*_metrics") into picard_multimetrics_metrics

	shell:

		// the temporary dictory and the metrics file
		tmp_dirname = "tmp"
		metrics_filename = sample + ".multimetrics"

		"""
		java -Xmx10g -Djava.io.tmpdir=${tmp_dirname} \
			-jar \$EBROOTPICARD/picard.jar CollectMultipleMetrics \
			VALIDATION_STRINGENCY=SILENT \
			INPUT=${bam} \
			OUTPUT=${metrics_filename} \
			PROGRAM=CollectAlignmentSummaryMetrics \
			R=${SEQ_FILEPATH} \
			TMP_DIR=${tmp_dirname}
		"""
}


///////////////////////////////////////////////////////////////////////////////
/* --                                                                     -- */
/* --                              RSEQC                                  -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////

process infer_experiment {

	// custom label
	tag { sample }

	// problem with slurm otherwise
	beforeScript "module purge"

	// modules
	module MODULE_RSEQC

	// HPC
	cpus 32
	executor "slurm"
	memory "6000"

	// output directory
	publishDir rseqc_dir.toString(),
		mode: PUBLISHDIR_MODE,
		overwrite: PUBLISHDIR_OVERWRITE

	input:
		set val(sample), file(bam) from picard_dupmarked_infer_experiment
	
	output:
		set val(sample), file("*.infer_experiment*") into infer_experiment
	
	shell:

		metrics_filename = sample + ".infer_experiment"

		"""
		infer_experiment.py \
			-i ${bam} \
			-r ${ANNOT_BED_FILEPATH} \
			> ${metrics_filename}
		"""
}

process junction_annotation {

	// custom label
	tag { sample }

	// problem with slurm otherwise
	beforeScript "module purge"

	// modules
	module MODULE_RSEQC

	// HPC
	cpus 32
	executor "slurm"
	memory "6000"

	// output directory
	publishDir rseqc_dir.toString(),
		mode: PUBLISHDIR_MODE,
		overwrite: PUBLISHDIR_OVERWRITE

	input:
		set val(sample), file(bam) from picard_dupmarked_junction_annotation

	output:
		set val(sample), file("*.junction_annotation*") into junction_annotation
	
	shell:

		metrics_filename = sample + ".junction_annotation"

		"""
		junction_annotation.py \
			-i ${bam} \
			-r ${ANNOT_BED_FILEPATH} \
			-o ${metrics_filename}
		"""
}

process junction_saturation {

	// custom label
	tag { sample }

	// problem with slurm otherwise
	beforeScript "module purge"

	// modules
	module MODULE_RSEQC

	// HPC
	cpus 32
	executor "slurm"
	memory "6000"

	// output directory
	publishDir rseqc_dir.toString(),
		mode: PUBLISHDIR_MODE,
		overwrite: PUBLISHDIR_OVERWRITE

	input:
		set val(sample), file(bam) from picard_dupmarked_junction_saturation

	output:
		set val(sample), file("*.junction_saturation*") into junction_saturation

	shell:

		metrics_filename = sample + ".junction_saturation"

		"""
		junction_saturation.py \
			-i ${bam} \
			-r ${ANNOT_BED_FILEPATH} \
			-o ${metrics_filename}
		"""
}

process mismatch_profile {

	// custom label
	tag { sample }

	// problem with slurm otherwise
	beforeScript "module purge"

	// modules
	module MODULE_RSEQC

	// HPC
	cpus 32
	executor "slurm"
	memory "6000"

	// output directory
	publishDir rseqc_dir.toString(),
		mode: PUBLISHDIR_MODE,
		overwrite: PUBLISHDIR_OVERWRITE

	input:
		set val(sample), file(bam) from picard_dupmarked_mismatch_profile

	output:
		set val(sample), file("*.mismatch_profile*") into mismatch_profile

	shell:

		metrics_filename = sample + ".mismatch_profile"

		"""
		mismatch_profile.py \
			-i ${bam} \
			-l ${ROUGH_READ_LENGTH} \
			-o ${metrics_filename}
		"""
}

process read_distribution {

	// custom label
	tag { sample }

	// problem with slurm otherwise
	beforeScript "module purge"

	// modules
	module MODULE_RSEQC

	// HPC
	cpus 32
	executor "slurm"
	memory "6000"

	// output directory
	publishDir rseqc_dir.toString(),
		mode: PUBLISHDIR_MODE,
		overwrite: PUBLISHDIR_OVERWRITE

	input:
		set val(sample), file(bam) from picard_dupmarked_read_distribution

	output:
		set val(sample), file("*.read_distribution*") into read_distribution
	
	shell:

		metrics_filename = sample + ".read_distribution"

		"""
		read_distribution.py \
			-i ${bam} \
			-r ${ANNOT_BED_FILEPATH} \
			> ${metrics_filename}
		"""
}

process transcript_integrity_number {

	// custom label
	tag { sample }

	// problem with slurm otherwise
	beforeScript "module purge"

	// modules
	module MODULE_RSEQC

	// HPC
	cpus 32
	executor "slurm"
	memory "6000"

	// output directory
	publishDir rseqc_dir.toString(),
		mode: PUBLISHDIR_MODE,
		overwrite: PUBLISHDIR_OVERWRITE

	input:
		set val(sample), file(bam),
			file(bai) from bai_picard_transcript_integrity_number

	output:
		set val(sample), file("*.{xls,txt}") into transcript_integrity_number

	shell:
		"""
		tin.py \
			-i ${bam} \
			-r ${ANNOT_BED_FILEPATH}
		"""
}


///////////////////////////////////////////////////////////////////////////////
/* --                                                                     -- */
/* --                             RNASEQC                                 -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////

process rnaseqc {

	// custom label
	tag { sample }

	// problem with slurm otherwise
	beforeScript "module purge"

	// modules
	module MODULE_RNASEQC

	// HPC
	cpus 32
	executor "slurm"
	memory "6000"

	// output directory
	publishDir rnaseqc_dir.toString(),
		mode: PUBLISHDIR_MODE,
		overwrite: PUBLISHDIR_OVERWRITE

	input:
		set val(sample), file(bam), file(bai) from bai_picard_rnaseqc

	output:
		set val(sample), file("*rnaseqc*") into rnaseqc

	shell:

		if (SINGLE_END) {

			metrics_filename = sample + ".rnaseqc"

			"""
			java -Xmx10g -jar \$EBROOTRNAMINSEQC/RNA-SeQC_v[0-9].[0-9].[0-9].jar \
				-d 1000000 \
				-rRNA ${ANNOT_RRNA_FILEPATH} \
				-r ${SEQ_FILEPATH} \
				-t ${ANNOT_RNASEQC_FILEPATH} \
				-o ${metrics_filename} \
				-singleEnd \
				-s "${sample}|${bam}|${sample}"
			"""

		} else {

			metrics_filename = sample + ".rnaseqc"

			"""
			java -Xmx10g -jar \$EBROOTRNAMINSEQC/RNA-SeQC_v[0-9].[0-9].[0-9].jar \
				-d 1000000 \
				-rRNA ${ANNOT_RRNA_FILEPATH} \
				-r ${SEQ_FILEPATH} \
				-t ${ANNOT_RNASEQC_FILEPATH} \
				-o ${metrics_filename} \
				-s "${sample}|${bam}|${sample}"
			"""
		}

}


///////////////////////////////////////////////////////////////////////////////
/* --                                                                     -- */
/* --                            MULTIQC                                  -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////

process multiqc {

	echo true

	// problem with slurm otherwise
	beforeScript "module purge"

	// module
	module MODULE_MULTIQC

	// HPC
	cpus 1
	executor "slurm"
	memory "6000"

	// output directory
	publishDir log_dir.toString(),
		mode: PUBLISHDIR_MODE,
		overwrite: PUBLISHDIR_OVERWRITE

	input:
		file("cutadapt/*") from cutadapt_log.map{x->x[1]}.collect()

		file("*.html") from fastq_screen_html.map{x->x[1]}.collect()
		file("*.txt") from fastq_screen_txt.map{x->x[1]}.collect()

		file("star/*") from bam_star_multiqc.map{x->x[1]}.collect()
		file("star/*") from transcript_bam_star.map{x->x[1]}.collect()
		file("star/*") from star_results_multiqc.map{x->x[1]}.collect()
		file("star/*") from star_stat.map{x->x[1]}.collect()

		//file("mbscreen/*") from bam_star_mbscreen.map{x->x[1]}.collect()
		//file("mbscreen/*") from transcript_bam_star_mbscreen.
		//	map{x->x[1]}.collect()
		//file("mbscreen/*") from results_star_mbscreen.map{x->x[1]}.collect()
		//file("mbscreen/*") from stat_star_mbscreen.map{x->x[1]}.collect()

		file("sort_index_star/*") from bam_star_sorted_multiqc.
			map{x->x[1]}.collect()

		file("picard/*") from picard_rg_multiqc.map{x->x[1]}.collect()
		file("picard/*") from picard_duplicate.map{x->x[1]}.collect()
		file("picard/*") from bai_picard_multiqc.map{x->x[1]}.collect()
		file("picard/*") from picard_complexity.map{x->x[1]}.collect()
		file("picard/*") from picard_rnaseqmetrics.map{x->x[1]}.collect()
		file("picard/*") from picard_multimetrics_pdf.map{x->x[1]}.collect()
		file("picard/*") from picard_multimetrics_metrics.map{x->x[1]}.collect()
		
		file("rseqc/*") from infer_experiment.map{x->x[1]}.collect()
		file("rseqc/*") from junction_annotation.map{x->x[1]}.collect()
		file("rseqc/*") from junction_saturation.map{x->x[1]}.collect()
		file("rseqc/*") from mismatch_profile.map{x->x[1]}.collect()
		file("rseqc/*") from read_distribution.map{x->x[1]}.collect()
		file("rseqc/*") from transcript_integrity_number.map{x->x[1]}.collect()

		file("rnaseqc/*") from rnaseqc.map{x->x[1]}.collect()

	output:
		file "multiqc_data" into multiqc_data
		file "multiqc_report.html" into multiqc_report

	shell:
		
		if (MULTIQC_CONFIG_FILEPATH) {

			"""
			multiqc -f --config ${MULTIQC_CONFIG_FILEPATH} ./
			"""

		} else {

			"""
			multiqc -f ./
			"""
		}
}


///////////////////////////////////////////////////////////////////////////////
/* --                                                                     -- */
/* --                           NOTIFICATION                              -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////

workflow.onComplete {

	// the bioinformatician
	def USER = System.getenv("USER")

	// might need to be changed but don't want to give name of the
	// bioinformatician as a parameter
	def recipient = USER + "@crick.ac.uk"

	// if it failed or not
	def subject = "[nexflow] the workflow has successfuly finished"
	if (!workflow.success) {
		subject = "[nextflow] PROBLEM: the workflow has failed"
	}

	// just the strict necessary, need to add more
	def body =
	"""
	Start time: ${workflow.start}
	End time: ${workflow.complete}
	Exit status: ${workflow.exitStatus}
	"""

	// if it failed... when it fails
	if (workflow.errorMessage) {
		body = body +
	"""
	Error message: ${workflow.errorMessage}
	"""
	}

	// if it failed... when it fails
	if (workflow.errorReport) {
		body = body +
	"""
	Report message: ${workflow.errorReport}
	"""
	}
	
	// there is only mutt installed on the node
	["mutt", "-s", subject, recipient].execute() << body
}

