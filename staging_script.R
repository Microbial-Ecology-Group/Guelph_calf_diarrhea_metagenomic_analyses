## Start with this staging file to set up your analysis.
# Source the utility functions file, which should be in the scripts folder with this file
source('scripts/meg_utility_functions.R')
source('scripts/load_libraries.R')

# Set working directory to the MEG_R_metagenomic_analysis folder and add your data to that folder
#setwd("")

# Set the output directory for graphs:
graph_output_dir = 'graphs'
# Set the output directory for statistics:
stats_output_dir = 'stats'
# In which column of the metadata file are the sample IDs stored?
sample_column_id = 'ID'


## AMR analysis
## The files you want to use for input to this (for the MEG group analyses)
## is the AMR_analytic_matrix.csv. So you should have pulled these files from the output of the nextflow pipeline
## and you are now performing this analysis on your local machine.

## For the AMR analysis, you will also need to download the megares_annotations.csv
## file from the MEGARes website; the annotation file must be from the same version
## of the database as the file you used in the AmrPlusPlus pipeline, i.e. the headers
## must match between the annotation file and the database file.
# Where is the metadata file stored on your machine?
amr_metadata_filepath = 'weese_metadata.csv'
# Name of the megares annotation file used for this project
megares_annotation_filename = 'data/amr/megares_annotations_v1.02.csv'

# Load the data, MEGARes annotations, and metadata
amr <- newMRexperiment(read.table('SamDedup_AMR_analytic_matrix.csv', header=T, row.names=1, sep=','))
amr <- newMRexperiment(round(MRcounts(amr),0))
amr_temp_metadata <- read.csv(amr_metadata_filepath, header=T)
amr_temp_metadata[, sample_column_id] <- make.names(amr_temp_metadata[, sample_column_id])

# Annotation for regular AMR++ analysis
annotations <- data.table(read.csv(megares_annotation_filename, header=T))
setkey(annotations, header)  # Data tables are SQL objects with optional primary keys





###################
## User Controls ##
###################
## Hopefully, this section should be the only code you need to modify.
## However, you can look into the code in further sections if you need
## to change other, more subtle variables in the exploratory or
## statistical functions.

# to use for EXPLORATORY analysis (NMDS, PCA, alpha rarefaction, barplots)

# The following is a list of analyses based on variables in
# your metadata.csv file that you want
# to use for EXPLORATORY analysis (NMDS, PCA, alpha rarefaction, barplots)
AMR_exploratory_analyses = list(
  # Analysis 1
  # Description:
  list(
    name = 'Farm',
    subsets = list(),
    exploratory_var = 'Farm',
    order= ''
  ),
  # Analysis 2
  # Description:
  list(
    name = 'Time',
    subsets = list(),
    exploratory_var = 'Time',
    order= c("Pre","Post")
  ),
  # Analysis 4
  # Description:
  list(
    name = 'Group',
    subsets = list(),
    exploratory_var = 'Group',
    order = c("Farm 1 Pre","Farm 1 Post","Farm 2 Pre", "Farm 2 Post")
  ),
  # Farm 1 Analysis 1
  # Description:
  list(
    name = 'Farm1_by_Time',
    subsets = list('Farm == Farm1'),
    exploratory_var = 'Time',
    order= c("Pre","Post")
  ),
  # Farm 2 Analysis 1
  # Description:
  list(
    name = 'Farm2_by_Time',
    subsets = list('Farm == Farm2'),
    exploratory_var = 'Time',
    order=  c("Pre","Post")
  )
)

# Each analyses you wish to perform should have its own list in the following
# statistical_analyses list.  A template is provided to get you started.
# Multiple analyses, subsets, and contrasts are valid, but only one random
# effect can be used per analysis.  The contrasts of interest must have their
# parent variable in the model matrix equation.  Contrasts are named by
# parent variable then child variable without a space inbetween, for example:
# PVar1Cvar1 where the model matrix equation is ~ 0 + Pvar1.
AMR_statistical_analyses = list(
  # Analysis 1
  # Description:
  list(
    name = 'Time',
    subsets = list(),
    model_matrix = '~ 0 + Time',
    contrasts = list('TimePre - TimePost'),
    random_effect = NA
  ),
  # Analysis 2
  # Description:
  list(
    name = 'Group',
    subsets = list(),
    model_matrix = '~ 0 + Group_label',
    contrasts = list('Group_labelFarm1_Pre-Group_labelFarm1_Post','Group_labelFarm1_Pre-Group_labelFarm2_Pre','Group_labelFarm1_Pre-Group_labelFarm2_Post','Group_labelFarm1_Post-Group_labelFarm2_Post','Group_labelFarm1_Post-Group_labelFarm2_Pre','Group_labelFarm2_Pre-Group_labelFarm2_Post'),
    random_effect = NA
  ),
  # Analysis 3
  # Description:
  list(
    name = 'Time_Farm1',
    subsets = list('Farm == Farm1'),
    model_matrix = '~ 0 + Time',
    contrasts = list('TimePre - TimePost'),
    random_effect = NA
  ),
  # Analysis 4
  # Description:
  list(
    name = 'Time_Farm2',
    subsets = list('Farm == Farm2'),
    model_matrix = '~ 0 + Time',
    contrasts = list('TimePre - TimePost'),
    random_effect = NA
  )
)


## Run the analysis

## Run the script that handles resistome data.
source('scripts/metagenomeSeq_AMR.R')

# After running this script, these are the useful objects that contain all the data aggregated to different levels
# The metagenomeSeq objects are contained in these lists "AMR_analytic_data"
# Melted counts are contained in these data.table objects "amr_melted_analytic"

## Run code to make some exploratory figures, zero inflated gaussian model, and output count matrices.

#source('scripts/print_AMR_figures.R')
