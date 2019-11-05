***********
Description
***********

HabitatCorrectedNetwork.py is a python module for calculating 
correlation matrices for microial correlation network analysis 
that account for habitat filtering effects.  

HabitatCorrectedNetwork.py is compatible with both Python2 and Python3.  
In addition to standard python libraries, it relies on the numpy, 
scipy, and statsmodels packages, which must also be installed.  

The tool takes in 2 files:
  - an ASV (or OTU) count/abundance data file
  - a sample data file
It saves the results to 3 files for each analysis 
  - a list of ASV (or OTU) names in the network 
  - a correlation matrix of r values filtered to include only the 
desired correlations based on input parameters
  - a corresponding matrix of p-values

In addition to performing habitat filtering correction on a specified 
habitat (sub group) variable, this tool also provides flexibility to 
customize the correlation analysis with the following options:
  - Division of a large data set into multiple smaller datasets for 
generation of individual networks
  - Centered log ratio transformation of the data to account for 
compositional effects (particularly important for low complexity data)
  - False discovery rate p-value adjustment
  - Spearman or Pearson correlation methods
  - Specification of a minimum ASV (or OTU) prevalence (proportion of 
samples in which an ASV (or OTU) must be detected to be included in the
analysis)
  - Specification of a minimum correlation value (r) and the option to 
use Spearman or Pearson correlation coefficients
  - Specification of a significance threshold (maximum p-value)
  - Selection of direction of correlations to keep (positive, negative, 
or both)

These options can be determined from the command line (default), with a 
configuration file, or interactively, depending on the setup method 
selected. 


****************************
Basic command line operation
****************************

python HabitatCorrectedNetwork.py -A <ASV data filename> 
-S <sample data filename> -out <output file prefix>

NOTE: With this basic operation above, the program assumes that all 
samples are from the same habitat. In order to perform habitat 
filtering correction, you also need to provide information on which 
sample variable identifies the habitats you are correcting for. You 
can do this by using the -hab flag, or by using the configuration file 
or custom interactive setup options (see optional parameters below)


**********************************************
Habitat Filtering Correction With Example Data
**********************************************

The example data are from rhizosphere samples from 2 maize accessions 
grown in 3 soil types. In this example we split the data into 2 
networks (1 for each maize accession) and do habitat filtering 
correction for the 3 soil types. These examples show 2 different ways 
to run the analysis with the same parameters, either by using a 
configuration file or vy entering parameters through the command line. 

python HabitatCorrectedNetwork.py -A example_ASV_table.txt -S example_sample_data.txt -out example.txt -s config -sfn example_config.cfg

python HabitatCorrectedNetwork.py -A example_ASV_table.txt -S example_sample_data.txt -out example.txt -split accession -hab soil -clr True -padj False -cm spearman -prev 0.5 -cor 0.75 -pval 0.01 -d pos


*******************
Required Parameters
*******************
-----------------------------------------------------------------------
Short flag 	Long flag
-----------------------------------------------------------------------
-A		--ASV_filename

Filename for a tab delimited text file containing ASV (or OTU) 
count/abundance data  
Any rarefaction or normalization of the data should be done prior to 
this analysis.

Format: 
First row lists ASV (or OTU) IDs. 
Subsequent rows start with the sample name followed by ASV (or OTU) 
counts/abundances for that sample
-----------------------------------------------------------------------
-S		--sample_filename	

Filename for a tab delimited text file containing sample variable data
This data is used to split samples into groups for individual network 
analysis if desired, and to define sub groups (habitats) for habitat 
filtering correction.

Format: 
First row lists column names for sample variables (must include 
sample_name as one of the columns). 
Subsequent rows give data for each sample.
-----------------------------------------------------------------------
-out		--output_prefix	
	
Prefix for the output file names
-----------------------------------------------------------------------


*******************
Optional Parameters
*******************
-----------------------------------------------------------------------
Short flag	Long flag		Default value	Possible values
-----------------------------------------------------------------------
-s		--setup			default		default, custom, 
							config
Method to use to set up the analysis parameters

default: Start out with the default parameters. Individual parameters 
can be changed with additional command line arguments. 
custom: Interactively input the individual parameters.
config: Use a configuration file to define the parameters. If this 
option is selected, a valid configuration file must be supplied.

NOTE: Any individual parameter can then be overridden with the command 
line flags below
-----------------------------------------------------------------------
-sfn		--setup_filename			valid setup 
							filename	
Filename of the setup configuration file
Requited if --setup is config, otherwise ignored
-----------------------------------------------------------------------
-split		--split_group		none		none or valid
 							sample variable 
							name
Sample variable to use to split the data into multiple networks
If not none, this must be a valid sample variable (column name) in the 
sample data file.
-----------------------------------------------------------------------
-hab		--habitat_group		none		none or vaild 
							sample variable 
							name
Sample variable to define habitats for habitat filtering correction  
If not none, this must be a valid sample variable (column name) in the 
sample data file.
-----------------------------------------------------------------------
-clr		--CLR			True		True, False

Whether or not to do a centered log ratio transformation of the data in 
order to reduce compositionality effects.  This is particularly an 
issue with low complexity microbial communities.
-----------------------------------------------------------------------
-padj		--pval_adjustment	False		True, False

Whether or not to do a false discovery rate correction for multiple 
comparisons
-----------------------------------------------------------------------
-cm		--corr_method		spearman	spearman, 
							pearson
Which correlation method to use, Spearman or Pearson
-----------------------------------------------------------------------
-prev		--min_prev		0.66		value between 0 
							and 1
Minimum prevalence for an ASV (or OTU) in the network, specified as the
proportion of samples in which an ASV (or OTU) must be detected to be 
included in the analysis
ASVs (or OTUs) that only appear in a small proportion of samples can be 
problematic for identifying correlations since their abundances are 
largely below detection limits.
-----------------------------------------------------------------------
-cor		--min_corr		0		value between 0 
							and 1
Minimum correlation r value to keep a correlation as a connection in 
the network
-----------------------------------------------------------------------
-pval		--max_pval		0.05		value between 0
							and 1
Maximum p-value to keep a correlation as a connection in the network
-----------------------------------------------------------------------
-d		--direction		pos		pos, neg, both

Direction of correlations to keep as connections in the network 
Positive correlations indicate co-occurrence while negative 
correlations indicate co-exclusion.
-----------------------------------------------------------------------

