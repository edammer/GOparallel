# GOparallel
Uses BaderLab's monthly updated .GMT formatted ontology gene lists for Fisher exact enrichment into differential or co-expression or other user-supplied gene lists.
GO.obo file download is required for redundant term pruning.
WGCNA modules can be co-clustered for relatedness of modules by GO-cellular component terms hit or depleted per signed Z score calculation.

Compatible with the Seyfried systems biology pipeline.

See Instructions in the R function wrapper code to set all needed variables in the global environment before the function call.

Sample input for WGCNA modules, or ANOVA-Tukey table + volcano settings from pipeline outputs are in <a href="https://github.com/edammer/GOparallel/blob/main/GOparallel-SampleOutput.zip">this .zip</a>, and can be loaded prior to GOparallel using the provided sample RData for either. These sample analysis inputs used the 8619 fully corrected, age, sex, and PMI-regressed protein abundance data from <a href="https://www.nature.com/articles/s41593-021-00999-y">Johnson ECB, et al, Nat Neurosci., 2022</a>, processed through the Seyfried Lab systems biology pipeline, of which GOparallel is a late step.

Note the current analysis output samples are using June 2022 ontologies, much newer compared to those in the publication, which was based on the GO-Elite Ensembl v62+ database.

Requires R packages: piano, WGCNA, curl, doParallel, rvest, ontologyIndex, NMF.


Wrapper code for function GOparallel():
```
## Preload WGCNA standard pipeline R ression's data to memory, if desired (example minimal RData given here)
## otherwise below code runs as a step late in the Seyfried systems biology pipeline,
## or following the pipeline's volcano code block, or with .csv input formatted as simple lists
## in columns, or the .csv may be kME table output from the global network plots pipeline.
###############################################################################################

setwd("c:/BaderLabGO/")

## Sample pipeline outputs required to run as inputs for input scenario #2 below  -- contains empty cleanDat (only ordered rownames matching net$colors are needed here)
load("./PipelineSample-MEGATMT-WGCNA-NatNeurosci2022.RData")

## Sample pipeline outputs required to run as inputs for input scenario #3 below
load("./PipelineSample-MEGATMT-ANOVA+Volcano-NatNeurosci2022.RData")


###################################################################################################################################
## One-Step GSA FET (R piano implementation) WITH USER PARAMETERS
##  - by Eric Dammer, Divya Nandakumar
## Nicholas Seyfried Lab Bioinformatics - for the lab - 06/10/2022 version 1.0
####################################################################################################################################
### The piano package can be installed from bioconductor, enabling ontology enrichment analysis on gene sets. 
### Documentation for this package is found https://rdrr.io/bioc/piano/ before installation.
### Other dependencies for code below include doParallel and stringr packages.
###
### .GMT standard tab-separated format databases can be downloaded from http://download.baderlab.org/EM_Genesets/current_release/
### where the 3 core GO ontology types are updated monthly, and additional curated pathways are collected from various sources,
### e.g. the Molecular Signatures (MSig) C2 database, Reactome, wikipathways, etc., which are updated less regularly.
###
### If GO redundant term removal enabled, ontologyIndex package and the full go.obo will be specified, should have been downloaded
### from http://current.geneontology.org/ontology/go.obo
###
####################################################################################################################################
##  *Functional Description*
##  - performs hypergeometric overlap with Fisher Exact Test for enrichment p<0.05
##    and 5 minimum genes per ontology
##  - currently outputs files that can be used as GO-Elite input, not read by this script...
##  - outputs Z score barplots for each input list (module, or DiffEx Protein list) across 6 ontology types in the larger GMT files
##    available from the Bader Lab website.
##
##  *INSTRUCTIONS*
##  - Download "[SPECIES]_GO_AllPathways_with_GO_iea_[DATE]_symbol.gmt" from http://download.baderlab.org/EM_Genesets/current_release/
##  - (optional) Download go.obo from http://current.geneontology.org/ontology/go.obo
##    (Enables removal of redundant terms with option removeRedundantGOterms=TRUE)
##  - Change parameters between the lines below, and execute GOparallel("") from function in source file GOparallel-FET.R.
##
##  - Input can be:
##    1) a .csv file with uniqueIDs (Symbol|otherID(s)...) or species-appropriate gene symbols in columns.
##       (Each will be considered as a list of DEPs or module members); specify fileName and have below two flags =FALSE.
##       alternately, the input can be a WGCNA kME table output from the Seyfried systems pipeline saved as .csv.
##    2) cleanDat, net$colors, and kMEdat from Seyfried systems biology pipeline loaded in memory              (modulesInMemory=TRUE)
##    3) ANOVA/T-Test statistics table with Seyfried volcano code block data structures loaded in memory       (ANOVAgroups=TRUE)
##
##  - Outputs are:
##    1) Z score barplots with min z ~1.96 (p<0.05) calculated from Fisher Exact Test ORA (one-tailed) p values. One barplot is
##       created for each input list (module, or up DEGs, or down DEGs for each ANOVA/T-test comparison in ANOVAout).
##    2) Tab-separated .txt table of all lists' signed Zscores (each list is a column) for all ~13000 ontologies.
##       Z scores are converted from p values as two-tailed p values, since each enrichment and depletion p value is calculated
##       to get signed Z scores (less than zero if depleted more than by chance).  I.e. the conversion to Z uses p/2.
##    3) Tab-separated .txt table of all lists' Enrichment p values
##    4) Tab-separated .txt table of all lists' Enrichment Benjamini-Hochberg corrected FDR values.
##    5) (optional) GO:Cellular Component ontologies with at least one highly significant list enrichment are coclustered
##       in a heatmap of signed Zscores across all lists. Works well with WGCNA module input.                   (cocluster=TRUE)
##
#####################################################################################################################################

dev.off()
options(stringsAsFactors=FALSE)

######################## EDIT THESE VARIABLES (USER PARAMETERS SET IN GLOBAL ENVIRONMENT) ############################################
#fileName <- "ENDO_MG_TWO_WAY_LIST_NTS_v02b_forGOelite.csv"                                            #Sample File 1 - has full human background
#fileName <- "ModuleAssignments_Jingting32TH_BOOTaspRegr_power8_MergeHeight0.07_PAMstageTRUE_ds2.csv"  #Sample File 2 - WGCNA kME table for (Dai, et al, 2019)
            #INPUT CSV FILE - in the filePath folder.
            #Can be formatted as Kme table from WGCNA pipeline, or
            #can be a CSV of columns, one symbol or UniqueID (Symbol|...) list per column, with the LIST NAMEs in row 1
            #in this case, the longest list is used as background or the "universe" for the FET contingencies
            #  For simple columnwise list input, DON'T FORGET TO PUT THE APPROPRIATE BACKGROUND LIST IN, OR RESULTS WILL BE UNRELIABLE.

filePath <- "c:/BaderLabGO/"   #gsub("//","/",outputfigs)
            #Folder that (may) contain the input file specified above, and which will contain the outFilename project Folder.

outFilename <- "MEGATMT"
            #SUBFOLDER WITH THIS NAME WILL BE CREATED, and .PDF + .csv file using the same name will be created within this folder.

outputGOeliteInputs=FALSE
            #If TRUE, GO Elite background file and module or list-specific input files will be created in the outFilename subfolder.

maxBarsPerOntology=5
            #Ontologies per ontology type, used for generating the PDF report; does not limit tabled output

GMTdatabaseFile="c:/BaderLabGO/current.gmt"   #e.g. "Human_GO_AllPathways_with_GO_iea_June_01_2022_symbol.gmt"
                                              # Current month release will be download if file does not exist.
            #drive:folder/to/filename of ontology database for the appropriate species (no conversion is performed)
            #BaderLab website links to their current monthly update of ontologies for Human, Mouse, and Rat, minimally
            #http://download.baderlab.org/EM_Genesets/current_release/
            #For more information, see documentation:  http://baderlab.org/GeneSets

panelDimensions=c(3,2)    #dimensions of the individual parblots within a page of the main barplot PDF output
pageDimensions=c(8.5,11)  #main barplot PDF output page dimensions, in inches

color=c("darkseagreen3","lightsteelblue1","lightpink4","goldenrod","darkorange","gold")
            #colors respectively for ontology Types:
            #"Biological Process","Molecular Function","Cellular Component","Reactome","WikiPathways","MSig.C2"
            #must be valid R colors

modulesInMemory=TRUE
            #uses cleanDat, net, and kMEdat from WGCNA systems biology pipeline already run, and these variables must be in memory
            #input fileName will be ignored
ANOVAgroups=FALSE
            #if true, modulesInMemory ignored. Volcano pipeline code should already have been run!
            #input fileName will be ignored

############ MUST HAVE AT LEAST 2 THREADS ENABLED TO RUN ############################################################################

parallelThreads=30

removeRedundantGOterms=TRUE
            #if true, the 3 GO ontology types are collapased into a minimal set of less redundant terms using the below OBO file
GO.OBOfile<-"C:/BaderLabGO/go.obo"
            #only used and needed if above flag to remove redundant GO terms is TRUE.
	    #Download from http://current.geneontology.org/ontology/go.obo will commence if needed.
            #Does not appear to be species specific, stores all GO term relations and is periodically updated.

cocluster=TRUE
            #If TRUE, output PDF of signed Zscore coclustering on GO:cellular component terms (useful for WGCNA modules)
	    #Note: Currently fails after bar plot output completes, if symbol lists are not valid R colors (modules)...

######################## END OF PARAMETER VARIABLES ###################################################################################

source("GOparallel-FET.R")
GOparallel("")  # parameters must be set in global environment as above; no defaults are specified for the function or passed to it.


#save.image(paste0("saved.image.GSA.GO.FET-",outFilename,"-completed.Rdata"))
