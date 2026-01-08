# GOparallel
Uses BaderLab's monthly updated .GMT formatted ontology gene lists for Fisher exact enrichment into differential or co-expression or other user-supplied gene lists.
GO.obo file download is required for redundant term pruning.
WGCNA modules can be co-clustered for relatedness of modules by GO-cellular component terms hit or depleted per signed Z score calculation.

Compatible with the Seyfried systems biology pipeline.

See Instructions below and use R wrapper code to set all needed variables in the global environment before the function call. For 2 input options, Seyfried analysis pipeline should already have been run, or outputs loaded into memory.

Sample input for <a href="https://github.com/edammer/GOparallel/raw/main/PipelineSample-MEGATMT-WGCNA-NatNeurosci2022.RData">WGCNA modules</a>, or <a href="https://github.com/edammer/GOparallel/raw/main/PipelineSample-MEGATMT-ANOVA%2BVolcano-NatNeurosci2022.RData">ANOVA-Tukey table + volcano settings</a> from pipeline outputs can be downloaded and loaded from .RData files prior to using this wrapper. Outputs will reproduce what is contained in <a href="https://github.com/edammer/GOparallel/blob/main/GOparallel-SampleOutput.zip">this .zip</a>. These sample analysis inputs use the 8619 fully corrected, age, sex, and PMI-regressed protein abundance data from <a href="https://www.nature.com/articles/s41593-021-00999-y">Johnson ECB, et al, Nat Neurosci., 2022</a>, processed through the Seyfried Lab systems biology pipeline, of which GOparallel is a late step.

Note the current analysis output samples in the above linked .zip use June 2022 ontologies, much newer compared to those in the publication, which was based on the GO-Elite Ensembl v62+ database.
## **Functional Description**
  - performs hypergeometric overlap with Fisher Exact Test for enrichment p<0.05
    and 5 minimum genes per ontology
  - currently outputs files that can be used as GO-Elite input, not read by this script...
  - outputs Z score barplots for each input list (module, or DiffEx Protein list) across 6 ontology types in the larger GMT files
    available from the Bader Lab website.

## INSTRUCTIONS
  - Download <a href="https://github.com/edammer/GOparallel/raw/main/GOparallel-FET.R">GOparallel-FET.R</a> source, change variables serving as parameters in the wrapper code below and run GOparallel().

  - Input can be:
    1) a .csv file with uniqueIDs (Symbol|otherID(s)...) or species-appropriate gene symbols in columns.
       (Each will be considered as a list of DEPs or module members); specify fileName and have below two flags =FALSE.
       alternately, the input can be a WGCNA kME table output from the Seyfried systems pipeline saved as .csv.
    2) cleanDat, net$colors (or NETcolors) variables from Seyfried systems biology pipeline loaded in memory     (modulesInMemory=TRUE)
    3) ANOVA/T-Test statistics table with <a href="https://github.com/edammer/parANOVA">Seyfried volcano code block (plotVolc function)</a> data structures loaded in memory       (ANOVAgroups=TRUE)

    #3 Uses ANOVAout data in memory or passed to GOparallel(), and previous volcano selection criteria for DEx gene products if the plotVolc() function was already run from the <a href="https://github.com/edammer/parANOVA">parANOVA.dex source</a>.

  - Outputs are:
    1) Z score barplots with min z ~1.96 (p<0.05) calculated from Fisher Exact Test ORA (one-tailed) p values. One barplot is
       created for each input list (module, or up DEGs, or down DEGs for each ANOVA/T-test comparison in ANOVAout).
    2) Tab-separated .txt table of all lists' signed Zscores (each list is a column) for all ~13000 ontologies.
       Z scores are converted from p values as two-tailed p values, since each enrichment and depletion p value is calculated
       to get signed Z scores (less than zero if depleted more than by chance).  I.e. the conversion to Z uses p/2.
    3) Tab-separated .txt table of all lists' Enrichment p values
    4) Tab-separated .txt table of all lists' Enrichment Benjamini-Hochberg corrected FDR values.
    5) (optional) GO:Cellular Component ontologies with at least one highly significant list enrichment are coclustered
       in a heatmap of signed Zscores across all lists. Works well with WGCNA module input.                   (cocluster=TRUE)
    6) GSEA-style bubble plots of enrichment significance and gene (coverage) ratio, one PDF for each input list
                                                                                                              (bubble=TRUE)
##


<b>Wrapper code for function GOparallel():</b>
<BR>(All variables set in global environment, but defaults are tried if none are specified).
```
###################################################################################################################################
## One-Step GSA FET (R piano implementation) WITH USER PARAMETERS
##  - by Eric Dammer, Divya Nandakumar
## Nicholas Seyfried Lab Bioinformatics - for the lab - 07/27/2022 version 1.04
###################################################################################################################################
## Preload WGCNA standard pipeline R ression's data to memory, if desired (example minimal RData given here)
## otherwise below code runs as a step late in the Seyfried systems biology pipeline,
## or following the pipeline's volcano code block, or with .csv input formatted as simple lists
## in columns, or the .csv may be kME table output from the global network plots pipeline.
###############################################################################################

dev.off()
options(stringsAsFactors=FALSE)

## Sample pipeline outputs required to run as inputs when modulesInMemory=TRUE and ANOVAgroups=FALSE below:
#  -- contains empty cleanDat (only ordered rownames matching net$colors are needed here)
load("c:/BaderLabGO/PipelineSample-MEGATMT-WGCNA-NatNeurosci2022.RData")

## Sample pipeline outputs required to run as inputs when ANOVAgroups=TRUE below:
load("c:/BaderLabGO/PipelineSample-MEGATMT-ANOVA+Volcano-NatNeurosci2022.RData")


######################## EDIT THESE VARIABLES (USER PARAMETERS SET IN GLOBAL ENVIRONMENT) ############################################
#inputFile <- "ENDO_MG_TWO_WAY_LIST_NTS_v02b_forGOelite.csv"                                            #Sample File 1 - has full human background
#inputFile <- "ModuleAssignments_Jingting32TH_BOOTaspRegr_power8_MergeHeight0.07_PAMstageTRUE_ds2.csv"  #Sample File 2 - WGCNA kME table for (Dai, et al, 2019)
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

GMTdatabaseFile="c:/BaderLabGO/current.gmt"   # e.g. "Human_GO_AllPathways_with_GO_iea_June_01_2022_symbol.gmt"
                                              # Current month release will be downloaded if file does not exist.
					      # **Specify a nonexistent file to always download the current database to this folder.**
					      # Database .GMT file will be saved to the specified folder with its original date-specific name.
            #path/to/filename of ontology database for the appropriate species (no conversion is performed)
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
            #inputFile will be ignored
ANOVAgroups=FALSE
            #if true, modulesInMemory ignored. Volcano pipeline code should already have been run!
            #inputFile will be ignored

############ MUST HAVE AT LEAST 2 THREADS ENABLED TO RUN ############################################################################

parallelThreads=30

removeRedundantGOterms="kappa"
            # if TRUE, ontologyIndex parent terms are kept as a minimal term set (less redundancy for 3 main ontology types);
            # if "kappa", graph-based clustering redundancy removal is performed in parallel, reducing similar ontologies
            #   within each of the 6 types to keep the one best enrichment scorer for the cluster of terms with kappa <0.30.
            # if FALSE, no redundancy removal is performed.

GO.OBOfile<-"C:/BaderLabGO/go.obo"
            #only used and needed if above flag to remove redundant GO terms is TRUE.
	    #Download from http://current.geneontology.org/ontology/go.obo will commence into the specified folder if the specified filename does not exist.
            #Does not appear to be species specific, stores all GO term relations and is periodically updated.

cocluster=TRUE
            #If TRUE, output PDF of signed Zscore coclustering on GO:cellular component terms (useful for WGCNA modules)

bubble=TRUE
           #if TRUE, GSEA style bubble plots for each input list's enriched ontologies for each of the 6 ontology types will be output to PDF

######################## END OF PARAMETER VARIABLES ###################################################################################

source("GOparallel-FET.R")
GOparallel()  # parameters are set in global environment as above; if not set, the function falls back to defaults and looks for all inputs available.
              # priority is given to modulesInMemory
	      
```
## Other Notes
Requires R packages: piano, WGCNA, curl, doParallel, rvest, ontologyIndex, igraph, NMF, stringr. Bubble plots also require ggplot2, dplyr, tibble, scales, grid, and gtable packages. The piano package (documentation <a href="https://rdrr.io/bioc/piano/">here</a>) can be installed from bioconductor, enabling ontology enrichment analysis on gene sets.

The function for Fisher exact test calculation used within GOparallel is derived from a similar function in the piano package, with its dependencies, but further calculates the p values of both enrichment and depletion so that we have signed Zscores available.
##
.GMT standard tab-separated format databases are downloaded after scraping available file links from http://download.baderlab.org/EM_Genesets/current_release/
for any of the available species there. According to <a href="http://baderlab.org/GeneSets">the site's documentation</a>, the 3 core GO ontology types are updated monthly, and additional curated pathways are collected from various sources, e.g. the Molecular Signatures (MSig) C2 database, Reactome, wikipathways, etc., which are updated less regularly.

In interactive sessions run from an R console, needed ontology database download will commence of [SPECIES]_GO_AllPathways_with_GO_iea_[DATE]_symbol.gmt from the BaderLab download website online folder structure after interactive user input to select the species. The interaction is skipped if the .GMT file specified in GMTdatabaseFile variable exists.  Note that repetitive bulk downloads of large files during working hours can trigger a ban of the downloader's IP.
  

If GO redundant term removal is enabled, ontologyIndex package will be loaded and the full go.obo file downloaded from <a href="http://current.geneontology.org/ontology/go.obo">http://current.geneontology.org/ontology/go.obo</a> if not already present as specified in the GO.OBOfile variable.

<b>Wrapper code for function GO.bubblePlot():</b>
<BR>(if post-hoc on output Zscore txt file from GOparallel)
```
########## PARAMETER SECTION ###############
ZinputFile="GSA-GO-FET_MEGATMT-Zscores.txt"
GMTfile="./Human_GO_AllPathways_noPFOCR_with_GO_iea_March_01_2025_symbol.gmt"
GO.OBOfile="./go.obo"
removeRedundantGOterms=TRUE

keep_ontology_types <- c("GOBP", "GOMF", "GOCC")  #, "REACTOME", "WikiPathways", "MSIGDB_C2") #, "MSIGDBHALLMARK" also useable, but not plotted in barplots by GOparallel()
#ontology_to_color <- c("GOBP"="darkseagreen3", "GOMF"="lightsteelblue1", "GOCC"="lightpink4",
#                       "REACTOME"="goldenrod", "WikiPathways"="darkorange", "MSIGDB_C2"="gold", "MSIGDBHALLMARK"="greenyellow")
    # default colors are used if not specified. Vector element names are required.

MAX_TERMS_PER_ONTOLOGY_TYPE <- 5
max_label_width <- 50          # Define maximum character width for wrapping ontology labels (change conservatively from 50, if at all)
#########################################

source("GOparallel-FET.R")  # includes GO.bubblePlot() function

GO.bubblePlot(ZinputFile=ZinputFile,
              GMTfile=GMTfile,
#              GO.OBOfile=GO.OBOfile,  # Will be downloaded if not specified or found in working directory.
              removeRedundantGOterms=TRUE,
              keep_ontology_types=keep_ontology_types,
#              ontology_to_color=ontology_to_color,  # Default colors are automatic.
              MAX_TERMS_PER_ONTOLOGY_TYPE=MAX_TERMS_PER_ONTOLOGY_TYPE,
              max_label_width=max_label_width)

```
