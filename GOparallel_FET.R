## Preload WGCNA standard pipeline R ression's data to memory, if desired (example shown here)
## otherwise below code runs as a step late in the Seyfried systems biology pipeline,
## or following the pipeline's volcano code block, or with .csv input formatted as simple lists
## in columns, or the .csv may be kME table output from the global network plots pipeline.
###############################################################################################
rootdir="Q:/MEGATMT/Sample/"
outputfigs<-outputtabs<-rootdir
setwd(rootdir)
load("../PipelineSample-MEGATMT-NatNeurosci2022.RData")  #Sample pipeline outputs required to run as inputs type 2 or 3.
rootdir="Q:/MEGATMT/Sample/"
outputfigs<-outputtabs<-rootdir




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
##  - Change parameters between the lines below, and execute with code sub-blocks 1-3 which follow.
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

######################## EDIT THESE VARIABLES (USER PARAMETERS) #####################################################################
#fileName <- "ENDO_MG_TWO_WAY_LIST_NTS_v02b_forGOelite.csv"                                            #Sample File 1 - has full human background
#fileName <- "ModuleAssignments_Jingting32TH_BOOTaspRegr_power8_MergeHeight0.07_PAMstageTRUE_ds2.csv"  #Sample File 2 - WGCNA kME table for (Dai, et al, 2019)
            #INPUT CSV FILE - in the filePath folder.
            #Can be formatted as Kme table from WGCNA pipeline, or
            #can be a CSV of columns, one symbol or UniqueID (Symbol|...) list per column, with the LIST NAMEs in row 1
            #in this case, the longest list is used as background or the "universe" for the FET contingencies
            #  For simple columnwise list input, DON'T FORGET TO PUT THE APPROPRIATE BACKGROUND LIST IN, OR RESULTS WILL BE UNRELIABLE.

filePath <- "Q:/MEGATMT/Sample/"   #gsub("//","/",outputfigs)
            #Folder that (may) contain the input file specified above, and which will contain the outFilename project Folder.

outFilename <- "MEGATMT"
            #SUBFOLDER WITH THIS NAME WILL BE CREATED, and .PDF + .csv file using the same name will be created within this folder.

outputGOeliteInputs=FALSE
            #If TRUE, GO Elite background file and module or list-specific input files will be created in the outFilename subfolder.

maxBarsPerOntology=5
            #Ontologies per ontology type, used for generating the PDF report; does not limit tabled output

GMTdatabaseFile="c:/BaderLabGO/Human_GO_AllPathways_with_GO_iea_June_01_2022_symbol.gmt"
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
	    #Download from from http://current.geneontology.org/ontology/go.obo
            #Does not appear to be species specific, stores all GO term relations and is periodically updated.

cocluster=TRUE
            #If TRUE, output PDF of signed Zscore coclustering on GO:cellular component terms (useful for WGCNA modules)

######################## REST OF THE CODE IS AUTOMATIC ##############################################################################


## Clean out spaces and escaped backslashes from folder paths (folder names with spaces should not be used on non-windows systems with this script)
filePath=paste0(paste( sapply(do.call(c,strsplit(filePath,"[/\\]")),function(x) { if (grepl(" ",x)) { gsub(x,paste0(substr(gsub(" ","",x),1,6),"~1"),x) } else { x } } ),collapse="/"),"/")
#pythonPath=paste0(paste( sapply(do.call(c,strsplit(pythonPath,"[/\\]")),function(x) { if (grepl(" ",x)) { gsub(x,paste0(substr(gsub(" ","",x),1,6),"~1"),x) } else { x } } ),collapse="/"),"/")
#GOeliteFolder=paste0(paste( sapply(do.call(c,strsplit(GOeliteFolder,"[/\\]")),function(x) { if (grepl(" ",x)) { gsub(x,paste0(substr(gsub(" ","",x),1,6),"~1"),x) } else { x } } ),collapse="/"),"/")


## The files we create here are input files for GO-Elite, text files with the gene list as the 1st column, a symbol identified (gene symbol, uniprot etc) as the 2nd column
## Different accepted inputs are given in the tutorial
## Commonly used symbols - Gene Symbol - Sy (example of input file below)
### GeneSymbol		SystemCode (Symbol format)
###	  GFAP		Sy
###	  APOE		Sy
## All input files are placed in one folder

## The background file is prepared similarly and is placed in a separate folder
## The initial part of the code prepares files for GO-Elite. This can be skipped if the files are being made manually as described above.
## The second part of the code runs GO-ELite either from R (using the system command) or can be run using the terminal (in mac)
## The second part requires GO-Elite to be installed and path to the GO-Elite installation site indicated following python
## The 3rd part of the code plots the results from the GO-Elite results folder. When using the GUI the 1st 2 parts can be skipped and only the 3rd part can be used for plotting

##-------------------------------##
## Preparing files for GO-Elite ##
## Takes in the module assignment file as input with 1st column having gene names, 2nd column having color assignments followed by kME values


dir.create(file.path(filePath, outFilename))

##1a. Organize input gene lists, of the significant up and down (p<0.05) proteins in the current cleanDat
### List building from ANOVA-defined categories ###
if (ANOVAgroups) {
  sigThresh=0.05;

  #colnames(ANOVAout)
  #********************
  ##This code relies on a pre-existing data frame ANOVAout already in the environment, processed through volcano output for selected pairwise comparisons of interest.
  #numComp=6 #number of pairwise comparisons for ANOVA+Tukey p value columns, which are followed by the same-order log2(mean difference) columns
  #********************

  ANOVAout$Symbol <- do.call("rbind", strsplit(as.character(rownames(ANOVAout)), "[|]"))[, 1]

  DEXlistsForGO<-list()
  iter=0
  for (i in testIndexMasterList) {
    iter=iter+1;
    j=paste0(gsub(" ",".",comparisonIDs$Comparison[iter]),".down")
    k=paste0(gsub(" ",".",comparisonIDs$Comparison[iter]),".up")
    if (length(intersect(i,flip))==1) {
      #flipped sign (all diffs >0 for down)
      DEXlistsForGO[[j]]<-ANOVAout$Symbol[which(ANOVAout[,i]<sigThresh & ANOVAout[,i+numComp]>0)]
      DEXlistsForGO[[k]]<-ANOVAout$Symbol[which(ANOVAout[,i]<sigThresh & ANOVAout[,i+numComp]<0)]
    } else {
      #do not flip sign (all diffs <0 for down)
      DEXlistsForGO[[j]]<-ANOVAout$Symbol[which(ANOVAout[,i]<sigThresh & ANOVAout[,i+numComp]<0)]
      DEXlistsForGO[[k]]<-ANOVAout$Symbol[which(ANOVAout[,i]<sigThresh & ANOVAout[,i+numComp]>0)]
    }
  }

  #write lists to GOElite input files, and also the background file
  for (i in names(DEXlistsForGO)) { 
    dfGO<-data.frame(GeneSymbol=DEXlistsForGO[[i]],SystemCode=rep("Sy",length(DEXlistsForGO[[i]])))
    if(outputGOeliteInputs) write.table(unique(dfGO),file=paste(filePath,outFilename,"/",i,".txt",sep=""),row.names=FALSE,col.names=TRUE,sep="\t", quote=FALSE)
  }
  #write background
  background <- unique(ANOVAout$Symbol)
  background <- cbind(background,rep("Sy",length=length(background)))
  colnames(background) <- c("GeneSymbol","SystemCode")
  if(outputGOeliteInputs) dir.create(file.path(paste0(filePath,outFilename),"background"))
  if(outputGOeliteInputs) write.table(background,paste0(filePath,outFilename,"/background/background.txt"),row.names=FALSE,col.names=TRUE,quote=FALSE,sep="\t")
  nModules=length(names(DEXlistsForGO))
  WGCNAinput=FALSE

} else { #NOT creating lists from ANOVA/volcano up & down groups


##1b. Organize input gene lists from the WGCNA modules, either from specified input file or in the current cleanDat, net, and kME table
### List-building from WGCNA-defined modules ###

##use data structures in memory if modulesInMemory=TRUE; otherwise read csv following the template which is written in an earlier R session and edited in excel to produce module membership table saved as .csv:
if (modulesInMemory) {
  modulesData <- cbind(rownames(cleanDat),net$colors,kMEdat)
  WGCNAinput=TRUE
} else {
  modulesData <- read.csv(paste(filePath,fileName,sep=""),header=TRUE, sep=",");
  # check if this is a WGCNA modules/kME table or simple list input
  if(length(na.omit(match("net.colors",colnames(modulesData))))>0) {
    WGCNAinput=TRUE
    # Remove first and last column if it contains sort information ("Original Order" and "color|kMEin")
    modulesData <- modulesData[,-c(1,ncol(modulesData))];
  } else {
    WGCNAinput=FALSE
  }
}

if (WGCNAinput) {
  library(WGCNA) #for labels2colors

  # Include column with Symbol (if it is gene symbol, if not use appropriate code as given in GO-Elite manual)
  modulesData$SystemCode <- rep("Sy",nrow(modulesData)) 

  # Assign Names of First columns, in case they are non standard
  colnames(modulesData)[1]<-"Unique.ID" #This should have Symbol|UniprotID
  colnames(modulesData)[2]<-"net.colors" #This should have colors

  #Split out symbols from UniprotIDs, keep symbols in column 1
  rownames(modulesData)<-modulesData$Unique.ID
  modulesData$Unique.ID<-do.call("rbind",strsplit(as.character(modulesData$Unique.ID), "[|]"))[,1]

  ## Creating background file for GO Elite analysis
  background <- unique(modulesData[,"Unique.ID"])
  background <- cbind(background,rep("Sy",length=length(background)))
  colnames(background) <- c("GeneSymbol","SystemCode")
  if(outputGOeliteInputs) dir.create(file.path(paste0(filePath,outFilename),"background"))
  if(outputGOeliteInputs) write.table(background,paste0(filePath,outFilename,"/background/background.txt"),row.names=FALSE,col.names=TRUE,quote=FALSE,sep="\t")

  # Separate into independent module txt files for analysis by GO-Elite (CREATE INPUT FILES)
  greySubtractor=if(length(which(modulesData$net.colors=="grey"))>0) { 1 } else { 0 } #remove grey from count of modules
  nModules <- length(unique(modulesData$net.colors))-greySubtractor
  moduleColors <- uniquemodcolors <- labels2colors(c(1:nModules)) 
  for (i in 1:length(moduleColors)) {
    moduleName <- moduleColors[i]
    ind <- which(colnames(modulesData) == gsub("kMEME","kME",paste("kME",moduleName,sep="")))
    moduleInfo <- modulesData[modulesData$net.colors == gsub("ME","",moduleName), c(1,ncol(modulesData),ind)]
    colnames(moduleInfo) <- c("GeneSymbol","SystemCode","kME")
    if (moduleName == "blue" | moduleName == "brown" | moduleName == "green" | moduleName == "cyan") { if(outputGOeliteInputs) write.table(moduleInfo,file=paste(filePath,outFilename,"/",moduleName,"_2_Module.txt",sep=""),row.names=FALSE,col.names=TRUE,sep="\t", quote=FALSE)
    } else {
      if(outputGOeliteInputs) write.table(unique(moduleInfo),file=paste(filePath,outFilename,"/",moduleName,"_Module.txt",sep=""),row.names=FALSE,col.names=TRUE,sep="\t", quote=FALSE)
    }
  }
} else { #input is not WGCNA kME table format

##1c. List building from the WGCNA modules from specified input file, which must be in a column-wise list format, and including longest such list as background.
  # We process the input file as simple lists by column in the CSV (largest list used as background)

  #reread the file to a list of gene symbol (or UniqueID) lists
  modulesData <- as.list(read.csv(paste(filePath,fileName, sep=""),sep=",", stringsAsFactors=FALSE,header=T)) 

  nModules <- length(names(modulesData))
  for (a in 1:nModules) {
    modulesData[[a]] <- unique(modulesData[[a]][modulesData[[a]] != ""])
    modulesData[[a]] <- modulesData[[a]][!is.na(modulesData[[a]])]
    modulesData[[a]] <- do.call("rbind",strsplit(as.character(modulesData[[a]]), "[|]"))[,1]
  }
  ## Creating background file for GO Elite analysis
  background <- modulesData[order(sapply(modulesData,length),decreasing=TRUE)][[1]]
  background <- unique(background)
  background <- cbind(background,rep("Sy",length=length(background)))
  colnames(background) <- c("GeneSymbol","SystemCode")
  if(outputGOeliteInputs) dir.create(file.path(paste0(filePath,outFilename),"background"))
  if(outputGOeliteInputs) write.table(background,paste0(filePath,outFilename,"/background/background.txt"),row.names=FALSE,col.names=TRUE,quote=FALSE,sep="\t")

  # Separate Symbol Lists into independent module txt files for analysis by GO-Elite (not performed by this script)  (CREATE INPUT FILES)
  modulesData[[ names(modulesData[order(sapply(modulesData,length),decreasing=TRUE)])[1] ]] <- NULL
  nModules = nModules -1 #no background
  listNames <- uniquemodcolors <- names(modulesData)
  for (i in listNames) {
    listName <- i
    listInfo <- cbind(modulesData[[listName]],rep("Sy",length=length(modulesData[[listName]])))
    colnames(listInfo) <- c("GeneSymbol","SystemCode")
    if(outputGOeliteInputs) write.table(unique(listInfo),file=paste(filePath,outFilename,"/",listName,".txt",sep=""),row.names=FALSE,col.names=TRUE,sep="\t", quote=FALSE)
  }
} #end if (WGCNAinput)
} #end else for if (ANOVAgroups)




##2. GSA FET (parallelized within R, must have parallelThreads>1 to work currently)
####----------------------- piano package and dependencies required ------------------------------------#####

library(piano)

## Adapted version of piano::runGSAhyper() function with depletion p value also calculated (for signed Z score if we will use it to cocluster by module, e.g.)
runGSAhyper.twoSided <- function(genes, pvalues, pcutoff, universe, gsc, gsSizeLim = c(1,Inf), adjMethod = "fdr") {
    if (length(gsSizeLim) != 2) 
        stop("argument gsSizeLim should be a vector of length 2")
    if (missing(genes)) {
        stop("argument genes is required")
    }
    else {
        genes <- as.vector(as.matrix(genes))
        if (!is(genes, "character")) 
            stop("argument genes should be a character vector")
        if (length(unique(genes)) != length(genes)) 
            stop("argument genes should contain no duplicated entries")
    }
    if (missing(pvalues)) {
        pvalues <- rep(0, length(genes))
    }
    else {
        pvalues <- as.vector(as.matrix(pvalues))
        if (!is(pvalues, "numeric")) 
            stop("argument pvalues should be a numeric vector")
        if (length(pvalues) != length(genes)) 
            stop("argument pvalues should be the same length as argument genes")
        if (max(pvalues) > 1 | min(pvalues) < 0) 
            stop("pvalues need to lie between 0 and 1")
    }
    if (missing(pcutoff)) {
        if (all(pvalues %in% c(0, 1))) {
            pcutoff <- 0
        }
        else {
            pcutoff <- 0.05
        }
    }
    else {
        if (length(pcutoff) != 1 & !is(pcutoff, "numeric")) 
            stop("argument pcutoff should be a numeric of length 1")
        if (max(pcutoff) > 1 | min(pcutoff) < 0) 
            stop("argument pcutoff needs to lie between 0 and 1")
    }
    if (missing(gsc)) {
        stop("argument gsc needs to be given")
    }
    else {
        if (!is(gsc, "GSC")) 
            stop("argument gsc should be of class GSC, as returned by the loadGSC function")
    }
    if (missing(universe)) {
        if (!all(pvalues == 0)) {
            universe <- genes
            message("Using all genes in argument genes as universe.")
        }
        else {
            universe <- unique(unlist(gsc$gsc))
            message("Using all genes present in argument gsc as universe.")
        }
    }
    else {
        if (!is(universe, "character")) 
            stop("argument universe should be a character vector")
        if (!all(pvalues == 0)) 
            stop("if universe is given, genes should be only the genes of interest, i.e. pvalues should all be set to 0.")
    }
    if (!all(unique(unlist(gsc$gsc)) %in% universe)) 
        warning("there are genes in gsc that are not in the universe, these will be removed before analysis")
    if (!all(genes %in% universe)) {
        warning("not all genes given by argument genes are present in universe, these will be added to universe")
        universe <- c(universe, genes[!genes %in% universe])
    }
    if (length(unique(universe)) != length(universe)) 
        stop("argument universe should contain no duplicated entries")
    tmp <- try(adjMethod <- match.arg(adjMethod, c("holm", 
        "hochberg", "hommel", "bonferroni", 
        "BH", "BY", "fdr", "none"), several.ok = FALSE), 
        silent = TRUE)
    if (is(tmp, "try-error")) {
        stop("argument adjMethod set to unknown method")
    }
    pvalues[pvalues == 0] <- -1e-10
    goi <- genes[pvalues < pcutoff]
    if (length(goi) < 1) 
        stop("no genes selected due to too strict pcutoff")
    bg <- universe[!universe %in% goi]
    gsc <- gsc$gsc
    delInd <- vector()
    for (i in 1:length(gsc)) {
        gs <- gsc[[i]]
        gs <- gs[gs %in% universe]
        if (length(gs) < gsSizeLim[1] | length(gs) > gsSizeLim[2]) 
            delInd <- c(delInd, i)
        gsc[[i]] <- gs
    }
    gsc <- gsc[!c(1:length(gsc)) %in% delInd]
    message(paste("Analyzing the overrepresentation of ", 
        length(goi), " genes of interest in ", length(gsc), 
        " gene sets, using a background of ", length(bg), 
        " non-interesting genes.", sep = ""))
    p <- p.depletion <- rep(NA, length(gsc))
    names(p) <- names(p.depletion) <- names(gsc)
    padj <- rep(NA, length(gsc))
    names(padj) <- names(gsc)
    contTabList <- list()
    resTab <- matrix(nrow = length(gsc), ncol = 7)
    colnames(resTab) <- c("Pvalue.Enrichment", "Adjusted.Enr.Pvalue", "Pvalue.Depletion",
        "Significant (in gene set)", "Non-significant (in gene set)", 
        "Significant (not in gene set)", "Non-significant (not in gene set)")
    rownames(resTab) <- names(gsc)
    for (i in 1:length(gsc)) {
        gs <- gsc[[i]]
        nogs <- universe[!universe %in% gs]
        ctab <- rbind(c(sum(goi %in% gs), sum(goi %in% nogs)), 
            c(sum(bg %in% gs), sum(bg %in% nogs)))
        p[i] <- fisher.test(ctab, alternative = "greater")$p.value
        p.depletion[i] <- fisher.test(ctab, alternative = "less")$p.value
        rownames(ctab) <- c("Significant", "Non-significant")
        colnames(ctab) <- c("Genes in gene set", "Genes not in gene set")
        contTabList[[i]] <- ctab
        resTab[i, ] <- c(p[i], NA, p.depletion[i], sum(goi %in% gs), sum(bg %in% 
            gs), sum(goi %in% nogs), sum(bg %in% nogs))
    }
    padj.greater <- p.adjust(p, method = adjMethod)
    resTab[, 2] <- padj.greater
    res <- list()
    res$pvalues.greater <- p
    res$p.adj.greater <- padj.greater
    res$pvalues.depletion <- p.depletion
    res$resTab <- resTab
    res$contingencyTable <- contTabList
    res$gsc <- gsc
    return(res)
}


## Set up parallel backend.
library("doParallel")
# stopCluster(clusterLocal)
# ##when running at Emory (requires RSA public key-ssh, & manual run of script to start server backend):  
# parallelThreads set in configuration section.
# #clusterLocal <- makeCluster(c(rep("haplotein.biochem.emory.edu",parallelThreads)), type = "SOCK", port=10191, user="edammer", rscript="/usr/bin/Rscript",rscript_args="OUT=/dev/null SNOWLIB=/usr/lib64/R/library",manual=FALSE)
# ##when running elsewhere:
clusterLocal <- makeCluster(c(rep("localhost",parallelThreads)),type="SOCK")
registerDoParallel(clusterLocal)


## Load GMT file
GSCfromGMT<-loadGSC(file=GMTdatabaseFile)


## Output piano package GSA FET output tables as list assembly
GSA.FET.outlist<-list()

#  colnames(modulesData)[3:(ncol(modulesData)-1)]
if (ANOVAgroups) uniquemodcolors=names(DEXlistsForGO)  #otherwise, already set above.

# parallelized to speed up.
#for (this.geneList in uniquemodcolors) {
GSA.FET.outlist <- foreach(this.geneList=uniquemodcolors) %dopar% {
#  this.geneList=uniquemodcolors[i]
  zeroToKeep.idx= if (WGCNAinput) { which( background[,"GeneSymbol"] %in% unique(modulesData[which(modulesData$net.colors==this.geneList),"Unique.ID"]) ) } else {
                                    if(ANOVAgroups) { which( background[,"GeneSymbol"] %in% DEXlistsForGO[[this.geneList]] ) } else {
                                            which( background[,"GeneSymbol"] %in% unique(modulesData[[this.geneList]]) ) }}  #Handles file-based input modulesData
  zeroToKeep=rep(1,nrow(background))
  zeroToKeep[zeroToKeep.idx]<- 0
  cat( paste0("\n",this.geneList, "... ") )  #" (n=",length(zeroToKeep.idx)," gene symbols) now processing: ") )

  thislist <- runGSAhyper.twoSided(genes=background[which(zeroToKeep==0),"GeneSymbol"],universe=background[,"GeneSymbol"],gsc=GSCfromGMT,gsSizeLim=c(5,Inf),adjMethod="BH")
  #above line runs in time, ~30 sec/list, or 50 sec/list for .twoSided
  #list [[this.color]][["pvalues.greater"]] is enrichment p value vector
  #list [[this.color]][["padj.greater"]] is FDR vector
  #list [[this.color]][["resTab"]] is same-ordered (rows) matrix of: p-value, FDR, ... with rownames equal to the ontology name%ontology type%OntologyID
  #list [[this.color]][["pvalues.less"]] is depletion p value vector (relevant for signed Z score calculation), only in customized function

  return(list(thislist[["resTab"]], thislist[["gsc"]]))
}

  # re-combine list elements from two outputs, over all uniquemodcolors
  GSA.FET.resTab.list           = do.call(list,lapply(GSA.FET.outlist,function(x){x[[1]]}))
  GSA.FET.genesByOntology.list  = do.call(list,lapply(GSA.FET.outlist,function(x){x[[2]]}))

  names(GSA.FET.resTab.list) <- names(GSA.FET.genesByOntology.list) <- uniquemodcolors


#Add signed Zscore, pull out ontologyType, ontology (description, title case)
GSA.FET.outSimple <- lapply(GSA.FET.resTab.list, function(x) { 
  ontology=stringr::str_to_title(gsub("\\%WIKIPATHWAYS_\\d*","", gsub("\\%WP_\\d*","", gsub("\\&(.*);","\\1",gsub("<\\sI>","",gsub("<I>","", gsub("(.*)\\%.*\\%.*","\\1",rownames(x))))))))
  ontologyType=gsub("^WP\\d*","WikiPathways", gsub(".*\\%(.*)\\%.*","\\1",rownames(x)))
  ZscoreSign=rep(1,nrow(x))
  ZscoreSign[ x[,"Pvalue.Depletion"] < x[,"Pvalue.Enrichment"] ] <- -1
  Zscore=apply(x, 1, function(p) qnorm(min(p["Pvalue.Enrichment"], p["Pvalue.Depletion"])/2, lower.tail=FALSE))
  out=as.data.frame(x)
  out$Zscore=Zscore*ZscoreSign
  out$ontologyType=ontologyType
  out$ontology=ontology
  out })

#head(GSA.FET.outSimple[[1]])


## Available Ontology Types
 data.frame(x=table( gsub("^WP\\d*","WikiPathways",gsub(".*\\%(.*)\\%.*","\\1",rownames(GSA.FET.resTab.list[[1]])))))
#                                                 x.Var1 x.Freq
#1                                                BIOCYC     99
#2                                                  GOBP   6432
#3                                                  GOCC    980
#4                                                  GOMF   2015
#5                                              HUMANCYC      2
#6                                                   IOB     33
#7                                             MSIGDB_C2    395
#8                                       PANTHER PATHWAY     92
#9  PATHWAY INTERACTION DATABASE NCI-NATURE CURATED DATA    198
#10                                             PATHWHIZ    320
#11                                             REACTOME    744
#12                      REACTOME DATABASE ID RELEASE 80    795
#13                                                SMPDB    371
#14                                         WikiPathways    529




##3. Output Report of Z-Score Barplots, processing all GSA FET output tables
############################# ----------------------Plotting for modules ------------------------#######################
######## this script plots the top 3 ontologies for biological process, mol function and cell component for each module


##color scheme for ontology type key/legend (can be changed in user parameters, editing the "color" vector)
ontologyTypes=c("Biological Process","Molecular Function","Cellular Component","Reactome","WikiPathways","MSIG.C2")

if(ANOVAgroups) {
  xlabels <- names(DEXlistsForGO)
  xlabels.frame <- data.frame(Colors=rep(NA,length(xlabels)),Labels=xlabels)
  uniquemodcolors <- names(DEXlistsForGO) #not set above
} else {
  xlabels <- uniquemodcolors #labels2colors(c(1:nModules))
  xlabels1 <- paste("M",seq(1:nModules),sep="")
  xlabels.frame <- as.data.frame(data.frame(Colors=xlabels,Labels=paste0(xlabels1," ",xlabels)))
}

if (removeRedundantGOterms) {
	library(ontologyIndex)
	ontology.index<-get_OBO(file=GO.OBOfile,extract_tags="everything")
	#below function minimal_set uses ancestors list to dereplicate; what about using synonyms list from our go.obo as ancestors?
}

setwd(paste0(filePath,outFilename,"/"))
redundancyRemovalTag=if(removeRedundantGOterms) { "-redundancyRemoved" } else { "" }
filenameFinal=paste0(outFilename,redundancyRemovalTag)

pdf(paste0("GSA-GO-FET_",filenameFinal,".pdf"),height=pageDimensions[2],width=pageDimensions[1])
op <- par(mfrow=panelDimensions,oma=c(0,0,3,0))
frame()
legend(x="topleft",legend = ontologyTypes, fill=color, title=" ",cex=2,horiz=F,xpd=T)
legend(x="topleft",legend = c(" "," "," "), title="Ontology Types",cex=2.5,horiz=F,xpd=T, bty='n', title.adj=1.4)
#frame()

minZ<-qnorm(0.05/2, lower.tail=FALSE)
summary <- list()
for(i in c(1:(length(uniquemodcolors)))){
	thismod=uniquemodcolors[i]	
	tmp=GSA.FET.outSimple[[thismod]]

	if (length(tmp[,2]) == 0) next
	tmp = tmp[,c(10,9,8,1)] ## Select GO-terms,GO-Type, Z-score,pValues (and previously, gene Lists)
	tmp1 = tmp[order(tmp$Zscore,decreasing=T),]
	tmp2 = tmp1[order(tmp1$ontologyType,decreasing=T),] #was tmp2

	if(removeRedundantGOterms) {
		#dim(tmp2[which((tmp2$ontologyType == "GOCC" | tmp2$ontologyType == "GOBP" | tmp2$ontologyType == "GOMF") & tmp2$Zscore>=minZ),])
		##[1] 124   4
		go.full.tmp2<-tmp2[which((tmp2$ontologyType == "GOCC" | tmp2$ontologyType == "GOBP" | tmp2$ontologyType == "GOMF") & tmp2$Zscore>=minZ),]
		go.full.tmp2$Term=gsub("^.*\\%GO..\\%(GO:\\d*)$","\\1",rownames(go.full.tmp2))

		go.minimal.terms<-minimal_set(ontology.index, terms=go.full.tmp2$Term)
		#go.pruned.terms<-prune_descendants(ontology.index, roots=?, terms=)
		go.minimal.tmp2<-go.full.tmp2[which(go.full.tmp2$Term %in% go.minimal.terms),]

		tmp2$Zscore[which(tmp2$ontologyType=="GOCC" | tmp2$ontologyType=="GOBP" | tmp2$ontologyType=="GOMF")] <- 0  #keeps in all GO ontologies, but Z scores zeroed
		tmp2$Zscore[match(rownames(go.minimal.tmp2),rownames(tmp2))] <-go.minimal.tmp2$Zscore                       #put back Z scores of non-redundant ontologies, to be possibly kept below
#		tmp2<-rbind(tmp2,go.minimal.tmp2[,-ncol(go.minimal.tmp2)])
		tmp2 = tmp2[order(tmp2$Zscore,decreasing=T),]
		tmp2 = tmp2[order(tmp2$ontologyType,decreasing=T),]
	}

	tmp3 = tmp2[tmp2$ontologyType == "GOBP",][c(1:maxBarsPerOntology),]
	tmp3 = rbind(tmp3,tmp2[tmp2$ontologyType == "GOMF",][c(1:maxBarsPerOntology),] )
	tmp3 = rbind(tmp3,tmp2[tmp2$ontologyType == "GOCC",][c(1:maxBarsPerOntology),] )
	tmp3 = rbind(tmp3,tmp2[tmp2$ontologyType == "REACTOME",][c(1:maxBarsPerOntology),] )
	tmp3 = rbind(tmp3,tmp2[tmp2$ontologyType == "WikiPathways",][c(1:maxBarsPerOntology),] )
	tmp3 = rbind(tmp3,tmp2[tmp2$ontologyType == "MSIGDB_C2",][c(1:maxBarsPerOntology),] )

	tmp3 <- na.omit(tmp3)
	tmp3 <- tmp3[which(tmp3$Zscore>=minZ),]
#	tmp3 <- tmp3[order(tmp3$Zscore,decreasing=T),] #added this row, if you want to mix ontology types and sort by Z.Score only
	tmp3 <- tmp3[rev(rownames(tmp3)),]

	summary[[i]] <- tmp3
	
	### To color bars by mol function, cell component or biological process
	for (j in 1:nrow(tmp3)){
		if (tmp3$ontologyType[j] == "GOMF"){
			tmp3$color[j] <- color[2]
		} else if (tmp3$ontologyType[j] == "GOCC"){
			tmp3$color[j] <- color[3]
		} else if (tmp3$ontologyType[j] == "GOBP"){
			tmp3$color[j] <- color[1]
		} else if (tmp3$ontologyType[j] == "REACTOME"){
			tmp3$color[j] <- color[4]
		} else if (tmp3$ontologyType[j] == "WikiPathways"){
			tmp3$color[j] <- color[5]
		} else if (tmp3$ontologyType[j] == "MSIGDB_C2"){
			tmp3$color[j] <- color[6]
		}

	# tmp3$color[j] <- uniquemodcolors[i] #module color for all bars, instead of different colors by ontology type
	} 

	if (tmp3$Zscore == F) next
	par(mar=c(4,15,4,3))
	xlim <- c(0,1.1*max(tmp3$Zscore))	
	moduleTitle <- xlabels.frame[i,"Labels"]
	xh <- barplot(tmp3$Zscore,horiz = TRUE,width =0.85,las=1,main=moduleTitle, xlim=xlim,col=tmp3$color,cex.axis=0.7,xlab="Z Score",cex.lab=0.9,cex.main=0.95,ylim=c(0,nrow(tmp3)+0.8))
	abline(v=minZ,col="red", cex.axis = 0.5)
	axis(2, at=xh, labels = tmp3$ontology, tick=FALSE, las =2, line =-0.5, cex.axis = 0.7) 
}

par(op) # Leaves the last plot
dev.off()



## Master Tables Generation for output to csv and for GO:CC Z score coclustering

GSA.FET.collapsed.outSimple.Zscore <- cbind(GSA.FET.outSimple[[1]][,c("ontology","ontologyType")], matrix(unlist(lapply(GSA.FET.outSimple, function(x) x[,"Zscore"] )),byrow=FALSE,ncol=length(uniquemodcolors)))
GSA.FET.collapsed.outSimple.Pvalue.Enrichment <- cbind(GSA.FET.outSimple[[1]][,c("ontology","ontologyType")], matrix(unlist(lapply(GSA.FET.outSimple, function(x) x[,"Pvalue.Enrichment"] )),byrow=FALSE,ncol=length(uniquemodcolors)))
GSA.FET.collapsed.outSimple.Enrichment.FDR.BH <- cbind(GSA.FET.outSimple[[1]][,c("ontology","ontologyType")], matrix(unlist(lapply(GSA.FET.outSimple, function(x) x[,"Adjusted.Enr.Pvalue"] )),byrow=FALSE,ncol=length(uniquemodcolors)))
colnames(GSA.FET.collapsed.outSimple.Zscore)[3:(length(uniquemodcolors)+2)] <- colnames(GSA.FET.collapsed.outSimple.Pvalue.Enrichment)[3:(length(uniquemodcolors)+2)] <- colnames(GSA.FET.collapsed.outSimple.Enrichment.FDR.BH)[3:(length(uniquemodcolors)+2)] <- names(GSA.FET.outSimple)

#write master tables
write.table(GSA.FET.collapsed.outSimple.Zscore,file=paste0(filePath,outFilename,"/GSA-GO-FET_",outFilename,"-Zscores.txt"),row.names=FALSE,col.names=TRUE,sep="\t", quote=FALSE)
write.table(GSA.FET.collapsed.outSimple.Pvalue.Enrichment,file=paste0(filePath,outFilename,"/GSA-GO-FET_",outFilename,"-Enr.Pvalues.txt"),row.names=FALSE,col.names=TRUE,sep="\t", quote=FALSE)
write.table(GSA.FET.collapsed.outSimple.Enrichment.FDR.BH,file=paste0(filePath,outFilename,"/GSA-GO-FET_",outFilename,"-Enr.FDR.BH.txt"),row.names=FALSE,col.names=TRUE,sep="\t", quote=FALSE)


if(cocluster) {
## For GO:CC Z score coclustering
GSA.FET.GOCC.Zscore <- GSA.FET.collapsed.outSimple.Zscore[which(GSA.FET.collapsed.outSimple.Zscore$ontologyType=="GOCC"),]
GSA.FET.GOCC.terms <- gsub("^.*\\%GO..\\%(GO:\\d*)$","\\1",rownames(GSA.FET.GOCC.Zscore))
minZ<-qnorm(1e-5/2, lower.tail=FALSE)
GSA.FET.GOCC.terms.minZreached <- apply(GSA.FET.GOCC.Zscore[,3:ncol(GSA.FET.GOCC.Zscore)],1,function(x) if (max(x)>=minZ) { TRUE } else { FALSE } )
GSA.FET.GOCC.minimal.terms <- minimal_set(ontology.index, terms=GSA.FET.GOCC.terms[which(GSA.FET.GOCC.terms.minZreached)])
GSA.FET.GOCC.Zscore.minimal.terms <- GSA.FET.GOCC.Zscore[which(GSA.FET.GOCC.terms %in% GSA.FET.GOCC.minimal.terms),]
dim(GSA.FET.GOCC.Zscore)
# 980 rows of GO:CC terms
length(which(GSA.FET.GOCC.terms.minZreached))
# 596 have at least 1 Z score above minZ set immediately above (with p val max 0.001 used to calc minZ immediately above...)
# 381 at p val max 0.00001
length(GSA.FET.GOCC.minimal.terms)
# 329 kept out of 980 (with p val max 0.001 used to calc minZ immediately above...)
# 191 kept out of 980 (with p val max 0.00001 used to calc minZ immediately above...)

matrixdata <- data <- t(as.matrix(GSA.FET.GOCC.Zscore.minimal.terms[,3:ncol(GSA.FET.GOCC.Zscore.minimal.terms)]))
data[matrixdata>4]<-4
data[matrixdata< -4] <- -4
 
#NMF - initial approach
library(WGCNA)
bw<-colorRampPalette(c("#0058CC", "white"))
wr<-colorRampPalette(c("white", "#CC3300"))
colvec<-c(bw(50),wr(50))

heatmapLegendColors<-list(Modules=sort(uniquemodcolors))

library(NMF) #for aheatmap
pdf(file=paste0("GO_cc_clustering_from_GSA_FET_Z-",filenameFinal,".pdf"),width=18.5,height=18)
  aheatmap(x=data, ## Numeric Matrix
         main="Co-clustering with manhattan distance function, ward metric",
#         annCol=metdat,  ## Color swatch and legend annotation of columns/samples and rows (not used)
         annRow=data.frame(Modules=as.numeric(factor(uniquemodcolors,levels=sort(uniquemodcolors)))),
         annColors=heatmapLegendColors,
         border=list(matrix = TRUE),
         scale="none",  ## row, column, or none
         distfun="manhattan",hclustfun="ward", ## Clustering options
         cexRow=0.8, ## Character sizes
         cexCol=0.8,
         col=colvec,   #c("white","black"), ## Color map scheme
         treeheight=80,
         Rowv=TRUE,Colv=TRUE) ## Cluster columns
dev.off()

} # end if (cocluster) 


setwd(rootdir)
#save.image(paste0("saved.image.GSA.GO.FET-",outFilename,"-completed.Rdata"))
