###################################################################################################################################
## One-Step GSA FET (R piano implementation) WITH USER PARAMETERS
##  - by Eric Dammer, Divya Nandakumar
## Nicholas Seyfried Lab Bioinformatics - for the lab -
##   01/07/2026 version 1.4         -- New: GO.bubblePlot() function for post-processing to GSEA-style bubble plots
##                                     from the Zscore GOparallel output; and added a GOparallel(..., bubble=TRUE) argument
##   08/20/2025 version 1.3         -- malloc warning/error for parallelization on some Macs fixed; added Cohen's kappa graph-based
##                                     pruning within ontology types with option removeRedundantGOterms="kappa"
##   03/07/2023 version 1.2         -- Outputs include semicolon-separated Genes.Hit goi[goi %in% gs] for each input list/gene set
##   10/18/2022 version 1.1         -- (inputFile variable previously fileName and other syntax changes; now runs on R v4.2.0)
##                                  -- cocluster now runs even if removeRedundantGOterms not set or FALSE.
###################################################################################################################################

GOparallel <- function(dummyVar="",env=.GlobalEnv) {
	minHitsPerOntology=5  # ontologies hit by fewer than this number of genes in an input list will not be plotted in Z score barplots.

	if(!exists("filePath")) { cat(paste0("- filePath not set. Using current working directory: ",getwd(),"\n")); filePath=getwd(); }
	## Clean out spaces and escaped backslashes from folder paths (folder names with spaces should not be used on non-windows systems with this script)
	#filePath=paste0(paste( sapply(do.call(c,strsplit(filePath,"[/\\]")),function(x) { if (grepl(" ",x)) { gsub(x,paste0(substr(gsub(" ","",x),1,6),"~1"),x) } else { x } } ),collapse="/"),"/")
	filePath=paste0(paste( sapply(do.call(as.vector,strsplit(filePath,"[/\\]")),function(x) { if (grepl(" ",x)) { gsub(" ","\ ",x) } else { x } } ),collapse="/"),"/")
	if(!dir.exists(filePath)) { cat(paste0("- filePath set to ",filePath," ...this path was not found. Using current working directory: ",getwd(),"\n")); filePath=getwd(); }
	
	if(!exists("GMTdatabaseFile")) { cat(paste0("- GMTdatabaseFile variable not specified. Current BaderLab .GMT database file will be downloaded to ",filePath,"\n")); GMTdatabaseFile=paste0(filePath,"nonexistent.file"); }
	GMTdatabaseFile=paste0(paste( sapply(do.call(as.vector,strsplit(GMTdatabaseFile,"[/\\]")),function(x) { if (grepl(" ",x)) { gsub(" ","\ ",x) } else { x } } ),collapse="/"),"")
	if(!exists("GO.OBOfile")) { cat(paste0("- go.obo file not specified; if needed, current working directory will be checked and if not present, will be downloaded...\n")); GO.OBOfile=paste0(getwd(),"/go.obo"); }
	GO.OBOfile=paste0(paste( sapply(do.call(as.vector,strsplit(GO.OBOfile,"[/\\]")),function(x) { if (grepl(" ",x)) { gsub(" ","\ ",x) } else { x } } ),collapse="/"),"")
	
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
	

	if (!exists("filePath")) { cat(paste0(" - filePath variable not specified. Input/Output will take place in the current working directory: ",getwd(),"/ ...\n")); filePath<-paste0(getwd(),"/"); }
	if (!exists("outFilename")) { cat(paste0("- outFilename variable not specified. Output files will be saved to: ",getwd(),"/GOparallel/ ...\n")); outFilename="GOparallel"; }
	if (!dir.exists(file.path(filePath, outFilename))) dir.create(file.path(filePath, outFilename))
	
	if(!file.exists(GMTdatabaseFile)) {
		if (interactive()) {
			suppressPackageStartupMessages(require(rvest,quietly=TRUE))
			species.links <- html_attr(html_nodes(read_html("http://download.baderlab.org/EM_Genesets/current_release/"), xpath="//a"), "href")
			species.links <- species.links[grepl("^[A-z].*\\/$",species.links)]
			cat("- GMT File not found: ", GMTdatabaseFile,"\n\n")
			print(data.frame(Species=species.links))
			input.idx <- readline(paste0("[INTERACTIVE]\nChoose one of the above species from http://download.baderlab.org/EM_Genesets/current_release/ [1-",length(species.links),"]: "))
			input.idx <- as.integer(input.idx)
		
			find.symbol.in.links <- html_attr(html_nodes(read_html(paste0("http://download.baderlab.org/EM_Genesets/current_release/",species.links[input.idx])), xpath="//a"), "href")
			find.symbol.in.links <- find.symbol.in.links[grepl("[Ss][Yy][Mm][Bb][Oo][Ll]\\/",find.symbol.in.links)]
			gmt.links <- html_attr(html_nodes(read_html(paste0("http://download.baderlab.org/EM_Genesets/current_release/",species.links[as.integer(input.idx)],find.symbol.in.links)), xpath="//a"), "href")
			file.candidates1<-which(grepl("*\\_GO\\_AllPathways\\_.*\\.[Gg][Mm][Tt]",gmt.links))    # main regEx filter for file to download
			file.candidates2<-which(grepl("*\\_noPFOCR.*\\.[Gg][Mm][Tt]",gmt.links))                # files after March 2024 are a subset, excluding PMID-linked lists
			file.candidates3<-which(grepl("*\\_with\\_GO\\_iea\\_.*\\.[Gg][Mm][Tt]",gmt.links))     # take files with automated ontologies
			this.file.idx=if(length(file.candidates2)>0) { intersect(file.candidates2, intersect(file.candidates1,file.candidates3)) } else { intersect(file.candidates1,file.candidates3) }
			if(length(this.file.idx)<1) stop(paste0("Web scraping of the Bader Lab Website could not find an expected GMT filename pattern match.\nDownload and specify a GMTdatabaseFile prior to running this function."))

			full.dl.file=gmt.links[this.file.idx[1] ]
		
			GMTtargetPath=gsub("\\/\\/","/", gsub("(.*\\/).*$","\\1",GMTdatabaseFile) )
			gmt.url<-paste0("http://download.baderlab.org/EM_Genesets/current_release/",species.links[input.idx],find.symbol.in.links,full.dl.file)
			if(file.exists(file.path(GMTtargetPath,full.dl.file))) {
				cat(paste0("- Found that the full current GMT file online matches a file name you already have:\n  ",full.dl.file," [skipping download]\n"))
				GMTdatabaseFile=paste0(GMTtargetPath,full.dl.file)
			} else {
				cat("Found full current GMT file online:  ",gmt.url,"\n")
				cat("Download this file to folder:  ",GMTtargetPath,"\n")
				input.dlYN <- readline("[Y/n]?")
				if(input.dlYN == "Y" | input.dlYN == "y" | input.dlYN == "") {
					suppressPackageStartupMessages(require(curl,quietly=TRUE))
					if (!dir.exists(GMTtargetPath)) dir.create(GMTtargetPath)
					curr.dir<-getwd()
					setwd(GMTtargetPath)
					cat("Downloading .gmt file for ",species.links[input.idx],"...\n")
					curl_download(url=gmt.url, destfile=full.dl.file, quiet = TRUE, mode = "w")
					setwd(curr.dir)
					cat("Using new downloaded .gmt file: ", paste0(GMTtargetPath,full.dl.file),"\n")
					GMTdatabaseFile=paste0(GMTtargetPath,full.dl.file)
				}
			}
		} else { stop(paste0("This is not an interactive session and required GMT file not found.\n",GMTdatabaseFile," must be downloaded interactively or prior to running this function.")) }
	}

	if(!exists("removeRedundantGOterms")) { cat("- removeRedundantGOterms not specified TRUE/FALSE, or 'kappa'. Removing them as the default, using go.obo and ontologyIndex package (you can also try removeRedundantGOterms='kappa').\n"); removeRedundantGOterms=TRUE; } else {
        if (removeRedundantGOterms=="Kappa" | removeRedundantGOterms=="KAPPA") removeRedundantGOterms="kappa"
        if (removeRedundantGOterms=="kappa") { cat("- removeRedundantGOterms='kappa'; will use Cohen's kappa<0.30 in igraph representation of gene symbol similarity to prune within ontology types.\n") } else {
	   if(removeRedundantGOterms) {
		if (!file.exists(GO.OBOfile)) {
			suppressPackageStartupMessages(require(curl,quietly=TRUE))
			OBOtargetPath=gsub("(.*\\/).*$","\\1",GO.OBOfile)
			if (!dir.exists(OBOtargetPath)) dir.create(OBOtargetPath)
			curr.dir<-getwd()
			setwd(OBOtargetPath)
			cat(paste0("- Downloading go.obo file for main GO term redundancy cleanup...\n...to location:  ",OBOtargetPath,"go.obo\n"))
			curl_download(url="http://current.geneontology.org/ontology/go.obo", destfile="go.obo", quiet = TRUE, mode = "w")
			setwd(curr.dir)
			cat("GO.OBOfile set to downloaded file: ", paste0(OBOtargetPath,"go.obo"),"\n")
			GO.OBOfile=paste0(OBOtargetPath,"go.obo")
		}
	     }
	   }
	}
	
	## Check what type of input the user wants. modulesInMemory/ANOVAgroups/file(lists as columns, or kME module membership table)
	if (!exists("ANOVAgroups") & !exists("modulesInMemory")) {
	  if(!exists("inputFile")) {
	    cat("- No input specified as modulesInMemory=TRUE, ANOVAgroups=TRUE, and no inputFile variable for input either.\n  Trying modulesInMemory=TRUE ...\n")
	    modulesInMemory=TRUE
	    ANOVAgroups=FALSE
	  } else {
	    if(!file.exists(file.path(filePath,inputFile))) stop(paste0("\ninputFile input specified but not found where expected, in: ",paste0(filePath,inputFile),"\nDid you mean to set ANOVAgroups or modulesInMemory=TRUE?\n\n"))
	    ANOVAgroups=FALSE
	    modulesInMemory=FALSE
	  }
	} else {
	  if(exists("ANOVAgroups")) {
	    if(is.logical(ANOVAgroups)) {
	      if(ANOVAgroups) {
	        if(exists("ANOVAout")) { cat("- Found ANOVAgroups=TRUE. Proceeding to process ANOVAout table from memory, using any selections and thresholds set during volcano plotting...\n"); modulesInMemory=FALSE; } else { if (!length(dummyVar)==1) { cat("- ANOVAout not in memory, trying to use input provided to this function (could be ANOVAout or CORout).\n"); ANOVAout=as.data.frame(dummyVar); modulesInMemory=FALSE; } else { stop("Variable ANOVAout not found or no input was provided.\nPlease run parANOVA.dex() function first, and save output to ANOVAout variable or pass its output to this function.\n\n") } }
	      } else {
	        if(exists("modulesInMemory")) if(is.logical(modulesInMemory)) if(!modulesInMemory) {  # both flags are FALSE
	          if(!exists("outFilename")) { stop("modulesInMemory=FALSE, ANOVAgroups=FALSE, and no outFilename variable for input either.\nOne of these must be used.\n") }
	        } else {
	          cat("- modulesInMemory=TRUE. We will use network module colors vector and symbols found in rownames of cleanDat table for gene lists to check for ontology enrichment.\n")
	        }
	      }
	    } else {  #ANOVAgroups not logical TRUE/FALSE
	      if(exists("modulesInMemory")) if(is.logical(modulesInMemory)) if(!modulesInMemory) {  # modulesInMemory is FALSE
	        if(!exists("inputFile")) { cat("- modulesInMemory=FALSE, ANOVAgroups not TRUE/FALSE, and no inputFile variable for input either.\nOne of these must be used.  Trying ANOVAgroups=TRUE ...\n"); ANOVAgroups=TRUE; 
	                                 if(exists("ANOVAout")) { cat("- Found ANOVAout. Proceeding to process ANOVAout table from memory, using any selections and thresholds set during volcano plotting...\n") } else { stop("\nANOVAout table not found in memory.\n\n") }
	        } else {  # inputFile variable set. Does the file exist?
	          if(file.exists(file.path(filePath,inputFile))) { cat("- found inputFile as set for input: ",paste0(filePath,inputFile),"\n"); ANOVAgroups=FALSE; } else { stop("\ninputFile as specified not found: ",paste0(filePath,inputFile),"\n\n"); }
	        }
	      } else {  # modulesInMemory is TRUE
	        cat("- modulesInMemory=TRUE. We will use network module colors vector and symbols found in rownames of cleanDat table for gene lists to check for ontology enrichment.\n")
	        ANOVAgroups=FALSE
	      }
	    }
	  } else { # ANOVAgroups does not exist, but modulesInMemory exists. Is it TRUE?
	    if(exists("modulesInMemory")) {
	      if(is.logical(modulesInMemory)) {
	        if(modulesInMemory) {
	          cat("- modulesInMemory=TRUE.\n  We will use network module colors vector and symbols found in rownames of cleanDat table for gene lists to check for ontology enrichment.\n  Note rownames of cleanDat table must contain gene symbol, and NETcolors or net$colors vector of module color assignments must be available.\n")
	          ANOVAgroups=FALSE
	        } else {
	          if(!exists("inputFile")) { cat("- modulesInMemory=FALSE, ANOVAgroups not TRUE/FALSE, and no inputFile variable for input either.\nOne of these must be used.  Trying ANOVAgroups=TRUE ...\n"); ANOVAgroups=TRUE; 
	                                   if(exists("ANOVAout")) { cat("- Found ANOVAout. Proceeding to process ANOVAout table from memory, using any selections and thresholds set during volcano plotting...\n") } else { stop("\nANOVAout table not found in memory.\n\n") }
	          } else {  # outFilenName variable set. Does the file exist?
	            if(file.exists(file.path(filePath,inputFile))) { cat("- found inputFile as set for input: ",paste0(filePath,inputFile),"\n"); ANOVAgroups=FALSE; } else { stop("\ninputFile as specified not found: ",paste0(filePath,inputFile),"\n\n"); }
	          }
	        }
	      } else {  #modulesInMemory not logical TRUE/FALSE
	        if(!exists("inputFile")) {
	          cat("- ANOVAgroups not set, modulesInMemory not TRUE/FALSE, and no inputFile variable for input either.\nOne of these must be used.  Trying modulesInMemory=TRUE ...\n")
	          modulesInMemory=TRUE;
	          ANOVAgroups=FALSE;
	        } else {  # inputFile variable set. Does the file exist?
	          if(file.exists(file.path(filePath,inputFile))) { cat("- found inputFile as set for input: ",paste0(filePath,inputFile),"\n"); ANOVAgroups=FALSE; modulesInMemory=FALSE; } else { stop("\ninputFile as specified not found: ",paste0(filePath,inputFile),"\n\n"); }
	        }
	      }
	    } #else { #both ANOVAgroups and modulesInMemory do not exist; handled first above.
	  }
	}
	          
	## Further checks if modulesInMemory=TRUE for nrow(cleanDat)==NETcolors; and choice of vector for NETcolors
	if(modulesInMemory) {
	  if (!exists("cleanDat")) stop("cleanDat variable must exist, rownames expected to hold gene symbols in the form of 'Symbol' or 'Symbol|...'.\n\n")
	  if (!exists("NETcolors")) {
	    if(exists("net")) if ("colors" %in% names(net)) {
	      NETcolors=net$colors
	    } else {
	      if ("NETcolors" %in% colnames(ANOVAout)) {
	        NETcolors=ANOVAout$NETcolors
	      } else {
	        NETcolors=c()
	      }
	    }
	  }
	  if (!length(NETcolors)==nrow(cleanDat)) { stop("\n\nNetwork color assignment vector not supplied or not of length in rows of cleanDat.\nEvery gene product needs a module (color) assignment when modulesInMemory=TRUE.\n\n") }
	}

	if (!exists("parallelThreads")) { cat("- parallelThreads variable not set. Attempting to run with 8 threads.\n"); parallelThreads=8; }
	if (!exists("outputGOeliteInputs")) outputGOeliteInputs=FALSE

	##1a. Organize input gene lists, of the significant up and down (p<0.05) proteins in the current cleanDat
	### List building from ANOVA-defined categories ###
	if(!exists("corVolc")) corVolc=FALSE
	if (ANOVAgroups) {

	  if(corVolc & exists("CORout")) { cat("- corVolc=TRUE so getting gene lists from significant correlations in statistics table stored in variable CORout.\n"); ANOVAout<-CORout; }
	  if(corVolc & !exists("CORout")) if (!length(dummyVar)==1) { cat("- corVolc=TRUE, but CORout not in memory, using input provided to this function.\n"); ANOVAout=as.data.frame(dummyVar); } else { stop("Variable CORout not found or no input was provided.\nPlease run trait.corStat() function first, and save output to CORout variable or pass its output to this function.\n\n") }

	  if (!exists("ANOVAout")) if (!length(dummyVar)==1) { cat("- ANOVAout not in memory, trying to use input provided to this function.\n"); ANOVAout=as.data.frame(dummyVar); } else { stop("Variable ANOVAout not found or no input was provided.\nPlease run parANOVA.dex() or trait.corStat() function first, and save output to ANOVAout variable or pass its output to this function.\n\n") }
	  if (!ncol(ANOVAout)>3) stop("\n\nInput or ANOVAout variable contents are not a data (frame) with at least 4 columns. It is not valid output from the parANOVA.dex() or trait.corStat() function.\n  Please run one of these functions first.\n\n")
	  
	  numberOfNonComparisonColumns=length(colnames(ANOVAout)) - if(!corVolc) { length(which(grepl("^diff ",colnames(ANOVAout))))*2 } else { length(which(grepl("^p ",colnames(ANOVAout))))*2 }
	  numComp <- (length(colnames(ANOVAout)) - numberOfNonComparisonColumns) / 2 # of columns separating comparisons from matched column of log2(diffs), i.e. # of comparisons

	  if (!exists("testIndexMasterList")) {
	    if(exists("selectComps")) { cat("- Volcano (plotVolc) selection of pairwise comparisons in ANOVAout may apply from the variable selectComps.\n"); testIndexMasterList=selectComps; } else { cat("- No comparison p value columns previously selected by running plotVolc(). Using ALL comparisons.\n"); testIndexMasterList="ALL"; }
	  }
	  if (max(testIndexMasterList)>numComp+2 | min(testIndexMasterList)<3) {
	    cat(" - Selected comparison p value columns numbers may not reference valid integer p value column indexes of ANOVAout (or CORout).\n   Output will be for all valid comparisons or correlation(s).\n")
	    testIndexMasterList="ALL"
	  }
	  if (testIndexMasterList[1]=="ALL" | testIndexMasterList[1]=="all" | testIndexMasterList[1]=="All") testIndexMasterList=c(3:(numComp+2))
	  testIndexMasterList=as.integer(testIndexMasterList)

	  if(corVolc & exists("flip")) { cat("- Significant groups were defined by trait corrrelation. Variable flip will be ignored so positive correlations remain positive.\n"); flip=c(); }
	  if(!exists("flip")) { cat("- No comparisons selected for flipping numerator and denominator. Variable flip=c().\n"); flip=c(); }

	  if(!exists("dexComps") | !exists("comparisonIDs")) {
	    cat("- Extracting names of groups being compared in ANOVAout table from p value column headers...\n")
	    dexComps <- list()
	    iter <- length(testIndexMasterList) + 1
	    comparisonIDs <- data.frame(dfVariable = rep(NA, length(testIndexMasterList)), Comparison = rep(NA, length(testIndexMasterList)))
	    for (i in testIndexMasterList) {
	      iter <- iter - 1
	      if(!corVolc) {
	        comparisonIDs[iter, ] <- as.vector(c(paste0("dexTargets.", gsub("-", ".", colnames(ANOVAout)[i])), paste0(as.character(gsub("-", " vs ", colnames(ANOVAout)[i])))))
	      } else { 
	        comparisonIDs[iter, ] <- as.vector(c(paste0("dexTargets.", gsub("^p ", "", colnames(ANOVAout)[i])), paste0(as.character(gsub("^(.*)[' '](.*)\\.(.*)$", "\\1 (\\2 in \\3 samples)", colnames(ANOVAout)[i+numComp])))))
	      }
	      dexComps[[comparisonIDs[iter, 1] ]] <- ANOVAout
	      if (!is.na(match(i, flip))) {
	        dexComps[[comparisonIDs[iter, 1] ]][, i + numComp] <- -1 * as.numeric(dexComps[[comparisonIDs[iter, 1] ]][, i + numComp])
	        comparisonIDs[iter, 2] <- gsub("(*.*) vs (*.*)", "\\2 vs \\1", comparisonIDs[iter, 2]) # flip label "vs" in comParisonIDs$Comparison[iter]
	      }
	    }
	  }
	  # comparisonIDs  # list element names and Logical comparisons for those retrievable Dex measurements in the list elements
	  # ls(dexComps)   # list elements are dataframes with the DEX entries for that comparison


	  ANOVAout$Symbol <- suppressWarnings(do.call("rbind", strsplit(as.character(rownames(ANOVAout)), "[|]"))[, 1])
	  if(length(which(grepl(";",ANOVAout$Symbol)))>0) {
	    cat("- *Found some gene symbols have semicolons! Splitting these and keeping only symbol *before* semicolon.\n")
	    ANOVAout$Symbol<-suppressWarnings(do.call("rbind",strsplit(as.character(ANOVAout$Symbol), "[;]"))[,1])
	  }

	  cat("\n")

	  if (!exists("sigThresh")) if(exists("sigVolcCutoff")) { sigThresh=sigVolcCutoff } else { sigThresh=0.05 }
	  print(paste0("...Applying a minimum p value cutoff of ",sigThresh," for ",if(corVolc) { "correlation statistics" } else { "ANOVA" }, " lists..."))

	  if(corVolc & exists("FCmin")) { cat("- Using correlation significance for defining groups. Variable FCmin will be ignored.\n"); FCmin=0; }
	  if (!exists("FCmin")) FCmin=0
	  cutoff=log2(1+FCmin)
	  # shows what your cutoff for log2(FC) calculates as
	  if(!corVolc) { print(paste0("...Applying a ", FCmin*100,"% minimum fold change threshold at + and - x=", signif(cutoff,2),"...")) }


	  DEXlistsForGO<-list()
	  iter=0
	  for (i in testIndexMasterList) {
	    iter=iter+1;
	    j=paste0(gsub(" ",".",comparisonIDs$Comparison[iter]),".down")
	    k=paste0(gsub(" ",".",comparisonIDs$Comparison[iter]),".up")
	    if (length(intersect(i,flip))==1) {
	      #flipped sign (all diffs >cutoff for down)
	      DEXlistsForGO[[j]]<-ANOVAout$Symbol[which(ANOVAout[,i]<sigThresh & ANOVAout[,i+numComp]> cutoff)]
	      DEXlistsForGO[[k]]<-ANOVAout$Symbol[which(ANOVAout[,i]<sigThresh & ANOVAout[,i+numComp]< -cutoff)]
	    } else {
	      #do not flip sign (all diffs < -cutoff for down)
	      DEXlistsForGO[[j]]<-ANOVAout$Symbol[which(ANOVAout[,i]<sigThresh & ANOVAout[,i+numComp]< -cutoff )]
	      DEXlistsForGO[[k]]<-ANOVAout$Symbol[which(ANOVAout[,i]<sigThresh & ANOVAout[,i+numComp]> cutoff)]
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
	    modulesData <- as.data.frame(cbind(rownames(cleanDat),NETcolors))
	    WGCNAinput=TRUE
	  } else {
	    modulesData <- read.csv(paste(filePath,inputFile,sep=""),header=TRUE, sep=",");
	    # check if this is a WGCNA modules/kME table containing net.colors, or NETcolors column; (otherwise it is assumed to be simple list input)
	    if("net.colors" %in% colnames(modulesData)) {
	      WGCNAinput=TRUE
	      # .csv column before colors is assumed to contain Unique IDs (Symbol|...) unless colors are in first column; then it's assumed to be in column following colors
	      NETcolors.idx=which(colnames(modulesData) %in% "net.colors")[1]
	      if(NETcolors.idx==1) { modulesData <- as.data.frame(modulesData[,c(NETcolors.idx+1,NETcolors.idx)]); } else { modulesData <- as.data.frame(modulesData[,c(NETcolors.idx-1,NETcolors.idx)]); }
	    } else {
	      if("NETcolors" %in% colnames(modulesData)) {
	        WGCNAinput=TRUE
	        # .csv column before colors is assumed to contain Unique IDs (Symbol|...) unless colors are in first column; then it's assumed to be in column following colors
	        NETcolors.idx=which(colnames(modulesData) %in% "NETcolors")[1]
	        if(NETcolors.idx==1) { modulesData <- as.data.frame(modulesData[,c(NETcolors.idx+1,NETcolors.idx)]); } else { modulesData <- as.data.frame(modulesData[,c(NETcolors.idx-1,NETcolors.idx)]); }
	      } else {
	        WGCNAinput=FALSE
	      }
	    }
	  }
	
	  if (WGCNAinput) {
	    suppressPackageStartupMessages(require(WGCNA,quietly=TRUE)) #for labels2colors
	  
	    # Include column with Symbol (if it is gene symbol, if not use appropriate code as given in GO-Elite manual)
	    modulesData$SystemCode <- rep("Sy",nrow(modulesData)) 
	  
	    # Assign Names of First columns, in case they are non standard
	    colnames(modulesData)[1]<-"Unique.ID" #This should have Symbol|UniprotID
	    colnames(modulesData)[2]<-"net.colors" #This should have colors
	    
	    #Split out symbols from UniprotIDs, keep symbols in column 1
	    rownames(modulesData)<-modulesData$Unique.ID
	    modulesData$Unique.ID<-suppressWarnings(do.call("rbind",strsplit(as.character(modulesData$Unique.ID), "[|]"))[,1])
	    if(length(which(grepl(";",modulesData$Unique.ID)))>0) {
	      cat("- *Found some gene symbols have semicolons! Splitting these and keeping only symbol *before* semicolon.\n")
	      modulesData$Unique.ID<-suppressWarnings(do.call("rbind",strsplit(as.character(modulesData$Unique.ID), "[;]"))[,1])
	    }
	    
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
	      colnames(moduleInfo) <- c("GeneSymbol","SystemCode")  #,"kME")
	      if (moduleName == "blue" | moduleName == "brown" | moduleName == "green" | moduleName == "cyan") { if(outputGOeliteInputs) write.table(moduleInfo,file=paste(filePath,outFilename,"/",moduleName,"_2_Module.txt",sep=""),row.names=FALSE,col.names=TRUE,sep="\t", quote=FALSE)
	      } else {
	        if(outputGOeliteInputs) write.table(unique(moduleInfo),file=paste(filePath,outFilename,"/",moduleName,"_Module.txt",sep=""),row.names=FALSE,col.names=TRUE,sep="\t", quote=FALSE)
	      }
	    }
	  } else { #input is not WGCNA kME table format 
	
	##1c. List building from the columns of a user-specified .csv input file, which must be in a column-wise list format, and including longest such list as background.
	    # We process the input file as simple lists by column in the CSV (largest list used as background)
	    
	    #reread the file to a list of gene symbol (or UniqueID) lists
	    modulesData <- as.list(read.csv(paste(filePath,inputFile, sep=""),sep=",", stringsAsFactors=FALSE,header=T)) 
	    
	    nModules <- length(names(modulesData))
	    semicolonsFound=FALSE
	    for (a in 1:nModules) {
	      modulesData[[a]] <- unique(modulesData[[a]][modulesData[[a]] != ""])
	      modulesData[[a]] <- modulesData[[a]][!is.na(modulesData[[a]])]
	      modulesData[[a]] <- suppressWarnings(do.call("rbind",strsplit(as.character(modulesData[[a]]), "[|]"))[,1])
	      if(length(which(grepl(";",modulesData[[a]])))>0) {
	        modulesData[[a]]<-suppressWarnings(do.call("rbind",strsplit(as.character(modulesData[[a]]), "[;]"))[,1])
	      }
	    }
	    if(semicolonsFound) cat("- *Found some gene symbols have semicolons! Splitting these and keeping only symbol *before* semicolon.\n")

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
	
	suppressPackageStartupMessages(require(piano,quietly=TRUE))
	
	## Adapted version of piano::runGSAhyper() function with depletion p value also calculated (for signed Z score if we will use it to cocluster by module, e.g.)
	runGSAhyper.twoSided <- function(genes, pvalues, pcutoff, universe, gsc, gsSizeLim = c(1,Inf), adjMethod = "fdr") {
	    if (length(gsSizeLim) != 2) 
	        stop("argument gsSizeLim should be a vector of length 2")
	    if (missing(genes)) {
	        stop("argument genes is required")
	    } else {
	        genes <- as.vector(as.matrix(genes))
	        if (!is(genes, "character")) 
	            stop("argument genes should be a character vector")
	        if (length(unique(genes)) != length(genes)) 
	            stop("argument genes should contain no duplicated entries")
	    }
	    if (missing(pvalues)) {
	        pvalues <- rep(0, length(genes))
	    } else {
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
	        } else {
	            pcutoff <- 0.05
	        }
	    } else {
	        if (length(pcutoff) != 1 & !is(pcutoff, "numeric")) 
	            stop("argument pcutoff should be a numeric of length 1")
	        if (max(pcutoff) > 1 | min(pcutoff) < 0) 
	            stop("argument pcutoff needs to lie between 0 and 1")
	    }
	    if (missing(gsc)) {
	        stop("argument gsc needs to be given")
#	    } else {
#	        if (!is(gsc, "GSC")) 
#	            stop("argument gsc should be of class GSC, as returned by the loadGSC function")  # disabled since the list we create is not of GSC class
	    }
	    if (missing(universe)) {
	        if (!all(pvalues == 0)) {
	            universe <- genes
	            message("Using all genes in argument genes as universe.")
	        } else {
	            universe <- unique(unlist(gsc$gsc))
	            message("Using all genes present in argument gsc as universe.")
	        }
	    } else {
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
	    if (length(goi) < 1) {
	        cat("\nrunGSEAhyper: no genes selected due to too strict pcutoff. (no genes of interest made an input list)\n")
	        res<-list()
	        res$resTab <- NA
	        res$gsc <- NA
	        return(res)
	    }
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
	    resTab <- matrix(nrow = length(gsc), ncol = 8)  #added 8th column to hold "Genes.Hit"
	    colnames(resTab) <- c("Pvalue.Enrichment", "Adjusted.Enr.Pvalue", "Pvalue.Depletion",
	        "Significant (in gene set)", "Non-significant (in gene set)", 
	        "Significant (not in gene set)", "Non-significant (not in gene set)", "Genes.Hit")
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
	            gs), sum(goi %in% nogs), sum(bg %in% nogs), paste0(goi[goi %in% gs],collapse=";"))  #*** added semicolon separated Genes.Hit to 8th column
	    }
	    padj.greater <- p.adjust(p, method = adjMethod)
	    resTab[, 2] <- padj.greater
	    res <- list()
	    res$pvalues.greater <- p
	    res$p.adj.greater <- padj.greater
	    res$pvalues.depletion <- p.depletion
	    res$resTab <- resTab   #*** includes Genes.Hit in 8th column.
	    res$contingencyTable <- contTabList
	    res$gsc <- gsc
	    return(res)
	}
	
	
	## Set up parallel backend.
	suppressPackageStartupMessages(require("doParallel",quietly=TRUE))
	clusterLocal <- makeCluster(c(rep("localhost",parallelThreads)) )  #,type="SOCK") <- may cause malloc error on macOS/non-windows."PSOCK" or leave off...
	registerDoParallel(clusterLocal)
	
	## Load GMT file; Clean UTF-8 characters (since Dec 2023); Write clean.GMT back out
	#GMT.df <- read.delim(GMTdatabaseFile, encoding = "utf-8",quote="", sep="\t",header=FALSE) 
	GMT.df <- readLines(con <- file(GMTdatabaseFile, encoding = "utf-8"))
        close(con)
        GMT.df <- unlist(sapply(GMT.df, function(x) iconv(gsub("^(PMC\\d*__.+?)\\\t(.*)$","\\1%PMC%\\2", 
                                                          gsub("\\\"","", gsub("\\x83\\x80.","-",x) )),
                                                          "utf-8","ASCII", "")))
	names(GMT.df)<-NULL
        GMT.df <- lapply(GMT.df, function(x) stringr::str_split_fixed(x, pattern="\t", n=Inf))

        # Create list object that is identical to a GSC class object, just not of this class, since not loaded by the loadGSC() function in piano package.
        GSCfromGMT<-list()
        GSCfromGMT[["addInfo"]]<-do.call(rbind, lapply(GMT.df, function(x) if(grepl("^PMC.*\\%PMC\\%",x[1])) { c(x[1],gsub("^(PMC.*)\\%PMC\\%.*$","\\1",x[1])) } else { x[c(1:2)] } ))
        GSCfromGMT[["gsc"]]<-lapply(GMT.df, function(x) if(grepl("^PMC.*\\%PMC\\%",x[1])) { x[c(2:length(x))][!x[c(2:length(x))]==""] } else { x[c(3:length(x))][!x[c(3:length(x))]==""] })
        names(GSCfromGMT$gsc)<-GSCfromGMT$addInfo[,1]
        
        # Time and memory overhead are too great to write and read back in a clean.GMT.  We process the provided .GMT with UTF-8 and inconsistencies every time this script is run.
        #write.table(GMT.df,file="clean.GMT",sep='\t',quote=FALSE, col.names=FALSE, row.names=FALSE)
        #GSCfromGMT<-loadGSC(file="clean.GMT")  # loadGSC(file=GMTdatabaseFile)
	
	## Be sure cluster nodes for parallel processing inherit needed variables from both .GlobalEnv and current function environment (error seen in R 4.2.1 in RStudio on Windows).
	if(!exists("DEXlistsForGO")) DEXlistsForGO<-list()
	parallel::clusterExport(cl=clusterLocal, list("ANOVAgroups","WGCNAinput","background","DEXlistsForGO","GSCfromGMT"), envir=environment())   ## avoid error during foreach below:  Error in { : task 1 failed - "object 'ANOVAgroups' not found"

	
	## Output piano package GSA FET output tables as list assembly
	GSA.FET.outlist<-list()
	
	#  colnames(modulesData)[3:(ncol(modulesData)-1)]
	if (ANOVAgroups) uniquemodcolors=names(DEXlistsForGO)  #otherwise, already set above.
	
	# parallelized to speed up.
	cat("\nRunning FET overlap statistics in parallel for ",length(uniquemodcolors)," symbol lists using up to ", parallelThreads," threads...\n\n")

	#for (this.geneList in uniquemodcolors) {
	GSA.FET.outlist <- foreach(this.geneList=uniquemodcolors) %dopar% {
	#  this.geneList=uniquemodcolors[i]
	  zeroToKeep.idx= if (WGCNAinput) { which( background[,"GeneSymbol"] %in% unique(modulesData[which(modulesData$net.colors==this.geneList),"Unique.ID"]) ) } else {
	                                    if(ANOVAgroups) { which( background[,"GeneSymbol"] %in% DEXlistsForGO[[this.geneList]] ) } else {
	                                            which( background[,"GeneSymbol"] %in% unique(modulesData[[this.geneList]]) ) }}  #Handles file-based input modulesData
	  zeroToKeep=rep(1,nrow(background))
	  zeroToKeep[zeroToKeep.idx]<- 0
	  cat( paste0("\n",this.geneList, "... ") )  #" (n=",length(zeroToKeep.idx)," gene symbols) now processing: ") )
	
	  thislist <- runGSAhyper.twoSided(genes=background[which(zeroToKeep==0),"GeneSymbol"],universe=background[,"GeneSymbol"],gsc=GSCfromGMT,gsSizeLim=c(minHitsPerOntology,Inf),adjMethod="BH",)
	  #above line runs in time, ~30 sec/list, or 50 sec/list for .twoSided
	  #list [[this.color]][["pvalues.greater"]] is enrichment p value vector
	  #list [[this.color]][["padj.greater"]] is FDR vector
	  #list [[this.color]][["resTab"]] is same-ordered (rows) matrix of: p-value, FDR, ... with rownames equal to the ontology name%ontology type%OntologyID
	  #list [[this.color]][["pvalues.less"]] is depletion p value vector (relevant for signed Z score calculation), only in customized function

	  return(list(thislist[["resTab"]], thislist[["gsc"]]))
	}
        stopCluster(clusterLocal)
	
	  # re-combine list elements from two outputs, over all uniquemodcolors
	  GSA.FET.resTab.list           = do.call(list,lapply(GSA.FET.outlist,function(x){x[[1]]}))
	  GSA.FET.genesByOntology.list  = do.call(list,lapply(GSA.FET.outlist,function(x){x[[2]]}))
	
	  names(GSA.FET.resTab.list) <- names(GSA.FET.genesByOntology.list) <- uniquemodcolors
	
	
	#Add signed Zscore, pull out ontologyType, ontology (description, title case)
	GSA.FET.outSimple <- lapply(GSA.FET.resTab.list, function(x) { 
	  if(is.na(x[1])) {    #*** occurs when "no genes selected due to too strict pcutoff" in GSEA-FET piano
	    NA
	  } else {
	    #PMC\\d*__F\\d* ontologyTypes for PMC gene sets added December 2023 -- collapse to ontologyType "PMC"
	    #rownames(x)=gsub("^(PMC\\d*__.*)\\t(.*)\\t(.*)",iconv("\\2%PMC%\\1\\t\\3", "latin1","ASCII", ""),rownames(x))
	    ontology=stringr::str_to_title(gsub("\\%WIKIPATHWAYS_\\d*","", gsub("\\%WP_\\d*","", gsub("\\&(.*);","\\1",gsub("<\\sI>","",gsub("<I>","", gsub("(.*)\\%.*\\%.*","\\1",rownames(x))))))))
	    ontologyType=gsub("^WP\\d*","WikiPathways", gsub(".*\\%(.*)\\%.*","\\1",rownames(x)))

	    #force all caps for ontologyType of GObp GOmf GOcc (changed in downloaded GMT files Sept 2022 and/or different in mouse GMT compared to human)
	    ontologyType=gsub("GObp","GOBP",ontologyType)
	    ontologyType=gsub("GOmf","GOMF",ontologyType)
	    ontologyType=gsub("GOcc","GOCC",ontologyType)

	    ZscoreSign=rep(1,nrow(x))
	    ZscoreSign[ as.numeric(x[,"Pvalue.Depletion"]) < as.numeric(x[,"Pvalue.Enrichment"]) ] <- -1
	    Zscore=apply(x, 1, function(p) qnorm(min(as.numeric(p["Pvalue.Enrichment"]), as.numeric(p["Pvalue.Depletion"]))/2, lower.tail=FALSE))
	    out=as.data.frame(x)
	    out$Zscore=Zscore*ZscoreSign
	    out$ontologyType=ontologyType
	    out$ontology=ontology
	    out
	  }
	})
	
	#head(GSA.FET.outSimple[[1]])
	
	
	## Available Ontology Types
#	 data.frame(x=table( gsub("^WP\\d*","WikiPathways",gsub(".*\\%(.*)\\%.*","\\1",rownames(GSA.FET.resTab.list[[1]])))))

##Circa 2022
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

##Dec2023
#	 data.frame(x=table( gsub("^PMC\\d*__.*","PMC", gsub("^WP\\d*","WikiPathways",gsub(".*\\%(.*)\\%.*","\\1",rownames(GSA.FET.resTab.list[[1]]))))))
#	 data.frame(x=table( GSA.FET.outSimple[[1]]$ontologyType ))

	#                                                 x.Var1 x.Freq
	#1                                                BIOCYC     94
	#2                                                  GOBP   6518
	#3                                                  GOCC    784
	#4                                                  GOMF   1704
	#5                                              HUMANCYC      3
	#6                                                   IOB     34
	#7                                             MSIGDB_C2    415
	#8                                        MSIGDBHALLMARK     50
	#9                                       PANTHER PATHWAY     89
	#10 PATHWAY INTERACTION DATABASE NCI-NATURE CURATED DATA    201
	#11                                             PATHWHIZ    260
	#12                                             REACTOME    843
	#13                      REACTOME DATABASE ID RELEASE 38    802
	#14                                                SMPDB    290
	#15                                         WikiPathways    632
	
##Jan2024
#	 data.frame(x=table( GSA.FET.outSimple[[1]]$ontologyType ))
# note: ontologies kept from input GMT is dependent on background overlap with 'universe' of genes in the GMT  -- this is for a smaller background of ~1100 symbols
	#                                                 x.Var1 x.Freq
	#1                                                BIOCYC      9
	#2                                                  GOBP   2426
	#3                                                  GOCC    335
	#4                                                  GOMF    573
	#5                                                   IOB      4
	#6                                             MSIGDB_C2     38
	#7                                        MSIGDBHALLMARK     35
	#8                                       PANTHER PATHWAY     13
	#9  PATHWAY INTERACTION DATABASE NCI-NATURE CURATED DATA     53
	#10                                             PATHWHIZ     31
	#11                                                  PMC    880
	#12                                             REACTOME    188
	#13                      REACTOME DATABASE ID RELEASE 65    191
	#14                                                SMPDB     37
	#15                                         WikiPathways    122

	
	##3. Output Report of Z-Score Barplots, processing all GSA FET output tables
	############################# ----------------------Plotting for modules ------------------------#######################
	######## this script plots the top 3 ontologies for biological process, mol function and cell component for each module
	
	
	##color scheme for ontology type key/legend (can be changed in user parameters, editing the "color" vector)
	ontologyTypes=c("Biological Process","Molecular Function","Cellular Component","Reactome","WikiPathways","MSIG.C2") # ,"PMC")
	
	if(ANOVAgroups) {
	  xlabels <- names(DEXlistsForGO)
	  xlabels.frame <- data.frame(Colors=rep(NA,length(xlabels)),Labels=xlabels)
	  uniquemodcolors <- names(DEXlistsForGO) #not set above
	} else {
	  xlabels <- uniquemodcolors #labels2colors(c(1:nModules))
	  xlabels1 <- paste("M",seq(1:nModules),sep="")
	  xlabels.frame <- as.data.frame(data.frame(Colors=xlabels,Labels=paste0(xlabels1," ",xlabels)))
	}
	
	if (!removeRedundantGOterms=="kappa") {
	  if(removeRedundantGOterms) {
		suppressPackageStartupMessages(require(ontologyIndex))
		ontology.index<-get_OBO(file=GO.OBOfile,extract_tags="everything")
		#below function minimal_set uses ancestors list to dereplicate; what about using synonyms list from our go.obo as ancestors?
	  }
	} else {  # removeRedundantGOterms=="kappa"

          # --- Cohen's kappa clustering on gene-set membership -------------------------
          kappa_cluster_ontologies <- function(gene_sets_present, kappa_cut = 0.30,
                                               method = c("hier", "graph")) {
            method <- match.arg(method, c("hier","graph"))
          
            # Normalize & validate
            terms <- names(gene_sets_present)
            if (is.null(terms)) stop("gene_sets_present must be a *named* list.")
            gene_sets_present <- lapply(gene_sets_present, function(x) unique(as.character(x)))
            n <- length(gene_sets_present)
          
            if (n == 0L) return(data.frame(term = character(0), cluster = integer(0)))
            if (n == 1L) return(data.frame(term = terms, cluster = 1L))
          
            # Universe & encodings
            U <- sort(unique(unlist(gene_sets_present)))
            N <- length(U)
            idx_map  <- setNames(seq_along(U), U)
            sets_idx <- lapply(gene_sets_present, function(gs) if (length(gs)) sort(idx_map[gs]) else integer(0))
            sizes    <- vapply(sets_idx, length, integer(1))
          
            # If nothing in universe, just give unique clusters
            if (N == 0L || all(sizes == 0L)) {
              return(data.frame(term = terms, cluster = seq_len(n)))
            }
          
            # Pairwise intersections (upper triangle)
            intersec_mat <- matrix(0L, n, n)
            for (i in seq_len(n - 1L)) {
              gi <- sets_idx[[i]]
              for (j in (i + 1L):n) {
                intersec_mat[i, j] <- length(intersect(gi, sets_idx[[j]]))
              }
            }
            intersec_mat <- intersec_mat + t(intersec_mat)
            diag(intersec_mat) <- sizes
          
            # Cohen’s kappa
            A  <- intersec_mat
            Bi <- matrix(sizes, n, n, byrow = FALSE) - A
            Cj <- matrix(sizes, n, n, byrow = TRUE)  - A
            D  <- N - (A + Bi + Cj)
          
            p0 <- (A + D) / max(N, 1L)
            pe <- (outer(sizes, sizes, "*") + outer(N - sizes, N - sizes, "*")) / max(N^2, 1L)
          
            denom <- pmax(1 - pe, .Machine$double.eps)
            K <- (p0 - pe) / denom
            diag(K) <- 1
            K[!is.finite(K)] <- 0
          
            if (method == "hier") {
              Dmat <- as.dist(pmax(0, 1 - K))
              hc <- hclust(Dmat, method = "average")
              cl <- cutree(hc, h = 1 - kappa_cut)
              data.frame(term = terms, cluster = as.integer(cl))
            } else {
              if (!requireNamespace("igraph", quietly = TRUE)) {
                stop("Please install 'igraph' for graph-based clustering, pruning via 'kappa'.")
              }
              adj <- (K >= kappa_cut)
              diag(adj) <- FALSE
              g <- igraph::graph_from_adjacency_matrix(adj, mode = "undirected", diag = FALSE)
              comp <- igraph::components(g)$membership
              data.frame(term = terms, cluster = as.integer(comp))
            }
          }
          
          # --- Prune a per-list table (tmp2) to 1 representative per cluster -----------
          # tmp2 must have columns: ontology, ontologyType, Zscore, Pvalue.Enrichment, Genes.Hit
          prune_ontologies_by_kappa <- function(tmp2, kappa_cut = 0.30, method = c("hier", "graph")) {
            method <- match.arg(method, c("hier","graph"))
            need <- c("ontology","ontologyType","Zscore","Pvalue.Enrichment","Genes.Hit")
            stopifnot(all(need %in% names(tmp2)))
          
            df <- tmp2
            df$ontology          <- as.character(df$ontology)
            df$ontologyType      <- as.character(df$ontologyType)
            df$Zscore            <- suppressWarnings(as.numeric(df$Zscore))
            df$Pvalue.Enrichment <- suppressWarnings(as.numeric(df$Pvalue.Enrichment))
          
            # Split “Genes.Hit” into symbol vectors
            split_genes <- function(x) {
              if (is.na(x) || !nzchar(x)) return(character(0))
              parts <- trimws(unlist(strsplit(x, ";", fixed = TRUE)))
              parts[nzchar(parts)]
            }
            genes_list <- lapply(df$Genes.Hit, split_genes)
            df$.n_genes <- vapply(genes_list, length, integer(1))
          
            # Cluster **within each ontologyType**; keep 1 best per cluster
            keep_blocks <- lapply(split(seq_len(nrow(df)), df$ontologyType), function(ii) {
              sub <- df[ii, , drop = FALSE]
              gs  <- lapply(sub$Genes.Hit, split_genes); names(gs) <- sub$ontology
          
              cl_df <- kappa_cluster_ontologies(gs, kappa_cut = kappa_cut, method = method)
              stopifnot(nrow(cl_df) == nrow(sub))
              sub$cluster <- cl_df$cluster
          
              reps <- do.call(rbind, lapply(split(seq_len(nrow(sub)), sub$cluster), function(jj) {
                block <- sub[jj, , drop = FALSE]
                # choose: most significant (lowest p), then fewest genes, then highest Z, then alpha
                o <- order(block$Pvalue.Enrichment, block$.n_genes, -block$Zscore, block$ontology, na.last = TRUE)
                block[o[1L], , drop = FALSE]
              }))
              reps
            })
          
            kept <- do.call(rbind, keep_blocks)
            kept <- kept[, c("ontology","ontologyType","Zscore","Pvalue.Enrichment","Genes.Hit")]
            rownames(kept) <- NULL
            kept
          }

          # Precompute kappa-pruned ontology tables for each list in parallel
          GSA.KAPPA <- NULL
          if (identical(removeRedundantGOterms, "kappa")) {
            clusterLocal <- makeCluster(c(rep("localhost",parallelThreads)) )
            registerDoParallel(clusterLocal)

            # Make sure workers see these
#            parallel::clusterExport(
#              cl = clusterLocal,
#              varlist = c("GSA.FET.outSimple", "minHitsPerOntology",
#                          "kappa_cluster_ontologies", "prune_ontologies_by_kappa"),
#              envir = environment()
#            )

            cat("\nPerforming Cohen's kappa-based graph clustering (k<0.30) for all output ontology lists in parallel using up to ", parallelThreads," threads...\n")
            # igraph is needed (kappa pruning method="graph")
            GSA.KAPPA <- foreach(thismod = uniquemodcolors,
                                 .packages = c("igraph"),
                                 .export = c("kappa_cluster_ontologies", "prune_ontologies_by_kappa"),
                                 .errorhandling = "remove") %dopar% {
              tmp <- GSA.FET.outSimple[[thismod]]
              # Skip empty modules (occur if no genes passed the cutoff)
              if (is.null(tmp) || (is.na(tmp[[1]][1]))) return(NULL)
          
              # Keep rows with enough hits
              idx <- which(tmp[, "Significant (in gene set)"] >= minHitsPerOntology)
              if (!length(idx)) return(NULL)
          
              # Reconstruct tmp2 (the trimmed 5-column data.frame you used before)
              tmp2 <- tmp[idx, c("ontology","ontologyType","Zscore","Pvalue.Enrichment","Genes.Hit")]
              tmp2 <- tmp2[order(tmp2$Zscore, decreasing = TRUE), ]
              tmp2 <- tmp2[order(tmp2$ontologyType, decreasing = TRUE), ]
          
              # Prune via kappa
              prune_ontologies_by_kappa(tmp2, kappa_cut = 0.30, method = "graph")
            }
            names(GSA.KAPPA) <- uniquemodcolors
          }

        }  # end: removeRedundantGOterms=="kappa"
	
	setwd(paste0(filePath,outFilename,"/"))
	redundancyRemovalTag=if(is.logical(removeRedundantGOterms)) { if(removeRedundantGOterms) { "-redundancyRemoved.OBOparent" } else { "" } } else if (removeRedundantGOterms=="kappa") { "-redundancyRemoved.Kbest" } else { "" }
	filenameFinal=paste0(outFilename,redundancyRemovalTag)

	if(!exists("pageDimensions")) pageDimensions=c(8.5,11)
	if(!exists("panelDimensions")) panelDimensions=c(3,2)
	if(!exists("color")) color=c("darkseagreen3","lightsteelblue1","lightpink4","goldenrod","darkorange","gold")
            #colors respectively for ontology Types:
            #"Biological Process","Molecular Function","Cellular Component","Reactome","WikiPathways","MSig.C2","PMC"
	if(!length(color)==6) color=c("darkseagreen3","lightsteelblue1","lightpink4","goldenrod","darkorange","gold")
	if(!exists("maxBarsPerOntology")) maxBarsPerOntology=5
	
	
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
#		cat(unlist(tmp))
		if (is.na(tmp[[1]][1])) {    #*** occurs when "no genes selected due to too strict pcutoff" in GSEA-FET piano -- error bypass; output empty frame for this comparison.
		  moduleTitle <- xlabels.frame[i,"Labels"]
		  frame();
		  mtext(paste0(moduleTitle,"\n-no genes made it into input"), adj=0.5, line=1, cex=0.85, font=2);
		  next;
		}
		filter.minHits.idx<-which(tmp[,"Significant (in gene set)"]>=minHitsPerOntology)
		if (length(tmp[,2]) == 0 | length(filter.minHits.idx) == 0)  { frame(); next; }
		#tmp = tmp[filter.minHits.idx,c(11,10,9,1, 8)] ## Select GO-terms,GO-Type,Z-score,pValues (and previously/again, semicolon-separated gene Lists(8))
                tmp2 <- tmp[filter.minHits.idx, c("ontology","ontologyType","Zscore","Pvalue.Enrichment","Genes.Hit")]
		tmp2 = tmp2[order(tmp2$Zscore,decreasing=T),]
		tmp2 = tmp2[order(tmp2$ontologyType,decreasing=T),]
	
		if(identical(removeRedundantGOterms, "kappa")) { 
                  tmp2 <- GSA.KAPPA[[thismod]]
                  if (is.null(tmp2) || !nrow(tmp2)) { frame(); next; }
		} else if (isTRUE(removeRedundantGOterms)) {
		  #dim(tmp2[which((tmp2$ontologyType == "GOCC" | tmp2$ontologyType == "GOBP" | tmp2$ontologyType == "GOMF") & tmp2$Zscore>=minZ),])
		  ##[1] 124   4
		  go.full.tmp2<-tmp2[which((tmp2$ontologyType == "GOCC" | tmp2$ontologyType == "GOBP" | tmp2$ontologyType == "GOMF") & tmp2$Zscore>=minZ),]
		  go.full.tmp2$Term=gsub("^.*\\%GO..\\%(GO:\\d*)$","\\1",rownames(go.full.tmp2))

		  go.minimal.terms<-minimal_set(ontology.index, terms=go.full.tmp2$Term)
		  #go.pruned.terms<-prune_descendants(ontology.index, roots=?, terms=)
		  go.minimal.tmp2<-go.full.tmp2[which(go.full.tmp2$Term %in% go.minimal.terms),]
	
		  tmp2$Zscore[which(tmp2$ontologyType=="GOCC" | tmp2$ontologyType=="GOBP" | tmp2$ontologyType=="GOMF")] <- 0  #keeps in all GO ontologies, but Z scores zeroed
		  tmp2$Zscore[match(rownames(go.minimal.tmp2),rownames(tmp2))] <-go.minimal.tmp2$Zscore                       #put back Z scores of non-redundant ontologies, to be possibly kept below
	#	  tmp2<-rbind(tmp2,go.minimal.tmp2[,-ncol(go.minimal.tmp2)])
		  tmp2 = tmp2[order(tmp2$Zscore,decreasing=T),]
		  tmp2 = tmp2[order(tmp2$ontologyType,decreasing=T),]
		}
		tmp3 = tmp2[tmp2$ontologyType == "GOBP",][c(1:maxBarsPerOntology),]
		tmp3 = rbind(tmp3,tmp2[tmp2$ontologyType == "GOMF",][c(1:maxBarsPerOntology),] )
		tmp3 = rbind(tmp3,tmp2[tmp2$ontologyType == "GOCC",][c(1:maxBarsPerOntology),] )
		tmp3 = rbind(tmp3,tmp2[tmp2$ontologyType == "REACTOME",][c(1:maxBarsPerOntology),] )
		tmp3 = rbind(tmp3,tmp2[tmp2$ontologyType == "WikiPathways",][c(1:maxBarsPerOntology),] )
		tmp3 = rbind(tmp3,tmp2[tmp2$ontologyType == "MSIGDB_C2",][c(1:maxBarsPerOntology),] )
#		tmp3 = rbind(tmp3,tmp2[tmp2$ontologyType == "PMC",][c(1:maxBarsPerOntology),] )   # PMCid(s) of publication-associated gene lists enriched in your input; Bader Lab added these entries to GMT files starting December 2023.
	
		tmp3 <- na.omit(tmp3)
		tmp3 <- tmp3[which(tmp3$Zscore>=minZ),]
	#	tmp3 <- tmp3[order(tmp3$Zscore,decreasing=T),] #added this row, if you want to mix ontology types and sort by Z.Score only
		tmp3 <- tmp3[rev(rownames(tmp3)),]
	
		summary[[i]] <- tmp3
		
		moduleTitle <- xlabels.frame[i,"Labels"]
		if (is.na(tmp3$ontologyType[1])) { frame(); mtext(paste0(moduleTitle,"\n-no terms hit"), adj=0.5, line=1, cex=0.85, font=2); next; }  #*** occurs when "no genes selected due to too strict pcutoff" in GSEA-FET piano -- error bypass; output empty frame for this comparison.
		
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
#			} else if (tmp3$ontologyType[j] == "PMC"){
#				tmp3$color[j] <- color[7]
			}
	
		# tmp3$color[j] <- uniquemodcolors[i] #module color for all bars, instead of different colors by ontology type
		} 
	
		if (tmp3$Zscore[1] == F) { frame(); mtext(paste0(moduleTitle,"\n-no terms hit"), adj=0.5, line=1, cex=0.85, font=2); next; }
		par(mar=c(4,15,4,3))
		xlim <- c(0,1.1*max(tmp3$Zscore))	
		xh <- barplot(tmp3$Zscore,horiz = TRUE,width =0.85,las=1,main=moduleTitle, xlim=xlim,col=tmp3$color,cex.axis=0.7,xlab="Z Score",cex.lab=0.9,cex.main=0.95,ylim=c(0,nrow(tmp3)+0.8))
		abline(v=minZ,col="red", cex.axis = 0.5)
		axis(2, at=xh, labels = tmp3$ontology, tick=FALSE, las =2, line =-0.5, cex.axis = 0.7) 
	}
	
	par(op) # Leaves the last plot
	dev.off()



	## Master Tables Generation for output to csv and for GO:CC Z score coclustering

	for (this.input in names(GSA.FET.outSimple)) {
	  if (length(GSA.FET.outSimple[[this.input]])==1) { if (is.na(GSA.FET.outSimple[[this.input]][1])) {
	    cat(paste0("x - ",this.input," list had no genes and is excluded form table output.\n"))
	    GSA.FET.outSimple[[this.input]]<-NULL
	    uniquemodcolors<-uniquemodcolors[-which(uniquemodcolors==this.input)]
	  }}
	}
	
	if(length(uniquemodcolors)<1) {
	  cat("x - no lists had any genes/were enriched in any terms below cutoffs. Skipping output of tables, etc.\n")
	} else {
	  GSA.FET.collapsed.outSimple.Zscore <- cbind(GSA.FET.outSimple[[1]][,c("ontology","ontologyType")], matrix(unlist(lapply(GSA.FET.outSimple, function(x) x[,"Zscore"] )),byrow=FALSE,ncol=length(uniquemodcolors)), matrix(unlist(lapply(GSA.FET.outSimple, function(x) x[,"Genes.Hit"] )),byrow=FALSE,ncol=length(uniquemodcolors)))
	  GSA.FET.collapsed.outSimple.Pvalue.Enrichment <- cbind(GSA.FET.outSimple[[1]][,c("ontology","ontologyType")], matrix(unlist(lapply(GSA.FET.outSimple, function(x) x[,"Pvalue.Enrichment"] )),byrow=FALSE,ncol=length(uniquemodcolors)), matrix(unlist(lapply(GSA.FET.outSimple, function(x) x[,"Genes.Hit"] )),byrow=FALSE,ncol=length(uniquemodcolors)))
	  GSA.FET.collapsed.outSimple.Enrichment.FDR.BH <- cbind(GSA.FET.outSimple[[1]][,c("ontology","ontologyType")], matrix(unlist(lapply(GSA.FET.outSimple, function(x) x[,"Adjusted.Enr.Pvalue"] )),byrow=FALSE,ncol=length(uniquemodcolors)), matrix(unlist(lapply(GSA.FET.outSimple, function(x) x[,"Genes.Hit"] )),byrow=FALSE,ncol=length(uniquemodcolors)))
	  colnames(GSA.FET.collapsed.outSimple.Zscore)[3:(length(uniquemodcolors)*2+2)] <- colnames(GSA.FET.collapsed.outSimple.Pvalue.Enrichment)[3:(length(uniquemodcolors)*2+2)] <- colnames(GSA.FET.collapsed.outSimple.Enrichment.FDR.BH)[3:(length(uniquemodcolors)*2+2)] <- c(names(GSA.FET.outSimple), paste0(names(GSA.FET.outSimple),"_Genes.Hit"))
	  
	  #write master tables
	  write.table(GSA.FET.collapsed.outSimple.Zscore,file=paste0(filePath,outFilename,"/GSA-GO-FET_",outFilename,"-Zscores.txt"),row.names=FALSE,col.names=TRUE,sep="\t", quote=FALSE)
	  write.table(GSA.FET.collapsed.outSimple.Pvalue.Enrichment,file=paste0(filePath,outFilename,"/GSA-GO-FET_",outFilename,"-Enr.Pvalues.txt"),row.names=FALSE,col.names=TRUE,sep="\t", quote=FALSE)
	  write.table(GSA.FET.collapsed.outSimple.Enrichment.FDR.BH,file=paste0(filePath,outFilename,"/GSA-GO-FET_",outFilename,"-Enr.FDR.BH.txt"),row.names=FALSE,col.names=TRUE,sep="\t", quote=FALSE)
	
	  if(!exists("cocluster")) cocluster=TRUE

	  if(cocluster) {
	  ## For GO:CC Z score coclustering
	  GSA.FET.GOCC.Zscore <- GSA.FET.collapsed.outSimple.Zscore[which(GSA.FET.collapsed.outSimple.Zscore$ontologyType=="GOCC"), 1:(length(uniquemodcolors)+2)]
	  temp.rownames=rownames(GSA.FET.GOCC.Zscore)
	  GSA.FET.GOCC.Zscore<-cbind(GSA.FET.GOCC.Zscore[,1:2], apply(GSA.FET.GOCC.Zscore[,3:ncol(GSA.FET.GOCC.Zscore)],2,as.numeric))
	  rownames(GSA.FET.GOCC.Zscore)<-temp.rownames
	  GSA.FET.GOCC.terms <- gsub("^.*\\%GO..\\%(GO:\\d*)$","\\1",rownames(GSA.FET.GOCC.Zscore))
	  minZ<-qnorm(1e-5/2, lower.tail=FALSE)
	  if(nrow(GSA.FET.GOCC.Zscore)>1) {
	    GSA.FET.GOCC.terms.minZreached <- apply(GSA.FET.GOCC.Zscore[,3:ncol(GSA.FET.GOCC.Zscore)],1,function(x) if (max(x)>=minZ) { TRUE } else { FALSE } )

	    if (removeRedundantGOterms=="kappa") {
              
              ## --- Build gene sets for significant GO:CC terms
              sig_idx   <- which(GSA.FET.GOCC.terms.minZreached)
              sig_terms <- GSA.FET.GOCC.terms[sig_idx]
              
              # Assemble a named list of gene symbols per GO term (using GSCfromGMT$gsc names like "...%GO:####")
              gene_sets_present <- setNames(vector("list", length(sig_terms)), sig_terms)
              for (i in seq_along(sig_terms)) {
                go_id <- sig_terms[i]
                hit   <- grep(paste0("%", go_id, "$"), names(GSCfromGMT$gsc))
                gene_sets_present[[i]] <- if (length(hit)) unique(GSCfromGMT$gsc[[hit[1] ]]) else character(0)
              }
              
              ## --- Kappa clustering of the significant terms
              # Expecting a data.frame with columns: term, cluster
              kc <- kappa_cluster_ontologies(gene_sets_present, kappa_cut = 0.30, method = "graph")
              
              ## --- Pick one representative per cluster: max positive Z anywhere in the row
              # Compute per-term max positive Z across all Z-score columns (3:ncol)
              row_maxZ <- apply(GSA.FET.GOCC.Zscore[, 3:ncol(GSA.FET.GOCC.Zscore), drop = FALSE],
                                1, function(z) suppressWarnings(max(as.numeric(z), na.rm = TRUE)))
              names(row_maxZ) <- GSA.FET.GOCC.terms
              
              # Restrict to clustered significant terms
              ok <- kc$term %in% names(row_maxZ)
              kc  <- kc[ok, , drop = FALSE]
              
              # For each cluster, keep the term with the largest maxZ
              split_terms <- split(kc$term, kc$cluster)
              keep_terms  <- unlist(lapply(split_terms, function(tt) tt[which.max(row_maxZ[tt])]), use.names = FALSE)
              
              # Update the logical filter to keep only the chosen representatives
              GSA.FET.GOCC.terms.minZreached[] <- FALSE
              GSA.FET.GOCC.terms.minZreached[ match(keep_terms, GSA.FET.GOCC.terms) ] <- TRUE
              
              # For downstream convenience, also set the minimal terms vector
              GSA.FET.GOCC.minimal.terms <- keep_terms
              
	    } else if (removeRedundantGOterms) {
	      GSA.FET.GOCC.minimal.terms <- minimal_set(ontology.index, terms=GSA.FET.GOCC.terms[which(GSA.FET.GOCC.terms.minZreached)])
	    } else {
	      GSA.FET.GOCC.minimal.terms <- GSA.FET.GOCC.terms[which(GSA.FET.GOCC.terms.minZreached)]
	      cat("- Note: removeRedundantGOterms=FALSE.  You can reduce cellular component terms in the coclustering by setting this variable to TRUE, or 'kappa'.\n\n")
	    }

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
            
            if(ncol(data)==0) cat(paste0("- No highly significant (Z>",signif(minZ,3),") Cellular Component ontologies found. Skipping GOCC Cluster Heatmap output.\n\n"))
            if(ncol(data)==1) cat(paste0("- Only one highly significant (Z>",signif(minZ,3),") Cellular Component ontologies found. Skipping GOCC Cluster Heatmap output.\n\n"))
            if(ncol(data)>1) {

	      data[matrixdata>4]<-4
              data[matrixdata< -4] <- -4
              
              #NMF - initial approach
              suppressPackageStartupMessages(require(WGCNA,quietly=TRUE))
              bw<-colorRampPalette(c("#0058CC", "white"))
              wr<-colorRampPalette(c("white", "#CC3300"))
              colvec<-c(bw(50),wr(50))

              colnames(data)<-gsub("\\%GOCC\\%"," | ",gsub("GOcc","GOCC",colnames(data)))

              if(!modulesInMemory) {
                uniquemodcolors=labels2colors(1:((ncol(GSA.FET.collapsed.outSimple.Zscore)-2)/2))  #variable reused for color annotation here; meaningless for non-WGCNA lists  #/2 because Genes.Hit double the columns.
                myRowAnnotation=data.frame(Lists=as.numeric(factor(uniquemodcolors,levels=sort(uniquemodcolors))))
                heatmapLegendColors<-list(Lists=factor(sort(uniquemodcolors)))
              } else {
                heatmapLegendColors<-list(Modules=sort(uniquemodcolors))
                myRowAnnotation=data.frame(Modules=as.numeric(factor(uniquemodcolors,levels=sort(uniquemodcolors))))
              }

              suppressPackageStartupMessages(require(NMF,quietly=TRUE))  # for aheatmap
              pdf(file=paste0("GO_cc_clustering_from_GSA_FET_Z-",filenameFinal,".pdf"),width=18.5,height=18,onefile=FALSE)
                aheatmap(x=data, ## Numeric Matrix
                       main="Co-clustering with manhattan distance function, ward metric",
              #         annCol=metdat,  ## Color swatch and legend annotation of columns/samples and rows (not used)
                       annRow=myRowAnnotation,
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
            } # end if(ncol(data)>1)
            } else { cat("- None or only one significant Cellular Component ontologies (rows) found. Skipping GOCC Cluster Heatmap output.\n\n") }  # end if(nrow(GSA.FET.GOCC.Zscore)>1)
          } # end if (cocluster) 
        } # end if(length(uniquemodcolors)<1)

	# Optional (if bubble=TRUE): save a GSEA-style bubble plot to a PDF file for each input list's enriched ontologies.
	if (exists("bubble")) if (bubble) {
          cat("- Plotting GSEA-style bubble plots...\n")
          
	  GO.bubblePlot(ZinputFile = paste0(filePath,outFilename,"/GSA-GO-FET_",outFilename,"-Zscores.txt"),
	                GMTfile = if(file.exists(paste0("../", GMTdatabaseFile))) { paste0("../", GMTdatabaseFile) } else {
	                 GMTdatabaseFile },  # direct full path to GMT required if not up one level in the folder structure from the working directory at the time GOparallel was called.
	                GO.OBOfile = if(file.exists(paste0("../",GO.OBOfile))) paste0("../",GO.OBOfile),  # function will fallback to redownload if not found.
	                keep_ontology_types = c("GOBP", "GOMF", "GOCC", "REACTOME", "WikiPathways", "MSIGDB_C2"),
	                ontology_to_color = c("GOBP"=color[1], "GOMF"=color[2], "GOCC"=color[3], "REACTOME"=color[4], "WikiPathways"=color[5], "MSIGDB_C2"=color[6]),
	                MAX_TERMS_PER_ONTOLOGY_TYPE = maxBarsPerOntology,
	                max_label_width = 50)
	}

        setwd(filePath)

}


## GSEA-style Bubble plot PDF output function - stand alone to run after GOparallel on Zscore txt output, or with bubble=TRUE argument within GOparallel()
GO.bubblePlot <- function(ZinputFile=NULL, GMTfile=NULL, GO.OBOfile="./go.obo", removeRedundantGOterms=TRUE, 
                          keep_ontology_types = c("GOBP", "GOMF", "GOCC"),
                          ontology_to_color = c("GOBP"="darkseagreen3", "GOMF"="lightsteelblue1", "GOCC"="lightpink4", "REACTOME"="goldenrod", "WikiPathways"="darkorange", "MSIGDB_C2"="gold", "MSIGDBHALLMARK"="greenyellow"),
                          MAX_TERMS_PER_ONTOLOGY_TYPE = 5,
                          max_label_width = 50) {

  # Load required libraries
  suppressPackageStartupMessages(require(ggplot2,quietly=TRUE))
  suppressPackageStartupMessages(require(dplyr,quietly=TRUE))
  suppressPackageStartupMessages(require(tibble,quietly=TRUE))
  suppressPackageStartupMessages(require(stringr,quietly=TRUE))
  suppressPackageStartupMessages(require(scales,quietly=TRUE)) # For rescale function
  suppressPackageStartupMessages(require(grid,quietly=TRUE)) # For grob manipulation - to take out one tick mark!
  suppressPackageStartupMessages(require(gtable,quietly=TRUE))
  
  minHitsPerOntology <- 5        # Ontologies hit by fewer than this number of genes will not be plotted; 5 is hard coded into GOparallel function and required to get the same ontologies for the 3 ontology types output by GOparallel
  keep_ontology_types <- intersect(keep_ontology_types, c("GOBP", "GOMF", "GOCC", "REACTOME", "WikiPathways", "MSIGDB_C2", "MSIGDBHALLMARK"))
  if(length(keep_ontology_types)<1) stop("Error:  Argument 'keep_ontology_types' must specify one or more supported ontologies ('GOBP', 'GOMF', 'GOCC', 'REACTOME', 'WikiPathways', 'MSIGDB_C2', 'MSIGDBHALLMARK')")
  
  if (is.null(ZinputFile)) stop("Error:  No GOparallel Zscore txt file was specified as input.") 
  if (!file.exists(ZinputFile)) stop("Error:  ZinputFile not found.")
   
  if (is.null(GMTfile)) stop("Error:  No GMT file (ontology database from Bader Lab website) specified. This file must be specified and match the version used to generate the ZinputFile.") 
  if (!file.exists(GMTfile)) stop("Error:  GMT file (ontology database from Bader Lab website) not found. This file must be specified and match the version used to generate the ZinputFile.") 
  
  
  # Check if removeRedundantGOterms is specified
  if(!exists("removeRedundantGOterms")) {
      cat("- removeRedundantGOterms not specified TRUE/FALSE. Removing them as the default, using go.obo and ontologyIndex package.\n")
      removeRedundantGOterms <- TRUE
  }
  if(removeRedundantGOterms) {
      if (!exists("GO.OBOfile")) GO.OBOfile <- "../go.obo"
      if (!file.exists(GO.OBOfile)) {
          suppressPackageStartupMessages(require(curl, quietly=TRUE))
          OBOtargetPath <- gsub("(.*\\/).*$", "\\1", GO.OBOfile)
          if (!dir.exists(OBOtargetPath)) dir.create(OBOtargetPath)
          curr.dir <- getwd()
          setwd(OBOtargetPath)
          cat(paste0("- Downloading go.obo file for main GO term redundancy cleanup...\n...to location:  ", OBOtargetPath, "go.obo\n"))
          curl_download(url="http://current.geneontology.org/ontology/go.obo", destfile="go.obo", quiet = TRUE, mode = "w")
          setwd(curr.dir)
          cat("GO.OBOfile set to downloaded file: ", paste0(OBOtargetPath, "go.obo"), "\n")
          GO.OBOfile <- paste0(OBOtargetPath, "go.obo")
      }
  }
  
  if (removeRedundantGOterms) {
      suppressPackageStartupMessages(require(ontologyIndex))
      ontology.index <- get_OBO(file=GO.OBOfile, extract_tags="everything")
  }
  
  # Read input files
  df1 <- read.table(ZinputFile, header=TRUE, sep="\t", stringsAsFactors=FALSE, quote = "\"")
  
  # Ensure the second column is named "OntologyType"
  colnames(df1)[2] <- "OntologyType"
  
  # Read GMT file
  lines <- readLines(GMTfile)
  df2_list <- strsplit(lines, "\t")
  
  # Parse the GMT file
  df2_parsed <- lapply(df2_list, function(x) {
      matches <- regmatches(x[1], regexec("^(.*)%(.+)%(.+)$", x[1]))
      if (length(matches[[1]]) >= 4) {
          OntologyType <- matches[[1]][3]
          ontology <- matches[[1]][2]
          term <- matches[[1]][4]
      } else {
          OntologyType <- NA
          ontology <- NA
          term <- NA
      }
      genes <- x[3:length(x)]
      list(OntologyType=OntologyType, ontology=ontology, term=term, genes=genes)
  })
  
  # Create data frame
  df2_df <- tibble(
      OntologyType = sapply(df2_parsed, function(x) x$OntologyType),
      ontology = sapply(df2_parsed, function(x) stringr::str_to_title(x$ontology)),
      term = sapply(df2_parsed, function(x) x$term),
      genes = lapply(df2_parsed, function(x) x$genes)
  )
  
  # Copy WP# pathway IDs to term column
  df2_df$term[which(grepl("^WP\\d+",df2_df$OntologyType))]<-df2_df$OntologyType[which(grepl("^WP\\d+",df2_df$OntologyType))]
  # Revert WP# to 'WikiPathways' in OntologyType
  df2_df$OntologyType[which(grepl("^WP\\d+",df2_df$OntologyType))]<-"WikiPathways"
  # Remove '%Wikipathways_YYYYMMDD' from terms
  df2_df$ontology[which(df2_df$OntologyType=="WikiPathways")]<-gsub("%Wikipathways_\\d+$","",df2_df$ontology[which(df2_df$OntologyType=="WikiPathways")])
  
  
  df2_df <- df2_df[df2_df$OntologyType %in% keep_ontology_types, ]
  
  # Get the column names of df1
  colnames_df1 <- colnames(df1)
  
  # Identify Z score columns (excluding "ontology" and "OntologyType", and columns ending with "_Genes.Hit")
  zscore_cols <- colnames_df1[!(colnames_df1 %in% c("ontology", "OntologyType")) & !grepl("_Genes\\.Hit$", colnames_df1)]
  
  # Identify Genes.Hit columns
  genes_hit_cols <- paste0(zscore_cols, "_Genes.Hit")
  
  # Check that the Genes.Hit columns exist
  if (!all(genes_hit_cols %in% colnames_df1)) {
      stop("Some Genes.Hit columns corresponding to Z score columns are missing in input_file1.txt")
  }
  
  # Loop over each Z score column
  for (i in seq_along(zscore_cols)) {
      zscore_col <- zscore_cols[i]
      genes_hit_col <- genes_hit_cols[i]
      plot_title <- zscore_col
  
      # Prepare the data
      df_current <- df1[, c("ontology", "OntologyType", zscore_col, genes_hit_col)]
      colnames(df_current) <- c("ontology", "OntologyType", "Zscore", "Genes.Hit")
  
      # Merge data frames
      df_merged <- merge(df_current, df2_df, by=c("OntologyType", "ontology"))
  
      # Compute number of genes hit
      df_merged$num_genes_hit <- sapply(df_merged$Genes.Hit, function(x) length(strsplit(x, ";")[[1]]))
  
      # Compute total number of genes in ontology
      df_merged$num_genes_total <- sapply(df_merged$genes, length)
  
      # Compute gene ratio
      df_merged$gene_ratio <- df_merged$num_genes_hit / df_merged$num_genes_total
  
      # Correct p-value calculation: one-tailed test using log.p to handle very small p-values
      # Compute log p-value
      log_pvalue <- pnorm(df_merged$Zscore, lower.tail = FALSE, log.p = TRUE)
      # Compute -log10(p-value) directly from log_pvalue
      df_merged$neg_log10_pvalue <- -log_pvalue / log(10)
  
      # Handle infinite values by setting a maximum threshold
      max_neg_log10_pvalue <- 50  # Adjust this threshold as needed
      df_merged$neg_log10_pvalue <- pmin(df_merged$neg_log10_pvalue, max_neg_log10_pvalue)
  
      # Ensure neg_log10_pvalue is at least 1.30
      df_merged$neg_log10_pvalue <- pmax(df_merged$neg_log10_pvalue, 1.30)
          
      # Filter for Zscore > 1.96
      minZ <- qnorm(0.05/2, lower.tail=FALSE)
      df_filtered <- df_merged[df_merged$Zscore > minZ, ]
  
      # Filter for minHitsPerOntology
      df_filtered <- df_filtered[df_filtered$num_genes_hit >= minHitsPerOntology, ]
      
      # Keep only specified OntologyTypes
      df_filtered <- df_filtered[df_filtered$OntologyType %in% keep_ontology_types, ]
      
      # If no significant ontologies found, skip plotting
      if (nrow(df_filtered) == 0) {
          message(paste("No significant ontologies found for", zscore_col))
          next
      }
      
      if(removeRedundantGOterms) {
          df_filtered.coreGO <- df_filtered[df_filtered$OntologyType %in% c("GOCC", "GOBP", "GOMF") & df_filtered$Zscore >= minZ, ]
          go.minimal.terms <- minimal_set(ontology.index, terms=df_filtered.coreGO$term)
          df_filtered.minimal.terms <- df_filtered.coreGO[df_filtered.coreGO$term %in% go.minimal.terms, ]
  
          df_filtered.noCoreGO <- df_filtered[!df_filtered$OntologyType %in% c("GOCC", "GOBP", "GOMF"), ]
          df_filtered <- rbind(df_filtered.noCoreGO, df_filtered.minimal.terms)
      }
  
      rename_map <- c(
        "GOBP"           = "BP",
        "GOMF"           = "MF",
        "GOCC"           = "CC",
        "REACTOME"       = "RE",
        "WikiPathways"   = "WP",
        "MSIGDB_C2"      = "MS",
        "MSIGDBHALLMARK" = "HM"
      )
  
      # Select top 5 (or fewer) ontologies per OntologyType
      df_top <- df_filtered %>%
        group_by(OntologyType) %>%
        arrange(desc(Zscore), .by_group = TRUE) %>%
        slice_head(n = MAX_TERMS_PER_ONTOLOGY_TYPE) %>%
        ungroup() %>%
        mutate(
          OntologyType = as.character(OntologyType),
          OntologyType = recode(
            OntologyType,
            !!!rename_map[keep_ontology_types]
          ),
          OntologyType = factor(OntologyType)   # optional
        )
  
  
      # Create ontology label
      df_top$ontology_label <- df_top$ontology #paste(df_top$OntologyType, df_top$ontology, sep=": ")
      
      # Reorder the ontology types so that CC follows MF
      df_top$OntologyType <- factor(df_top$OntologyType, levels = c(if ("GOBP" %in% keep_ontology_types) "BP",
                                                                    if ("GOMF" %in% keep_ontology_types) "MF",
                                                                    if ("GOCC" %in% keep_ontology_types) "CC",
                                                                    if ("REACTOME" %in% keep_ontology_types) "RE",
                                                                    if ("WikiPathways" %in% keep_ontology_types) "WP",
                                                                    if ("MSIGDB_C2" %in% keep_ontology_types) "MS",
                                                                    if ("MSIGDBHALLMARK" %in% keep_ontology_types) "HM"))
  
      # Order ontology_label factor levels after ordering OntologyType
      df_top <- df_top %>%
          arrange(OntologyType, desc(Zscore))
      
      # Handle duplicates in ontology labels if they exist
      if(length(which(duplicated(df_top$ontology_label)))>0) df_top$ontology_label[which(duplicated(df_top$ontology_label))]<-paste0(".",df_top$ontology_label[which(duplicated(df_top$ontology_label))])
      if(length(which(duplicated(df_top$ontology_label)))>0) df_top$ontology_label[which(duplicated(df_top$ontology_label))]<-paste0("..",df_top$ontology_label[which(duplicated(df_top$ontology_label))])
      
      df_top$ontology_label <- factor(df_top$ontology_label, levels=rev(unique(df_top$ontology_label)))
      
      # Map OntologyType to specified colors
      ontology_colors <- c()
      if ("GOBP" %in% keep_ontology_types) ontology_colors <- c(ontology_colors, BP = ontology_to_color[["GOBP"]])
      if ("GOMF" %in% keep_ontology_types) ontology_colors <- c(ontology_colors, "MF" = ontology_to_color[["GOMF"]])
      if ("GOCC" %in% keep_ontology_types) ontology_colors <- c(ontology_colors, "CC" = ontology_to_color[["GOCC"]])
      if ("REACTOME" %in% keep_ontology_types) ontology_colors=c(ontology_colors, "RE" = ontology_to_color[["REACTOME"]])
      if ("WikiPathways" %in% keep_ontology_types) ontology_colors=c(ontology_colors, "WP" = ontology_to_color[["WikiPathways"]])
      if ("MSIGDB_C2" %in% keep_ontology_types) ontology_colors=c(ontology_colors, "MS" = ontology_to_color[["MSIGDB_C2"]])
      if ("MSIGDBHALLMARK" %in% keep_ontology_types) ontology_colors=c(ontology_colors, "HM" = ontology_to_color[["MSIGDBHALLMARK"]])
      
      # Wrap the ontology labels
      df_top$ontology_label <- stringr::str_wrap(df_top$ontology_label, width = max_label_width)
      
      # Order ontology_label factor levels (after wrapping)
      df_top$ontology_label <- as.character(df_top$ontology_label)
      df_top$ontology_label <- factor(df_top$ontology_label, levels = rev(unique(df_top$ontology_label)))
      
      dummy_label <- paste(rep("_", ceiling(max_label_width*0.85)), collapse = "")
  
      # Adjust the breaks for gene_ratio
      min_gene_ratio <- min(df_top$gene_ratio)
      max_gene_ratio <- max(df_top$gene_ratio)
      
      # Generate breaks for gene_ratio size scale
      breaks_gene_ratio <- seq(min_gene_ratio, max_gene_ratio, length.out = 5)
      breaks_gene_ratio <- unique(round(breaks_gene_ratio, digits = 4))  # Round to 4 decimals to avoid duplicates
      
      # Convert breaks to percentages with 1 decimal place
      labels_gene_ratio <- sprintf("%.1f%%", breaks_gene_ratio * 100)
  
      # Create a dummy row
      dummy_row <- df_top[1, ]  # Copy structure
      dummy_row$OntologyType <- NA
      dummy_row$ontology <- NA
      dummy_row$Zscore <- NA
      dummy_row$neg_log10_pvalue <- NA  #min(df_top$neg_log10_pvalue)
      dummy_row$gene_ratio <- NA  #min(df_top$gene_ratio)
      dummy_row$Genes.Hit <- NA
      dummy_row$num_genes_hit <- NA
      dummy_row$num_genes_total <- NA
      dummy_row$genes <- NA
      dummy_row$ontology_label <- dummy_label
      dummy_row$y_position <- nrow(df_top)+1
      dummy_row$txtColor="#FFFFFF00"
      
      df_top$y_position<-as.numeric(rownames(df_top))  #c(1:nrow(df_top))
      df_top$txtColor<-"black"
      
      # Add the dummy row to df_top
      df_top <- rbind(df_top, dummy_row)
  #    df_top$OntologyType.label <- paste0(df_top$OntologyType,".",df_top$ontology_label)
  #    df_top$OntologyType.label <- factor(df_top$OntologyType.label, levels=rev(unique(df_top$OntologyType.label)))
      
      # Create a named vector of labels, (the dummy label is set to a long set-length string)
      y_labels <- setNames(as.character(df_top$ontology_label), df_top$ontology_label)
      #y_labels[dummy_label] <- ""
      
      #all_black_last_white=c(rep("black",length(y_labels)-1),"#FFFFFF00")
  
      # Get plot height: 342 px = 4.75 in; 252 px = 15.6x;  1x=16.15px
      this.plotHeight=4.75*(16.15*0.6 +(nrow(df_top)-1)*16.15 +90)/342
      # Plot
      pdf(paste0("bubblePlot_", zscore_col, ".pdf"), width=8.8, height=this.plotHeight, bg = "transparent")
      
      # Create the plot and store it in a variable
      p <- ggplot(df_top, aes(x = neg_log10_pvalue, y = ontology_label, size = gene_ratio, fill = OntologyType)) +
          geom_point(shape = 21, color = "black", stroke = 0.5) +
          scale_size_continuous(
              name = "Gene Ratio",
              breaks = breaks_gene_ratio,
              labels = labels_gene_ratio,
              guide = guide_legend(override.aes = list(fill = "grey70", color = "black")),
              range = c(3, 8)
          ) +
          scale_fill_manual(
              values = ontology_colors,
              guide = "none"
          ) +
          xlab(expression(bold(-log[10](italic(p))))) +
          ylab("Ontology") +
          scale_y_discrete(labels = y_labels) +
  #        ggtitle(bquote(bold(.(plot_title)))) +
          theme_bw(base_size = 14, base_family = "Helvetica") +
          theme(
              panel.background = element_rect(fill = "transparent", color = NA),
              plot.background = element_rect(fill = "transparent", color = NA),
              axis.text.y = element_text(size = 9.4, hjust = 1, vjust=0.5, colour=df_top$txtColor),
              axis.title.y = element_blank(),
  #            axis.ticks.length.y = c(rep(unit(3,"points"),nrow(df_top)-1),unit(0,"points")),
              panel.grid.major.x = element_blank(),
              panel.grid.minor = element_blank(),
              panel.grid.major.y = element_line(color = "grey", linetype = 2, size = 0.25),
              legend.background = element_rect(fill = "transparent", colour = NA),
              legend.key = element_rect(fill = "transparent", colour = NA),
              plot.margin = unit(c(1, 4, 1, 1), "lines"),  # Adjust left margin
              panel.border = element_rect(fill = NA)
          )
      
      # Compute x-axis limits: data range plus 5% on each extreme
      x_min_data <- min(df_top$neg_log10_pvalue, na.rm = TRUE)
      x_max_data <- max(df_top$neg_log10_pvalue, na.rm = TRUE)
      x_range <- x_max_data - x_min_data
      
      # Handle the case when x_range is zero
      if (x_range == 0) {
          # Set a small range around the data point
          x_min <- x_min_data - 1
          x_max <- x_max_data + 1
          x_range <- x_max - x_min
          x_min_hbar <- x_min - 0.1
      } else {
          x_min <- x_min_data - 0.05 * x_range
          x_max <- x_max_data + 0.05 * x_range
          x_min_hbar <- x_min_data - (x_range * 0.05) * 2
      }
  
      # Add coordinate limits and prevent clipping
      p <- p + coord_cartesian(xlim = c(x_min, x_max), clip = "off")
      
      # Determine the position of the vertical bar outside the plot area
      offset <- x_range * 0.055
      bar_width <- x_range * 0.05 * 2
      
      # Ensure offset and bar_width have minimum values
      offset <- max(offset, 0.1)
      bar_width <- max(bar_width, 0.2)
      
      bar_xmin <- x_max + offset
      bar_xmax <- bar_xmin + bar_width
      text_x <- bar_xmin + (bar_xmax - bar_xmin) / 2
      
      # Create label data with y positions
      df_top <- df_top %>%
          mutate(y_position = as.numeric(ontology_label))  # concatenated Type.label ensures ontologies named identically across ontology types are uniquely numbered!
      
      # Get the total number of ontologies
      total_terms <- length(unique(df_top$ontology_label))
  
      # Group labels by OntologyType to get y_min and y_max for each group
      grouped_labels <- df_top %>%
          group_by(OntologyType) %>%
          summarise(
              y_min = min(y_position) - 0.5,
              y_max = max(y_position) + 0.5
          ) %>%
          arrange(y_min)
      
      # Adjust y_min and y_max to ensure bars are contiguous and cover the entire y-axis
      grouped_labels$y_min[1] <- 0.5
      grouped_labels$y_max[nrow(grouped_labels)] <- total_terms + 0.5
      
      if (nrow(grouped_labels) > 1) {
          for (ii in 2:nrow(grouped_labels)) {
              grouped_labels$y_min[ii] <- grouped_labels$y_max[ii - 1]
          }
      }
      
      # Determine the position of the vertical bar outside the plot area
      offset <- x_range * 0.055
      bar_xmin <- x_max + offset
      bar_xmax <- bar_xmin + (x_range * 0.05) * 2
      text_x <- bar_xmin + (bar_xmax - bar_xmin) / 2
  
      grouped_labels1<-grouped_labels
      grouped_labels1[1,"y_min"]<-floor(grouped_labels1[1,"y_min"])
      grouped_labels1[nrow(grouped_labels1),"y_max"]<-ceiling(grouped_labels1[nrow(grouped_labels1),"y_max"])
  
      title_box1<-grouped_labels1[nrow(grouped_labels1),]    
      grouped_labels1<-grouped_labels1[-nrow(grouped_labels1),]    
  
      # Add geom_rect for the vertical bars
      p <- p + geom_rect(
          data = grouped_labels1,
          aes(
              xmin = bar_xmin,
              xmax = bar_xmax,
              ymin = y_min,
              ymax = y_max,
              fill = OntologyType
          ),
          color = "black",
          size = 0.5,
          inherit.aes = FALSE
      )
      
      # Add geom_text for the vertical rotated bold text
      p <- p + geom_text(
          data = grouped_labels,
          aes(
              x = text_x,
              y = (y_min + y_max) / 2,
              label = gsub("GO", "", OntologyType)
          ),
          angle = 270,
          fontface = "bold",
          color = "black",
          size = 4,
          inherit.aes = FALSE
      )
  
      # Add horizontal title bar
      p <- p + geom_rect(
          data = title_box1,
          aes(
              xmin = x_min_hbar,
              xmax = bar_xmin,
              ymin = y_min+(y_max-y_min)*0.1,
              ymax = y_max,
              fill = "grey50"
          ),
          color = "black",
          size = 0,
          inherit.aes = FALSE
      )
      # Add horizontal title bar at the top
      title_bar_height <- 0.5  # Adjust as needed
      title_y_min <- total_terms + 0.5
      title_y_max <- title_y_min + title_bar_height
      title_text1<-data.frame(plot_title=plot_title, x=x_min_hbar, y=max(df_top$y_position))
      # Add geom_text for the plot title within the horizontal bar
      p <- p + geom_text(
          data=title_text1,
          aes(
              x = x,
              y = y,
              label = paste0("  ",plot_title)
          ),
          angle = 0,
          fontface = "bold",
          color = "white",
          size = 5.2,
          hjust = 0,
          vjust = -0.1,
          inherit.aes = FALSE
      )
  
      
      # Adjust the left margin to allocate fixed space for y-axis labels
      p <- p + theme(
          plot.margin = unit(c(1, 4, 1, 0.5), "lines"),  # Top, right, bottom, left
          axis.text.y = element_text(size = 10, hjust = 1, vjust = 0.5)
      )
      
  #    # Explicitly print the plot
  #    print(p)
      
  # Build the ggplot object
  g <- ggplot_build(p)
  gt <- ggplot_gtable(g)
  
  # Find the index of the y-axis in the gtable layout
  index <- which(gt$layout$name == "axis-l")
  
  # Extract the y-axis grob
  axis.grob <- gt$grobs[[index]]
  
  ## Extract the ticks grob
  ##ticks.grob <- axis.grob$children$axis$children[[which(axis.grob$children$axis$childrenOrder == "ticks")]]
  #
  ## Get the number of ticks
  #n.ticks <- length(ticks.grob)
  #
  ## Set the length of the last tick to zero
  #ticks.grob[n.ticks] <- 2 #ticks.grob$x0[n.ticks]
  #
  ## Assign the modified ticks.grob back to the axis.grob
  #axis.grob$children$axis$grobs[[1]][[4]] <- ticks.grob
  tick.parent.grob.idx<-which(sapply(axis.grob$children$axis$grobs, function(x) inherits(x, "polyline")))
  tick.attr.idx<-which(sapply(axis.grob$children$axis$grobs[[1]], function(x) inherits(x, "gpar")))
  axis.grob$children$axis$grobs[[tick.parent.grob.idx]][[tick.attr.idx]]$col <- c(rep(axis.grob$children$axis$grobs[[tick.parent.grob.idx]][[tick.attr.idx]]$col[1],nrow(df_top)-1),"transparent")
  
  # Assign the modified axis.grob back to the gtable
  gt$grobs[[index]] <- axis.grob
      
  grid.draw(gt)
      
      dev.off()
  }
}
