# GOparallel
Uses BaderLab's monthly updated .GMT formatted ontology gene lists for Fisher exact enrichment into differential or co-expression or other user-supplied gene lists.
GO.obo file download is required for redundant term pruning.
WGCNA modules can be co-clustered for relatedness of modules by GO-cellular component terms hit or depleted per signed Z score calculation.

Compatible with the Seyfried systems biology pipeline.

See Instructions in the R code before and in the configuration section (top).

Sample output for WGCNA modules, ANOVA-Tukey table from pipeline outputs are in <a href="https://github.com/edammer/GOparallel/GOparallel-SampleOutput.zip">this .zip</a>,
using the provided sample pipeline input RData for either, using the 8619 protein data from <a href="https://www.nature.com/articles/s41593-021-00999-y">Johnson ECB, et al, Nat Neurosci., 2022</a>.

Note the current analysis output samples are using June 2022 ontologies, much newer compared to those in the publication, which was based on the GO-Elite Ensembl v62+ database.
