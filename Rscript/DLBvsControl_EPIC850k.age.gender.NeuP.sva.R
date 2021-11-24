################################################################################
# DNA methylation analysis with RnBeads
# For project Stewart_DLB
# ------------------------------------------------------------------------------
#  Vanilla analysis of the Stewart DLBvsHealth EPIC 850K dataset
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# created: 2021-01-21
# 
# http://rnbeads.mpi-inf.mpg.de/
################################################################################
# (0) Preliminaries
################################################################################
# load the package
library(RnBeads)
library(RnBeads.hg19)

# define the directory structure
dataDir <- file.path(getwd(), "data")
resultDir <- file.path(getwd(), "results")

# dataset and file locations
datasetDir <- file.path(dataDir, "DBLvsHealth")
idatDir <- file.path(datasetDir, "idat")
sampleSheet <- file.path(datasetDir, "sample_annotation_DLBonly.age.gender.NeuN.csv")
reportDir <- file.path(resultDir, "DLBvsControl_age_gender_NeuN_sva")
################################################################################
# (1) Set analysis options
################################################################################
rnb.options(
	filtering.sex.chromosomes.removal = TRUE,
	identifiers.column                = "Sample_Name"
)
# optionally disable some parts of the analysis to reduce runtime
rnb.options(
	exploratory = TRUE,
	exploratory.correlation.qc        = FALSE,
	exploratory.intersample           = TRUE,
	exploratory.region.profiles       = c("genes","promoters"),
	exploratory.gene.symbols = c("APOE","MAOA","MAOB","DRD2","ANK1","MAPT","PRKN","SNCA","NBL1"),
	exploratory.clustering            = "top",
	exploratory.clustering.top.sites  = 1000,
	region.types                      = c("tiling"),
	filtering.cross.reactive = TRUE,
	normalization.background.method = "methylumi.noob",
	inference = TRUE,
	inference.age.prediction = TRUE,
	inference.age.column = "Age",
	inference.targets.sva = c("Sample_Group"),
	inference.sva.num.method ="be",
	differential.report.sites         = TRUE,
	differential.comparison.columns   = c("Sample_Group"), 
	export.to.csv = TRUE,
	export.to.bed = FALSE,
	export.to.trackhub = NULL,
	differential.adjustment.celltype = FALSE,
	differential.adjustment.sva = TRUE,
	min.group.size = 10,
	covariate.adjustment.columns = c("Age","Gender","NeuN_neg","NeuN_pos"),
	differential.variability = FALSE,
	differential.enrichment.go = TRUE,
	differential.enrichment.lola = FALSE
)

num.cores <- 30
parallel.setup(num.cores)

################################################################################
# (2) Run the analysis
################################################################################
rnb.run.analysis(
	dir.reports=reportDir,
	sample.sheet=sampleSheet,
	data.dir=idatDir,
	data.type="idat.dir",
	initialize.reports = TRUE
)

parallel.teardown()

################################################################################
# Link to finished analysis
################################################################################
# see the results at:
# http:???
