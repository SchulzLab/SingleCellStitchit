# SingleCellStitchit
Learning enhancer-gene interactions from single cell data

Recently, sequencing of RNA or DNA nucleotides with single-cell resolution has become possible, allowing to capture tissue heterogeneities as well as transcriptional dynamics.
To understand gene regulation mechanisms, it is essential to fathom out the role of enhancers, distal regulatory elements (REMs) that regulate one or several genes, in this process.
Integrative analysis of single-cell epigenetic and transcriptomic data can be used to gain insights into distal gene-expression regulation in specific phenotypes.
For this purpose, we introduce the sc-STITCHIT pipeline that is an extension of the STITCHIT algorithm.
STITCHIT utilizes both epigenetic and transciptomic information for the identification of regulatory elements of a specific target gene.
Our pipeline enables the application of STITCHIT to single-cell ATAC-seq data (chromatin accessibility) and single-cell RNA-seq data (transcriptome).
Since extremely small amounts of, for instance, mRNA molecules are present in an individual cell, single-cell sequencing experiments inherently generate highly sparse data.
Hence, extracting biologically meaningful information from single-cell sequencing data proves difficult.
This problem is addressed by summarizing the epigenetic and transcriptomic signal of individual cells into so-called metacells based on the similarity of their gene activity measurements.

## Installation STITCHIT: https://github.com/SchulzLab/STITCHIT

To build the STITCHIT cmake, a C++11 compiler, and the boost library must be available on your system.
We have tested STITCHIT with Linux operating systems.
Make sure you are in the projects main directory. Then, generate a build folder by typing: _mkdir build_
Change into that folder: _cd build_
To run cmake use: _cmake .._
To finally build the project use: _make_
You can speed up the build process by using the -j option, specifying how many cores should be used for building
To execute all tests use: make test
Tests can also be executed individually. They are located in the folder build/test

## Installation bam_merge-master:

To build the STITCHIT cmake, a C++17 compiler, and the boost library must be available on your system.
We have tested STITCHIT with Linux operating systems.
Make sure you are in the projects main directory. Then, generate a build folder by typing: _mkdir build_
Change into that folder: _cd build_
To run cmake use: _cmake .._
To finally build the project use: _make_

**Requirements: deeptools 3.5.1
                nextflow 21.04.3**

## Pipeline execution:
                nextflow main.nf -c nextlflow.config

## Set parameters in nextflow.config:

#### MetaCellaR:

MetaCellaR is a cell summarization strategy, in where similar cells (similarity based on their gene expression measurements)
are aggregated together to tackle the sparsity problem of single-cell data.

    params.singlecells = <path to Seurat object or csv file holding integrated scRNA-seq and scATAC-seq data>
    params.rna = <slot where the RNA counts are stored in the Seurat object, e.g. 'assays$RNA@data>
    params.atac = <slot where the RNA counts are stored in the Seurat object, 'assays$peaks@counts'>
    params.ctype = <csv file: 1st col: cell names provided in the csv file, and 2nd column the cell types, Seurat: slot cell type annotations for the single cell gene                        expression data>
    params.clusters = <number of clusters for the k-medoids clustering, or null for default: (k <- #cells in a cell type / 30)>
    params.mout = <directory name output folder>
    params.assay = <'meta.data$identity'>
    params.umapflag = <Project the RNA expression counts to a UMAP space? 'T' or 'F'>
    params.mccells = <e = expected number of cells per cluster, default = 30>
    params.treshold = <ATAC capping threshold, null -> default = 3*e>
    params.umapdims = <number of dimensions for UMAP reduction, default = 20> // (Default)
    params.reductions = <the cell embedding of the integration of interest , default = 'umap'>

For more detailed documentation see https://github.com/SchulzLab/MetaCellaR.

#### Gene Expression Discretization:

Gene expression discretization using Gaussian mixture modeling on the estimated signal density curve for to determine the expression/inexpression thresholds.

    params.vector = <either txt file holding gene (1st col: gene symbol, 2nd col: gene length) lengths or 'null' for automated gene length query if RNA counts should                           be normalized -> gene lengths are needed for TPM normalization>
    params.normalize = <should the RNA counts be TPM and DESeq2 normalized? true or false>
    params.genome = <choose genome for automated gene length query: 'h' for human, 'mm' for mouse>
    params.plot = <plot signal density curve and expression/inexpression threshold? true or false>
    params.binary = <binary classification? true; three classes? false >
    params.dout = <name file discretized gene expression matrix>
    params.tpmout = <name file normalized gene expression matrix>
    params.glout = <name file gene length vector>

Detailed explanations can be found on https://github.com/LrRmpf/fastged.

#### Bam Merge Master:

Generate one BAM file for each metacell holding all corresponding reads of the aggregated individual cells (MetaCellaR).

    Flag options: -X (10x), -D (non 10x ATAC)
    Start position of cell barcode in read name for non-10x data required
    params.bamfile = <Corresponding BAM file to Seurat object or csv file holding integrated scRNA-seq and scATAC-seq data>
    params.bamflag = '-X'
    params.pos = null
    params.outdir = <name output directory for metacell BAM files>

#### DeepTools bamcoverage: 

Generate bigWig files

    env.bamcoverageFolder = <name output folder>
    params.norm = <ATAC normalization, e.g. 'RPKM'>
    params.binsize = < bin size coverage track, Default 50>

All parameters are listed on https://deeptools.readthedocs.io/en/develop/content/tools/bamCoverage.html

####  STITCHIT
##### STITCHIT segmentation parameters

    STITCHIT segmentation parameters:  
    params.bw = <same path as bamcoverageFolder!!!> (Path to big wig files holding epigenetic signal that should be used for the segmentation)
    params.gtf = <path genome annotation file used to find the genomic position of the target gene (.gtf)>
    params.chromsizes = <txt file holding chromosome sizes of target organism>
    params.window = <extension up and downstream of the target gene, default=5000>
    params.cores = <number of cores used within STITCHIT>
    params.pvalue = <p-value threshold used to select a regulatory element, default=0.05>
    params.resolution = <resolution used to merge the initial data to cut runtime, default=10>
    params.stitchitout = <name output directory>
    params.searchspace = <maximum size of the entire considered search space, default=1100000>
    params.segmentsize = <size of the segments used within STITCHIT, default=2000>

##### Nested STITCHIT:

    params.folds = <number folds Monte-Carlo Cross Validation>

##### STITCHIT regression parameters

    params.response= <name of the response variable, "Expression"  >
    params.fixedAlpha = <fixed value for the alpha parameter in elastic net regulatisation, do not perform a grid search; grid search: null>
    params.alpha = <stepsize to optimisation alpha parameter in elastic net regularisation, default = 0.05>
    params.testsize =<proportion test data, default=0.2>
    params.regularization = <regularisation 'L' for Lasso, 'R' for Ridge, and 'E' for Elastic net (default 'E')>
    params.innerCV = <umber of folds for inner cross-validation, default=6>
    params.outerCV = <number of iterations of outer cross-validation to determine test error, default=3>
    params.constraint = <constraint on the coefficent sign, enter 'N' for negative and 'P' for positive constraint>
    params.performance = <flag indiciating whether the performance of the model should be assessed, default='TRUE'>
    params.seed = <random seed used for random number generation, default = null (random)>
    params.leaveOneOutCV = <flag indicating whether a leave one out cross-validation should be used, default='FALSE'>
    params.asRData = <store feature coefficients as RData files, default='FALSE'>
    params.randomise = <randomise the feature matrix,  default='FALSE'>
    params.logResponse = <flag indicating whether the response variable should be log transformed, default='TRUE'>
    params.ftest = <flag indicating whether partial F-test should be computed to assess the significance of each feature, default='FALSE'>
    params.coefP = <p-value threshold for model coefficient, default= 1 -> all OLS coefs will be returned)>

A detailed description of the parameters can be found here: https://github.com/SchulzLab/STITCHIT.

##### Flag STITCHIT: nested or normal version:
env.mode = <Default: 'normal', Options: ['normal', 'nested']>

