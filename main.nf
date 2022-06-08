#! /usr/bin/env nextflow

/*
 * Author: Laura Rumpf
 * 
 * Pipeline Enhancer-Gene Prediction based on single-cell data
 *
 */

println "$metacellarFolder"
f1 = file("$metacellarFolder")
f1.mkdirs()

f2 = file("$discretizationFolder")
f2.mkdirs()

f3 = file("$bam_mergeFolder")
f3.mkdirs()

f4 = file("$bamcoverageFolder")
f4.mkdirs()

f5 = file("$StitchItFolder")
f5.mkdirs()

f6 = file("$RegressionFolder")


// process definitions

process metacellar {

    publishDir metacellarFolder, mode: 'copy'

    input:
      path sc from params.singlecells

    output:
      path "${params.mout}/cellSummarized_*_sum.csv" into counts_ch //contains sc-RNA and sc-ATAC summarized counts into preprocessing step to use only metacells that contain both sc-atac and sc-rna 
      path "${params.mout}/RNA_cell2metacell_info_*.csv" // into rna_mc_ch 
      path "${params.mout}/ATAC_cell2metacell_info_*.csv" into atac_ch //input preprocessing step
      path "${params.mout}/*_clustered.RData"
      path "${params.mout}/cellSummarized_*_means.csv"
      path "${params.mout}/*.pdf"
      path 'plot_umap_seurat_*.pdf'		      
	
    script:

    if ( ! (params.clusters || params.treshold) ) {
	"""
	Rscript '${projectDir}/bin/MetaCellaR.R' -file $sc -RNA '$params.rna' -celltype '$params.ctype' -output '$params.mout' -umap '$params.umapflag' -assay '$params.assay' -ATAC '$params.atac' -e '$params.mccells' -d '$params.umapdims' -reduction '$params.reductions'
  	"""
    }
	
    else  if (params.clusters && !params.treshold ) {
	"""
	Rscript '${projectDir}/bin/MetaCellaR.R' -file $sc -RNA '$params.rna' -celltype '$params.ctype' -output '$params.mout' -k '$params.clusters' -umap '$params.umapflag' -assay '$params.assay' -ATAC '$params.atac' -e '$params.mccells' -d '$params.umapdims'  -reduction '$params.reductions'
        """
	}
    else  if (!params.clusters && params.treshold ) {
        """
        Rscript '${projectDir}/bin/MetaCellaR.R' -file $sc -RNA '$params.rna' -celltype '$params.ctype' -output '$params.mout' -umap '$params.umapflag' -assay '$params.assay' -ATAC '$params.atac' -e '$params.mccells' -t '$params.treshold' -d '$params.umapdims'  -reduction '$params.reductions'
        """	
	}
    else {
	"""
	Rscript '${projectDir}/bin/MetaCellaR.R' -file $sc -RNA '$params.rna' -celltype '$params.ctype' -output '$params.mout' -k '$params.clusters' -umap '$params.umapflag' -assay '$params.assay' -ATAC '$params.atac' -e '$params.mccells' -t '$params.treshold' -d '$params.umapdims'  -reduction '$params.reductions'
	"""
	}	

}

atac_ch.view()
counts_ch
	.collect()
	.view()


process treshold_metacell_counts {

        publishDir "${metacellarFolder}/${params.mout}", mode: 'copy'
	echo true

        input:
        val counts from counts_ch //sc-rna and sc-atac count files as list 				.
	path atac_cell2metacell from atac_ch
       
        output:
        path "cellSummarized_*_sum_RNA_final_metacells.csv" into rna_counts_ch //input discretization
	path "ATAC_cell2metacell_info_*_final_metacells.csv"into atac_counts_ch //input bam_merge
	path "cellSummarized_*_sum_ATAC_final_metacells.csv"
	 
        script:
	
        c1 = counts[0]
			.getSimpleName()
			.toString()
	println c1
		
	if (c1.contains("ATAC")) {
		atac_counts = counts[0]
		rna_counts = counts[1]
	}
	else {
		atac_counts = counts[1]
                rna_counts = counts[0]
	}
	
	
        atac_name = atac_cell2metacell.getSimpleName()
	println atac_name

	rna_name = rna_counts.getSimpleName()
	atac_counts_name = atac_counts.getSimpleName()
 

        """
        #!/usr/bin/env Rscript

        require(tidyverse)
	
	# count data -> rna counts input discretization
        rna.counts <- read.csv("${rna_counts}", header = TRUE, check.names=FALSE)
        atac.counts <- read.csv("${atac_counts}", header = TRUE, check.names=FALSE)
	
	# atac cell2metacell info -> input bam_merge
	atac.cell2metacell <- read.csv("${atac_cell2metacell}", header = TRUE, check.names=FALSE)
	
	# consistent naming
	atac.cell2metacell[,1] <- sub("_2", "", atac.cell2metacell[,1])
	#atac.cell2metacell[,1] <- sub("_1", "", atac.cell2metacell[,1])
	#atac.cell2metacell[,1] <- sub("atac_", "", atac.cell2metacell[,1])
	#atac.cell2metacell[,2] <- gsub(" ", "_", atac.cell2metacell[,2])
	#atac.cell2metacell[,2] <- gsub("-", "_", atac.cell2metacell[,2])
	
	 if (is.character(rna.counts[1,1])){
                rownames(rna.counts) <- rna.counts[,1]
                rna.counts <- rna.counts[,-1]
        }

	 if (is.character(atac.counts[1,1])){
                rownames(atac.counts) <- atac.counts[,1]
                atac.counts <- atac.counts[,-1]
        }
	
	# treshold number of reads in atac metacells = 200000
	print(ncol(atac.counts))
	atac.counts <- atac.counts[,colSums(atac.counts) > 200000]
	print(ncol(atac.counts))
	
	 #colnames(rna.counts) <- gsub("[.]", "_", colnames(rna.counts))
	 #colnames(rna.counts) <- gsub(" ", "_", colnames(rna.counts))
	 #colnames(rna.counts) <- gsub("-", "_", colnames(rna.counts))

        # extract unique metacell_ids from sc-rna counts = column names
        rna.mc.ids <- colnames(rna.counts)
	print("rna:")
	print(rna.mc.ids)
	
	# extract unique metacell_ids from sc-atac counts = column names -> filtered
	atac.mc.ids <- colnames(atac.counts)
	#atac.mc.ids <- unique(atac.cell2metacell[,2])
	print("atac:")
	print(atac.mc.ids)
	
	print("atac_cell2metacell: ")
	print(unique(atac.cell2metacell[,2]))
	
	# union	
	mc.ids <- unique(c(atac.mc.ids,rna.mc.ids))
	
	# intersection: metacells containing sc-rna and sc-atac 
	mc.ids <- intersect(atac.mc.ids,rna.mc.ids)
	
	# MetaCellaR -> ATAC cells are assigned to closest RNA metacells
	# Use ATAC metacells would be enough
	# In case MetaCellaR algo changes: use intersection 
 
	print("mc:")
	print(mc.ids)

	# valid atac counts 
        atac.counts <- atac.counts[,colnames(atac.counts) %in% mc.ids]
	print("valid atac metacells:")
        print(colnames(atac.counts))
	
	# keep only samples (metacells) that also contain atac counts -> rna counts for discretization
	rna.counts <- rna.counts[,colnames(rna.counts) %in% mc.ids]
	print("valid rna metacells:")
	print(colnames(rna.counts))

        # keep only rows that contain mc ids containing sc-rna and sc-atac and more than 200000 atac reads -> bam_merge input
        atac.cell2metacell <- atac.cell2metacell[atac.cell2metacell[,2] %in% mc.ids,]
	print("valid cell2mc atac:")
	print(unique(atac.cell2metacell[,2]))

        write.csv(atac.cell2metacell,"${atac_name}_final_metacells.csv", row.names=FALSE)
	write.csv(rna.counts, "${rna_name}_RNA_final_metacells.csv", row.names=TRUE)
	write.csv(atac.counts, "${atac_counts_name}_ATAC_final_metacells.csv", row.names=TRUE)

        """


}

rna_counts_ch.view()
atac_counts_ch.view()


process discretization {

    publishDir discretizationFolder, mode: 'copy'

    input:
      path counts from rna_counts_ch

    output:
      path params.dout into discrete_rna_ch1, discrete_rna_ch2 //input StitchIt
      path params.tpmout into ensembl_counts_ch
      path params.glout
      path "Plots/Sample_*.png" optional true
      path "deleted_genes.log" optional true

    script:
	
        if(!(params.normalize || params.plot || params.binary)){
        """
        Rscript '${projectDir}/bin/discretization.R' -f $counts -o '$params.dout' 

        """
	}	
	else if(params.normalize && !params.plot && !params.binary){
		if(params.vector != null){
    		"""
    		Rscript '${projectDir}/bin/discretization.R' -f $counts -n -g '$params.genome' -v '$params.vector' -o '$params.dout' --tpmout '$params.tpmout' --lout '$params.glout' 
    
   		"""
		}
		else{
		 """
                Rscript '${projectDir}/bin/discretization.R' -f $counts -n -g '$params.genome' -o '$params.dout' --tpmout '$params.tpmout' --lout '$params.glout'

                """
		}
    	} 
	else if(!params.normalize && params.plot && !params.binary){
        """
        Rscript '${projectDir}/bin/discretization.R' -f $counts -p -o '$params.dout'  

        """
        }
 	else if(!params.normalize && !params.plot && params.binary){
        """
        Rscript '${projectDir}/bin/discretization.R' -f $counts -b -o '$params.dout'

        """
        } 
        else if(!params.normalize && params.plot && params.binary){
        """
        Rscript '${projectDir}/bin/discretization.R' -f $counts -p -b -o '$params.dout'

        """
        }
 	else if(params.normalize && params.plot && params.binary){

		if(params.vector != null){
       		"""
        	Rscript '${projectDir}/bin/discretization.R' -f $counts -n -g '$params.genome' -v '$params.vector' -p -b -o '$params.dout' --tpmout '$params.tpmout' --lout '$params.glout'

        	"""
		} 
		else {
		"""
                Rscript '${projectDir}/bin/discretization.R' -f $counts -n -g '$params.genome' -p -b -o '$params.dout' --tpmout '$params.tpmout' --lout '$params.glout'

                """	
		}
	}

 	else if(params.normalize && !params.plot && params.binary){
		if(params.vector != null){
        	"""
        	Rscript '${projectDir}/bin/discretization.R' -f $counts -n -g '$params.genome' -v '$params.vector' -b -o '$params.dout' --tpmout '$params.tpmout' --lout '$params.glout'

        	"""
		} 
		else {
		"""
		Rscript '${projectDir}/bin/discretization.R' -f $counts -n -g '$params.genome' -b -o '$params.dout' --tpmout '$params.tpmout' --lout '$params.glout'

                """
		}
        }
	else if(params.normalize && params.plot && !params.binary){
		if(params.vector != null){
        	"""
        	Rscript '${projectDir}/bin/discretization.R' -f $counts -n -g '$params.genome' -v '$params.vector' -p -o '$params.dout' --tpmout '$params.tpmout' --lout '$params.glout'

        	"""
		} 
		else {
		"""
                Rscript '${projectDir}/bin/discretization.R' -f $counts -n -g '$params.genome' -p -o '$params.dout' --tpmout '$params.tpmout' --lout '$params.glout'

                """

		}
        }

}


process bam_merge {

    publishDir bam_mergeFolder, mode: 'copy' 

    input:
      path atac_cell2metacell from atac_counts_ch
      path bam from params.bamfile

    output:
      path "${params.outdir}/output_*.bam" into merged_atac_ch //input Deeptools

    script:
    if(params.bamflag == '-D')
    	"""
    	"${bamMergeDir}" $params.bamflag $bam $atac_cell2metacell $params.pos $params.outdir
    	"""
    else 
    	"""
	 "${bamMergeDir}" $params.bamflag $bam $atac_cell2metacell $params.outdir
    	"""
     
}




process preprocessing {

	publishDir "${bam_mergeFolder}/${params.outdir}"

	input:
        path bam from merged_atac_ch.flatten()
	
	output:
        tuple file("output_*_sorted.bam"), file("output_*_sorted.bam.bai") into index_ch //input bam_coverage
	// path "output_*_sorted.bam" into sorted_ch 

	script:

	filename = bam.getSimpleName()
        println filename

        """
	samtools sort -o ${filename}_sorted.bam $bam
        samtools index ${filename}_sorted.bam
        """

}


process bam_coverage {
				  
    publishDir bamcoverageFolder, mode: 'copy'
    // publishDir "${bam_mergeFolder}/${params.outdir}"

    input:
    // path bam_file from sort_ch2.flatten()
    // path bam_index from index_ch
    // tuple val(mc_id), file(indexed_bam) from paired_bam_ch
    tuple sorted_bam, indexed_bam from index_ch
	      

    output: 
      path "*.bw" into bigwig_ch //input StitchIt
      // path bamcoverageFolder into bigwig_ch
	
    script:

     /* index_sorted_bam_ch = Channel
                                .fromFilePairs("${bam_mergeFolder}/${params.outdir}/output_*_sorted{.bam,.bam.bai}")
                                .toList()*/
       
         
    //fname = bam_file.getSimpleName()
    // def bam_file = indexed_bam.first()
    mc_id = sorted_bam.getSimpleName() - "output_" - "_sorted" 
     

    """
    bamCoverage -b $sorted_bam -o ${mc_id}.bw --normalizeUsing $params.norm --binSize $params.binsize
    """
}


//def geneList = [];
		
process get_genes {

	input:
	   path dge from discrete_rna_ch1
					 

	output:
	   val geneList into genes_ch
	
	 script:
	   file_path = discrete_rna_ch1
                        .evaluate()
 
	   geneList = [];
  	   String line
	   f = file(file_path)
	   //f = file(discrete_rna_ch1.view())
	   //f = file(dge.toString())
           //f = file("${discretizationFolder}/dge.txt")	    
           int lineCount = 0
	   String gene
           myReader = f.newReader()
                while ((line = myReader.readLine()) != null) {
                //skip header
                        if (lineCount==0) {
                                lineCount++
                                continue
                        }

                        String[] words = line.split("\t");
                        gene = words[0]
                        geneList << gene
                        lineCount++;
                }

	    myReader.close() 
 	     
	   """
	   echo $geneList  	
	   """  	    
}

genes_ch.view()

//Default mode == normal
process stitch_it {
	
        //publishDir StitchItFolder
	echo true
	errorStrategy 'ignore'
	
	input:
	   path bigwig from bigwig_ch.collect() 
	   path dge from discrete_rna_ch2
	   path counts from ensembl_counts_ch
	   each item from genes_ch.flatten()

	output:
	   //optional true
	   //path "Segmentation_*.txt" into segmentation_ch
	   val item into segmentation_ch
				  
	script:
	
	if ("$mode" == "normal"){
		"""
		"${stitchDir}"  -b '$bamcoverageFolder' -a '$params.gtf' -d '$dge' -o '$counts' -s '$params.chromsizes' -w $params.window -c $params.cores -p $params.pvalue -g '$item' -f '$StitchItFolder' -z $params.resolution -r $params.searchspace -t $params.segmentsize	
		"""
	}
        else if (mode == 'nested') {
                """
                python ${projectDir}/bin/NestedSTITCHIT.py '$bamcoverageFolder' '$params.gtf' '$dge' '$counts' '$StitchItFolder' '$params.chromsizes' '$item' 'params.folds' '$RegressionFolder'
                """
	}
        else {
                 error "Invalid regression mode: ${mode}"
	}
		
	

}


process regression {

	echo true
		
	input:	
	val gene_id from segmentation_ch
					.collect()
					.view()


when:
	mode == 'normal'
	
	script:
		
	if ( ! (params.fixedAlpha || params.seed) ) {
		"""
		Rscript '${projectDir}/bin/Two_level_learning.R' --dataDir="${StitchItFolder}" --outDir="${RegressionFolder}" --response=$params.response --cores=$params.cores --alpha=$params.alpha --testsize=$params.testsize --regularisation=$params.regularization --innerCV=$params.innerCV --outerCV=$params.outerCV --constraint=$params.constraint --performance=$params.performance --leaveOneOutCV=$params.leaveOneOutCV --asRData=$params.asRData --randomise=$params.randomise --logResponse=$params.logResponse --ftest=$params.ftest --coefP=$params.coefP 
		   
		"""
	}
	else if ( !params.fixedAlpha &&  params.seed ) {
		"""
		Rscript '${projectDir}/bin/Two_level_learning.R' --dataDir="${StitchItFolder}" --outDir="${RegressionFolder}" --response=$params.response --cores=$params.cores --alpha=$params.alpha --testsize=$params.testsize --regularisation=$params.regularization --innerCV=$params.innerCV --outerCV=$params.outerCV --constraint=$params.constraint --performance=$params.performance --seed=$params.seed --leaveOneOutCV=$params.leaveOneOutCV --asRData=$params.asRData --randomise=$params.randomise --logResponse=$params.logResponse --ftest=$params.ftest --coefP=$params.coefP	

		"""
	}
	else if (params.fixedAlpha &&  !params.seed) {
		"""
        	Rscript '${projectDir}/bin/Two_level_learning.R' --dataDir="${StitchItFolder}" --outDir="${RegressionFolder}" --response=$params.response --cores=$params.cores --fixedAlpha=$params.fixedAlpha --alpha=$params.alpha --testsize=$params.testsize --regularisation=$params.regularization --innerCV=$params.innerCV --outerCV=$params.outerCV --constraint=$params.constraint --performance=$params.performance --leaveOneOutCV=$params.leaveOneOutCV --asRData=$params.asRData --randomise=$params.randomise --logResponse=$params.logResponse --ftest=$params.ftest --coefP=$params.coefP

        	"""
	}
	else {
		"""
		 Rscript '${projectDir}/bin/Two_level_learning.R' --dataDir="${StitchItFolder}" --outDir="${RegressionFolder}" --response=$params.response --cores=$params.cores --fixedAlpha=$params.fixedAlpha --alpha=$params.alpha --testsize=$params.testsize --regularisation=$params.regularization --innerCV=$params.innerCV --outerCV=$params.outerCV --constraint=$params.constraint --performance=$params.performance --seed=$params.seed --leaveOneOutCV=$params.leaveOneOutCV --asRData=$params.asRData --randomise=$params.randomise --logResponse=$params.logResponse --ftest=$params.ftest --coefP=$params.coefP
		"""
	}   
    
}
