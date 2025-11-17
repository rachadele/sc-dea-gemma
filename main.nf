process LOAD_PREFERRED_CTA {
    tag "$experiment"

     input:
        val experiment

    output :
		val experiment, emit: cta_done
        path "message.txt"

	script:
	def experiment_name = experiment.tokenize('.')[0]
	def cta_file = "${params.cta_file_prefix}/${experiment_name}/${experiment}_${params.level}_cell_type.tsv"

	"""
	gemma-cli loadSingleCellData -loadCta -e ${experiment} \\
				-replaceCta \\
			   -ctaFile ${cta_file} -preferredCta \\
			   -ctaName "sc-pipeline-${params.version}-${params.level}" \\
			   -ctaProtocol "sc-pipeline-${params.version}" \\
			   --data-type NULL \\
			   -ignoreSamplesLackingData 2> "message.txt"
	"""
}



process AGGREGATE_DATA {
	tag "$experiment"

	input:
	val experiment

	output:
	path "message.txt"
	val experiment, emit: agg_done


	script:
	"""
	gemma-cli aggregateSingleCellData --make-preferred -e ${experiment} 2> "message.txt"
	"""
}


process RUN_DEA {
	tag "$experiment"

	input:
	val experiment

	output:
	path "message.txt"
	val experiment, emit: dea_done

	script:
	factors = params.experiment_factors[experiment].join(",")
	"""
	gemma-cli diffExAnalyze -subset cell_type -e ${experiment} \\
							--ignore-failing-subsets \\
							--factors ${factors} \\
							--no-files \\
							 2> "message.txt"
	"""

}


process GET_DEA_RESULTS {
	tag "$experiment"
	conda "/home/rschwartz/anaconda3/envs/gemmapy"
	publishDir "${params.outdir}/dea_results_files/${experiment}", mode: 'copy'

	input:
	val experiment

	output:
	tuple val(experiment), path("**tsv")

	script:
	"""
	python3 $projectDir/bin/get_dea_results.py \\
		--experiment ${experiment} \\
		--username ${params.GEMMA_USERNAME} \\
		--password ${params.GEMMA_PASSWORD}
	"""

}

process GET_DEA_META {
	tag "$experiment"
	conda "/home/rschwartz/anaconda3/envs/gemmapy"
	publishDir "${params.outdir}/meta_files", mode: 'copy'

	input:
	val experiment

	output:
	tuple val(experiment), path("**tsv")

	script:
	"""
	python $projectDir/bin/get_dea_meta.py \\
		--experiment ${experiment} \\
		--username ${params.GEMMA_USERNAME} \\
		--password ${params.GEMMA_PASSWORD}
	"""

}

process MAP_DEA_META {
	tag "$experiment"
	conda "/home/rschwartz/anaconda3/envs/scanpyenv"
	publishDir "${params.outdir}/dea_meta_mapped/${experiment}", mode: 'copy'

	input:
	tuple val(experiment), val(result_id), path(dea_results), path(dea_meta)

	output:
	tuple val(experiment), val(result_id), path("**tsv"), emit: mapped_contrasts

	script:
	"""
	python3 $projectDir/bin/map_contrast_meta.py \\
		--experiment ${experiment} \\
		--dea_results ${dea_results} \\
		--dea_meta ${dea_meta} \\
		--result_id ${result_id}
	"""
}


process PLOT_PVAL_DISTS {
    tag "$experiment"
    conda "/home/rschwartz/anaconda3/envs/scanpyenv"
    publishDir "${params.outdir}/pvalue_dists/${experiment}/${result_id}/", mode: 'copy'

    input:
    tuple val(experiment), val(result_id), path(dea_results_mapped)

    output:
    path "**.png"

    script:
    """
    python $projectDir/bin/plot_pval_dists.py \
        --dea_results_mapped ${dea_results_mapped}
    """
}


workflow {

	experiments = Channel.fromList(params.experiments)

	// run preferred cta loading
	LOAD_PREFERRED_CTA(experiments)
	
	

	// aggregate data
	AGGREGATE_DATA(LOAD_PREFERRED_CTA.out.cta_done)


	// run dea
	RUN_DEA(AGGREGATE_DATA.out.agg_done)
	

	//only run these after DEA is done
	dea_results_files = GET_DEA_RESULTS(RUN_DEA.out.dea_done)
	dea_meta_files = GET_DEA_META(RUN_DEA.out.dea_done)


	dea_results_files.flatMap { experiment, dea_results ->
		dea_results
		.collect { dea_result_file ->
			def result_id = dea_result_file.getName().split("_")[1].split("\\.")[0]
			[experiment, result_id, dea_result_file]
		}
	}.set { dea_results_ch }

	dea_results_ch.combine(dea_meta_files, by: 0)
	.set { dea_meta_combined_ch }

	

	MAP_DEA_META(dea_meta_combined_ch)
	// view mapped contrasts
	//MAP_DEA_META.out.mapped_contrasts.view()

	MAP_DEA_META.out.mapped_contrasts.flatMap { experiment, result_id, mapped_tsv_list ->
		// mapped_tsv_list may not be a list if there's only one contrast
		// check if no_res in filename
		if ( mapped_tsv_list instanceof Path ) {
			mapped_tsv_list = [ mapped_tsv_list ]
		}
		else {
			mapped_tsv_list = mapped_tsv_list
		}
		// Only keep mapped_tsvs that do NOT have 'no_res' in the filename
		mapped_tsv_list.findAll { mapped_tsv ->
			!mapped_tsv.getName().contains('no_res')
		}.collect { mapped_tsv ->
			[experiment, result_id, mapped_tsv]
		}
	}.set { mapped_tsvs }
	//view
	//mapped_tsvs.view()
		


   	PLOT_PVAL_DISTS(mapped_tsvs)
}