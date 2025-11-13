

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
	gemma-cli loadSingleCellData -loadCta -e ${experiment} \
			   -ctaFile ${cta_file} -preferredCta \
			   -ctaName "sc-pipeline-${params.version}-${params.level}" \
			   -ctaProtocol "sc-pipeline-${params.version}" \
			   --data-type NULL \
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
	"message.txt"
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
	conda "/home/rschwartz/anaconda3/envs/r4.3"
	publishDir "${params.outdir}/dea_results_files/${experiment}", mode: 'copy'

	input:
	val experiment

	output:
	tuple val(experiment), path("**tsv")

	script:
	"""
	Rscript $projectDir/bin/get_dea_results.R \\
		--experiment ${experiment} \\
		--username ${params.GEMMA_USERNAME} \\
		--password ${params.GEMMA_PASSWORD}
	"""

}

process GET_DEA_META {
	tag "$experiment"
	conda "/home/rschwartz/anaconda3/envs/r4.3"
	publishDir "${params.outdir}/meta_files", mode: 'copy'

	input:
	val experiment

	output:
	tuple val(experiment), path("**tsv")

	script:
	"""
	Rscript $projectDir/bin/get_dea_meta.R \\
		--experiment ${experiment} \\
		--username ${params.GEMMA_USERNAME} \\
		--password ${params.GEMMA_PASSWORD}
	"""

}


workflow {

	experiments = Channel.fromList(params.experiments)

	// run preferred cta loading
	LOAD_PREFERRED_CTA(experiments)
	
	

	// aggregate data
	AGGREGATE_DATA(	LOAD_PREFERRED_CTA.out.cta_done)


	// run dea
	RUN_DEA(AGGREGATE_DATA.out.agg_done)
	

	//only run these after DEA is done
	dea_results_files = GET_DEA_RESULTS(RUN_DEA.out.dea_done)
	dea_meta_files = GET_DEA_META(RUN_DEA.out.dea_done)


	// gt experiment from filenames

	//dea_results_files.flatten()
		//.map { file -> 
			//def fullname = file.getName()
			//def experiment = fullname.split("_")[0]
			//def resultSet = fullname.split("_")[1]
			//[experiment, resultSet, file]
		//}
	//.set { dea_results_ch}

	//dea_meta_files.flatten()
		//.map { file -> 
			//def fullname = file.getName()
			//def experiment = fullname.split("_")[0]
			//[experiment, file]
		//}
	
	//.set { dea_meta_ch }


	//// group tuples by experiment
	//dea_results_ch
		//.groupTuple(by: [0,1])
		//.set { grouped_contrast_ch }

	//grouped_contrast_ch.combine(dea_meta_ch, by: 0)
		//.set { final_grouped_ch }

	//final_grouped_ch.view()




}