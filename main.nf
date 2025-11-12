

process LOAD_PREFERRED_CTA {
    tag "$experiment"

     input:
        tuple val(experiment)

    output :
        path "message.txt"

    script:
	cta_file = "${params.cta_file_prefix}/${experiment}_${level}_cell_type.tsv"

    """
   	gemma-cli loadSingleCellData -loadCta -e ${experiment} \\
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


	script:
	"""
	gemma-cli aggregateSingleCellData --make-preferred -e ${experiment}
	"""
}



process GET_DEA {
	tag "$experiment"
	conda "/home/rschwartz/anaconda3/envs/r4.3"
	publishDir "${params.outdir}/dea_results_file/$", mode: 'copy'

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

	dea_results_files = GET_DEA(experiments)
	dea_meta_files = GET_DEA_META(experiments)


	// gt experiment from filenames

	//dea_results_files.flatten()
		//.map { file -> 
			//def fullname = file.getName()
			//def experiment = fullname.split("_")[0]
			//def contrast = fullname.split("_")[1]
			//[experiment, contrast, file]
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