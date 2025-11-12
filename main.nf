process WRANGLE_META {
	conda "/home/rschwartz/anaconda3/envs/r4.3"

	publishDir "${params.outdir}/meta", mode: 'copy'


	input:
		path(dea_meta)
	output:
		path "**.tsv", emit: wrangled_meta_files
	script:
		"""
		Rscript $projectDir/bin/wrangle_meta.R --dea_meta ${dea_meta}	
		"""
}

process WRANGLE_DEA {
	conda "/home/rschwartz/anaconda3/envs/r4.3"
	publishDir "${params.outdir}/dea_results", mode: 'copy'
	
	input:
		path dea_results
	output:
		path "**.tsv", emit: wrangled_dea_files
	script:
		"""
		Rscript $projectDir/bin/wrangle_dea.R --dea_results ${dea_results}	
		"""
}

process MAP_CONTRAST_META {
	conda "/home/rschwartz/anaconda3/envs/scanpyenv"
	
	input:
		tuple val(experiment), val(contrast), path(dea_results_file), path(dea_meta_file)
	output:
		path "${experiment}_mapped_dea_results.tsv"
	script:
		"""
		"""
}


workflow {

	dea_meta = Channel.fromPath(params.dea_meta)
	dea_results = Channel.fromPath(params.dea_results)

	WRANGLE_META(dea_meta)
	WRANGLE_DEA(dea_results)

	dea_results_files = WRANGLE_DEA.out.wrangled_dea_files
	dea_meta_files = WRANGLE_META.out.wrangled_meta_files

	// gt experiment from filenames

	dea_results_files.flatten()
		.map { file -> 
			def fullname = file.getName()
			def experiment = fullname.split("_")[0]
			def contrast = fullname.split("_")[1]
			[experiment, contrast, file]
		}
	.set { dea_results_ch}

	dea_meta_files.flatten()
		.map { file -> 
			def fullname = file.getName()
			def experiment = fullname.split("_")[0]
			[experiment, file]
		}
	
	.set { dea_meta_ch }


	// group tuples by experiment
	dea_results_ch
		.groupTuple(by: [0,1])
		.set { grouped_contrast_ch }

	grouped_contrast_ch.combine(dea_meta_ch, by: 0)
		.set { final_grouped_ch }

	final_grouped_ch.view()




}