process split_matrix {
    label 'bigmem'

    publishDir "${params.outdir}/${prefix}", pattern: "${prefix}.binary_matrix.npy"
    publishDir "${params.outdir}/${prefix}", pattern: "${prefix}.dhs_masterlist.bed"

    input:
        tuple val(prefix), val(ag_ids)

    output:
        tuple val(prefix), path("${prefix}_ag_ids.txt"), path(binary_matrix_subset), path(dhs_masterlist_subset)

    script:
    ag_ids_str = ag_ids.join(",")
    binary_matrix_subset = "${prefix}.binary_matrix.npy"
    dhs_masterlist_subset = "${prefix}.dhs_masterlist.bed"
    """
    echo ${ag_ids_str} > ${prefix}_ag_ids.txt

    ${moduleDir}/bin/split_matrix.py ${params.samples_order} ${params.binary_matrix} ${params.dhs_masterlist} ${prefix}_ag_ids.txt

    mv prefix_binary_matrix.npy ${binary_matrix_subset}
    mv dhs_masterlist_subset.bed ${dhs_masterlist_subset}
    """
}

process permutations {
 
    errorStrategy 'ignore'
   
    input:
        tuple val(prefix), val(k), val(signal_percentile), path(ag_ids), path(binary_matrix_subset), path(dhs_masterlist_subset)
    
    output:
        tuple val(prefix), val(signal_percentile), val(k), path(name)
    
    script:
    name = "${prefix}.${k}.${signal_percentile}.Rdata"
    """
    Rscript ${moduleDir}/bin/code_num_new_DHSs_nf.R ${binary_matrix_subset} ${dhs_masterlist_subset} ${k} ${signal_percentile} ${name}
    """   
}


process transform_to_saturation_curve {

    publishDir "${params.outdir}/${prefix}/pdfs", pattern: "*saturation_curve.pdf"
    publishDir "${params.outdir}/${prefix}/added_DHSs_matrices", pattern: "${prefix}.${signal_percentile}.RData"

    input:
        tuple val(prefix), val(signal_percentile), val(k), path(rdata_files), path(binary_matrix_subset)
    
    output:
        tuple val(prefix), val(signal_percentile), val(k), path(saturation_curve), path(percent_data)
    
    script:
    saturation_curve = "${prefix}.${signal_percentile}.saturation_curve.pdf"
    percent_data = "${prefix}.${signal_percentile}.RData"
    """

    Rscript ${moduleDir}/bin/code_parse_num_new_DHSs_nf.R "${k}" ${signal_percentile} ${saturation_curve} ${percent_data} ${binary_matrix_subset} ${prefix}
    mv *pdf ${saturation_curve}

    echo "Done"   
    """
}

process quantile_curve {


    publishDir "${params.outdir}/${prefix}/pdfs", pattern: "*quantile_curves.pdf"


    input:
   	tuple val(prefix), path(percent_data), path(binary_matrix_subset)

    output:
	tuple val(prefix), path(quantile_curves)

    script:

    quantile_curves = "${prefix}.quantile_curves.pdf"

    """
	
    	Rscript ${moduleDir}/bin/code_parse_quantiles_nf.R  > ${quantile_curves}

    """
	
}

workflow {
    percentiles = Channel.of(0..9)
        | map(it -> it * 10)
    
    metadata = Channel.fromPath(params.samples_file)
        | splitCsv(header: true, sep: "\t")
        | map(row -> tuple(row.prefix, row.ag_id))
        | groupTuple() // prefix, ag_ids
    
    n_counts = metadata
	| flatMap { it -> 
            def prefix = it[0]
            def k_values = (2..(it[1].size() - 2)) // Define range from 2 to number of samples
            k_values.collect { k -> tuple(prefix, k) } // Create individual (prefix, k) pairs
        }
	
    matrices = metadata
        | split_matrix //prefix, ag_ids_file, binary_matrix_subset, dhs_masterlist_subset
    
    binary_matrix_subset = matrices.map { it[2] } // Extract the third element from each tuple in matrices

    output = n_counts
        | combine(percentiles) // Combine each row with each percentile
	| combine(matrices, by: 0)
	| permutations
        | groupTuple(by: [0,1])
        | combine(binary_matrix_subset)
	| transform_to_saturation_curve
        | groupTuple(by: [0])
        | map(it -> tuple(it[0], it[4]))
        | combine(binary_matrix_subset)
        | quantile_curve
}

//| groupTuple(by: [0, 1])
//| transform_to_plotting_data

