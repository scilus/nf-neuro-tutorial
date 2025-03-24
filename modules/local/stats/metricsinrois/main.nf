process STATS_METRICSINROI {
    tag "$meta.id"
    label 'process_single'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://scil.usherbrooke.ca/containers/scilus_latest.sif':
        'scilus/scilus:latest' }"

    input:
    tuple val(meta), path(t1), path(metrics)

    output:
    tuple val(meta), path("*.json")                         , emit: stats
    tuple val(meta), path("*mask_wm.nii.gz")                , emit: wm_mask
    tuple val(meta), path("*mask_gm.nii.gz")                , emit: gm_mask
    tuple val(meta), path("*mask_csf.nii.gz")               , emit: csf_mask
    tuple val(meta), path("*map_wm.nii.gz")                 , emit: wm_map
    tuple val(meta), path("*map_gm.nii.gz")                 , emit: gm_map
    tuple val(meta), path("*map_csf.nii.gz")                , emit: csf_map

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def suffix = task.ext.first_suffix ? "${task.ext.first_suffix}_stats" : "stats"
    def bin = task.ext.bin ? "--bin " : ""
    def normalize_weights = task.ext.normalize_weights ? "--normalize_weights " : ""
    """
    export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=1
    export OMP_NUM_THREADS=1
    export OPENBLAS_NUM_THREADS=1

    fast -t 1 -n 3\
        -H 0.1 -I 4 -l 20.0 -g -o t1.nii.gz $t1
    scil_volume_math.py convert t1_seg_2.nii.gz ${prefix}__mask_wm.nii.gz --data_type uint8
    scil_volume_math.py convert t1_seg_1.nii.gz ${prefix}__mask_gm.nii.gz --data_type uint8
    scil_volume_math.py convert t1_seg_0.nii.gz ${prefix}__mask_csf.nii.gz --data_type uint8
    mv t1_pve_2.nii.gz ${prefix}__map_wm.nii.gz
    mv t1_pve_1.nii.gz ${prefix}__map_gm.nii.gz
    mv t1_pve_0.nii.gz ${prefix}__map_csf.nii.gz

    scil_volume_stats_in_ROI.py ${prefix}__mask*.nii.gz  \
        --metrics $metrics \
        --sort_keys \
        $bin $normalize_weights > ${prefix}__${suffix}.json

    """
}