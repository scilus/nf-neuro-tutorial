process EXTRACT_FSRECONALL_PARCELLATION {
    tag "$meta.id"
    label 'process_single'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://scil.usherbrooke.ca/containers/scilus_2.0.2.sif':
        'scilus/scilus:latest' }"

    input:
        tuple val(meta), path(folder)
    output:
        tuple val(meta), path("*__aparc_aseg.nii.gz"), emit: aparc_aseg
        tuple val(meta), path("*__wmparc.nii.gz"), emit: wmparc

    // Note. In dsl1, we used an additional option:   -parallel -openmp $params.nb_threads.
    // Removed here.
    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mrconvert $folder/mri/aparc+aseg.mgz ${prefix}__aparc_aseg.nii.gz
    mrconvert $folder/mri/wmparc.mgz ${prefix}__wmparc.nii.gz
    """
}
