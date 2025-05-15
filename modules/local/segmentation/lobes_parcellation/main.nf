process GENERATE_LOBES_PARCELLATION {
    tag "$meta.id"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://scil.usherbrooke.ca/containers/scilus_2.0.2.sif':
        'scilus/scilus:2.0.2' }"

    input:
    tuple val(meta), path(fs_folder)

    output:
    tuple val(meta), path("atlas_lobes.nii.gz")            , emit: labels
    tuple val(meta), path("atlas_lobes_dilate.nii.gz")     , emit: labels_dilate

    script:
    """
    export MPLCONFIGDIR=./
    mrconvert ${fs_folder}/mri/rawavg.mgz rawavg.nii.gz
    mrconvert ${fs_folder}/mri/wmparc.mgz wmparc.nii.gz
    mrconvert ${fs_folder}/mri/brainmask.mgz brain_mask.nii.gz

    scil_volume_reshape_to_reference.py wmparc.nii.gz rawavg.nii.gz \
        wmparc.nii.gz --interpolation nearest -f 
    scil_volume_math.py convert wmparc.nii.gz wmparc.nii.gz --data_type uint16 -f
    
    scil_volume_math.py lower_threshold brain_mask.nii.gz 0.001 brain_mask.nii.gz \
        --data_type uint8 -f
    scil_volume_math.py dilation brain_mask.nii.gz 1 brain_mask.nii.gz -f
    scil_volume_reshape_to_reference.py brain_mask.nii.gz rawavg.nii.gz \
        brain_mask.nii.gz --interpolation nearest \
        --interpolation nearest --keep_dtype -f 

    scil_labels_combine.py atlas_lobes.nii.gz \
        --volume_ids wmparc.nii.gz 1003 1012 1014 1017 1018 1019 1020 1024 1027 1028 1032 \
        --volume_ids wmparc.nii.gz 1008 1022 1025 1029 1031 \
        --volume_ids wmparc.nii.gz 1005 1011 1013 1021 \
        --volume_ids wmparc.nii.gz 1001 1006 1007 1009 1015 1015 1030 1033 \
        --volume_ids wmparc.nii.gz 1002 1010 1023 1026 \
        --volume_ids wmparc.nii.gz 8 \
        --volume_ids wmparc.nii.gz 10 11 12 13 17 18 26 28 \
        --volume_ids wmparc.nii.gz 2003 2012 2014 2017 2018 2019 2020 2024 2027 2028 2032 \
        --volume_ids wmparc.nii.gz 2008 2022 2025 2029 2031 \
        --volume_ids wmparc.nii.gz 2005 2011 2013 2021 \
        --volume_ids wmparc.nii.gz 2001 2006 2007 2009 2015 2015 2030 2033 \
        --volume_ids wmparc.nii.gz 2002 2010 2023 2026 \
        --volume_ids wmparc.nii.gz 49 50 51 52 53 54 58 60 \
        --volume_ids wmparc.nii.gz 47 \
        --volume_ids wmparc.nii.gz 16 --merge
    scil_labels_dilate.py atlas_lobes.nii.gz atlas_lobes_dilate.nii.gz \
        --distance 2 --labels_to_dilate 1 2 3 4 5 6 8 9 10 11 12 14 15 \
        --mask brain_mask.nii.gz
    """
}
