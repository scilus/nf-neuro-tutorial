#!/usr/bin/env nextflow

include { TRACTOFLOW } from './subworkflows/nf-neuro/tractoflow/main'
include { SEGMENTATION_FSRECONALL } from './modules/nf-neuro/segmentation/fsreconall/main'
include { REGISTRATION_ANTSAPPLYTRANSFORMS as TRANSFORM_LABELS } from './modules/nf-neuro/registration/antsapplytransforms/main'
include { CONNECTIVITY_DECOMPOSE } from './modules/nf-neuro/connectivity/decompose/main'


workflow get_data {
    main:
        if ( !params.input ) {
            log.info "You must provide an input directory containing all images using:"
            log.info ""
            log.info "    --input=/path/to/[input]   Input directory containing your subjects"
            log.info "                        |"
            log.info "                        ├-- S1"
            log.info "                        |    ├-- *dwi.nii.gz"
            log.info "                        |    ├-- *dwi.bval"
            log.info "                        |    ├-- *dwi.bvec"
            log.info "                        |    └-- *t1.nii.gz"
            log.info "                        └-- S2"
            log.info "                             ├-- *dwi.nii.gz"
            log.info "                             ├-- *dwi.bval"
            log.info "                             ├-- *dwi.bvec"
            log.info "                             └-- *t1.nii.gz"
            log.info ""
            error "Please resubmit your command with the previous file structure."
        }

        input = file(params.input)
        // ** Loading DWI files. ** //
        dwi_channel = Channel.fromFilePairs("$input/**/*dwi.{nii.gz,bval,bvec}", size: 3, flat: true)
            { it.parent.parent.parent.name + "_" + it.parent.parent.name} // Set the subject filename as subjectID + '_' + session.
            .map{ sid, bvals, bvecs, dwi -> [ [id: sid], dwi, bvals, bvecs ] } // Reordering the inputs.
        // ** Loading T1 file. ** //
        t1_channel = Channel.fromFilePairs("$input/**/*t1.nii.gz", size: 1, flat: true)
            { it.parent.parent.parent.name + "_" + it.parent.parent.name } // Set the subject filename as subjectID + '_' + session.
            .map{ sid, t1 -> [ [id: sid], t1 ] }
    emit:
        dwi = dwi_channel
        anat = t1_channel
}

workflow {
    inputs = get_data()

    TRACTOFLOW(
        inputs.dwi, // channel : [required] dwi, bval, bvec
        inputs.anat, // channel : [required] t1
        Channel.empty(),
        Channel.empty(),
        Channel.empty(),
        Channel.empty(),
        Channel.empty(),
        Channel.empty(),
        Channel.empty(),
        Channel.empty(),
        Channel.empty()
    )
}