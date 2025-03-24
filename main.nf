#!/usr/bin/env nextflow

include { PREPROC_T1 } from './subworkflows/nf-neuro/preproc_t1/main'
include { STATS_METRICSINROI } from './modules/local/stats/metricsinrois/main'
include { PREPROC_DIFF } from './subworkflows/local/preproc_diff/main'


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
            log.info "                             ├-- *bval"
            log.info "                             ├-- *bvec"
            log.info "                             └-- *t1.nii.gz"
            log.info ""
            error "Please resubmit your command with the previous file structure."
        }

        input = file(params.input)
        // ** Loading DWI files. ** //
        dwi_channel = Channel.fromFilePairs("$input/**/**/dwi/*dwi.{nii.gz,bval,bvec}", size: 3, flat: true)
            { it.parent.parent.parent.name + "_" + it.parent.parent.name} // Set the subject filename as subjectID + '_' + session.
            .map{ sid, bvals, bvecs, dwi -> [ [id: sid], dwi, bvals, bvecs ] } // Reordering the inputs.
        // ** Loading T1 file. ** //
        t1_channel = Channel.fromFilePairs("$input/**/**/anat/*T1w.nii.gz", size: 1, flat: true)
            { it.parent.parent.parent.name + "_" + it.parent.parent.name } // Set the subject filename as subjectID + '_' + session.
            .map{ sid, t1 -> [ [id: sid], t1 ] }
    emit:
        dwi = dwi_channel
        anat = t1_channel
}

workflow {
    inputs = get_data()

    //Processing DWI
    PREPROC_DIFF( inputs.dwi )

    // Preprocessing T1 images
    //inputs.anat.view()

    PREPROC_T1(
        inputs.anat,
        Channel.empty(),
        Channel.empty(),
        Channel.empty(),
        Channel.empty(),
        Channel.empty(),
        Channel.empty()
    )

    // Extract FA value
     input_extract_metric = PREPROC_T1.out.image_bet
             .join(PREPROC_DIFF.out.fa)
             .map{ it }

    STATS_METRICSINROI( input_extract_metric )
}