#!/usr/bin/env nextflow

//include { RECONST_DTIMETRICS } from 'modules/nf-neuro/reconst/dtimetrics/main'
workflow get_data {
    main:
        if ( !params.input ) {
            log.info "You must provide an input directory containing all files using:"
            log.info ""
            log.info "    --input=/path/to/[input]   Input directory containing the file needed"
            log.info "                        |"
            log.info "                        └-- Input"
            log.info "                             └-- participants.*"
            log.info ""
            error "Please resubmit your command with the previous file structure."
        }

        input = file(params.input)
        // ** Loading all files. ** //
        participants_channel = Channel.fromFilePairs("$input/participants.*", size: 2, flat: true)

    emit:
        participants = participants_channel
}

workflow {
    // ** Now call your input workflow to fetch your files ** //
    data = get_data()
    data.participants.view()
}