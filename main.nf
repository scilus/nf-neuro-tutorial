#!/usr/bin/env nextflow

 include { RECONST_DTIMETRICS } from './modules/nf-neuro/reconst/dtimetrics/main'
 include { RECONST_DTIMETRICS as RECONST_DTI_DENOISED } from './modules/nf-neuro/reconst/dtimetrics/main'
 include { DENOISING_MPPCA } from './modules/nf-neuro/denoising/mppca/main'
 include { PREPROC_T1 } from './subworkflows/nf-neuro/preproc_t1/main'
 include { STATS_METRICSINROI } from './modules/local/stats/metricsinrois/main'



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
        dwi = dwi_channel // 
        anat = t1_channel
}

// workflow {
//     // ** Now call your input workflow to fetch your files ** //
//     data = get_data()
//     data.dwi.view() // Contains your DWI data: [meta, dwi, bval, bvec]
//     data.anat.view() // Contains your anatomical data (T1 in this case): [meta, t1]
// }
// execute nextflow run 

workflow {
    // 1. Create input channel for dwi and bvec_bval inputs (Branch 1)
    inputs = get_data()

    ch_dwi_bvalbvec = inputs.dwi
        .multiMap { meta, dwi, bval, bvec ->
            dwi:            [ meta, dwi ]
            bvs_files:      [ meta, bval, bvec ]
            dwi_bval_bvec:  [ meta, dwi, bval, bvec ]
        }
    // ch_dwi_bvalbvec.view()
    // execute nextflow run 

    // 2. Create DTI channel (Branch 2)
    dti_ch = ch_dwi_bvalbvec.dwi_bval_bvec
            .map{ it + [[]] }
    //dti_ch.view()

    // 3. Bind DTI module
    RECONST_DTIMETRICS( dti_ch )

    // 4. add parameters and publish stiff to nextflow.config
    // execute nextflow run 

    // 5. Install and add DENOISING_MPPCA module (branch 3)
    ch_denoise_dwi = ch_dwi_bvalbvec.dwi
                    .map{ it + [[]] }
    DENOISING_MPPCA( ch_denoise_dwi )
    ch_dwi_denoised = DENOISING_MPPCA.out.image

    // 6. Join output from module to input or other
    dti_denoised_ch = ch_dwi_denoised
            .join(ch_dwi_bvalbvec.bvs_files)
            .map{ it + [[]] }
    //dti_ch.view()
    // execute nextflow run 

    // 7. Bind DTI module
    RECONST_DTI_DENOISED( dti_denoised_ch )
    // execute nextflow run 


    // 8. Install and Add subworkflow for processing T1 (PREPROC_T1) (branch 4)
    // 8.1 Install swf
    // 8.2 Bind preproc
    // 8.3 Add element to netxflow.config

    PREPROC_T1(
        inputs.anat,
        Channel.empty(),
        Channel.empty(),
        Channel.empty(),
        Channel.empty(),
        Channel.empty(),
        Channel.empty()
    )
    

    // 9. Channel to join dwi_denoised from DENOISING_MPPCA and t1_final from PREPROC_T1 (Branch 5)
    ch_extract_metric = PREPROC_T1.out.image_bet
                .join(RECONST_DTIMETRICS.out.fa)
                .map{ it }
    // ch_extract_metric.view()
    // execute nextflow run

    // 10. Add custom local module : EXTRACT_METRIC module
    // 10.1 create local folder in modules
    // 10.2 create main.nf
    // 10.3 Bind local module
    // 10.4 Add required stuff to nextflow.config
    STATS_METRICSINROI( ch_extract_metric )

    // 11. Add custom local module : EXTRACT_METRIC module (branch 6, optional)
    // publish, SaveAs, ... in nextflow.config

}