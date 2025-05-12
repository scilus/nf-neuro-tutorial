#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { TRACTOFLOW } from './subworkflows/nf-neuro/tractoflow/main'
include { REGISTRATION_ANTS } from './modules/nf-neuro/registration/ants/main'
include { SEGMENTATION_FSRECONALL } from './modules/nf-neuro/segmentation/fsreconall/main'
include { REGISTRATION_ANTSAPPLYTRANSFORMS as TRANSFORM_LABELS } from './modules/nf-neuro/registration/antsapplytransforms/main'
include { CONNECTIVITY_DECOMPOSE } from './modules/nf-neuro/connectivity/decompose/main'
include { BETCROP_SYNTHBET} from './modules/nf-neuro/betcrop/synthbet/main'

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
            { it.parent.parent.name + "_" + it.parent.name} // Set the subject filename as subjectID + '_' + session.
            .map{ sid, bvals, bvecs, dwi -> [ [id: sid], dwi, bvals, bvecs ] } // Reordering the inputs.
        // ** Loading T1 file. ** //
        t1_channel = Channel.fromFilePairs("$input/**/*t1.nii.gz", size: 1, flat: true)
            { it.parent.parent.name + "_" + it.parent.name } // Set the subject filename as subjectID + '_' + session.
            .map{ sid, t1 -> [ [id: sid], t1 ] }

         // ** Fetch license file ** //
        ch_fs_license = params.fs_license
            ? Channel.fromPath(params.fs_license, checkIfExists: true, followLinks: true)
            : Channel.empty().ifEmpty { error "No license file path provided. Please specify the path using --fs_license parameter." }
    
        ch_mni_template = params.mni_template
            ? Channel.fromPath(params.mni_template, checkIfExists: true, followLinks: true)
            : Channel.empty().ifEmpty { error "No template file path provided. Please specify the path using --mni_template parameter." }
    
    emit:
        dwi = dwi_channel
        anat = t1_channel
        fs_license = ch_fs_license
        mni_template = ch_mni_template
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

    ch_anat_license = inputs.anat
        .combine(inputs.fs_license)

    SEGMENTATION_FSRECONALL(ch_anat_license)

    ch_bet_t1 = inputs.anat
        .combine(Channel.empty())
    BETCROP_SYNTHBET(ch_bet_t1)

    ch_register = BETCROP_SYNTHBET.out.bet_image
        .combine(inputs.mni_template)
    REGISTRATION_ANTS ( ch_register )
    // ch_versions = REGISTRATION_ANTS.out.versions.first()

        // ch_transforms = ANATTODWI.out.warp
        //     .join(ANATTODWI.out.affine)
        // ch_peaks = RECONST_FODF.out.peaks
        // ch_fodf = RECONST_FODF.out.fodf
        // ch_dwi_bval_bvec = ch_processed_dwi
        //     .join(ch_processed_bval)
        //     .join(ch_processed_bvec)
        // ch_anat = ANATTODWI.out.t1_warped
        // ch_metrics = RECONST_DTIMETRICS.out.fa
        //     .join(RECONST_DTIMETRICS.out.md)
        //     .join(RECONST_DTIMETRICS.out.ad)
        //     .join(RECONST_DTIMETRICS.out.rd)
        //     .join(RECONST_DTIMETRICS.out.mode)
        //     .join(RECONST_FODF.out.afd_total)
        //     .join(RECONST_FODF.out.nufo)
        //     .map{ meta, fa, md, ad, rd, mode, afd_total, nufo ->
        //         tuple(meta, [ fa, md, ad, rd, mode, afd_total, nufo ])}

        // //
        // // MODULE : Run AntsApplyTransforms.
        // //
        // ch_labels = ch_labels.branch {
        //     reg: it.size() > 2
        //         return [it[0], it[2]]
        //     notreg: it.size() < 3
        //         return [it[0], it[1]]
        // }

        // ch_antsapply = ch_labels.notreg
        //     .join(ch_anat)
        //     .join(ch_transforms)

        // TRANSFORM_LABELS ( ch_antsapply )
        // ch_versions = ch_versions.mix(TRANSFORM_LABELS.out.versions.first())
        // // ch_multiqc_files = ch_multiqc_files.mix(TRANSFORM_LABELS.out.zip.collect{it[1]})
        // ch_nifti_files_to_transform = ch_nifti_files_to_transform
        //     .mix(TRANSFORM_LABELS.out.warped_image)

        // //
        // // MODULE: Run DECOMPOSE.
        // //
        // ch_decompose = ch_trk
        //     .join(ch_labels.reg, remainder: true)
        //     .map { id, trk, reg_label ->
        //         reg_label ? [id, trk, reg_label] : [id, trk, null]
        //     }
        //     .join(TRANSFORM_LABELS.out.warped_image.map { id, warped -> [id, warped] }, remainder: true)
        //     .map { id, trk, reg_label, warped_label ->
        //         def label = reg_label ?: warped_label
        //         [id, trk, label]
        //     }

        // CONNECTIVITY_DECOMPOSE ( ch_decompose )
        // ch_versions = ch_versions.mix(CONNECTIVITY_DECOMPOSE.out.versions.first())


}