#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { TRACTOFLOW } from './subworkflows/nf-neuro/tractoflow/main'

include { SEGMENTATION_FSRECONALL } from './modules/nf-neuro/segmentation/fsreconall/main'
include { GENERATE_LOBES_PARCELLATION } from './modules/local/segmentation/lobes_parcellation/main'

include { BETCROP_ANTSBET} from './modules/nf-neuro/betcrop/antsbet/main'
include { REGISTRATION_ANTS as REGISTER_TO_DWI } from './modules/nf-neuro/registration/ants/main'
include { REGISTRATION_ANTS as REGISTER_TO_MNI} from './modules/nf-neuro/registration/ants/main'

include { REGISTRATION_ANTSAPPLYTRANSFORMS as TRANSFORM_LABELS_TO_DWI} from './modules/nf-neuro/registration/antsapplytransforms/main'
include { REGISTRATION_ANTSAPPLYTRANSFORMS as TRANSFORM_LABELS_TO_MNI } from './modules/nf-neuro/registration/antsapplytransforms/main'
include { REGISTRATION_ANTSAPPLYTRANSFORMS as TRANSFORM_IMAGE_T1_MNI } from './modules/nf-neuro/registration/antsapplytransforms/main'
include { REGISTRATION_ANTSAPPLYTRANSFORMS as TRANSFORM_IMAGE_NUFO_MNI } from './modules/nf-neuro/registration/antsapplytransforms/main'
include { REGISTRATION_ANTSAPPLYTRANSFORMS as TRANSFORM_IMAGE_AFD_MNI } from './modules/nf-neuro/registration/antsapplytransforms/main'
include { REGISTRATION_ANTSAPPLYTRANSFORMS as TRANSFORM_IMAGE_FA_MNI } from './modules/nf-neuro/registration/antsapplytransforms/main'
include { REGISTRATION_ANTSAPPLYTRANSFORMS as TRANSFORM_IMAGE_MD_MNI } from './modules/nf-neuro/registration/antsapplytransforms/main'
include { REGISTRATION_ANTSAPPLYTRANSFORMS as TRANSFORM_MASK_WM_MNI } from './modules/nf-neuro/registration/antsapplytransforms/main'

include { TRACTOGRAM_MATH as CONCATENATE} from './modules/local/tractogram/math/main'
include { REGISTRATION_TRACTOGRAM as TRANSFORM_TRACTOGRAM_MNI} from './modules/nf-neuro/registration/tractogram/main'

include { CONNECTIVITY_DECOMPOSE } from './modules/nf-neuro/connectivity/decompose/main'
include { GENERATE_JUNCTION_SIGNATURES } from './modules/local/connectivity/signatures/main'

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
    
        ch_ants_template = params.ants_template
            ? Channel.fromPath(params.ants_template, checkIfExists: true, followLinks: true)
            : Channel.empty().ifEmpty { error "No ANTS template file path provided. Please specify the path using --ants_template parameter." }
        ch_ants_probability_map = params.ants_probability_map
            ? Channel.fromPath(params.ants_probability_map, checkIfExists: true, followLinks: true)
            : Channel.empty().ifEmpty { error "No ANTS probability map file path provided. Please specify the path using --ants_probability_map parameter." }

        ch_all_signatures = params.all_signatures
            ? Channel.fromPath(params.all_signatures, checkIfExists: true, followLinks: true)
            : Channel.empty().ifEmpty { error "No all signatures file path provided. Please specify the path using --all_signatures parameter." }

        ch_signatures_mapping = params.signatures_mapping
            ? Channel.fromPath(params.signatures_mapping, checkIfExists: true, followLinks: true)
            : Channel.empty().ifEmpty { error "No signatures mapping file path provided. Please specify the path using --signatures_mapping parameter." }

    emit:
        dwi = dwi_channel
        anat = t1_channel
        fs_license = ch_fs_license
        mni_template = ch_mni_template
        ants_template = ch_ants_template
        ants_probability_map = ch_ants_probability_map
        all_signatures = ch_all_signatures
        signatures_mapping = ch_signatures_mapping
}

workflow {
    inputs = get_data()

    // Somehow Tractoflow needs to have a template with a meta: use combine
    ch_template = inputs.anat
        .combine(inputs.ants_template)
        .map { id, t1, template ->
            [id, template]
        }
    ch_probability_map = inputs.anat
        .combine(inputs.ants_probability_map)
        .map { id, t1, probability_map ->
            [id, probability_map]
        }

    TRACTOFLOW(
        inputs.dwi, // channel : [required] dwi, bval, bvec
        inputs.anat, // channel : [required] t1
        Channel.empty(),
        Channel.empty(),
        Channel.empty(),
        Channel.empty(),
        Channel.empty(),
        Channel.empty(),
        ch_template, // Provided template (one per project)
        ch_probability_map, // Provided probability map (one per project)
        Channel.empty()
    )

    // Run Freesurfer on the T1w image
    ch_anat_license = inputs.anat
        .combine(inputs.fs_license)
    SEGMENTATION_FSRECONALL(ch_anat_license)
    
    ch_labels = SEGMENTATION_FSRECONALL.out.recon_all_out_folder
    GENERATE_LOBES_PARCELLATION(ch_labels)

    // Since Tractoflow does not output a registration (affine+warp)
    // from T1w to DWI, we do it again using ANTS.
    // First we skull strip the T1w image using ANTS.
    ch_bet_t1 = inputs.anat
        .combine(inputs.ants_template)
        .combine(inputs.ants_probability_map)
        .map { id, t1, template, probability_map ->
            [id, t1, template, probability_map, [], []]
        }
    BETCROP_ANTSBET(ch_bet_t1)

    ch_register_to_dwi = BETCROP_ANTSBET.out.t1
        .join(TRACTOFLOW.out.t1)
        .map { id, t1, template ->
            [id, template, t1, []]
        }
    REGISTER_TO_DWI ( ch_register_to_dwi )

    // Labels (lobes) can now be transformed to DWI space
    ch_transforms_to_dwi = REGISTER_TO_DWI.out.warp
        .join(REGISTER_TO_DWI.out.affine)
    ch_ants_apply_labels_to_dwi = GENERATE_LOBES_PARCELLATION.out.labels_dilate
        .join(TRACTOFLOW.out.dti_fa)
        .join(ch_transforms_to_dwi)
    TRANSFORM_LABELS_TO_DWI ( ch_ants_apply_labels_to_dwi )
    
    // We will do all processing in MNI space, we will use the skull stripped
    // T1w in DWI space as input for registration to MNI space.
    ch_register = TRACTOFLOW.out.t1
        .combine(inputs.mni_template)
        .map { id, t1, template ->
            [id, template, t1, []]
        }
    REGISTER_TO_MNI ( ch_register )

    // We need to transform the labels and all relevant images to MNI space.
    ch_transform_to_mni = REGISTER_TO_MNI.out.warp
        .join(REGISTER_TO_MNI.out.affine)
    ch_ants_apply_labels_mni = TRANSFORM_LABELS_TO_DWI.out.warped_image
        .combine(inputs.mni_template)
        .join(ch_transform_to_mni)
    TRANSFORM_LABELS_TO_MNI ( ch_ants_apply_labels_mni )

    ch_ants_apply_image_t1_mni = TRACTOFLOW.out.t1
        .combine(inputs.mni_template)
        .join(ch_transform_to_mni)
    TRANSFORM_IMAGE_T1_MNI ( ch_ants_apply_image_t1_mni )

    ch_ants_apply_image_fa_mni = TRACTOFLOW.out.dti_fa
        .combine(inputs.mni_template)
        .join(ch_transform_to_mni)
    TRANSFORM_IMAGE_FA_MNI ( ch_ants_apply_image_fa_mni )

    ch_ants_apply_image_md_mni = TRACTOFLOW.out.dti_md
        .combine(inputs.mni_template)
        .join(ch_transform_to_mni)
    TRANSFORM_IMAGE_MD_MNI ( ch_ants_apply_image_md_mni )

    ch_ants_apply_image_nufo_mni = TRACTOFLOW.out.nufo
        .combine(inputs.mni_template)
        .join(ch_transform_to_mni)
    TRANSFORM_IMAGE_NUFO_MNI ( ch_ants_apply_image_nufo_mni )

    ch_ants_apply_image_afd_mni = TRACTOFLOW.out.afd_total
        .combine(inputs.mni_template)
        .join(ch_transform_to_mni)
    TRANSFORM_IMAGE_AFD_MNI ( ch_ants_apply_image_afd_mni )

    ch_ants_apply_mask_wm_mni = TRACTOFLOW.out.wm_mask
        .combine(inputs.mni_template)
        .join(ch_transform_to_mni)
    TRANSFORM_MASK_WM_MNI ( ch_ants_apply_mask_wm_mni )
    
    ch_concatenate = TRACTOFLOW.out.pft_tractogram
        .concat(TRACTOFLOW.out.local_tractogram)
        .groupTuple()
        .join(TRACTOFLOW.out.dti_fa)

    CONCATENATE ( ch_concatenate )

    // We need to transform the tractogram to MNI space.
    ch_ants_apply_tractogram_to_mni = CONCATENATE.out.trk
        .combine(inputs.mni_template)
        .join(REGISTER_TO_MNI.out.affine)
        .join(TRACTOFLOW.out.dti_fa)
        .join(REGISTER_TO_MNI.out.inverse_warp)
        .map { id, tractogram, template, affine, fa, warp ->
            [id, template, affine, tractogram, fa, warp]
        }
    TRANSFORM_TRACTOGRAM_MNI ( ch_ants_apply_tractogram_to_mni )

    ch_decompose = TRANSFORM_TRACTOGRAM_MNI.out.warped_tractogram
        .join(TRANSFORM_LABELS_TO_MNI.out.warped_image)

    CONNECTIVITY_DECOMPOSE ( ch_decompose )

    ch_junction_signatures = CONNECTIVITY_DECOMPOSE.out.hdf5
        .join(TRANSFORM_MASK_WM_MNI.out.warped_image)
        .join(TRANSFORM_IMAGE_NUFO_MNI.out.warped_image)
        .combine(inputs.all_signatures)
        .combine(inputs.signatures_mapping)
    GENERATE_JUNCTION_SIGNATURES ( ch_junction_signatures )
}