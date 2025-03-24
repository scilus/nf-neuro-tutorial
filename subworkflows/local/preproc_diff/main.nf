include { RECONST_DTIMETRICS } from '../../../modules/nf-neuro/reconst/dtimetrics/main'
include { DENOISING_MPPCA } from '../../../modules/nf-neuro/denoising/mppca/main'

workflow PREPROC_DIFF {

    take:
        ch_dwi           // channel: [ val(meta), dwi, bval, bvec ]

    main:
        ch_multiqc_files = Channel.empty()

        // ** Denoise DWI ** //
        if (params.preproc_dwi_run_denoising) {
            ch_dwi_bvalbvec = ch_dwi
                .multiMap { meta, dwi, bval, bvec ->
                    dwi:    [ meta, dwi ]
                    bvs_files: [ meta, bval, bvec ]
                }

            ch_denoise_dwi = ch_dwi_bvalbvec.dwi
                .map{ it + [[]] }

            DENOISING_MPPCA ( ch_denoise_dwi )

            // Fetch specific output
            ch_dwi = DENOISING_MPPCA.out.image
                .join(ch_dwi_bvalbvec.bvs_files)
        }

        // Input DTI update with DWI denoised output
        input_dti = ch_dwi.map{ it + [[]] }

        // DTI-derived metrics
        RECONST_DTIMETRICS( input_dti )
        ch_multiqc_files = ch_multiqc_files.mix(RECONST_DTIMETRICS.out.mqc)

    emit:
        dwi                 = ch_dwi_bvalbvec.dwi           // channel: [ val(meta), dwi-raw ]
        dwi_denoised        = DENOISING_MPPCA.out.image     // channel: [ val(meta), dwi-after-mppca ]
        bvs_files           = ch_dwi_bvalbvec.bvs_files     // channel: [ val(meta), bval, bvec ]
        fa                  = RECONST_DTIMETRICS.out.fa     // channel: [ val(meta), fa ]
        md                  = RECONST_DTIMETRICS.out.md     // channel: [ val(meta), md ]
        mqc                 = ch_multiqc_files              // channel: [ val(meta), mqc ]

}