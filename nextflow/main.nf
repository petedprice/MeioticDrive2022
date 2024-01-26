nextflow.enable.dsl=2

include { alevin_fry_index } from './modules/alevin_fry_index.nf'
include { alevin_fry_quant } from './modules/alevin_fry_quant.nf'

include { cellranger_mkref } from './modules/cellranger_mkref.nf'
include { cellranger_count } from './modules/cellranger_count.nf'
include { seurat_filter } from './modules/seurat/seurat_filter.nf'
include { seurat_doublet } from './modules/seurat/seurat_doublet.nf'

include { seurat_SCT_integrate } from './modules/seurat/seurat_SCT_integrate.nf'
include { seurat_SCT_integrate_noX } from './modules/seurat/seurat_SCT_integrate_noX.nf'

include { seurat_markers } from './modules/seurat/seurat_markers.nf'
include { seurat_markers_noX } from './modules/seurat/seurat_markers_noX.nf'


include { sc_var_call } from './modules/sc_var_call.nf'
include { R_var_filtering } from './modules/R_var_filtering.nf'
include { contig_names } from './modules/contig_names.nf'
include { split_bam } from './modules/split_bam.nf'

include { ws_var_call_cleanup.nf } from './modules/ws_varcall/ws_var_call_cleanup.nf'
include { ws_var_call_splitNCR.nf } from './modules/ws_varcall/ws_var_call_splitNCR.nf'
include { ws_var_call_HC } from './modules/ws_varcall/ws_var_call_HC.nf'
include { ws_var_call_VF_stringent } from './modules/ws_varcall/ws_var_call_VF_stringent.nf'
include { ws_var_call_VF_relaxed } from './modules/ws_varcall/ws_var_call_VF_relaxed.nf'

workflow {
    //Channels species name and reference name
    species_ch=Channel.fromPath(params.metadata)
	.splitCsv()
	.map {row -> tuple(row[1], row[2])}
	.unique()
	

    //Channels sample name and species name
    samples_ch = Channel
	.fromPath(params.metadata)
	.splitCsv()
	.map {row -> tuple(row[1], row[0])}

    species_mt_ch=Channel.fromPath(params.metadata)
        .splitCsv()
        .map {row -> tuple(row[1], row[2], row[5])}
        .unique()

    if (params.alevin_quant == 'FALSE'){
	    ref_made=cellranger_mkref(species_ch).view()
	    counted=cellranger_count(ref_made.combine(samples_ch, by: 0))
    } else if (params.alevin_quant == 'TRUE'){
    	    ref_made_alevin=alevin_fry_index(species_ch)
	    counted=alevin_fry_quant(ref_made_alevin.combine(samples_ch, by: 0).view())
    }


    seurat_filtered=seurat_filter(counted.combine(species_mt_ch, by: 0))


    seurat_doubleted=seurat_doublet(seurat_filtered)
	.collect()

    seurat_integrated=seurat_SCT_integrate(seurat_doubleted)

    seurat_marked=seurat_markers(seurat_integrated)


    if (params.sc_var_call == 'TRUE'){
	cns=contig_names(species_ch)
        	.transpose()
		.view()

        splitted=split_bam(counted.combine(cns, by: 0))
	//var_called=sc_var_call(splitted.combine(species_ch, by: 0))
	//var_filtered=R_var_filtering(var_called)
	cleaned=ws_var_call_cleanup(counted)
	splitNCRed=ws_var_call_splitNCR(cleaned.combine(species_ch, by: 0))
	ws_VCed=ws_var_call_HC(splitNCRed)
	VFs=ws_var_call_VF_stringent(ws_VCed)
	VFr=ws_var_call_VF_relaxed(ws_VCed)
	var_called_sc=sc_var_call(VFr)


        

    }









/*    

    pre_noX=seurat_doubleted
    pre_noX.view()
    seurat_integrated_nX=seurat_SCT_integrate_noX(pre_noX)
    seurat_marked_nX=seurat_markers_noX(seurat_integrated_nX)

 */
}




