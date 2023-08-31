"""
List of sources for reference data
"""

import dataclasses
import os
from typing import Callable


PROJECT = os.environ['PROJECT']


def gcs_rsync(src: str, dst: str) -> str:
    assert src.startswith('gs://')
    return f'gsutil -u {PROJECT} -m rsync -d -r {src} {dst}'


def gcs_cp_r(src: str, dst: str) -> str:
    assert src.startswith('gs://')
    return f'gcloud --billing-project {PROJECT} storage cp -r {src} {dst}'


def curl(src: str, dst: str) -> str:
    assert src.startswith('https://')
    return (
        f'curl -L {src} -o tmp && '
        f'gcloud --billing-project {PROJECT} storage cp -r tmp {dst}'
    )


@dataclasses.dataclass
class Source:
    """
    Specifies one source location to pull data from (either a GCS bucket or an
    HTTP URL), along with an optional map of keys/files to expand this source as
    a section in the finalised config. E.g.
    """

    name: str  # name of the resource/section
    dst: str  # destination suffix, to be appended to PREFIX
    src: str | None = None  # fully qualified source URL
    files: dict[str, str] | None = None  # map of other suffixes appended to `dst`
    transfer_cmd: Callable[[str, str], str] | None = None

    def is_folder(self) -> bool:
        return self.files or (
            self.dst.endswith('.ht')
            or self.dst.endswith('.mt')
            or self.dst.endswith('.vds')
        )


# Genome build. Only GRCh38 is currently supported.
GENOME_BUILD = 'GRCh38'

SOURCES = [
    Source(
        'vep_mount',
        # Folder with uncompressed VEP tarballs for mounting with cloudfuse.
        # No `src` field: the process of building it is described in `vep/README.md`.
        # Hopefully to be deprecated once VEP for Hail Query is finalised:
        # https://github.com/hail-is/hail/pull/12428)
        dst='vep/105.0/mount',
    ),
    Source(
        'liftover_38_to_37',
        # Liftover chain file to translate from GRCh38 to GRCh37 coordinates
        src='gs://hail-common/references/grch38_to_grch37.over.chain.gz',
        dst='liftover/grch38_to_grch37.over.chain.gz',
        transfer_cmd=gcs_rsync,
    ),
    Source(
        'somalier_sites',
        # Site list for somalier fingerprinting
        src='https://github.com/brentp/somalier/files/3412456/sites.hg38.vcf.gz',
        dst='somalier/sites.hg38.vcf.gz',
        transfer_cmd=curl,
    ),
    Source(
        'broad',
        # The Broad hg38 reference bundle
        src='gs://gcp-public-data--broad-references/hg38/v0',
        dst='hg38/v0',
        transfer_cmd=gcs_rsync,
        files=dict(
            dragmap_prefix='dragen_reference',
            ref_fasta='dragen_reference/Homo_sapiens_assembly38_masked.fasta',
            # Primary contigs BED file
            noalt_bed='sv-resources/resources/v1/primary_contigs_plus_mito.bed.gz',
            # Calling intervals lists
            genome_calling_interval_lists='wgs_calling_regions.hg38.interval_list',
            exome_calling_interval_lists='exome_calling_regions.v1.interval_list',
            genome_evaluation_interval_lists='wgs_evaluation_regions.hg38.interval_list',
            exome_evaluation_interval_lists='exome_evaluation_regions.v1.interval_list',
            genome_coverage_interval_list='wgs_coverage_regions.hg38.interval_list',
            unpadded_intervals_file='hg38.even.handcurated.20k.intervals',  # for SV
            # VQSR
            dbsnp_vcf='Homo_sapiens_assembly38.dbsnp138.vcf',
            dbsnp_vcf_index='Homo_sapiens_assembly38.dbsnp138.vcf.idx',
            hapmap_vcf='hapmap_3.3.hg38.vcf.gz',
            hapmap_vcf_index='hapmap_3.3.hg38.vcf.gz.tbi',
            omni_vcf='1000G_omni2.5.hg38.vcf.gz',
            omni_vcf_index='1000G_omni2.5.hg38.vcf.gz.tbi',
            one_thousand_genomes_vcf='1000G_phase1.snps.high_confidence.hg38.vcf.gz',
            one_thousand_genomes_vcf_index='1000G_phase1.snps.high_confidence.hg38.vcf.gz.tbi',
            mills_vcf='Mills_and_1000G_gold_standard.indels.hg38.vcf.gz',
            mills_vcf_index='Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi',
            axiom_poly_vcf='Axiom_Exome_Plus.genotypes.all_populations.poly.hg38.vcf.gz',
            axiom_poly_vcf_index='Axiom_Exome_Plus.genotypes.all_populations.poly.hg38.vcf.gz.tbi',
            # Genome contamination check
            genome_contam_ud='contamination-resources/1000g/1000g.phase3.100k.b38.vcf.gz.dat.UD',
            genome_contam_bed='contamination-resources/1000g/1000g.phase3.100k.b38.vcf.gz.dat.bed',
            genome_contam_mu='contamination-resources/1000g/1000g.phase3.100k.b38.vcf.gz.dat.mu',
            # Exome contamination check
            exome_contam_ud='contamination-resources/1000g/whole_exome_illumina_coding_v1.Homo_sapiens_assembly38.1000g.contam.UD',
            exome_contam_bed='contamination-resources/1000g/whole_exome_illumina_coding_v1.Homo_sapiens_assembly38.1000g.contam.bed',
            exome_contam_mu='contamination-resources/1000g/whole_exome_illumina_coding_v1.Homo_sapiens_assembly38.1000g.contam.mu',
            # shifted from the gatk-sv-resources-public section
            wham_include_list_bed_file = "sv-resources/resources/v1/wham_whitelist.bed",
            primary_contigs_list = "sv-resources/resources/v1/primary_contigs.list",
            primary_contigs_fai = "sv-resources/resources/v1/contig.fai",
            manta_region_bed = "sv-resources/resources/v1/primary_contigs_plus_mito.bed.gz",
            manta_region_bed_index = "sv-resources/resources/v1/primary_contigs_plus_mito.bed.gz.tbi",
            genome_file = "sv-resources/resources/v1/hg38.genome",
            wgd_scoring_mask = "sv-resources/resources/v1/wgd_scoring_mask.hg38.gnomad_v3.bed",
            allosomal_contigs = "sv-resources/resources/v1/allosome.fai",
            inclusion_bed = "sv-resources/resources/v1/hg38_primary_contigs.bed",
            autosome_file = "sv-resources/resources/v1/autosome.fai",
            allosome_file = "sv-resources/resources/v1/allosome.fai",
            cnmops_exclude_list = "sv-resources/resources/v1/GRCh38_Nmask.bed",
            cytoband = "sv-resources/resources/v1/cytobands_hg38.bed.gz",
            mei_bed = "sv-resources/resources/v1/mei_hg38.bed.gz",
            ped_file = "sv-resources/ref-panel/1KG/v1/ped/1kg_ref_panel_v1.ped",
            clean_vcf = "sv-resources/ref-panel/1KG/v1/calls/ref_panel_1kg_v1.cleaned.vcf.gz",
            ref_panel_bincov_matrix = "sv-resources/ref-panel/1KG/v1/merged_evidence/ref_panel_1kg_v1.bincov.bed.gz",
        ),
    ),
    Source(
        'gatk_sv',
        # The Broad resources for running the GATK-SV workflow
        src='gs://gatk-sv-resources-public/hg38/v0/sv-resources',
        dst='gatk-sv/hg38/v0/sv-resources',
        transfer_cmd=gcs_rsync,
        files=dict(
            preprocessed_intervals='resources/v1/preprocessed_intervals.interval_list',
            melt_standard_vcf_header='resources/v1/melt_standard_vcf_header.txt',
            contig_ploidy_priors='resources/v1/hg38.contig_ploidy_priors_homo_sapiens.tsv',
            rmsk='resources/v1/hg38.randomForest_blacklist.withRepMask.bed.gz',
            segdups='resources/v1/hg38.SD_gaps_Cen_Tel_Heter_Satellite_lumpy.blacklist.sorted.merged.bed.gz',
            seed_cutoffs='resources/v1/seed_cutoff.txt',
            pesr_exclude_list='resources/v1/PESR.encode.peri_all.repeats.delly.hg38.blacklist.sorted.bed.gz',
            depth_exclude_list='resources/v1/depth_blacklist.sorted.bed.gz',
            bin_exclude='resources/v1/bin_exclude.hg38.gatkcov.bed.gz',
            empty_file='resources/v1/empty.file',
            protein_coding_gtf='resources/v1/MANE.GRCh38.v1.2.ensembl_genomic.gtf',
            # ref panel
            qc_definitions='ref-panel/1KG/v2/single_sample.qc_definitions.tsv',
            contig_ploidy_model_tar='ref-panel/1KG/v2/gcnv/ref_panel_1kg_v2-contig-ploidy-model.tar.gz',
            model_tar_tmpl='ref-panel/1KG/v2/gcnv/model_files/ref_panel_1kg_v2-gcnv-model-shard-{shard}.tar.gz',
            ref_panel_PE_file_tmpl='ref-panel/tws_SVEvidence/pe/{sample}.pe.txt.gz',
            ref_panel_SR_file_tmpl='ref-panel/tws_SVEvidence/sr/{sample}.sr.txt.gz',
            ref_panel_SD_file_tmpl='ref-panel/tws_SVEvidence/sd/{sample}.sd.txt.gz',
            recalibrate_gq_repeatmasker='resources/v1/ucsc-genome-tracks/hg38-RepeatMasker.bed.gz',
            recalibrate_gq_segmental_dups='resources/v1/ucsc-genome-tracks/hg38-Segmental-Dups.bed.gz',
            recalibrate_gq_simple_reps='resources/v1/ucsc-genome-tracks/hg38-Simple-Repeats.bed.gz',
            recalibrate_gq_umap_s100='resources/v1/ucsc-genome-tracks/hg38_umap_s100.bed.gz',
            recalibrate_gq_umap_s24='resources/v1/ucsc-genome-tracks/hg38_umap_s24.bed.gz',
            recalibrate_gq_repeatmasker_index='resources/v1/ucsc-genome-tracks/hg38-RepeatMasker.bed.gz.tbi',
            recalibrate_gq_segmental_dups_index='resources/v1/ucsc-genome-tracks/hg38-Segmental-Dups.bed.gz.tbi',
            recalibrate_gq_simple_reps_index='resources/v1/ucsc-genome-tracks/hg38-Simple-Repeats.bed.gz.tbi',
            recalibrate_gq_umap_s100_index='resources/v1/ucsc-genome-tracks/hg38_umap_s100.bed.gz.tbi',
            recalibrate_gq_umap_s24_index='resources/v1/ucsc-genome-tracks/hg38_umap_s24.bed.gz.tbi',
        ),
    ),
    Source(
        'gnomad',
        # The Broad resources for running the GnomAD QC pipeline
        src='gs://gcp-public-data--gnomad/resources/grch38',
        dst='gnomad/v0',
        transfer_cmd=gcs_cp_r,
        files=dict(
            tel_and_cent_ht='telomeres_and_centromeres/hg38.telomeresAndMergedCentromeres.ht',
            lcr_intervals_ht='lcr_intervals/LCRFromHengHg38.ht',
            seg_dup_intervals_ht='seg_dup_intervals/GRCh38_segdups.ht',
            clinvar_ht='clinvar/clinvar_20190923.ht',
            hapmap_ht='hapmap/hapmap_3.3.hg38.ht',
            kgp_omni_ht='kgp/1000G_omni2.5.hg38.ht',
            kgp_hc_ht='kgp/1000G_phase1.snps.high_confidence.hg38.ht',
            mills_ht='mills/Mills_and_1000G_gold_standard.indels.hg38.ht',
            predetermined_qc_variants='sample_qc/pre_ld_pruning_qc_variants.ht',
        ),
    ),
    Source(
        'seqr_combined_reference_data',
        # The Broad resources for annotation for the Seqr Loader
        src='gs://seqr-reference-data/GRCh38/all_reference_data/combined_reference_data_grch38.ht',
        dst='seqr/v0/combined_reference_data_grch38.ht',
        transfer_cmd=gcs_cp_r,
    ),
    Source(
        'seqr_clinvar',
        # The Broad resources for annotation for the Seqr Loader
        src='gs://seqr-reference-data/GRCh38/clinvar/clinvar.GRCh38.ht',
        dst='seqr/v0/clinvar.GRCh38.ht',
        transfer_cmd=gcs_cp_r,
    ),
    Source(
        'syndip',
        # The Broad resources for running variant calling validation on common benchmark datasets
        src='gs://gcp-public-data--gnomad/resources/grch38/syndip',
        dst='validation/syndip',
        transfer_cmd=gcs_rsync,
        files=dict(
            truth_vcf='full.38.20180222.vcf.gz',
            regions_bed='syndip.b38_20180222.bed',
            truth_mt='syndip.b38_20180222.mt',
            regions_ht='syndip_b38_20180222_hc_regions.ht',
        ),
    ),
    Source(
        'na12878',
        # The Broad resources for running variant calling validation on common benchmark datasets
        src='gs://gcp-public-data--gnomad/resources/grch38/na12878',
        dst='validation/na12878',
        transfer_cmd=gcs_rsync,
        files=dict(
            truth_vcf='HG001_GRCh38_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_PGandRTGphasetransfer.vcf.gz',
            regions_bed='HG001_GRCh38_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_nosomaticdel_noCENorHET7.bed',
            truth_mt='HG001_GRCh38_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_PGandRTGphasetransfer.mt',
            regions_ht='HG001_GRCh38_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_nosomaticdel_noCENorHET7_hc_regions.ht',
        ),
    ),
    # validation related content
    Source(
        'HG001_NA12878',
        dst='validation/HG001_NA12878',
        files=dict(
            vcf='HG001_GRCh38_1_22_v4.2.1_benchmark.vcf.gz',
            index='HG001_GRCh38_1_22_v4.2.1_benchmark.vcf.gz.tbi',
            bed='HG001_GRCh38_1_22_v4.2.1_benchmark.bed'
        )
    ),
    Source(
        'SYNDIP',
        dst='validation/SYNDIP',
        files=dict(
            vcf='syndip_truth.vcf.gz',
            index='syndip_truth.vcf.gz.tbi',
            bed='syndip.b38_20180222.bed'
        )
    ),
    Source(
        'HG002_NA24385',
        dst='validation/HG002_NA24385',
        files=dict(
            vcf='HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz',
            index='HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz.tbi',
            bed='HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed'
        )
    ),
    Source(
        'HG003_NA24149',
        dst='validation/HG003_NA24149',
        files=dict(
            vcf='HG003_GRCh38_1_22.vcf.gz',
            index='HG003_GRCh38_1_22.vcf.gz.tbi',
            bed='HG003_GRCh38_1_22.bed'
        )
    ),
    Source(
        'HG004_NA24143',
        dst='validation/HG004_NA24143',
        files=dict(
            vcf='HG004_GRCh38_1_22.vcf.gz',
            index='HG004_GRCh38_1_22.vcf.gz.tbi',
            bed='HG004_GRCh38_1_22.bed'
        )
    ),
    Source(
        'VCGS_NA12878',
        dst='validation/VCGS_NA12878',
        files=dict(
            vcf='twist_exome_benchmark_truth.vcf.gz',
            index='twist_exome_benchmark_truth.vcf.gz.tbi',
            bed='Twist_Exome_Core_Covered_Targets_hg38.bed'
        )
    ),
    Source('stratification', dst='validation/stratification'),
    Source('refgenome_sdf', dst='validation/masked_reference_sdf'),
    Source(
        'gnomad_mito',
        # The Broad resources for running the gnomAD mitochondrial pipeline.
        # Contains two versions of the mito genome + bwa indexes. One wt and one
        # wtih the linearisation point "shifted" by 8000nt
        src='gs://gcp-public-data--broad-references/hg38/v0/chrM/',
        dst='hg38/v0/chrM',
        transfer_cmd=gcs_rsync,
        files=dict(
            dict='Homo_sapiens_assembly38.chrM.dict',
            fasta='Homo_sapiens_assembly38.chrM.fasta',
            shifted_dict='Homo_sapiens_assembly38.chrM.shifted_by_8000_bases.dict',
            shifted_fasta='Homo_sapiens_assembly38.chrM.shifted_by_8000_bases.fasta',
            shift_back_chain='ShiftBack.chain',
            shifted_control_region_interval='control_region_shifted.chrM.interval_list',
            non_control_region_interval='non_control_region.chrM.interval_list',
            blacklist_sites='blacklist_sites.hg38.chrM.bed',
        ),
    ),

]
