"""
List of sources for reference data
"""

import dataclasses
from shlex import quote
from typing import Protocol

CANONICAL_CHROMOSOMES = [f'chr{x}' for x in list(range(1, 23)) + ['X', 'Y']]


class SyncCommandProtocol(Protocol):
    def __call__(self, src: str, dst: str, project: str) -> str: ...


def quote_command(cmd: list[str]) -> str:
    return ' '.join(map(quote, cmd))


def gcs_rsync(src: str, dst: str, project: str) -> str:
    """
    defines a gcs rsync function
    -u sets the billing project
    -d for deleting files in the destination that are not in the source
    -r for recursive
    """
    assert src.startswith('gs://')
    c = ['gcloud', '--billing-project', project, 'storage', 'rsync', '-r', src, dst]
    return quote_command(c)


def gcs_rsync_no_billing_project(src: str, dst: str, project: str) -> str:
    """
    defines a gcs rsync function, without setting the billing project
    within CI the billing project attempts to reset after ~60 mins, killing transfers
    only use this with public (i.e. not requester-pays) buckets
    -r for recursive
    """
    print(f'ignoring {project} - attempting standard transfer')
    assert src.startswith('gs://')
    c = ['gcloud', 'storage', 'rsync', '-r', src, dst]
    return quote_command(c)


def gcs_cp_single(src: str, dst: str, project: str) -> str:
    """
    defines a single-file gcs copy function
    """
    assert src.startswith('gs://')
    c = ['gcloud', '--billing-project', project, 'storage', 'cp', src, dst]
    return quote_command(c)


def gcs_cp_r(src: str, dst: str, project: str) -> str:
    """
    defines a recursive gcs copy function
    """
    assert src.startswith('gs://')
    c = ['gcloud', '--billing-project', project, 'storage', 'cp', '-r', src, dst]
    return quote_command(c)


def curl(src: str, dst: str, project: str) -> str:
    """
    defines a curl & recursive copy upload function
    """
    assert src.startswith('https://')
    return f'curl -L {quote(src)} | gcloud --billing-project {quote(project)} storage cp - {quote(dst)}'


def curl_with_user_agent(src: str, dst: str, project: str) -> str:
    """
    defines a curl & recursive copy upload function, with a user agent
    """
    assert src.startswith('https://')
    return f'curl -A "Mozilla/5.0" -L {quote(src)} | gcloud --billing-project {quote(project)} storage cp - {quote(dst)}'


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
    transfer_cmd: SyncCommandProtocol | None = None

    def is_folder(self) -> bool:
        """simple folder check using known extensions"""
        if self.files:
            return True

        return (
            self.dst.endswith('.ht')
            or self.dst.endswith('.mt')
            or self.dst.endswith('.vds')
        )


# Genome build. Only GRCh38 is currently supported.
GENOME_BUILD = 'GRCh38'

SOURCES = [
    Source(
        'vep_105_mount',
        # Folder with uncompressed VEP tarballs for mounting with cloudfuse.
        # No `src` field: the process of building it is described in `vep/README.md`.
        # Hopefully to be deprecated once VEP for Hail Query is finalised:
        # https://github.com/hail-is/hail/pull/12428)
        dst='vep/105.0/mount',
    ),
    Source(
        # Folder with uncompressed VEP 110 tarballs for mounting with cloudfuse.
        # see documentation in vep/README.md
        'vep_110_mount',
        dst='vep/110/mount',
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
            wham_include_list_bed_file='sv-resources/resources/v1/wham_whitelist.bed',
            primary_contigs_list='sv-resources/resources/v1/primary_contigs.list',
            primary_contigs_fai='sv-resources/resources/v1/contig.fai',
            manta_region_bed='sv-resources/resources/v1/primary_contigs_plus_mito.bed.gz',
            manta_region_bed_index='sv-resources/resources/v1/primary_contigs_plus_mito.bed.gz.tbi',
            genome_file='sv-resources/resources/v1/hg38.genome',
            wgd_scoring_mask='sv-resources/resources/v1/wgd_scoring_mask.hg38.gnomad_v3.bed',
            allosomal_contigs='sv-resources/resources/v1/allosome.fai',
            inclusion_bed='sv-resources/resources/v1/hg38_primary_contigs.bed',
            autosome_file='sv-resources/resources/v1/autosome.fai',
            allosome_file='sv-resources/resources/v1/allosome.fai',
            cnmops_exclude_list='sv-resources/resources/v1/GRCh38_Nmask.bed',
            cytoband='sv-resources/resources/v1/cytobands_hg38.bed.gz',
            mei_bed='sv-resources/resources/v1/mei_hg38.bed.gz',
            ped_file='sv-resources/ref-panel/1KG/v1/ped/1kg_ref_panel_v1.ped',
            clean_vcf='sv-resources/ref-panel/1KG/v1/calls/ref_panel_1kg_v1.cleaned.vcf.gz',
            ref_panel_bincov_matrix='sv-resources/ref-panel/1KG/v1/merged_evidence/ref_panel_1kg_v1.bincov.bed.gz',
        ),
    ),
    Source(
        'gatk_sv',
        # The Broad resources for running the GATK-SV workflow
        src='gs://gatk-sv-resources-public/hg38/v0/sv-resources',
        dst='gatk-sv/hg38/v0/sv-resources',
        transfer_cmd=gcs_rsync,
        files=dict(
            clustering_config_part1='resources/v1/clustering_config.part_one.tsv',
            clustering_config_part2='resources/v1/clustering_config.part_two.tsv',
            stratification_config_part1='resources/v1/stratify_config.part_one.tsv',
            stratification_config_part2='resources/v1/stratify_config.part_two.tsv',
            clustering_track_sr='resources/v1/hg38.SimpRep.sorted.pad_100.merged.bed',
            clustering_track_sd='resources/v1/hg38.SegDup.sorted.merged.bed',
            clustering_track_rm='resources/v1/hg38.RM.sorted.merged.bed',
            hervk_reference='resources/v1/HERVK.sorted.bed.gz',
            line1_reference='resources/v1/LINE1.sorted.bed.gz',
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
            tel_and_cent_bed='telomeres_and_centromeres/hg38.telomeresAndMergedCentromeres.bed',
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
        'gnomad_sv',
        # Reference data related to gnomAD V4 SV
        src='gs://gatk-sv-resources-public/gnomad_AF/gnomad_v4_SV.Freq.tsv.gz',
        dst='gnomad_v4_SV.Freq.tsv.gz',
        transfer_cmd=gcs_cp_single,
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
        'igv_org_genomes',
        # The reference data used by Seqr to display reads with IGV.js
        src='https://s3.amazonaws.com/igv.org.genomes/hg38/',
        dst='igv_org_genomes/hg38',
        transfer_cmd=curl_with_user_agent,
        files=dict(
            cytoBandIdeo='annotations/cytoBandIdeo.txt.gz',
            hg38_alias='hg38_alias.tab',
            regGeneSorted='refGene.sorted.txt.gz',
        )
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
            bed='HG001_GRCh38_1_22_v4.2.1_benchmark.bed',
        ),
    ),
    Source(
        'SYNDIP',
        dst='validation/SYNDIP',
        files=dict(
            vcf='syndip_truth.vcf.gz',
            index='syndip_truth.vcf.gz.tbi',
            bed='syndip.b38_20180222.bed',
        ),
    ),
    Source(
        'HG002_NA24385',
        dst='validation/HG002_NA24385',
        files=dict(
            vcf='HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz',
            index='HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz.tbi',
            bed='HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed',
        ),
    ),
    Source(
        'HG003_NA24149',
        dst='validation/HG003_NA24149',
        files=dict(
            vcf='HG003_GRCh38_1_22.vcf.gz',
            index='HG003_GRCh38_1_22.vcf.gz.tbi',
            bed='HG003_GRCh38_1_22.bed',
        ),
    ),
    Source(
        'HG004_NA24143',
        dst='validation/HG004_NA24143',
        files=dict(
            vcf='HG004_GRCh38_1_22.vcf.gz',
            index='HG004_GRCh38_1_22.vcf.gz.tbi',
            bed='HG004_GRCh38_1_22.bed',
        ),
    ),
    Source(
        'VCGS_NA12878',
        dst='validation/VCGS_NA12878',
        files=dict(
            vcf='twist_exome_benchmark_truth.vcf.gz',
            index='twist_exome_benchmark_truth.vcf.gz.tbi',
            bed='Twist_Exome_Core_Covered_Targets_hg38.bed',
        ),
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
    Source(
        'star',
        # References for STAR
        dst='star',
        files=dict(ref_dir='2.7.10b/hg38', gtf='hg38/hg38.gtf', fasta='hg38/hg38.fa'),
    ),
    Source(
        'ancestry',
        # representative sites table for PCA
        # generated using the HGDP+1KG dataset
        # generated using the HGDP+1KG dataset + filtered by 10k10k vqsr
        # generated using the HGDP+1KG dataset + filtered by 10k10k vqsr + gnomad vsqr
        dst='ancestry',
        files=dict(
            sites_table='pruned_variants.ht',
            tenk10K_sites_table='tenk10K_pruned_variants.ht',
            gnomad_sites_table='gnomad_pruned_variants.ht',
        ),
    ),
    Source(
        'exomiser_core',
        # The Broad resources for running Exomiser (Default)
        src='gs://gcp-public-data--broad-references/hg38/v0/exomiser/2302_hg38',
        dst='exomiser/core',
        transfer_cmd=gcs_rsync,
        files=dict(
            clinvar_whitelist='2302_hg38_clinvar_whitelist.tsv.gz',
            clinvar_index='2302_hg38_clinvar_whitelist.tsv.gz.tbi',
            genome_h2='2302_hg38_genome.h2.db',
            ensembl_transcripts='2302_hg38_transcripts_ensembl.ser',
            refseq_transcripts='2302_hg38_transcripts_refseq.ser',
            ucsc_transcripts='2302_hg38_transcripts_ucsc.ser',
            variants='2302_hg38_variants.mv.db',
        ),
    ),
    Source(
        'exomiser_cadd',
        # The Broad resources for running Exomiser (CADD)
        src='gs://gcp-public-data--broad-references/hg38/v0/CADD/1.6',
        dst='exomiser/cadd',
        transfer_cmd=gcs_rsync,
        files=dict(
            indel_tsv='gnomad.genomes.r3.0.indel.tsv.gz',
            indel_index='gnomad.genomes.r3.0.indel.tsv.gz.tbi',
            snv_tsv='whole_genome_SNVs.tsv.gz',
            snv_index='whole_genome_SNVs.tsv.gz.tbi',
        ),
    ),
    Source(
        'exomiser_phenotype',
        # The Broad resources for running Exomiser (Phenotype)
        src='gs://gcp-public-data--broad-references/hg38/v0/exomiser/2302_phenotype',
        dst='exomiser/phenotype',
        transfer_cmd=gcs_rsync,
        files=dict(
            pheno_db='2302_phenotype.h2.db',
            hpo_obo='hp.obo',
            rw_string='rw_string_10.mv',
            phenix='phenix',
            phenix_tar='phenix.tar.gz',
        ),
    ),
    Source(
        'exomiser_remm',
        # The Broad resources for running Exomiser (REMM)
        src='gs://gcp-public-data--broad-references/hg38/v0/ReMM/v0.3.1',
        dst='exomiser/remm',
        files=dict(
            remm_tsv='ReMM.v0.3.1.post1.hg38.tsv.gz',
            remm_index='ReMM.v0.3.1.post1.hg38.tsv.gz.tbi',
        ),
    ),
    Source(
        # monarch annotations for running Exomiser 14.X+
        'exomiser_2402_pheno',
        src='https://data.monarchinitiative.org/exomiser/latest/2402_phenotype.zip',
        dst='exomiser_2402/phenotype',
        transfer_cmd=curl,
    ),
    Source(
        # monarch annotations for running Exomiser 14.X+
        'exomiser_2402_core',
        src='https://data.monarchinitiative.org/exomiser/latest/2402_hg38.zip',
        dst='exomiser_2402/core',
        transfer_cmd=curl,
    ),
    Source(
        'hg38_telomeres_and_centromeres_intervals',
        # gnomAD v3 hg38 coordinates for telomeres and centromeres converted to interval_list
        # Created with the convert_bed_to_interval_list_file.py script
        dst='hg38/v0',
        files=dict(
            interval_list='hg38.telomeresAndMergedCentromeres.interval_list'
        ),
    ),
    Source(
        'gnomad_4.1_vcfs',
        src='gs://gcp-public-data--gnomad/release/4.1/vcf/genomes',
        dst='gnomad/v4.1/vcfs',
        files={contig: f'gnomad.genomes.v4.1.sites.{contig}.vcf.bgz' for contig in CANONICAL_CHROMOSOMES},
        transfer_cmd=gcs_rsync,
    ),
    Source(
        # Mito frequencies are only available on 3.1, but we don't need/want the rest of the 3.1 files
        # this will copy the single MT file, then we'll have to get/generate the index separately
        'gnomad_3.1_mito',
        src='gs://gcp-public-data--gnomad/release/3.1/vcf/genomes/gnomad.genomes.v3.1.sites.chrM.vcf.bgz',
        dst='gnomad/v3.1/vcf/gnomad.genomes.r3.1.sites.chrM.vcf.bgz',
        transfer_cmd=gcs_cp_single,
    ),
    Source(
        'gnomad_4.1_ht',
        src='gs://gcp-public-data--gnomad/release/4.1/ht/genomes/gnomad.genomes.v4.1.sites.ht',
        dst='gnomad/v4.1/ht/gnomad.genomes.v4.1.sites.ht',
        transfer_cmd=gcs_rsync,
    ),
    Source(
        'gnomad_4.1_joint_vcfs',
        src='gs://gcp-public-data--gnomad/release/4.1/vcf/joint',
        dst='gnomad/v4.1/joint/vcfs',
        files={contig: f'gnomad.joint.v4.1.sites.{contig}.vcf.bgz' for contig in CANONICAL_CHROMOSOMES},
        transfer_cmd=gcs_rsync_no_billing_project,
    ),
    Source(
        'gnomad_4.1_joint_ht',
        src='gs://gcp-public-data--gnomad/release/4.1/ht/joint/gnomad.joint.v4.1.sites.ht',
        dst='gnomad/v4.1/joint/ht/gnomad.joint.v4.1.sites.ht',
        transfer_cmd=gcs_rsync_no_billing_project,
    ),
    Source(
        'alphamissense',
        # alphamissense raw data, processed HT, and compressed HT
        dst='alphamissense',
        files=dict(
            raw_tsv='alphamissense_38.tsv.gz',
            ht='alphamissense_38.ht',
            ht_tar='alphamissense_38.ht.tar.gz',
        ),
    ),
    Source(
        'ensembl_113',
        # ensembl GFF3, and derived BED files
        dst='ensembl_113',
        files=dict(
            gff3='GRCh38.gff3.gz',
            bed='GRCh38.bed',  # contains a column with each Gene's name/ID
            merged_bed='merged_GRCh38.bed',  # simplified regions, lacks per-gene data
            # copeid from https://ftp.ensembl.org/pub/release-113/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
            # decompressed locally and uploaded to GCP
            unmasked_reference='Homo_sapiens.GRCh38.dna.primary_assembly.fa'
        ),
    ),
    Source(
        'mane_1.4',
        # MANE v.1.4 digest and raw data
        dst='mane_1.4',
        files=dict(
            summary='mane_1.4.summary.txt.gz',  # raw data from MANE
            json='mane_1.4.json',  # parsed into a per-transcript lookup
        ),
    ),
    Source(
        'exome_probesets',
        # exome probset defintions (bed file and interval_list) format
        # downloaded from UCSC with the download_ucsc_exomes.py script
        dst='exome-probesets/hg38',
        files=dict(
            twist_refseq_exome_panel_target_regions_bed='Twist_Exome_RefSeq_targets_hg38.bed',
            twist_refseq_exome_panel_target_regions_interval_list='Twist_Exome_RefSeq_targets_hg38.interval_list',
            twist_exome_2_0_panel_target_regions_bed='TwistExome21.bed',
            twist_exome_2_0_panel_target_regions_interval_list='TwistExome21.interval_list',
            twist_bioscience_core_exome_panel_target_regions_bed='Twist_Exome_Target_hg38.bed',
            twist_bioscience_core_exome_panel_target_regions_interval_list='Twist_Exome_Target_hg38.interval_list',
            twist_comprehensive_exome_panel_target_regions_bed='Twist_ComprehensiveExome_targets_hg38.bed',
            twist_comprehensive_exome_panel_target_regions_interval_list='Twist_ComprehensiveExome_targets_hg38.interval_list',
            agilent_sureselect_all_exon_v7_target_regions_bed='S31285117_Regions.bed',
            agilent_sureselect_all_exon_v7_target_regions_interval_list='S31285117_Regions.interval_list',
            agilent_sureselect_all_exon_v7_covered_by_probes_bed='S31285117_Covered.bed',
            agilent_sureselect_all_exon_v7_covered_by_probes_interval_list='S31285117_Covered.interval_list',
            agilent_sureselect_all_exon_v6_utr_r2_covered_by_probes_bed='S07604624_Covered.bed',
            agilent_sureselect_all_exon_v6_utr_r2_covered_by_probes_interval_list='S07604624_Covered.interval_list',
            agilent_sureselect_all_exon_v6_cosmic_r2_target_regions_bed='S07604715_Regions.bed',
            agilent_sureselect_all_exon_v6_cosmic_r2_target_regions_interval_list='S07604715_Regions.interval_list',
            agilent_sureselect_all_exon_v6_cosmic_r2_covered_by_probes_bed='S07604715_Covered.bed',
            agilent_sureselect_all_exon_v6_cosmic_r2_covered_by_probes_interval_list='S07604715_Covered.interval_list',
            agilent_sureselect_all_exon_v6_r2_target_regions_bed='S07604514_Regions.bed',
            agilent_sureselect_all_exon_v6_r2_target_regions_interval_list='S07604514_Regions.interval_list',
            agilent_sureselect_all_exon_v6_r2_covered_by_probes_bed='S07604514_Covered.bed',
            agilent_sureselect_all_exon_v6_r2_covered_by_probes_interval_list='S07604514_Covered.interval_list',
            agilent_sureselect_all_exon_v6_utr_r2_target_regions_bed='S07604624_Regions.bed',
            agilent_sureselect_all_exon_v6_utr_r2_target_regions_interval_list='S07604624_Regions.interval_list',
            agilent_sureselect_all_exon_v5_utrs_target_regions_bed='S04380219_Regions.bed',
            agilent_sureselect_all_exon_v5_utrs_target_regions_interval_list='S04380219_Regions.interval_list',
            agilent_sureselect_all_exon_v5_utrs_covered_by_probes_bed='S04380219_Covered.bed',
            agilent_sureselect_all_exon_v5_utrs_covered_by_probes_interval_list='S04380219_Covered.interval_list',
            agilent_sureselect_all_exon_v4_utrs_target_regions_bed='S04380110_Regions.bed',
            agilent_sureselect_all_exon_v4_utrs_target_regions_interval_list='S04380110_Regions.interval_list',
            agilent_sureselect_all_exon_v4_utrs_covered_by_probes_bed='S04380110_Covered.bed',
            agilent_sureselect_all_exon_v4_utrs_covered_by_probes_interval_list='S04380110_Covered.interval_list',
            agilent_sureselect_focused_exome_target_regions_bed='S07084713_Regions.bed',
            agilent_sureselect_focused_exome_target_regions_interval_list='S07084713_Regions.interval_list',
            agilent_sureselect_focused_exome_covered_by_probes_bed='S07084713_Covered.bed',
            agilent_sureselect_focused_exome_covered_by_probes_interval_list='S07084713_Covered.interval_list',
            agilent_sureselect_clinical_research_exome_v2_target_regions_bed='S30409818_Regions.bed',
            agilent_sureselect_clinical_research_exome_v2_target_regions_interval_list='S30409818_Regions.interval_list',
            agilent_sureselect_clinical_research_exome_v2_covered_by_probes_bed='S30409818_Covered.bed',
            agilent_sureselect_clinical_research_exome_v2_covered_by_probes_interval_list='S30409818_Covered.interval_list',
            roche_seqcap_ez_medexome_mito_empirical_target_regions_bed='SeqCap_EZ_MedExomePlusMito_hg38_empirical_targets.bed',
            roche_seqcap_ez_medexome_mito_empirical_target_regions_interval_list='SeqCap_EZ_MedExomePlusMito_hg38_empirical_targets.interval_list',
            roche_seqcap_ez_medexome_mito_capture_probe_footprint_bed='SeqCap_EZ_MedExomePlusMito_hg38_capture_targets.bed',
            roche_seqcap_ez_medexome_mito_capture_probe_footprint_interval_list='SeqCap_EZ_MedExomePlusMito_hg38_capture_targets.interval_list',
            roche_seqcap_ez_medexome_empirical_target_regions_bed='SeqCap_EZ_MedExome_hg38_empirical_targets.bed',
            roche_seqcap_ez_medexome_empirical_target_regions_interval_list='SeqCap_EZ_MedExome_hg38_empirical_targets.interval_list',
            roche_seqcap_ez_medexome_capture_probe_footprint_bed='SeqCap_EZ_MedExome_hg38_capture_targets.bed',
            roche_seqcap_ez_medexome_capture_probe_footprint_interval_list='SeqCap_EZ_MedExome_hg38_capture_targets.interval_list',
            roche_kapa_hyperexome_primary_target_regions_bed='KAPA_HyperExome_hg38_primary_targets.bed',
            roche_kapa_hyperexome_primary_target_regions_interval_list='KAPA_HyperExome_hg38_primary_targets.interval_list',
            roche_kapa_hyperexome_capture_probe_footprint_bed='KAPA_HyperExome_hg38_capture_targets.bed',
            roche_kapa_hyperexome_capture_probe_footprint_interval_list='KAPA_HyperExome_hg38_capture_targets.interval_list',
            idt_xgen_exome_research_panel_v2_target_regions_bed='xgen-exome-research-panel-v2-targets-hg38.bed',
            idt_xgen_exome_research_panel_v2_target_regions_interval_list='xgen-exome-research-panel-v2-targets-hg38.interval_list',
            idt_xgen_exome_research_panel_v2_probes_bed='xgen-exome-research-panel-v2-probes-hg38.bed',
            idt_xgen_exome_research_panel_v2_probes_interval_list='xgen-exome-research-panel-v2-probes-hg38.interval_list',
            idt_xgen_exome_research_panel_v1_target_regions_bed='xgen-exome-research-panel-targets-hg38.bed',
            idt_xgen_exome_research_panel_v1_target_regions_interval_list='xgen-exome-research-panel-targets-hg38.interval_list',
            idt_xgen_exome_research_panel_v1_probes_bed='xgen-exome-research-panel-probes-hg38.bed',
            idt_xgen_exome_research_panel_v1_probes_interval_list='xgen-exome-research-panel-probes-hg38.interval_list',
            ),
    ),
]
