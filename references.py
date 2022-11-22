"""
List of sources for reference data
"""

import dataclasses


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
    ),
    Source(
        'somalier_sites',
        # Site list for somalier fingerprinting
        src='https://github.com/brentp/somalier/files/3412456/sites.hg38.vcf.gz',
        dst='somalier/sites.hg38.vcf.gz',
    ),
    Source(
        'broad',
        src='gs://gcp-public-data--broad-references/hg38/v0',
        dst='hg38/v0',
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
        )
    ),
    Source(
        'gatk_sv',
        src='gs://gatk-sv-resources-public/hg38/v0/sv-resources',
        dst='hg38/v0/sv-resources',
        files=dict(
            wham_include_list_bed_file='resources/v1/wham_whitelist.bed',
            primary_contigs_list='resources/v1/primary_contigs.list',
            primary_contigs_fai='resources/v1/contig.fai',
            preprocessed_intervals='resources/v1/preprocessed_intervals.interval_list',
            manta_region_bed='resources/v1/primary_contigs_plus_mito.bed.gz',
            melt_standard_vcf_header='resources/v1/melt_standard_vcf_header.txt',
            genome_file='resources/v1/hg38.genome',
            wgd_scoring_mask='resources/v1/wgd_scoring_mask.hg38.gnomad_v3.bed',
            allosomal_contigs='resources/v1/allosome.fai',
            contig_ploidy_priors='resources/v1/hg38.contig_ploidy_priors_homo_sapiens.tsv',
            inclusion_bed='resources/v1/hg38_primary_contigs.bed',
            autosome_file='resources/v1/autosome.fai',
            allosome_file='resources/v1/allosome.fai',
            cnmops_exclude_list='resources/v1/GRCh38_Nmask.bed',
            cytoband='resources/v1/cytobands_hg38.bed.gz',
            mei_bed='resources/v1/mei_hg38.bed.gz',
            rmsk='resources/v1/randomForest_blacklist.withRepMask.bed.gz',
            segdups='resources/v1/hg38.SD_gaps_Cen_Tel_Heter_Satellite_lumpy.blacklist.sorted.merged.bed.gz',
            seed_cutoffs='resources/v1/seed_cutoff.txt',
            pesr_exclude_list='resources/v1/PESR.encode.peri_all.repeats.delly.hg38.blacklist.sorted.bed.gz',
            depth_exclude_list='resources/v1/depth_blacklist.sorted.bed.gz',
            bin_exclude='resources/v1/bin_exclude.hg38.gatkcov.bed.gz',
            empty_file='resources/v1/empty.file',
            protein_coding_gtf='resources/v1/MANE.GRCh38.v0.95.select_ensembl_genomic.gtf',
            # ref panel
            ped_file='ref-panel/1KG/v1/ped/1kg_ref_panel_v1.ped',
            clean_vcf='ref-panel/1KG/v1/calls/ref_panel_1kg_v1.cleaned.vcf.gz',
            ref_panel_bincov_matrix='ref-panel/1KG/v1/merged_evidence/ref_panel_1kg_v1.bincov.bed.gz',
            qc_definitions='ref-panel/1KG/v2/single_sample.qc_definitions.tsv',
            contig_ploidy_model_tar='ref-panel/1KG/v2/gcnv/ref_panel_1kg_v2-contig-ploidy-model.tar.gz',
            model_tar_tmpl='ref-panel/1KG/v2/gcnv/model_files/ref_panel_1kg_v2-gcnv-model-shard-{shard}.tar.gz',
            ref_panel_PE_file_tmpl='ref-panel/tws_SVEvidence/pe/{sample}.pe.txt.gz',
            ref_panel_SR_file_tmpl='ref-panel/tws_SVEvidence/sr/{sample}.sr.txt.gz',
            ref_panel_SD_file_tmpl='ref-panel/tws_SVEvidence/sd/{sample}.sd.txt.gz',
        )
    ),
    Source(
        'gnomad',
        src='gs://gnomad-public-requester-pays/resources/grch38',
        dst='gnomad',
        files=dict(
            tel_and_cent_ht='telomeres_and_centromeres/hg38.telomeresAndMergedCentromeres.ht',
            lcr_intervals_ht='lcr_intervals/LCRFromHengHg38.ht',
            seg_dup_intervals_ht='seg_dup_intervals/GRCh38_segdups.ht',
            clinvar_ht='clinvar/clinvar_20190923.ht',
            hapmap_ht='hapmap/hapmap_3.3.hg38.ht',
            kgp_omni_ht='kgp/1000G_omni2.5.hg38.ht',
            kgp_hc_ht='kgp/1000G_phase1.snps.high_confidence.hg38.ht',
            mills_ht='mills/Mills_and_1000G_gold_standard.indels.hg38.ht',
        )
    ),
    Source(
        'seqr_combined_reference_data',
        src='gs://seqr-reference-data/GRCh38/all_reference_data/combined_reference_data_grch38.ht',
        dst='seqr/combined_reference_data_grch38.ht',
    ),
    Source(
        'seqr_clinvar',
        src='gs://seqr-reference-data/GRCh38/clinvar/clinvar.GRCh38.2022-09-17.ht',
        dst='seqr/clinvar.GRCh38.2022-09-17.ht',
    ),
    Source(
        'syndip',
        src='gs://gcp-public-data--gnomad/resources/grch38/syndip',
        dst='validation/syndip',
        files=dict(
            truth_vcf='full.38.20180222.vcf.gz',
            regions_bed='syndip.b38_20180222.bed',
            truth_mt='syndip.b38_20180222.mt',
            regions_ht='syndip_b38_20180222_hc_regions.ht',
        )
    ),
    Source(
        'na12878',
        src='gs://gcp-public-data--gnomad/resources/grch38/na12878',
        dst='validation/na12878',
        files=dict(
            truth_vcf='HG001_GRCh38_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_PGandRTGphasetransfer.vcf.gz',
            regions_bed='HG001_GRCh38_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_nosomaticdel_noCENorHET7.bed',
            truth_mt='HG001_GRCh38_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_PGandRTGphasetransfer.mt',
            regions_ht='HG001_GRCh38_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_nosomaticdel_noCENorHET7_hc_regions.ht',
        )
    )
]
