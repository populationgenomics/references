#!/usr/bin/env python3

"""
One-off scan of the masked SpliceAI VCFs for repeated (variant, gene) records.

One Batch job per file and contig streams the contig through bcftools and reports, per
gene symbol with any repeated record, the number of repeated keys, the position range,
the most copies seen, how many (key, score field) pairs differ between copies, and the
largest spread of one score field across copies with the record it occurs at. Nothing is
written to a bucket; the gather job prints the combined table. Each job localises its
whole VCF, since the bcftools image has no gcloud to mint a token for htslib's gs:// reads.
"""

from cpg_utils.config import reference_path
from cpg_utils.hail_batch import get_batch, image_path

CONTIGS = [str(c) for c in [*range(1, 23), 'X', 'Y']]

# Records at one POS are consecutive, so only one position's records are held at a time.
AWK = r"""
BEGIN { FS = OFS = "\t"; split("DS_AG DS_AL DS_DG DS_DL", name, " ") }
$1 != pos { flush(); pos = $1 }
{
    split($4, a, "|")
    k = $2 SUBSEP $3 SUBSEP a[2]
    n[k]++
    for (i = 1; i <= 4; i++) {
        v = a[i + 2] + 0
        if (!((k, i) in lo) || v < lo[k, i]) lo[k, i] = v
        if (!((k, i) in hi) || v > hi[k, i]) hi[k, i] = v
    }
}
function flush(  k, p, sym, i, delta) {
    for (k in n) if (n[k] > 1) {
        split(k, p, SUBSEP); sym = p[3]
        c[sym]++
        if (!(sym in first)) first[sym] = pos
        last[sym] = pos
        if (n[k] > copies[sym]) copies[sym] = n[k]
        for (i = 1; i <= 4; i++) {
            delta = hi[k, i] - lo[k, i]
            if (delta > 0) differing[sym]++
            if (delta > maxdelta[sym] + 0) {
                maxdelta[sym] = delta; where[sym] = pos ":" p[1] ">" p[2] " " name[i] " " lo[k, i] "-" hi[k, i]
            }
        }
    }
    delete n; delete lo; delete hi
}
END {
    flush()
    for (sym in c) print kind, contig, sym, c[sym], first[sym], last[sym], copies[sym], differing[sym] + 0, maxdelta[sym] + 0, where[sym]
}
"""


def main():
    b = get_batch('Scan SpliceAI VCFs for repeated (variant, gene) records')
    outputs = []
    for kind, storage in [('snv', '40Gi'), ('indel', '90Gi')]:
        vcf = b.read_input_group(
            **{
                'vcf.gz': reference_path(f'spliceai_resources/splice_ai_{kind}s'),
                'vcf.gz.tbi': reference_path(
                    f'spliceai_resources/splice_ai_{kind}s_index'
                ),
            }
        )
        for contig in CONTIGS:
            j = b.new_job(f'scan {kind} chr{contig}')
            j.image(image_path('bcftools_120'))
            j.cpu(2)
            j.storage(storage)
            j.command(
                'set -o pipefail && '
                f"bcftools query -r {contig} -f '%POS\\t%REF\\t%ALT\\t%INFO/SpliceAI\\n' "
                f"{vcf['vcf.gz']} | awk -v kind={kind} -v contig={contig} '{AWK}' > {j.ofile}"
                f' && cat {j.ofile}'
            )
            outputs.append(j.ofile)
    gather = b.new_job('gather')
    gather.image(image_path('bcftools_120'))
    gather.command(
        'echo -e "kind\\tcontig\\tsymbol\\trepeated_keys\\tfirst_pos\\tlast_pos\\tmax_copies\\tdiffering_scores"; '
        + 'cat ' + ' '.join(outputs)
    )
    b.run(wait=False)


if __name__ == '__main__':
    main()
