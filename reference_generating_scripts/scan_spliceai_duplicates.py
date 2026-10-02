#!/usr/bin/env python3

"""
One-off scan of the masked SpliceAI VCFs for repeated (variant, gene) records.

One Batch job per file and contig streams the contig through bcftools and reports, per
gene symbol with any repeated record, the number of repeated keys, the position range,
the most copies seen and how many repeated keys have differing scores. Nothing is
written to a bucket; the gather job prints the combined table.
"""

from cpg_utils.config import reference_path
from cpg_utils.hail_batch import (
    authenticate_cloud_credentials_in_job,
    get_batch,
    image_path,
)

CONTIGS = [str(c) for c in [*range(1, 23), 'X', 'Y']]

# Records at one POS are consecutive, so only one position's records are held at a time.
AWK = r"""
BEGIN { FS = OFS = "\t" }
$1 != pos { flush(); pos = $1 }
{
    split($4, a, "|")
    k = $2 SUBSEP $3 SUBSEP a[2]
    s = a[3] "|" a[4] "|" a[5] "|" a[6]
    n[k]++
    if (!(k in f)) f[k] = s; else if (f[k] != s) d[k] = 1
}
function flush(  k, p, sym) {
    for (k in n) if (n[k] > 1) {
        split(k, p, SUBSEP); sym = p[3]
        c[sym]++
        if (!(sym in lo)) lo[sym] = pos
        hi[sym] = pos
        if (n[k] > mx[sym]) mx[sym] = n[k]
        if (k in d) dd[sym]++
    }
    delete n; delete f; delete d
}
END {
    flush()
    for (sym in c) print kind, contig, sym, c[sym], lo[sym], hi[sym], mx[sym], dd[sym] + 0
}
"""


def main():
    b = get_batch('Scan SpliceAI VCFs for repeated (variant, gene) records')
    outputs = []
    for kind in ['snv', 'indel']:
        path = reference_path(f'spliceai_resources/splice_ai_{kind}s')
        for contig in CONTIGS:
            j = b.new_job(f'scan {kind} chr{contig}')
            j.image(image_path('cpg_workflows'))
            j.cpu(2)
            authenticate_cloud_credentials_in_job(j)
            j.command(
                'export GCS_OAUTH_TOKEN=$(gcloud auth print-access-token) && '
                f"bcftools query -r {contig} -f '%POS\\t%REF\\t%ALT\\t%INFO/SpliceAI\\n' {path}"
                f" | awk -v kind={kind} -v contig={contig} '{AWK}' > {j.ofile}"
                f' && cat {j.ofile}'
            )
            outputs.append(j.ofile)
    gather = b.new_job('gather')
    gather.image(image_path('cpg_workflows'))
    gather.command(
        'echo -e "kind\\tcontig\\tsymbol\\trepeated_keys\\tfirst_pos\\tlast_pos\\tmax_copies\\tdiffering_scores"; '
        + 'cat ' + ' '.join(outputs)
    )
    b.run(wait=False)


if __name__ == '__main__':
    main()
