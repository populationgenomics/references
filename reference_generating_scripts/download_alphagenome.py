"""Download AlphaGenome atlas files and upload to GCS.

The DeepMind download page requires a per-link session cookie.
For each link, click it in your browser, open DevTools > Network,
right-click the request > "Copy as cURL", and extract the Cookie
header value.

    python data_management_scripts/download_alphagenome.py \
        --cookie-avi 'NID=...; SID=...' \
        --cookie-splicing 'NID=...; SID=...' \
        --cookie-importances 'NID=...; SID=...'
"""

import argparse

from cpg_utils import hail_batch

DEST = 'gs://cpg-common-main/references/alphagenome'

DOWNLOADS = [
    ('avi_scores_snvs_tabix.zip', 'cookie_avi'),
    ('combined_splicing_snvs_tabix.zip', 'cookie_splicing'),
    ('avi_feature_importances_snvs_tabix.zip', 'cookie_importances'),
]

BASE_URL = 'https://deepmind.google.com/science/alphagenome/_/download/atlas'

CURL_HEADERS = (
    "-H 'User-Agent: Mozilla/5.0 (Macintosh; Intel Mac OS X 10.15; rv:153.0) Gecko/20100101 Firefox/153.0' "
    "-H 'Accept: text/html,application/xhtml+xml,application/xml;q=0.9,*/*;q=0.8' "
    "-H 'Referer: https://deepmind.google.com/'"
)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--cookie-avi', required=True)
    parser.add_argument('--cookie-splicing', required=True)
    parser.add_argument('--cookie-importances', required=True)
    args = parser.parse_args()

    batch = hail_batch.get_batch(name='Download AlphaGenome atlas files')

    for filename, cookie_attr in DOWNLOADS:
        cookie = getattr(args, cookie_attr)
        url = f'{BASE_URL}/{filename}'

        j = batch.new_bash_job(f'download-{filename}')
        j.storage('500Gi')

        j.command(
            f'curl -L -f {CURL_HEADERS} '
            f'-H \'Cookie: {cookie}\' '
            f'-o {j.output} \'{url}\''
        )
        batch.write_output(j.output, f'{DEST}/{filename}')

    batch.run(wait=False)


if __name__ == '__main__':
    main()
