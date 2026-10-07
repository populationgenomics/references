"""
Tests for download.fetch against a local HTTP server (no internet).

From the repo root:

    uv run --no-project --with pytest pytest gwas_sumstats/test_download.py
"""

import hashlib
import threading
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer

import pytest

import download

BODY = bytes(range(256)) * 4  # 1024 bytes


class Handler(BaseHTTPRequestHandler):
    """
    /whole      the full body with its Content-Length
    /short      Content-Length 1024, then the connection closes after 500 bytes
    /chunked    chunked encoding, closed before the terminating chunk
    /flaky      short on the first request, whole after that
    """

    requests = 0

    def do_GET(self):  # noqa: N802
        Handler.requests += 1
        self.send_response(200)
        if self.path == '/chunked':
            self.send_header('Transfer-Encoding', 'chunked')
            self.end_headers()
            self.wfile.write(b'%x\r\n%s\r\n' % (500, BODY[:500]))
        else:
            self.send_header('Content-Length', str(len(BODY)))
            self.end_headers()
            short = self.path == '/short' or (
                self.path == '/flaky' and Handler.requests == 1
            )
            self.wfile.write(BODY[:500] if short else BODY)
        self.close_connection = True

    def log_message(self, *args):
        pass


@pytest.fixture
def server(monkeypatch):
    monkeypatch.setattr(download.time, 'sleep', lambda seconds: None)
    Handler.requests = 0
    httpd = ThreadingHTTPServer(('127.0.0.1', 0), Handler)
    threading.Thread(target=httpd.serve_forever, daemon=True).start()
    yield f'http://127.0.0.1:{httpd.server_address[1]}'
    httpd.shutdown()


def test_whole_body_is_accepted(server, tmp_path):
    md5, size = download.fetch(f'{server}/whole', tmp_path / 'out')
    assert (md5, size) == (hashlib.md5(BODY).hexdigest(), len(BODY))  # noqa: S324


@pytest.mark.parametrize('path', ['/short', '/chunked'])
def test_truncated_body_is_retried_then_refused(server, tmp_path, path):
    with pytest.raises(OSError, match='bytes|IncompleteRead|incomplete'):
        download.fetch(f'{server}{path}', tmp_path / 'out')
    assert Handler.requests == download.RETRIES
    assert not (tmp_path / 'out').exists()


def test_truncated_once_then_whole_succeeds(server, tmp_path):
    _, size = download.fetch(f'{server}/flaky', tmp_path / 'out')
    assert size == len(BODY)
    assert Handler.requests == 2


def test_md5_mismatch_is_retried_then_refused(server, tmp_path):
    with pytest.raises(OSError, match='MD5'):
        download.fetch(f'{server}/whole', tmp_path / 'out', expected_md5='0' * 32)
    assert Handler.requests == download.RETRIES
    assert not (tmp_path / 'out').exists()


def test_original_name_survives_a_build_correction():
    row = {
        'file_id': '2017_Wheeler_PLoSMed_HbA1c_SAS_GCST007951',
        'source_kind': 'catalog_original',
        'source_build': 'GRCh37',
        'suffix': '.txt.gz',
    }
    before = download.original_name(row)
    assert before == '2017_Wheeler_PLoSMed_HbA1c_SAS_GCST007951_original.txt.gz'
    assert download.original_name(row | {'source_build': 'NCBI36'}) == before


@pytest.mark.parametrize(
    'body, expected',
    [
        (b'h\n1\n2\n', 'h\n1\n2\n'),  # complete file: every line kept
        (b'h\n1\n2', 'h\n1\n'),  # cut mid-line: only the partial line dropped
    ],
)
def test_fetch_head_drops_only_a_partial_last_line(monkeypatch, body, expected):
    import io

    monkeypatch.setattr(
        download, 'open_url', lambda url, headers=None: io.BytesIO(body)
    )
    assert download.fetch_head('u') == expected
