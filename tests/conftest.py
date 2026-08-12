import netrc
import os
import pytest
from pathlib import Path


def _has_credentials(host: str) -> bool:
    """Return True if ~/.netrc contains an entry for *host*."""
    usrhome = os.path.expanduser("~")
    for fname in (".netrc", "_netrc"):
        netrc_path = os.path.join(usrhome, fname)
        if os.path.isfile(netrc_path):
            try:
                auth = netrc.netrc(netrc_path).authenticators(host)
                return auth is not None
            except Exception:
                pass
    return False


def pytest_collection_modifyitems(config, items):
    """Skip tests marked need_credentials when no credentials file is found."""
    skip_marker = pytest.mark.skip(
        reason="requires credentials (netrc); run with -m need_credentials to force"
    )
    for item in items:
        if item.get_closest_marker("need_credentials"):
            if not _has_credentials("nrt.cmems-du.eu"):
                item.add_marker(skip_marker)


@pytest.fixture
def test_data():
    return Path(__file__).parent.resolve().joinpath('data')

