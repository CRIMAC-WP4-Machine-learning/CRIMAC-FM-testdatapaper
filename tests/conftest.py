import pytest
from helpers import run


# List existing datasets
result = run(["list"])
_TEST_DATA_SETS = result.stdout.strip().splitlines()[-1].split()


@pytest.fixture(
    params=_TEST_DATA_SETS,
    ids=_TEST_DATA_SETS,
)
def dataset_id(request):
    """Return one available test dataset."""
    return request.param

