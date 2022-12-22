import os

import pytest


def pytest_addoption(parser):
    parser.addoption("--answer_dir", help="Directory where answers are stored.")
    parser.addoption("--check_dir", help="Directory where spectrum checks are stored.")
    parser.addoption(
        "--answer_store",
        action="store_true",
        help="Generate new answers, but don't test.",
    )


@pytest.fixture()
def answer_store(request):
    return request.config.getoption("--answer_store")


@pytest.fixture()
def answer_dir(request):
    ad = request.config.getoption("--answer_dir")
    if ad is None:
        pytest.skip("No answer directory")
    ad = os.path.abspath(ad)
    if not os.path.exists(ad):
        os.makedirs(ad)
    return ad


@pytest.fixture()
def check_dir(request):
    cd = request.config.getoption("--check_dir")
    if cd is not None:
        cd = os.path.abspath(cd)
        if not os.path.exists(cd):
            os.makedirs(cd)
        cd = os.path.join(cd, f"{request.node.name}.dat")
    return cd
