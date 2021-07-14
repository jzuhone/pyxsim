import pytest
import os


def pytest_addoption(parser):
    parser.addoption("--answer_dir",
                     help="Directory where answers are stored.")
    parser.addoption("--answer_store", action="store_true",
                     help="Generate new answers, but don't test.")


@pytest.fixture()
def answer_store(request):
    return request.config.getoption('--answer_store')


@pytest.fixture()
def answer_dir(request):
    ad = os.path.abspath(request.config.getoption('--answer_dir'))
    if not os.path.exists(ad):
        os.makedirs(ad)
    return ad
