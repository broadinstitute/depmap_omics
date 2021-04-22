from depmapomics.tests.config import VIRTUAL_RELEASE, REFERENCE_RELEASE

def pytest_sessionstart(session):
    print('REFERENCE_RELEASE: ', REFERENCE_RELEASE)
    print('VIRTUAL_RELEASE: ', VIRTUAL_RELEASE)
