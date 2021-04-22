from depmapomics.tests.config import VIRTUAL_RELEASE, REFERENCE_RELEASE

def pytest_sessionstart(session):
    """
    Called after the Session object has been created and
    before performing collection and entering the run test loop.
    """
    print('REFERENCE_RELEASE: ', REFERENCE_RELEASE)
    print('VIRTUAL_RELEASE: ', VIRTUAL_RELEASE)
