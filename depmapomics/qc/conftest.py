from depmapomics import constants
from depmapomics import env_config
from depmapomics.qc.config import NEW_constants.RELEASE, PREV_constants.RELEASE

def pytest_sessionstart(session):
    print('REFERENCE_constants.RELEASE: {}.{}'.format(PREV_constants.RELEASE['name'], PREV_constants.RELEASE['version']))
    print('constants.VIRTUAL_constants.RELEASE: {}.{}'.format(NEW_constants.RELEASE['name'], NEW_constants.RELEASE['version']))
