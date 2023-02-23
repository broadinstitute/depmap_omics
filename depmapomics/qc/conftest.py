from depmapomics.qc.config import NEW_RELEASE, PREV_RELEASE


def pytest_sessionstart(session):
    print(
        "REFERENCE_RELEASE: {}.{}".format(PREV_RELEASE["name"], PREV_RELEASE["version"])
    )
    print("VIRTUAL_RELEASE: {}.{}".format(NEW_RELEASE["name"], NEW_RELEASE["version"]))
