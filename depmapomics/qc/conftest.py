from depmapomics.qc.config import NEW_RELEASE, PREV_RELEASE, PORTALS


def pytest_sessionstart(session):
    for portal in PORTALS:
        print(portal, ":")
        print(
            "REFERENCE_RELEASE: {}.{}".format(
                PREV_RELEASE[portal]["name"], PREV_RELEASE[portal]["version"]
            )
        )
        print(
            "VIRTUAL_RELEASE: {}.{}".format(
                NEW_RELEASE[portal]["name"], NEW_RELEASE[portal]["version"]
            )
        )
