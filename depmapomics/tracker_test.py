import pandas as pd
import numpy as np
from depmapomics.config import *
from depmapomics import tracker as track
from depmapomics import expressions


def test_expression_outputs():
    trackerobj = track.initTracker()
    testset_name = "test_postprocessing"

    updated_tracker = expressions._CCLEPostProcessing(
        samplesetname=SAMPLESETNAME,
        trackerobj=trackerobj,
        samplesetToLoad=testset_name,
        dry_run=True,
    )

    return True
