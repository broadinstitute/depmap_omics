import os
from depmapomics.config_global import *

if os.getenv("DEPMAP_ENV") == "PROD":
    from depmapomics.config_prod import *
else:
    from depmapomics.config_dev import *
