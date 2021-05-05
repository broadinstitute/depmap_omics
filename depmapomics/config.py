CACHE_PATH = '~/.depmapomics/'
TMP_PATH = '/tmp/'
ENSEMBL_SERVER_V = "http://nov2020.archive.ensembl.org/biomart"

RNASEQC_THRESHOLDS_LOWQUAL = {'minmapping': 0.85, 'minendmapping': 0.75, 'minefficiency': 0.75,
                              'maxendmismatch': 0.02, 'maxmismatch': 0.02, 'minhighqual': 0.8,
                              'minexon': 0.7, "maxambiguous": 0.05, "maxsplits": 0.1,
                              "maxalt": 0.2, "maxchim": 0.05, "minreads": 20000000,
                              "minlength": 80, "maxgenes": 35000, "mingenes": 12000}


RNASEQC_THRESHOLDS_FAILED = {'minmapping': 0.7, 'minendmapping': 0.66, 'minefficiency': 0.6,
                             'maxendmismatch': 0.02, 'maxmismatch': 0.02, 'minhighqual': 0.7,
                             'minexon': 0.66, "maxambiguous": 0.1, "maxsplits": 0.1,
                             "maxalt": 0.5, "maxchim": 0.2, "minreads": 20000000,
                             "minlength": 80, "maxgenes": 35000, "mingenes": 10000}
