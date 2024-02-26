from taigapy import TaigaClient

# there are some issues in the older versions of omics data on virtual that this flag deals with
# some columns got renamed in the mutation file and some data was uploaded as tsv instead of csv
LEGACY_PATCH_FLAGS = {"rename_column": True, "tsv2csv": False}

# release ids on taiga
tc = TaigaClient()


def taiga_latest_version(dataset_name):
    return int(tc.get_dataset_metadata(dataset_name)["versions"][-1]["name"])


def taiga_latest_path(dataset_name):
    return {
        "name": dataset_name,
        "version": taiga_latest_version(dataset_name),
    }


# from depmapomics.qc.create_dataset import update_tentative_virtual
# update_tentative_virtual()
TENTATIVE_VIRTUAL = taiga_latest_path("tentative-virtual-d84e")

VIRTUAL_RELEASES = {
    "21Q4v2": {
        "internal": taiga_latest_path("internal-21q4v2-403b"),
        "ibm": taiga_latest_path("ibm-21q4v2-2d92"),
        "public": taiga_latest_path("public-21q4v2-103d"),
        "dmc": taiga_latest_path("dmc-21q4v2-0e7c"),
    },
    "21Q3": {
        "internal": {"name": "internal-21q3-fe4c", "version": 12},
        "ibm": {"name": "ibm-21q3-179f", "version": 8},
        "dmc": {"name": "dmc-21q3-482c", "version": 7},
        "public": {"name": "public-21q3-bf1e", "version": 7},
    },
    "21Q4": {
        "internal": taiga_latest_path("internal-21q4-ac0a"),
        "ibm": taiga_latest_path("ibm-21q4-4e18"),
        "dmc": taiga_latest_path("dmc-21q4-5725"),
        "public": taiga_latest_path("public-21q4-a0d6"),
    },
    "22Q1": {
        "internal": taiga_latest_path("internal-22q1-1778"),
        "ibm": taiga_latest_path("ibm-22q1-cce1"),
        "dmc": taiga_latest_path("dmc-22q1-d00a"),
        "public": taiga_latest_path("public-22q1-305b"),
    },
    "22Q2": {
        "internal": taiga_latest_path("internal-22q2-097a"),
        "ibm": taiga_latest_path("ibm-22q2-a71e"),
        "dmc": taiga_latest_path("dmc-22q2-5e51"),
        "public": taiga_latest_path("public-22q2-de04"),
    },
    "22Q2dryrun": {
        "internal": taiga_latest_path("internal-22q2dryrun-5ffb"),
        "ibm": taiga_latest_path("ibm-22q2dryrun-ddbe"),
        "dmc": taiga_latest_path("dmc-22q2dryrun-0f8d"),
        "public": taiga_latest_path("public-22q2dryrun-e259"),
    },
    "22Q2renamed": {
        "internal": taiga_latest_path("internal-22q2-08ee"),
    },
    "22Q4": {
        "internal": taiga_latest_path("internal-22q4-56d4"),
        "dmc": taiga_latest_path("dmc-22q4-a73a"),
        "public": taiga_latest_path("public-22q4-6837"),
    },
    "23Q2": {
        "internal": taiga_latest_path("internal-23q2-1e49"),
        "dmc": taiga_latest_path("dmc-23q2-d87c"),
        "public": taiga_latest_path("public-23q2-19de"),
    },
    "23Q4": {
        "internal": taiga_latest_path("internal-23q4-ac2b"),
        "dmc": taiga_latest_path("dmc-23q4-32b7"),
        "public": taiga_latest_path("public-23q4-d816"),
    },
}  # release ids on taiga

PORTALS = ["dmc", "public", "internal"]  # used for 'bookkeeping' markers
PORTAL = "internal"  # used for 'not bookkeeping' markers
PREV_QUARTER = "23Q2"
NEW_QUARTER = "23Q4"

PREV_RELEASE = VIRTUAL_RELEASES[PREV_QUARTER][PORTAL]
NEW_RELEASE = VIRTUAL_RELEASES[NEW_QUARTER][PORTAL]
# NEW_RELEASE = TENTATIVE_VIRTUAL
PORTALS = [PORTAL]

LINES_TO_DROP = {"DNA": set([]), "RNA": set([])}
LINES_TO_RELEASE = {}
LINES_TO_RELEASE["public"] = set([])
LINES_TO_RELEASE["dmc"] = set([])
LINES_TO_RELEASE["internal"] = set([])
IGNORE_FAILED_TO_RELEASE = []


# these are the columns that if merged with an older release (assuming that old data was not altered),
# should uniquely identify each row of the file to find equal values in each column
FUSIONS_FILTERED_MERGE_COLS = [
    "ModelID",
    "LeftGene",
    "RightGene",
    "LeftBreakpoint",
    "RightBreakpoint",
    "SpliceType",
]
FUSIONS_UNFILTERED_MERGE_COLS = [
    "ProfileID",
    "LeftGene",
    "RightGene",
    "LeftBreakpoint",
    "RightBreakpoint",
    "SpliceType",
]
SEGMENT_CN_MERGE_COLS = ["ProfileID", "Chromosome", "Start", "End"]
MUTATIONS_MERGE_COLS = [
    "ModelID",
    "Chromosome",
    "Start_position",
    "End_position",
    "Alternate_Allele",
]

# if there are new files that were added in the previous release, add them here
FILES_RELEASED_BEFORE = [
    "OmicsExpressionProteinCodingGenesTPMLogp1",
    "OmicsExpressionGenesExpectedCountProfile"
    "OmicsExpressionTranscriptsExpectedCountProfile",
    "OmicsExpressionGeneSetEnrichment",
    "OmicsExpressionGeneSetEnrichmentProfile",
    "OmicsFusionFiltered",
    "OmicsFusionUnfilteredProfile",
    "OmicsCNGene",
    "OmicsCNSegmentsProfile",
    "OmicsSomaticMutations",
    "OmicsSomaticMutationsProfile",
    "OmicsSomaticMutationsMatrixDamaging",
    "OmicsSomaticMutationsMatrixHotspot",
]

# correlation thresholds above which we consider two releases as 'similar'
CORRELATION_THRESHOLDS = {"OmicsCNGene": 0.8, "all_expressions": 0.98}

SKIP_ARXSPAN_COMPARISON = (
    True  # set to False if you want to test whether some arxspans were added/removed
)

PLOTS_OUTPUT_FILENAME_PREFIX = "/tmp/plots_"  # location/prefix for saving output plots

# all the file attributes
FILE_ATTRIBUTES = [
    # all model level
    {
        "file": "OmicsExpressionProteinCodingGenesTPMLogp1",
        "ismatrix": True,
        "hasNA": False,
        "gene_id": "entrez",
        "omicssource": "RNA",
        "id": "ModelID",
    },
    {
        "file": "OmicsExpressionGenesExpectedCountProfile",
        "ismatrix": True,
        "hasNA": False,
        "gene_id": "ensg",
        "omicssource": "RNA",
        "id": "ProfileID",
    },
    {
        "file": "OmicsExpressionTranscriptsExpectedCountProfile",
        "ismatrix": True,
        "hasNA": False,
        "gene_id": "enst",
        "omicssource": "RNA",
        "id": "ProfileID",
    },
    {
        "file": "OmicsFusionFiltered",
        "ismatrix": False,
        "omicssource": "RNA",
        "merge_cols": FUSIONS_FILTERED_MERGE_COLS,
        "expected_changed_cols": ["CCLE_count"],
        "id": "ModelID",
    },
    {
        "file": "OmicsFusionUnfilteredProfile",
        "ismatrix": False,
        "omicssource": "RNA",
        "merge_cols": FUSIONS_UNFILTERED_MERGE_COLS,
        "expected_changed_cols": ["CCLE_count"],
        "id": "ProfileID",
    },
    {
        "file": "OmicsExpressionGeneSetEnrichment",
        "ismatrix": True,
        "hasNA": False,
        "omicssource": "RNA",
        "gene_id": None,
        "id": "ModelID",
    },
    {
        "file": "OmicsExpressionGeneSetEnrichmentProfile",
        "ismatrix": True,
        "hasNA": False,
        "omicssource": "RNA",
        "gene_id": None,
        "id": "ProfileID",
    },
    {
        "file": "OmicsCNGene",
        "ismatrix": True,
        "hasNA": True,
        "gene_id": "entrez",
        "omicssource": "DNA",
        "id": "ModelID",
    },
    {
        "file": "OmicsCNSegmentsProfile",
        "ismatrix": False,
        "omicssource": "DNA",
        "merge_cols": SEGMENT_CN_MERGE_COLS,
        "expected_changed_cols": [],
        "id": "ProfileID",
    },
    {
        "file": "OmicsSomaticMutations",
        "ismatrix": False,
        "omicssource": "DNA",
        "merge_cols": MUTATIONS_MERGE_COLS,
        "expected_changed_cols": [],
        "id": "ModelID",
    },
    {
        "file": "OmicsSomaticMutationsProfile",
        "ismatrix": False,
        "omicssource": "DNA",
        "merge_cols": MUTATIONS_MERGE_COLS,
        "expected_changed_cols": [],
        "id": "ProfileID",
    },
    {
        "file": "OmicsSomaticMutationsMatrixDamaging",
        "ismatrix": True,
        "hasNA": False,
        "gene_id": "entrez",
        "omicssource": "DNA",
        "id": "ModelID",
    },
    {
        "file": "OmicsSomaticMutationsMatrixHotspot",
        "ismatrix": True,
        "hasNA": False,
        "gene_id": "entrez",
        "omicssource": "DNA",
        "id": "ModelID",
    },
]

# comment/uncomment to use all/subset of files for testing
FILE_ATTRIBUTES = [
    x
    for x in FILE_ATTRIBUTES
    if (x["file"] in ["OmicsCNSegmentsProfile", "OmicsCNGene"])
]
# FILE_ATTRIBUTES = [x for x in FILE_ATTRIBUTES if (x["file"] in ["CCLE_mutations"])]
# FILE_ATTRIBUTES = [
#     x for x in FILE_ATTRIBUTES if (x["omicssource"] in ["RNA"]) and x["ismatrix"]
# ]
# FILE_ATTRIBUTES = [
#     x
#     for x in FILE_ATTRIBUTES
#     if (x["file"] in ["OmicsFusionFiltered", "OmicsFusionUnfilteredProfile"])
# ]
# FILE_ATTRIBUTES = [x for x in FILE_ATTRIBUTES if (x["file"] == "OmicsFusionFiltered")]

# the following information is used to create a tentative virtual
MUTATIONS_TAIGA_ID = "mutations-latest-ed72"
FUSIONS_TAIGA_ID = "fusions-95c9"
EXPRESSION_TAIGA_ID = "expression-d035"
CN_TAIGA_ID = "cn-achilles-version-06ca"
# CN_TAIGA_ID = 'cn-latest-d8d4'

TAIGA_IDS_LATEST = {
    MUTATIONS_TAIGA_ID: [
        ("OmicsSomaticMutationsProfile", "somaticMutations_profile"),
    ],
    FUSIONS_TAIGA_ID: [
        ("OmicsFusionUnfilteredProfile", "fusions_unfiltered_profile"),
    ],
    EXPRESSION_TAIGA_ID: [
        ("OmicsExpressionGenesExpectedCountProfile", "genes_expectedCount_profile"),
        (
            "OmicsExpressionTranscriptsExpectedCountProfile",
            "transcripts_expectedCount_profile",
        ),
        ("OmicsExpressionGeneSetEnrichmentProfile", "gene_set_enrichment_profile"),
    ],
    CN_TAIGA_ID: [
        ("OmicsCNSegmentsProfile", "merged_segments_profile"),
    ],
}
