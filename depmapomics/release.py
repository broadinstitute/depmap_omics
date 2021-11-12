# release.py
import os
import yaml
from genepy.utils import helper as h
import pandas as pd
from depmapomics.config import *
from gsheets import Sheets
from taigapy import TaigaClient

tc = TaigaClient()

virtual = {}
INFO = {}


def clean():
    global virtual
    global INFO
    virtual = {}
    INFO = {}


def run(release=RELEASE, notebooks=RUN_NOTEBOOKS):
    cmd = ""
    for val in notebooks:
        cmd += "jupyter nbconvert --to notebook --execute " + val + " & "
    cmd = cmd[:-3]
    os.system(cmd)
    make(release=release)
    save(release, notebooks=notebooks)


def findLatestVersion(dataset, approved_only=True):
    highest = 0
    latest_version = 0
    data = tc.get_dataset_metadata(dataset)
    for val in data["versions"]:
        if val["state"] == "approved" or not approved_only:
            if int(val["name"]) > highest:
                highest = int(val["name"])
                latest_version = highest
    if latest_version == 0:
        raise ValueError("could not find a version")
    return data["permanames"][0] + "." + str(latest_version)


def save(release=RELEASE, notebooks=RUN_NOTEBOOKS):
    for val in notebooks:
        rn = val.replace("ipynb", "html")
        os.system(
            "jupyter nbconvert --to html "
            + val
            + " && mv "
            + rn
            + " release_notebooks/"
            + release
            + "/"
        )

    os.system(
        "cd ../ccle_processing && git add . && git commit -m "
        + '"depmap omics '
        + release
        + ' final" && git push'
    )


def make(
    release=RELEASE,
    my_id=MY_ID,
    mystorage_id=MYSTORAGE_ID,
    potential_list_url=POTENTIAL_LIST,
    prev_virtual=PREV_VIRTUAL,
    eternal_dataset=TAIGA_ETERNAL,
    onlyRelease=["internal", "ibm", "dmc", "public"],
    eternal_include=[
        # CN
        "CCLE_gene_cn",
        "CCLE_segment_cn",
        # Mutations
        "CCLE_mutations",
        "CCLE_mutations_bool_damaging",
        "CCLE_mutations_bool_nonconserving",
        "CCLE_mutations_bool_otherconserving",
        "CCLE_mutations_bool_hotspot",
        # Expression
        "CCLE_expression_full",
        "CCLE_RNAseq_transcripts",
        "CCLE_RNAseq_reads",
        "CCLE_expression",
        "CCLE_expression_proteincoding_genes_expected_count",
        "CCLE_expression_transcripts_expected_count",
        # Fusions
        "CCLE_fusions_unfiltered",
        "CCLE_fusions",
    ],
    updateReadmes=True,
    doCN=True,
    doMut=True,
    doRNA=True,
    doFusions=True,
    changes=CHANGES,
):

    new = {}
    global virtual
    global INFO

    save_output = "temp/" + release

    sheets = Sheets.from_files(my_id, mystorage_id)
    gsheets = sheets.get(potential_list_url).sheets[0].to_frame()
    # making the virtual dataset
    print("loading new samples")
    new["internal"] = set(
        [i for i in gsheets["Internal"].values.tolist() if str(i) != "nan"]
    )
    new["dmc"] = set([i for i in gsheets["DMC"].values.tolist() if str(i) != "nan"])
    new["ibm"] = set([i for i in gsheets["IBM"].values.tolist() if str(i) != "nan"])
    new["public"] = set(
        [i for i in gsheets["Public"].values.tolist() if str(i) != "nan"]
    )

    new["dmc"] = new["dmc"] | new["public"]

    new["ibm"] = new["ibm"] | new["dmc"]

    new["internal"] = new["internal"] | new["ibm"]

    # getting what was released before
    print("loading prev samples")
    prevmut = {}
    prevrna = {}
    prevcn = {}
    prevwes = {}
    prev = {}
    for val in ["internal", "dmc", "ibm", "public"]:
        print(val)
        prevmut[val] = set(
            tc.get(name=prev_virtual[val], file="CCLE_mutations").DepMap_ID
        )
        prevrna[val] = set(tc.get(name=prev_virtual[val], file="CCLE_expression").index)
        prevcn[val] = set(
            tc.get(name=prev_virtual[val], file="CCLE_segment_cn").DepMap_ID
        )
        prev[val] = prevmut[val] | prevrna[val] | prevcn[val]
        prevwes[val] = prevmut[val] | prevcn[val]
    prevmut["dmc"] = prevmut["dmc"] | prevmut["public"]

    prevrna["dmc"] = prevrna["dmc"] | prevrna["public"]
    prevcn["dmc"] = prevcn["dmc"] | prevcn["public"]
    prev["dmc"] = prev["dmc"] | prev["public"]
    prevwes["dmc"] = prevwes["dmc"] | prevwes["public"]

    prevmut["ibm"] = prevmut["ibm"] | prevmut["dmc"]
    prevrna["ibm"] = prevrna["ibm"] | prevrna["dmc"]
    prevcn["ibm"] = prevcn["ibm"] | prevcn["dmc"]
    prev["ibm"] = prev["ibm"] | prev["dmc"]
    prevwes["ibm"] = prevwes["ibm"] | prevwes["dmc"]

    prevmut["internal"] = prevmut["internal"] | prevmut["ibm"]
    prevrna["internal"] = prevrna["internal"] | prevrna["ibm"]
    prevcn["internal"] = prevcn["internal"] | prevcn["ibm"]
    prev["internal"] = prev["internal"] | prev["ibm"]
    prevwes["internal"] = prevwes["internal"] | prevwes["ibm"]

    print("in cn but not mut")
    print(prevwes["internal"] - prevmut["internal"])
    print("in mut but not cn")
    print(prevwes["internal"] - prevcn["internal"])
    print("in rna but no wes/wgs")
    print(prev["internal"] - prevwes["internal"])
    print("in wes/wgs but not rna")
    print(prev["internal"] - prevrna["internal"])

    if updateReadmes:
        print("updating the readmes")
        # managing the readmes
        os.system("cd ../depmap-release-readmes && git pull - -no-commit")
        os.system(
            "cd ../depmap-release-readmes/ && python3 make_new_release.py "
            + release
            + " && git add . && git commit -m "
            + release
            + " && git push "
        )

        # update the reamdes
        folder = "../depmap-release-readmes/" + release
        a = os.listdir(folder)
        print(a)
        for val in a:
            with open(folder + val, "r+") as file:
                # The FullLoader parameter handles the conversion from YAML
                # scalar values to Python the dictionary format
                dict_file = yaml.load(file, Loader=yaml.FullLoader)

                dict_file["extra_sections"][0]["description"] = changes
                yaml.dump(dict_file, file)
        # save the readmes
        # and save them
        os.system(
            "cd ../depmap-release-readmes && cp release-"
            + release
            + "/* ../ccle_processing/readmes/ && git add . && git commit -m "
            + '"Omics: updating readmes to new release" && git push'
        )

    # doing CN
    segmentcn = pd.read_csv("temp/all_" + release + "_segment.csv")
    if doCN:
        genecn = pd.read_csv("temp/all_" + release + "_gene_cn.csv", index_col=0)
        # virtual= {}

    print("looking at dnaseq changes")
    for val in onlyRelease:
        print("_________________________________________________")
        print(val)
        print("removed")
        removed = set(prevcn[val]) - set(segmentcn.DepMap_ID)
        print(removed)
        missing = set(new[val]) - set(segmentcn.DepMap_ID)
        blacklist = set(segmentcn.DepMap_ID) - (prevcn[val] | set(new[val]))
        print("missing")
        print(missing)
        newlines = set(new[val])
        print("blacklist")
        print(len(blacklist), blacklist)

        INFO[val] = (
            "# "
            + val
            + """ dataset:

## DNAseq Omics:

NEW LINES:
"""
            + str(newlines)
            + """

BLACKLIST:
"""
            + str(blacklist)
            + """

MISSING:
"""
            + str(missing)
            + """

REMOVED:
"""
            + str(removed)
        )

        if doCN:
            print("doing CN upload")
            ## for segment removing first blacklisted, then embargoed, to create two datasets
            print(len(segmentcn))
            a = segmentcn[~segmentcn.DepMap_ID.isin(blacklist)]
            print(len(segmentcn) - len(a))
            a.to_csv("temp/all_merged_segments.csv", index=False)
            print(len(genecn))
            a = genecn[~genecn.index.isin(blacklist)]
            print(len(genecn) - len(a))
            a.to_csv("temp/all_merged_genes_cn.csv")

            upload_files = [
                {
                    "path": "temp/all_merged_genes_cn.csv",
                    "name": "CCLE_gene_cn",
                    "format": "NumericMatrixCSV",
                    "encoding": "utf-8",
                },
                {
                    "path": "temp/all_merged_segments.csv",
                    "name": "CCLE_segment_cn",
                    "format": "TableCSV",
                    "encoding": "utf-8",
                },
            ]
            if updateReadmes:
                upload_files += (
                    [
                        {
                            "path": "readmes/" + val + "-" + release + ".yaml",
                            "name": "README",
                            "format": "Raw",
                            "encoding": "utf-8",
                        }
                    ],
                )
            # Add to Taiga
            # 1.(replace create_dataset() by update_dataset() if the dataset already exists)
            # 2.remove folder_id=virtual_folder
            # 3.add a changes_description="something",
            # 4.add add_all_existing_files=True,
            # 5.replace val+"_"+release by the virtual dataset id
            virtual[val] = tc.update_dataset(
                virtual[val],  # val+"_"+release,
                dataset_description=release
                + " relese of the DepMap dataset for the DepMap"
                + " Portal. Please look at the README file for additional information about this dataset. ",
                upload_files=upload_files,
                add_all_existing_files=True,
                changes_description="rerunning the CN upload script",
            )  # folder_id=virtual_folder)

    # mutations
    if doMut:
        print("doing mutations")
        mutations = pd.read_csv(
            "temp/all_somatic_mutations_withlegacy_" + release + "_depmapversion.csv"
        )
        damaging = pd.read_csv(
            "temp/all_somatic_mutations_boolmatrix_fordepmap_damaging.csv", index_col=0
        )
        othercons = pd.read_csv(
            "temp/all_somatic_mutations_boolmatrix_fordepmap_othercons.csv", index_col=0
        )
        othernoncons = pd.read_csv(
            "temp/all_somatic_mutations_boolmatrix_fordepmap_othernoncons.csv",
            index_col=0,
        )
        hotspot = pd.read_csv(
            "temp/all_somatic_mutations_boolmatrix_fordepmap_hotspot.csv", index_col=0
        )
        try:
            mutations = mutations.drop(columns="Unnamed: 0")
        except:
            print("no Unnamed col")
        blacklist = set()
        for val in onlyRelease:
            print("_________________________________________________")
            print(val)
            missing = set(new[val]) - set(mutations.DepMap_ID)
            print("not present")
            removed = set(prev[val]) - set(mutations.DepMap_ID)
            print(removed)
            print("removed")
            removed = set(prevmut[val]) - set(mutations.DepMap_ID)
            print(removed)
            blacklist = (
                set(mutations.DepMap_ID) - (prevmut[val] | set(new[val]))
            ) | blacklist
            print("missing")
            print(missing)
            newlines = set(new[val])
            print("blacklist")
            print(blacklist)

            # adding files
            a = mutations[~mutations.DepMap_ID.isin(blacklist)]
            print(len(mutations) - len(a))
            a.to_csv("temp/all_somatic_mutations_withlegacy.csv", index=False)
            a = damaging[~damaging.index.isin(blacklist)]
            print(len(damaging) - len(a))
            a.to_csv("temp/all_somatic_mutations_boolmatrix_fordepmap_damaging.csv")
            a = othercons[~othercons.index.isin(blacklist)]
            print(len(othercons) - len(a))
            a.to_csv("temp/all_somatic_mutations_boolmatrix_fordepmap_othercons.csv")
            a = othernoncons[~othernoncons.index.isin(blacklist)]
            print(len(othernoncons) - len(a))
            a.to_csv("temp/all_somatic_mutations_boolmatrix_fordepmap_othernoncons.csv")
            a = hotspot[~hotspot.index.isin(blacklist)]
            print(len(hotspot) - len(a))
            a.to_csv("temp/all_somatic_mutations_boolmatrix_fordepmap_hotspot.csv")
            os.popen("cp readmes/" + val + "-" + release + ".txt temp/README").read()

            # updating on taiga
            tc.update_dataset(
                dataset_id=virtual[val],
                changes_description="adding mutations",
                upload_files=[
                    {
                        "path": "temp/all_somatic_mutations_withlegacy.csv",
                        "name": "CCLE_mutations",
                        "format": "TableCSV",
                        "encoding": "utf-8",
                    },
                    {
                        "path": "temp/all_somatic_mutations_boolmatrix_fordepmap_damaging.csv",
                        "name": "CCLE_mutations_bool_damaging",
                        "format": "NumericMatrixCSV",
                        "encoding": "utf-8",
                    },
                    {
                        "path": "temp/all_somatic_mutations_boolmatrix_fordepmap_othernoncons.csv",
                        "name": "CCLE_mutations_bool_nonconserving",
                        "format": "NumericMatrixCSV",
                        "encoding": "utf-8",
                    },
                    {
                        "path": "temp/all_somatic_mutations_boolmatrix_fordepmap_othercons.csv",
                        "name": "CCLE_mutations_bool_otherconserving",
                        "format": "NumericMatrixCSV",
                        "encoding": "utf-8",
                    },
                    {
                        "path": "temp/all_somatic_mutations_boolmatrix_fordepmap_hotspot.csv",
                        "name": "CCLE_mutations_bool_hotspot",
                        "format": "NumericMatrixCSV",
                        "encoding": "utf-8",
                    },
                ],
                add_all_existing_files=True,
            )

    # expression
    genes_tpm = pd.read_csv(
        "temp/expression_" + release + "_genes_tpm_logp1.csv", index_col=0
    )
    if doRNA:
        transcripts_tpm = pd.read_csv(
            "temp/expression_" + release + "_transcripts_tpm_logp1.csv", index_col=0
        )
        genes_expected_count = pd.read_csv(
            "temp/expression_" + release + "_genes_expected_count.csv", index_col=0
        )
        proteincoding_genes_expected_count = pd.read_csv(
            "temp/expression_" + release + "_proteincoding_genes_expected_count.csv",
            index_col=0,
        )
        proteincoding_genes_tpm = pd.read_csv(
            "temp/expression_" + release + "_proteincoding_genes_tpm_logp1.csv",
            index_col=0,
        )
        transcripts_expected_count = pd.read_csv(
            "temp/expression_" + release + "_transcripts_expected_count.csv",
            index_col=0,
        )
        enrichments = pd.read_csv("temp/gene_sets_" + release + "_all.csv", index_col=0)
    print("loading rna changes")
    rnafailed = h.fileToList(save_output + "rna_failed_samples.txt")

    blacklist = set()
    for val in onlyRelease:
        print("_________________________________________________")
        print(val)
        print("not present")
        removed = set(prev[val]) - set(genes_tpm.index)
        print(removed)
        print("removed for QC reasons")
        print(rnafailed)
        print("removed")
        removed = set(prevrna[val]) - set(genes_tpm.index)
        # print(removed - set(rename.keys()))
        missing = set(new[val]) - set(genes_tpm.index)
        blacklist = (set(genes_tpm.index) - (prevrna[val] | set(new[val]))) | blacklist
        print("missing")
        print(missing)
        newlines = set(new[val])
        print("blacklist")
        print(len(blacklist), blacklist)

        INFO[val] += (
            """

## RNAseq Omics:

NEW LINES:
"""
            + str(newlines)
            + """

BLACKLIST:
"""
            + str(blacklist)
            + """

MISSING:
"""
            + str(missing)
            + """

REMOVED:
"""
            + str(removed)
            + """

REMOVED FOR QC REASONS:
"""
            + str(rnafailed)
        )

        if doRNA:
            print("doing expression upload")
            ## removing first blacklisted, then embargoed, to create two datasets
            print(len(genes_expected_count))
            a = genes_expected_count[~genes_expected_count.index.isin(blacklist)]
            print(len(a))
            a.to_csv("temp/expression_genes_expected_count.csv")
            print(len(genes_tpm))
            a = genes_tpm[~genes_tpm.index.isin(blacklist)]
            print(len(a))
            a.to_csv("temp/expression_genes_tpm.csv")
            print(len(proteincoding_genes_tpm))
            a = proteincoding_genes_tpm[~proteincoding_genes_tpm.index.isin(blacklist)]
            print(len(a))
            a.to_csv("temp/expression_proteincoding_genes_tpm.csv")
            print(len(transcripts_tpm))
            a = transcripts_tpm[~transcripts_tpm.index.isin(blacklist)]
            print(len(a))
            a.to_csv("temp/expression_transcripts_tpm.csv")
            print(len(proteincoding_genes_expected_count))
            a = proteincoding_genes_expected_count[
                ~proteincoding_genes_expected_count.index.isin(blacklist)
            ]
            print(len(a))
            a.to_csv("temp/expression_proteincoding_genes_expected_count.csv")
            print(len(transcripts_expected_count))
            a = transcripts_expected_count[
                ~transcripts_expected_count.index.isin(blacklist)
            ]
            print(len(a))
            a.to_csv("temp/expression_transcripts_expected_count.csv")
            print(len(enrichments))
            a = enrichments[~enrichments.index.isin(blacklist)]
            print(len(a))
            a.to_csv("temp/enrichments_ssGSEA.csv")

            # adding to taiga
            tc.update_dataset(
                virtual[val],
                changes_description="adding expression",
                upload_files=[
                    {
                        "path": "temp/expression_genes_expected_count.csv",
                        "name": "CCLE_RNAseq_reads",
                        "format": "NumericMatrixCSV",
                        "encoding": "utf-8",
                    },
                    {
                        "path": "temp/expression_transcripts_tpm.csv",
                        "name": "CCLE_RNAseq_transcripts",
                        "format": "NumericMatrixCSV",
                        "encoding": "utf-8",
                    },
                    {
                        "path": "temp/expression_genes_tpm.csv",
                        "name": "CCLE_expression_full",
                        "format": "NumericMatrixCSV",
                        "encoding": "utf-8",
                    },
                    {
                        "path": "temp/expression_proteincoding_genes_tpm.csv",
                        "name": "CCLE_expression",
                        "format": "NumericMatrixCSV",
                        "encoding": "utf-8",
                    },
                    {
                        "path": "temp/expression_proteincoding_genes_expected_count.csv",
                        "name": "CCLE_expression_proteincoding_genes_expected_count",
                        "format": "NumericMatrixCSV",
                        "encoding": "utf-8",
                    },
                    {
                        "path": "temp/expression_transcripts_expected_count.csv",
                        "name": "CCLE_expression_transcripts_expected_count",
                        "format": "NumericMatrixCSV",
                        "encoding": "utf-8",
                    },
                    {
                        "path": "temp/enrichments_ssGSEA.csv",
                        "name": "CCLE_ssGSEA",
                        "format": "NumericMatrixCSV",
                        "encoding": "utf-8",
                    },
                ],
                upload_async=False,
                add_all_existing_files=True,
            )

    # fusions
    if doFusions:
        print("doing fusions")
        fusions = pd.read_csv("temp/fusions_" + release + ".csv")
        filtered = pd.read_csv("temp/filtered_fusions_" + release + ".csv")

        blacklist = set()
        for val in onlyRelease:
            print("_________________________________________________")
            print(val)
            print("not present")
            removed = set(prev[val]) - set(fusions.DepMap_ID)
            print(removed)
            print("removed for QC reasons")
            print(rnafailed)
            print("removed")
            removed = set(prevrna[val]) - set(fusions.DepMap_ID)
            print(removed)
            missing = set(new[val]) - set(fusions.DepMap_ID)
            blacklist = (
                set(fusions.DepMap_ID) - (prevrna[val] | set(new[val]))
            ) | blacklist
            print("missing")
            print(missing)
            newlines = set(new[val])
            print("blacklist")
            print(len(blacklist), blacklist)
            ## removing first blacklisted, then embargoed, to create two datasets
            print(len(fusions))
            a = fusions[~fusions.DepMap_ID.isin(blacklist)]
            print(len(a))
            a.to_csv("temp/fusions.csv", index=False)
            print(len(filtered))
            a = filtered[~filtered.DepMap_ID.isin(blacklist)]
            print(len(a))
            a.to_csv("temp/filtered_fusions.csv", index=False)

            # uploading to taiga
            tc.update_dataset(
                virtual[val],
                changes_description="adding fusions",
                upload_files=[
                    {
                        "path": "temp/fusions.csv",
                        "name": "CCLE_fusions_unfiltered",
                        "format": "TableCSV",
                        "encoding": "utf-8",
                    },
                    {
                        "path": "temp/filtered_fusions.csv",
                        "name": "CCLE_fusions",
                        "format": "TableCSV",
                        "encoding": "utf-8",
                    },
                ],
                dataset_description=INFO[val],
                add_all_existing_files=True,
            )

    # updating eternal
    print("updating eternal")
    latest_version = findLatestVersion(virtual["internal"])
    tc.update_dataset(
        eternal_dataset,
        changes_description="new " + release + " omics dataset.",
        add_taiga_ids=[
            {"taiga_id": latest_version + "/" + file, "name": file}
            for file in eternal_include
        ],
        add_all_existing_files=True,
    )

    print("think of cleanning your temp/ folder once the release is fully validated...")

