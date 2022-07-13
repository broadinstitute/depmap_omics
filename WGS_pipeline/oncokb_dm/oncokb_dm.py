import sys
from cravat import BaseAnnotator
from cravat import InvalidData
import sqlite3
import os
import requests
import json


class CravatAnnotator(BaseAnnotator):
    def setup(self):
        dir_path = os.path.dirname(os.path.realpath(__file__))
        datafile_path = os.path.join(dir_path, "data", "oncokb.txt")
        if os.path.exists(datafile_path) == False:
            datafile_path = os.path.join(dir_path, "data", "token.txt")
        if os.path.exists(datafile_path) == False:
            datafile_path = None
        if datafile_path is not None:
            with open(datafile_path) as f:
                for l in f:
                    self.token = l.strip()
        else:
            self.token = None

        datafile_path = os.path.join(dir_path, "data", "oncokb_annotated_genes.txt")
        with open(datafile_path) as f:
            self.oncokb_genes = f.read().splitlines()

        self.output_columns = self.conf["output_columns"]
        self.make_json_colnames()

    def make_json_colnames(self):
        self.json_colnames = []
        for col in self.output_columns:
            if "table" in col and col["table"] == True:
                self.json_colnames.append(col["name"])

    def fill_empty_output(self, output_dict):
        for output_col in self.conf["output_columns"]:
            col_name = output_col["name"]
            if col_name not in output_dict:
                output_dict[col_name] = ""
        return output_dict

    def annotate(self, input_data):
        self.uid = input_data["uid"]
        chrom = input_data["chrom"]
        pos = input_data["pos"]
        ref = input_data["ref_base"]
        alt = input_data["alt_base"]
        hugo = input_data["hugo"]
        return {"chrom": chrom, "pos": pos, "ref": ref, "alt": alt, "hugo": hugo}

    def _get_hgvs_g(self, chrom, pos, ref, alt):
        hgvs_g = f"{chrom}:g."
        hgvs_g += self._get_hgvs_nuc(pos, ref, alt)
        return hgvs_g

    def _get_hgvs_nuc(self, pos, ref, alt):
        hgvs_nuc = ""
        start = pos
        end = str(int(pos) + len(ref) - 1)
        if ref == "-":
            hgvs_nuc += "%s_%sins%s" % (str(int(start) - 1), start, alt)
        elif alt == "-":
            if len(ref) == 1:
                hgvs_nuc += "%sdel" % start
            else:
                hgvs_nuc += "%s_%sdel" % (start, end)
        else:
            if len(ref) == 1 and len(alt) == 1:
                hgvs_nuc += "%s%s>%s" % (start, ref, alt)
            else:
                hgvs_nuc += "%s_%sdelins%s" % (start, end, alt)
        return hgvs_nuc

    def handle_jsondata(self, output_dict):
        for colname in self.json_colnames:
            json_data = output_dict.get(colname, None)
            if json_data is not None:
                json_data = json.dumps(json_data)
            output_dict[colname] = json_data
        return output_dict

    def postprocess(self):
        if self.token is None:
            print("no token found!!")
            return
        batch = []
        headers = {
            "accept": "application/json",
            "Authorization": "Bearer " + self.token,
            "Content-Type": "application/json",
        }
        max_lnum = self.uid - 1
        uids = []
        keys = {}
        count = 0
        reqd = True
        for lnum, line, input_data, secondary_data in self._get_input():
            try:
                chrom = input_data["chrom"]
                pos = str(input_data["pos"])
                ref = input_data["ref_base"]
                alt = input_data["alt_base"]
                uid = input_data["uid"]
                if "hugo" in input_data and input_data["hugo"] is not None:
                    if input_data["hugo"] in self.oncokb_genes:
                        uids.append(uid)
                        count = count + 1
                        hgvs_g = self._get_hgvs_g(chrom, pos, ref, alt)
                        batch.append(hgvs_g)
                        reqd = False
                datas = ""
                if (count % 100 == 0 and lnum != 0 or lnum == max_lnum) and not reqd:
                    data = "[ "
                    for b in batch:
                        data = (
                            data
                            + '{ "evidenceTypes": ["MUTATION_SUMMARY", "TUMOR_TYPE_SUMMARY", "PROGNOSTIC_SUMMARY", "DIAGNOSTIC_SUMMARY", "ONCOGENIC", "MUTATION_EFFECT", "PROGNOSTIC_IMPLICATION", "DIAGNOSTIC_IMPLICATION", "STANDARD_THERAPEUTIC_IMPLICATIONS_FOR_DRUG_SENSITIVITY", "STANDARD_THERAPEUTIC_IMPLICATIONS_FOR_DRUG_RESISTANCE"], "hgvsg":'
                            + '"'
                            + b
                            + '"'
                            + ', "id": "", "referenceGenome": "GRCh38"}, '
                        )
                    data = data[:-2] + "]"
                    response = requests.post(
                        "https://www.oncokb.org/api/v1/annotate/mutations/byHGVSg",
                        headers=headers,
                        data=data,
                    )
                    datas = response.json()
                    batch = []
                    reqd = True

                if len(datas) > 0:
                    do = True
                    for n, x in enumerate(datas):
                        oncogenic = x["oncogenic"]
                        mutaff = x["mutationEffect"]

                        knownEffect = mutaff["knownEffect"]
                        citations = mutaff["citations"]
                        if len(citations["pmids"]) < 1:
                            pmids = ""
                        pmids = "; ".join(citations["pmids"])
                        highestSensitiveLevel = x["highestSensitiveLevel"]
                        if highestSensitiveLevel:
                            highestSensitiveLevel = highestSensitiveLevel.replace(
                                "LEVEL_", ""
                            )
                        highestResistanceLevel = x["highestResistanceLevel"]
                        if highestResistanceLevel:
                            highestResistanceLevel = highestResistanceLevel.replace(
                                "LEVEL_R", ""
                            )
                        highestDiagnosticImplicationLevel = x[
                            "highestDiagnosticImplicationLevel"
                        ]
                        if highestDiagnosticImplicationLevel:
                            highestDiagnosticImplicationLevel = highestDiagnosticImplicationLevel.replace(
                                "LEVEL_Dx", ""
                            )
                        highestPrognosticImplicationLevel = x[
                            "highestPrognosticImplicationLevel"
                        ]
                        if highestPrognosticImplicationLevel:
                            highestPrognosticImplicationLevel = highestPrognosticImplicationLevel.replace(
                                "LEVEL_Px", ""
                            )
                        hotspot = x["hotspot"]
                        variantSummary = x["variantSummary"]
                        tumorSummary = x["tumorTypeSummary"]
                        diagnosticImplications = x["diagnosticImplications"]
                        precomp_data = []
                        for i in range(len(diagnosticImplications)):
                            dd = diagnosticImplications[i]
                            level = dd["levelOfEvidence"].replace("LEVEL_Dx", "")
                            diagnostic_data = [
                                level,
                                dd["tumorType"]["mainType"]["name"],
                                dd["tumorType"]["mainType"]["tumorForm"],
                                dd["pmids"],
                            ]
                            precomp_data.append([{"diagnostic_data": diagnostic_data}])
                        treatments = x["treatments"]
                        for i in range(len(treatments)):
                            tt = treatments[i]
                            drugs = tt["drugs"]
                            code = [x["ncitCode"] for x in drugs]
                            drugname = [x["drugName"] for x in drugs]
                            approved_indications = tt["approvedIndications"]
                            treatment_pmids = tuple(tt["pmids"])
                            levelAssociatedCancerType = tt["levelAssociatedCancerType"]
                            levelAssociatedCancerType_level = levelAssociatedCancerType[
                                "level"
                            ]
                            cancer_name = levelAssociatedCancerType["name"]
                            cancer_tissue = levelAssociatedCancerType["tissue"]
                            cancer_tumor = levelAssociatedCancerType["tumorForm"]
                            treatment_data = [
                                code,
                                drugname,
                                approved_indications,
                                treatment_pmids,
                                levelAssociatedCancerType_level,
                                cancer_name,
                                cancer_tissue,
                                cancer_tumor,
                            ]
                            precomp_data.append([{"treatment_data": treatment_data}])
                        prognosticImplications = x["prognosticImplications"]
                        for p in range(len(prognosticImplications)):
                            prog = prognosticImplications[p]
                            progLevelOfEvidence = prog["levelOfEvidence"].replace(
                                "LEVEL_Px", ""
                            )
                            progTumorType = prog["tumorType"]["mainType"]["name"]
                            progTumorForm = prog["tumorType"]["mainType"]["tumorForm"]
                            progTissue = prog["tumorType"]["tissue"]
                            progPmids = prog["pmids"]
                            prognostic_data = [
                                progLevelOfEvidence,
                                progTumorType,
                                progTumorForm,
                                progTissue,
                                progPmids,
                            ]
                            precomp_data.append([{"prognostic_data": prognostic_data}])
                        if len(precomp_data) < 1:
                            precomp_data = None
                        output_dict = {
                            "oncogenic": oncogenic if oncogenic != "Unknown" else None,
                            "knownEffect": knownEffect
                            if knownEffect != "Unknown"
                            else None,
                            "pmids": pmids,
                            "highestSensitiveLevel": highestSensitiveLevel,
                            "highestResistanceLevel": highestResistanceLevel,
                            "highestDiagnosticImplicationLevel": highestDiagnosticImplicationLevel,
                            "highestPrognosticImplicationLevel": highestPrognosticImplicationLevel,
                            "hotspot": hotspot if hotspot else None,
                            "tumorSummary": tumorSummary,
                            "variantSummary": None
                            if variantSummary.endswith("is unknown.")
                            | variantSummary.endswith("likely oncogenic.")
                            else variantSummary,
                            "all": precomp_data,
                        }
                        output_dict = self.handle_jsondata(output_dict)
                        if n in keys:
                            keys[n] = keys[n] + 100
                            output_dict["uid"] = uids[keys[n]]
                        else:
                            output_dict["uid"] = uids[n]
                            keys[n] = n
                        output_dict = self.fill_empty_output(output_dict)
                        self.output_writer.write_data(output_dict)
            except Exception as e:
                self._log_runtime_exception(lnum, line, input_data, e)

    def cleanup(self):
        pass


if __name__ == "__main__":
    annotator = CravatAnnotator(sys.argv)
    annotator.run()
