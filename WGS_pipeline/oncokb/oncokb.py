import sys
from cravat import BaseAnnotator
import os


class CravatAnnotator(BaseAnnotator):
    def setup(self):
        """
        Set up data sources.
        Cravat will automatically make a connection to
        data/example_annotator.sqlite using the sqlite3 python module. The
        sqlite3.Connection object is stored as self.dbconn, and the
        sqlite3.Cursor object is stored as self.cursor.
        """
        try:
            import pandas as pd
        except ModuleNotFoundError:
            res = os.system("pip install pandas")
            if res != 0:
                print("pip install pandas failed")
                sys.exit(1)
            import pandas as pd
        dir_path = os.path.dirname(os.path.realpath(__file__))
        datafile_path = os.path.join(dir_path, "data", "oncokb.csv")
        self.oncokb = pd.read_csv(datafile_path)
        self.oncokb["Hotspot"] = self.oncokb["Hotspot"].astype(str)
        self.oncokbset = {
            val["Chromosome"]
            + ":"
            + str(val["Position"])
            + ":"
            + val["Alt"]: [
                val["ProteinChange"],
                val["Oncogenic"],
                val["MutationEffect"],
                val["Hotspot"],
            ]
            for _, val in self.oncokb.iterrows()
        }

    def annotate(self, input_data, secondary_data=None):
        """
        The annotator parent class will call annotate for each line of the
        input file. It takes one positional argument, input_data, and one
        keyword argument, secondary_data.

        input_data is a dictionary containing the data from the current input
        line. The keys depend on what what file is used as the input, which can
        be changed in the module_name.yml file.
        Variant level includes the following keys:
            ('uid', 'chrom', 'pos', 'ref_base', 'alt_base')
        Variant level crx files expand the key set to include:
            ('hugo', 'transcript','so','all_mappings')
        Gene level files include
            ('hugo', 'num_variants', 'so', 'all_so')

        secondary_data is used to allow an annotator to access the output of
        other annotators. It is described in more detail in the CRAVAT
        documentation.

        annotate should return a dictionary with keys matching the column names
        defined in example_annotator.yml. Extra column names will be ignored,
        and absent column names will be filled with None. Check your output
        carefully to ensure that your data is ending up where you intend.
        """
        out = {}
        out["protein_change"] = ""
        out["oncogenic"] = ""
        out["mutation_effect"] = ""
        out["hotspot"] = ""

        chrom = input_data["chrom"]
        pos = input_data["pos"]
        alt = input_data["alt_base"]

        key = str(chrom) + ":" + str(pos) + ":" + alt
        if key in self.oncokbset:
            print("found oncokb oncogenic:" + key)
            out["protein_change"] = self.oncokbset[key][0]
            out["oncogenic"] = self.oncokbset[key][1]
            out["mutation_effect"] = self.oncokbset[key][2]
            out["hotspot"] = self.oncokbset[key][3]
        return out

    def cleanup(self):
        """
        cleanup is called after every input line has been processed. Use it to
        close database connections and file handlers. Automatically opened
        database connections are also automatically closed.
        """
        pass


if __name__ == "__main__":
    annotator = CravatAnnotator(sys.argv)
    annotator.run()
