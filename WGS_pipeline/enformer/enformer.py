import sys
from cravat import BaseAnnotator
from cravat import InvalidData
import sqlite3
import os

SEQ_LENGTH = 393_216

ENFORMER_LINEAGE_CELL_MAPPING = {"": ""}


class CravatAnnotator(BaseAnnotator):
    def setup(self):
        """
        Set up data sources.
        Cravat will automatically make a connection to
        data/example_annotator.sqlite using the sqlite3 python module. The
        sqlite3.Connection object is stored as self.dbconn, and the
        sqlite3.Cursor object is stored as self.cursor.
        """
        self.run_it = True
        try:
            import tensorflow as tf
            import tensorflow_hub as hub
            from kipoi_interpret.importance_scores.ism import Mutation
            from kipoiseq import Interval
            from . import enformer_utils as eu
            import pandas as pd
        except ModuleNotFoundError:
            # install it
            res = os.system(
                "pip install tensorflow tensorflow-hub kipoi_interpret kipoiseq pandas"
            )
            if res != 0:
                print(
                    "Failed to install some modules: tensorflow tensorflow-hub kipoi_interpret kipoiseq pandas"
                )
                self.run_it = False
                sys.exit(1)
            import tensorflow as tf
            import tensorflow_hub as hub
            from kipoi_interpret.importance_scores.ism import Mutation
            from kipoiseq import Interval
            import pandas as pd

            self.Interval = Interval

        # use kipoi to make variant level stuff / load the fasta sequences

        self.fasta = pyfaidx.Fasta(fasta_file)
        self._chromosome_sizes = {k: len(v) for k, v in self.fasta.items()}

        kipoiseq.transforms.functional.one_hot_dna(sequence).astype(np.float32)
        # load the kipoi model
        enformer = hub.Module("https://tfhub.dev/deepmind/enformer/1")
        # model = kipoi.get_model("DeepBind/Homo_sapiens/TF/D00765.001_ChIP-seq_GATA1")

        # if kipoi, use the variant effect predictor tool, else make your own scoring function

        # parametrize the model (define what to predict, where to look at, ...)
        # Numpy array [batch_size, SEQ_LENGTH, 4] one hot encoded in order 'ACGT'. The
        # `one_hot_encode` function is available in `enformer.py` and outputs can be
        # stacked to form a batch.
        inputs = tf.zeros((1, SEQ_LENGTH, 4), dtype=tf.float32)
        predictions = enformer.predict_on_batch(inputs)
        predictions["human"].shape  # [batch_size, 896, 5313]

    def extract(self, interval, **kwargs) -> str:
        # Truncate interval if it extends beyond the chromosome lengths.
        chromosome_length = self._chromosome_sizes[interval.chrom]
        trimmed_interval = self.Interval(
            interval.chrom,
            max(interval.start, 0),
            min(interval.end, chromosome_length),
        )
        # pyfaidx wants a 1-based interval
        sequence = str(
            self.fasta.get_seq(
                trimmed_interval.chrom,
                trimmed_interval.start + 1,
                trimmed_interval.stop,
            ).seq
        ).upper()
        # Fill truncated values with N's.
        pad_upstream = "N" * max(-interval.start, 0)
        pad_downstream = "N" * max(interval.end - chromosome_length, 0)
        return pad_upstream + sequence + pad_downstream

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
        out = {val: "" for val in where_to_look}
        target_interval = kipoiseq.Interval("chr11", 35_082_742, 35_197_430)  # @param

        kipoiseq.dataclasses.Variant(chrom=chrom, pos=pos, ref=ref, alt=alt, id=id)
        # predict variant level stuff by adding one variant only to the ref

        # save variant to list of variant in current seq

        # predict variant level stuff by adding all variant in the sequence
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
