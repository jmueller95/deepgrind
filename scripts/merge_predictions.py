import multiprocessing as mp
import pandas as pd
import numpy as np
from functools import partial
from pyteomics import mgf
import utils


def main():
    rt_prediction_df = pd.read_csv(snakemake.input['rt_prediction_output'], sep="\t", index_col=False)
    with mgf.read(snakemake.input['msms_prediction_output'], read_ions=True, convert_arrays=1,
                  index_by_scans=False) as msms_prediction_mgf:

        # Split indices into as many chunks as we have cores
        index_splits = np.array_split(np.asarray(list(msms_prediction_mgf.index.keys())), mp.cpu_count())
        # Merge RT with MSMS Chunk-Wise in parallel, then concatenate the results
        with mp.Pool(mp.cpu_count()) as pool:
            prediction_mgf_aslist = np.hstack(
                pool.map(partial(utils.add_rt_to_spectra, msms_prediction_mgf, rt_prediction_df),
                         index_splits))
        mgf.write(spectra=prediction_mgf_aslist, output=snakemake.output[0],
                  fragment_format="%.6f %.7f %s", header=msms_prediction_mgf.header, write_charges=False, write_ions=True)


if __name__ == '__main__':
    main()
