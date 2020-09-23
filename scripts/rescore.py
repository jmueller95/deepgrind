import multiprocessing as mp
import pandas as pd
import numpy as np
from functools import partial
from sklearn import preprocessing
from pyteomics import mgf
import utils
from math import factorial


def main():
    comet_df = pd.read_csv(snakemake.input['pin'], sep="\t", index_col=False)
    experimental_mgf = mgf.read(snakemake.input['experimental_spectra'], read_ions=False, convert_arrays=1,
                                index_by_scans=True)
    experimental_irt_present = "irt" in experimental_mgf.get_by_index(0)['params']
    if not experimental_irt_present:
        print("MGF File {} does not contain iRTs. Omitting RTDiff Rescoring.".format(
            snakemake.input['experimental_spectra']))
    if snakemake.params['score_function'].lower() in ["dp", "dot product", "dotproduct"]:
        score_function = "dot product"
    elif snakemake.params['score_function'].lower() in ["hs", "hyper score", "hyperscore"]:
        score_function = "hyperscore"
    else:
        raise Exception("Scoring function must either be 'dot product' or 'hyperscore'")

    # Split scans into as many chunks as we have processors available
    scan_splits = np.array_split(np.asarray(list(experimental_mgf.index.keys())), mp.cpu_count())
    # Read in the file with the predictions
    with mgf.read(snakemake.input['predicted_spectra'], read_ions=True, convert_arrays=1,
                  index_by_scans=False) as prediction_mgf:
        predicted_irt_present = "irt" in prediction_mgf.get_by_index(0)['params']
        if not predicted_irt_present:
            print("Prediction does not contain RT Predictions. Omitting RTDiff Rescoring.")
        print("Rescoring PSMs of file {} using {}{}...".format(snakemake.input['experimental_spectra'],
                                                               score_function,
                                                               " and RTDifference" if experimental_irt_present and predicted_irt_present else ""))
        # Run Rescoring on each chunk separately
        with mp.Pool(mp.cpu_count()) as pool:
            comet_rescored = pd.concat(pool.map(partial(rescore_psms,
                                                        experimental_mgf,
                                                        comet_df,
                                                        prediction_mgf,
                                                        {"msms": score_function,
                                                         "rt": "rtdiff" if experimental_irt_present and predicted_irt_present
                                                         else None},
                                                        snakemake.params['binsize']
                                                        ),
                                                scan_splits))

    comet_rescored.to_csv(snakemake.output[0], sep="\t", header=True, index=False, float_format='%.8f')

def rescore_psms(experimental_spectra, psm_df, predicted_spectra, rescoring_functions, binsize, scan_list):
    # We store the rescored PSMs in a list of Data Frames, which we then concatenate in the end
    rescored_psms_list = []
    # Iterate over scan numbers of mgf file
    for scan_nr in scan_list:
        # Get the PSMs for this scan
        psms_of_scan = pd.DataFrame(psm_df.loc[psm_df.ScanNr == int(scan_nr)])
        # If there were no hits for this spectrum, skip (mgf->comet is injective, but not surjective)
        if psms_of_scan.empty:
            continue
        # Get the actual charge from the one-hot encoding
        psms_of_scan['Charge'] = psms_of_scan.apply(lambda row: row.Charge1 + 2 * row.Charge2 + 3 * row.Charge3, axis=1)
        # Get the experimental spectrum
        experimental_spectrum = experimental_spectra.get_spectrum(scan_nr)

        # Get predictions for these peptides (Each prediction is formatted like this: TITLE=<Peptide>|<Modifications>|<Charge> )
        def get_predicted_spectrum(psm):
            psm_title = "{}|{}".format("|".join(utils.find_modifications(psm.Peptide[2:-2], style="pdeep") or [""]),
                                       psm.Charge)
            try:
                return predicted_spectra.get_by_id(psm_title)
            except KeyError:
                print("Did not find prediction for peptide: {}".format(psm.Peptide))
                return np.nan

        psms_of_scan['Predicted Spectrum'] = psms_of_scan.apply(get_predicted_spectrum, axis=1)

        # Drop rows where the prediction is nan
        psms_of_scan = psms_of_scan[~psms_of_scan['Predicted Spectrum'].isna()]
        # If no predictions could be found for this spectrum, skip
        if psms_of_scan.empty:
            continue
        # Perform the actual rescoring
        if rescoring_functions['msms'] == "dot product":
            psms_of_scan['DotProduct_rescored'] = dot_product(experimental_spectrum,
                                                              psms_of_scan['Predicted Spectrum'].values, binsize)
        elif rescoring_functions['msms'] == "hyperscore":
            psms_of_scan['Hyperscore_rescored'] = hyperscore(experimental_spectrum,
                                                             psms_of_scan['Predicted Spectrum'].values, binsize)
        else:
            raise Exception("Scoring function must either be dot product or hyperscore")
        if rescoring_functions['rt'] == "rtdiff":
            psms_of_scan['RTDifference_rescored'] = rt_difference(experimental_spectrum,
                                                                  psms_of_scan['Predicted Spectrum'].values)

        rescored_psms_list.append(psms_of_scan)

    # Merge all the rescored psms into a single dataframe
    rescored_psms_df = pd.concat(rescored_psms_list)

    # The columns of the new scores must be placed before the last two columns of the original dataframe, otherwise Percolator will not see them
    # So remember the column names to change their order later
    old_psm_columns = list(psm_df.columns)
    new_psm_columns = list(col for col in rescored_psms_df.columns if "rescored" in col)
    # Change the column order so Percolator can use the rescoring column(s)
    rescored_psms_df = rescored_psms_df[
        old_psm_columns[:old_psm_columns.index("Proteins") - 1] + new_psm_columns + old_psm_columns[
                                                                                    old_psm_columns.index(
                                                                                        "Proteins") - 1:]]
    # Rename the "Unnamed..." columns to empty strings (they contain additional proteins)
    unnamed_colnames = rescored_psms_df.columns[rescored_psms_df.columns.get_loc("Proteins") + 1:]
    rescored_psms_df.rename(columns={colname: "" for colname in unnamed_colnames}, inplace=True)
    return rescored_psms_df


def calc_dot_product_and_ions(experimental_spectrum, predicted_spectrum, binsize):
    """
    Worker function for both 'dot product' and 'hyperscore'.
    Dot Product is basically calculated like numpy.dot(), but by iterating over the m/z arrays instead of binning.
    This way, the number of matched b and y ions can be reported.
    Also, the intensities of peaks which match multiple times are divided equally among their matches.
    I could not achieve this with binning, which the following example illustrates:
    Consider predicted peaks at 99.995 mz and 100.010 mz, an experimental peak at 100.000 mz and bin size of 0.1 mz.
    Binning would put the first predicted peak into the bin starting at 99.9
    and both the second predicted peak and the experimentel peak into the bin starting at 100.0.
    Thus, only the second predicted peak would be recorded as a match to the experimental peak,
    even though the absolute distance of the first predicted peak to the experimental peak is smaller.
    """

    # Make sure predicted intensities are sorted by ascending m/zs
    predicted_ion_order = predicted_spectrum['m/z array'].argsort()
    predicted_mzs_sorted = predicted_spectrum['m/z array'][predicted_ion_order]
    predicted_intensities_sorted = predicted_spectrum['intensity array'][predicted_ion_order]
    predicted_ions_sorted = predicted_spectrum['ion array'][predicted_ion_order]

    # Ions are stored in sets to avoid counting an ion multiple times
    matched_b = set()
    matched_y = set()
    # Each experimental peak is assigned a 'bucket' of intensities of predicted ions that have matched it
    exp_buckets = [[] for _ in range(len(experimental_spectrum['m/z array']))]
    pred_index = 0
    exp_index = 0

    while pred_index < len(predicted_mzs_sorted) and exp_index < len(experimental_spectrum['m/z array']):
        # Test if the two peaks lie close enough to one another to be counted
        if abs(predicted_mzs_sorted[pred_index] - experimental_spectrum['m/z array'][exp_index]) < binsize and predicted_intensities_sorted[pred_index]>0:
            # Add ion to set of matched b or y ions
            if predicted_ions_sorted[pred_index].startswith('y'):
                # The splitting makes sure that neither fragments with different charge states nor neutral loss peaks are counted extra
                matched_y.add(predicted_ions_sorted[pred_index].split("+")[0].split("-")[0])
            else:
                matched_b.add(predicted_ions_sorted[pred_index].split("+")[0].split("-")[0])
            # Test if there are other experimental peaks which match this predicted one
            matched_exps = []
            k = exp_index
            while pred_index < len(predicted_mzs_sorted) and k < len(experimental_spectrum['m/z array']) and \
                    abs(predicted_mzs_sorted[pred_index] - experimental_spectrum['m/z array'][k]) < binsize:
                matched_exps.append(k)
                k += 1
            # Divide the intensity of the predicted peak among all experimental peaks that it has matched and add it to their bins
            for index in matched_exps:
                exp_buckets[index].append(predicted_intensities_sorted[pred_index] / len(matched_exps))
            # Increment predicted index, but not experimental (because the next predicted peak might still match it)
            pred_index += 1
        # If the two peaks are not close to one another, increment the index of the peak with smaller m/z and try again
        elif predicted_mzs_sorted[pred_index] > experimental_spectrum['m/z array'][exp_index]:
            exp_index += 1
        else:
            pred_index += 1
    # Now each experimental peak has a bucket of intensities of all predicted peaks that have matched it
    # To compute the dot product:
    # sum up those buckets and...
    # ...multiply these sums with the intensity of the respective experimental peak...
    # ...divided by the length of each bucket (i.e. distribute the experimental intensities equally among all matched predicted peaks)
    dot_product = sum([sum(bucket) * intensity / len(bucket)
                       for bucket, intensity in zip(exp_buckets, experimental_spectrum['intensity array']) if bucket])

    # Normalize the Dot Product by the vector norm of the two intensity vectors
    dot_product /= (
            np.linalg.norm(experimental_spectrum['intensity array']) * np.linalg.norm(predicted_intensities_sorted))
    return dot_product, len(matched_b), len(matched_y)


def dot_product(experimental_spectrum, predicted_spectra, binsize):
    dp_scores = [calc_dot_product_and_ions(experimental_spectrum, pred_spec, binsize)[0] for pred_spec in
                 predicted_spectra]
    return preprocessing.normalize(np.reshape(dp_scores, (1, -1)))[0]


def hyperscore(experimental_spectrum, predicted_spectra, binsize):
    dps_with_ions = [calc_dot_product_and_ions(experimental_spectrum, pred_spec, binsize) for pred_spec in
                     predicted_spectra]
    hyperscores = [dp * factorial(n_b) * factorial(n_y)
                   for dp, n_b, n_y in dps_with_ions]
    return preprocessing.normalize(np.reshape(hyperscores, (1, -1)))[0]


def rt_difference(experimental_spectrum, pred_spectra):
    """Return the absolute differences between the experimental RT and the predicted RTs
    """
    return [abs(float(experimental_spectrum['params']['irt']) - float(pred_spectrum['params']['irt'])) for pred_spectrum
            in pred_spectra]


if __name__ == '__main__':
    main()
