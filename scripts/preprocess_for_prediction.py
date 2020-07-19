import pandas as pd
import utils

def check_msms_model_name(converter):
    def wrapper(*args, **kwargs):
        if kwargs['style'] not in ["pdeep", "prosit"]:
            raise Exception("MSMS model must be 'pdeep' or 'prosit'")
        converter(*args, **kwargs)
    return wrapper


@check_msms_model_name
def _convert_for_msms(comet_df, style, output):
    if style == "prosit":
        res = pd.DataFrame(
            {"modified_sequence": comet_df.Peptide.apply(
                lambda pep: utils.find_modifications(pep[2:-2], style="prosit")).values,
             "collision_energy": snakemake.params['collision_energy'],
             "precursor_charge": comet_df.apply(lambda row: row.Charge1 + 2 * row.Charge2 + 3 * row.Charge3, axis=1)})
        res.dropna(inplace=True)
        res.to_csv(output, sep=",", header=True, index=False)
    else:
        res = pd.DataFrame(
            comet_df.Peptide.apply(lambda pep: utils.find_modifications(pep[2:-2], style="pdeep")).to_list(),
            columns=["peptide", "modification"])
        # The charge is one-hot encoded in the comet df, so we can resolve this into 1,2 or 3 by multiplying 1,2 and 3
        # with the entries of Charge1, Charge2 and Charge3
        res["charge"] = comet_df.apply(lambda row: row.Charge1 + 2 * row.Charge2 + 3 * row.Charge3, axis=1)
        res.dropna(inplace=True)
        res.to_csv(output, sep="\t", header=True, index=False)


@check_msms_model_name
def _convert_for_rt(comet_df, style, output):
    if style == "prosit":
        res = pd.DataFrame(
            {"modified_sequence": comet_df.Peptide.apply(lambda pep: utils.find_modifications(pep[2:-2], style="prosit")).values})
        res.dropna(inplace=True)
        res.to_csv(output, sep=",", header=True, index=False)
    else:
        raise Exception("Not implemented. Right now, the only accepted RT Model is 'prosit'.")


def main():
    # Parse the input file:
    comet_df = pd.read_csv(snakemake.input[0], sep="\t", header=0,
                           usecols=["Peptide", "Charge1", "Charge2", "Charge3"],
                           index_col=False)
    # Determine if MSMS and RT prediction will be performed jointly or separately
    if "msms_model" in dict(snakemake.params) and "rt_model" in dict(snakemake.params):
        _convert_for_msms(comet_df, style=snakemake.params['msms_model'].lower(),
                          output=snakemake.output['msms_prediction_input'])
        _convert_for_rt(comet_df, style=snakemake.params['rt_model'].lower(),
                        output=snakemake.output['rt_prediction_input'])
    else:
        # If only one model was supplied, the prediction will be joint
        # Only convert the input for msms in that case
        _convert_for_msms(comet_df, style=snakemake.params['model'].lower(),
                          output=snakemake.output['prediction_input'])


if __name__ == '__main__':
    main()
