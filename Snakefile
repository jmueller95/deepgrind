configfile: "config.yaml"

#Define some global variables that will be used in multiple rules
import os

#The 'basename' is the filename stripped of its directory and file extension
input_mgf_basename = os.path.splitext(os.path.basename(config['mgf']))[0]
percolator_output = "{}/{}_percolator_peptides.tsv".format(
    config['output_dir'] if "output_dir" in dict(config) else "deepgrind_out",
    input_mgf_basename)

rule all:
    input:
         percolator_output

rule run_comet:
    input:
         config['mgf']
    output:
          temp("tmp/{}.pin".format(input_mgf_basename))
    params:
          config['comet-parameterfile']
    shell:
         """
         OUTPUT={output}
         comet {input} -N${{OUTPUT%.*}} -P'{params}'
         #The Comet output files don't have the same number of columns in each line. The pandas module will have a hard time parsing them during rescoring later.
         #To avoid this, count the maximal number of columns in the file and add the appropriate number of '\\t' to the first line - Pandas will then be able to parse the file properly.
         ADD_COLUMNS=$(tr -d -c '\\t\\n' < '{output}' | awk '{{ if(length>max) {{max=length}}}} END {{print (max-24)}}'; )
         ADD_TABS=$( printf "%${{ADD_COLUMNS}}s" | sed 's/ /\\t/g' )
         sed -E -i "1s/Proteins/Proteins$ADD_TABS/" "{output}" 
         """

rule preprocess_for_prediction_joint:
    input:
         "tmp/{}.pin".format(input_mgf_basename),
    output:
          prediction_input=temp(
              "tmp/prediction_input.csv" if config['msms-model'] == config['rt-model'] else "culdesac")
    params:
          model=config['msms-model'],
          collision_energy=config['collision-energy']
    script:
          "scripts/preprocess_for_prediction.py"

rule preprocess_for_prediction_separate:
    input:
         "tmp/{}.pin".format(input_mgf_basename),
    output:
          msms_prediction_input=temp("tmp/msms_prediction_input.csv") if config['msms-model'] != config[
              'rt-model'] else "culdesac1",
          rt_prediction_input=temp("tmp/rt_prediction_input.csv") if config['msms-model'] != config[
              'rt-model'] else "culdesac2"
    params:
          msms_model=config['msms-model'],
          rt_model=config['rt-model'],
          collision_energy=config['collision-energy']
    script:
          "scripts/preprocess_for_prediction.py"

rule predict_joint:
    input:
         "tmp/prediction_input.csv"
    output:
          #If the user does NOT want neutral losses, add "NL" suffix to the output so snakemake uses the 'remove_nl_predictions' rule
          temp("tmp/prediction_output.mgf") if config["predict-neutral-losses"].lower() in ["yes",
                                                                                            "true"] else temp(
              "tmp/prediction_output_NL.mgf")
    params:
          model_address=config['msms-model-address']
    shell:
         "curl -o {output} -F 'peptides=@{input}' '{params.model_address}'"

rule predict_separate:
    input:
         msms_prediction_input="tmp/msms_prediction_input.csv",
         rt_prediction_input="tmp/rt_prediction_input.csv"
    output:
          msms_prediction_output=temp("tmp/msms_prediction_output.mgf"),
          rt_prediction_output=temp("tmp/rt_prediction_output.csv")
    params:
          msms_model_address=config['msms-model-address'],
          rt_model_address=config['rt-model-address']
    shell:"""
    curl -o {output.msms_prediction_output} -F 'peptides=@{input.msms_prediction_input}' '{params.msms_model_address}'&
    curl -o {output.rt_prediction_output} -F 'peptides=@{input.rt_prediction_input}' '{params.rt_model_address}'
    wait
    """

rule merge_predictions:
    input:
         msms_prediction_output="tmp/msms_prediction_output.mgf",
         rt_prediction_output="tmp/rt_prediction_output.csv"
    output:
          #If the user does NOT want neutral losses, add "NL" suffix to the output so snakemake uses the 'remove_nl_predictions' rule
          temp("tmp/prediction_output.mgf") if config["predict-neutral-losses"].lower() in ["yes", "true"] else temp(
              "tmp/prediction_output_NL.mgf")
    script:
          "scripts/merge_predictions.py"

rule remove_nl:
    input:
         "tmp/prediction_output_NL.mgf"
    output:
          temp("tmp/prediction_output.mgf")
    shell:
         "sed '/H2O+/d; /NH3+/d '  <{input} > {output}"

rule rescore:
    input:
         experimental_spectra=config["mgf"],
         predicted_spectra="tmp/prediction_output.mgf",
         pin="tmp/{}.pin".format(input_mgf_basename)
    output:
          temp("tmp/{}_rescored.pin".format(input_mgf_basename))
    params:
          score_function=config['score-function'],
          binsize=config['binsize']
    script:
          "scripts/rescore.py"

rule run_percolator:
    input:
         "tmp/{}_rescored.pin".format(input_mgf_basename)
    output:
          percolator_output
    shell:
         "percolator {input} --post-processing-tdc --results-peptides {output} --protein-decoy-pattern '_rev_' >/dev/null; "
         