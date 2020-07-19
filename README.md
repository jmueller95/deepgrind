# DeepGrind
DeepGrind is a _snakemake_ workflow that identifies peptides from LC/MS-MS data. It combines the well-established peptidomics 
identification tools __Comet__ (http://comet-ms.sourceforge.net/) and __Percolator__ (https://github.com/percolator/percolator) with the 
Deep Learning approaches offered by __pDeep__ (https://github.com/pFindStudio/pDeep) and __Prosit__ (https://github.com/kusterlab/prosit).
These tools can predict fragmentation spectra (and, in the case of Prosit, chromatographic retention times) of arbitrary peptides 
which can then be used to rescore the PSMs found by Comet. The new scores subsequently allow Percolator to separate 
true from false discoveries with more confidence, which results in a higher number of identifications and smaller False Discovery Rates.
  
## Installation 
### Requirements
In order to run DeepGrind, you need to have [Conda](https://docs.conda.io/projects/conda/en/latest/index.html) installed. 
All other packages required (including Comet, Percolator, and Snakemake) will be installed by Conda on-the-fly.
Furthermore, you need a __Prosit__ instance that is able to return predicted spectra in Mascot Generic Format.
The vanilla version of Prosit (as of July 2020) does not do that, but you can e.g. use my fork of Prosit (https://github.com/jmueller95/prosit).
If you also want to use __pDeep__ for MSMS prediction you'll need a server instance of pDeep, which you can obtain from https://github.com/jmueller95/pDeep.

### Installation Steps
1. Clone this repository.
2. Create the Conda environment for the workflow by running `conda env create -f envs/deepgrind_env.yml`. 
This will also install Comet and Percolator.  
3. Activate the environment: `conda activate deepgrind_env`.


## Usage
The input spectra need to be provided in Mascot Generic Format (http://www.matrixscience.com/help/data_file_help.html).
DeepGrind expects each spectrum to contain a local parameter `IRT` with the indexed (standardized) retention time.
If you do not possess indexed retention times for your spectra, the workflow will simply skip the RT-based rescoring step 
(it is not possible to use the `RTINSECONDS` parameter for rescoring).
The output is simply the tab-delimited PSM list that Percolator usually outputs (using the `-r` flag as described 
[here](https://github.com/percolator/percolator/wiki/Command-line-options#file-output-options)). 

The workflow parameters are given in the form of a configuration file. See `config_template.yaml` for an example and explanation
on each parameter. Create your own config file, put your custom parameter values in it and save it to `config.yaml`.

Then, start the workflow by `cd`ing into the directory containing the `Snakefile` and calling `snakemake --cores all`.
Instead of `all`, you can also put an integer specifying the number of processor cores you want to make available to the workflow.
