#!/bin/bash
VERSION=0.1
RELEASE_DATE="2020-07-06"

"
PARAM INFOS FROM OLD PYTHON FILE:
	'INPUT', type=str,
	                    help='INPUT must be either a comma-separated list of MGF files or a dictionary containing MGF files.\n'
	                         'In the latter case, all files ending on '.mgf' will be used as input.')
	'-m', '--msms_model', required=True, type=str.lower, choices=['prosit', 'pdeep'],
	                    help='Name of fragmentation prediction model to be used for the rescoring.\n'
	                         'Must be one of 'Prosit' or 'pDeep' (case insensitive).')
	'-r', '--rt_model', required=True, type=str.lower, default='prosit', choices=['prosit'],
	                    help='Name of retention time prediction model to be used for the rescoring.\n'
	                         'Right now, defaults to 'Prosit', but it is planned to also allow 'pDeep'\n'
	                         'once pDeep3 gets released, which will include a retention time model.')
	'-x', '--msms_model_address', required=True, type=str,
	                    help='URL of the server to request the MSMS prediction from.')
	'-y', '--rt_model_address', required=False, type=str,
	                    help='URL of the server to request the retention time prediction from.\n'
	                         'Required only if the model under the --msms-model-address does not already return retention times.')
	'-s', '--scoring_function', required=True, type=str.lower,
	                    choices=['dp', 'dotproduct', 'hs', 'hyperscore'],
	                    help='Scoring Function to be used for rescoring the experimental with the\n'
	                         'predicted fragmentation spectra. Must be one of:\n'
	                         ''dp' / 'dotproduct' --> Use Dot Product (Recommended for IT Spectra)\n'
	                         ''hs' / 'hyperscore' --> Use Hyperscore (Recommended for FT Spectra))\n')
	'-b', '--binsize', required=True, type=float,
	                    help='Binsize in m/z to be used in the rescoring. Recommended values are:\n'
	                         '0.005 [FT Spectra with Precursor Charges 2,3]\n'
	                         '0.01 [FT Spectra with Precursor Charge 1]\n'
	                         '0.1 [IT Spectra]')
	'-f', '--fragmentation_mode', required=True, type=str.lower, default='cid',
	                    choices=['cid', 'hcd'],
	                    help='Fragmentation mode used to acquire the MS2 spectra. Must be one of:\n'
	                         'CID [Default] 	[Collision Energy 35 will be used]\n'
	                         'HCD		[Collision Energy 25 will be used]')
	'-a', '--analyzer', required=True, type=str.lower, choices=['ft', 'it'],
	                    help='Type of the Mass Analyzer used to acquire the MS2 spectra. Must be one of:\n'
	                         ''FT' --> Fourier Transform Mass Mass Spectrometer\n'
	                         ''IT' --> Ion Trap Mass Spectrometer')
	'-o', '--output_directory', required=False, type=str, default='.',
	                    help='Path to the directory where the final percolator output will be sured.\n'
	                         'Please make sure you have read/write access to this directory.\n'
	                         'Default: Current Directory ('.')')
	'-n', '--neutral_losses', required=False, type=int, default=1,
	                    choices=[1, 0],
	                    help='Include or Exclude Neutral Loss Predictions. Must be one of:\n'
	                         ''1' -> Include\n'
	                         ''0' -> Exclude\n'
	                         'Default: 1')
	'-e', '--search_engine', required=False, type=str,
	                    help='Path to an alternative search engine to be used instead of Comet.\n'
	                         'Must return a tab-separated file of PSMs containing the following\n'
	                         'mandatory columns: <TODO>')
	'-p', '--psms', required=False, type=str,
	                    help='Path to a file (or a comma-separated list of files) of Search Engine\n'
	                         'Results. If provided, these well be used for rescoring and the Comet\n'
	                         'search will be skipped. The files must fulfil the same requirements as\n'
	                         'the files returned by an alternative search engine (see above).')
	'-l', '--library', required=False, type=str,
	                    help='Path to an SQLite file containing libraries of predicted msms spectra\n'
	                         'and retention times. If provided, only the predictions not present in\n'
	                         'the library will be predicted and added to the library afterwards.\n'
	                         'This can be used to speed up the pipeline.'

"


show_help()
{
	printf "Usage: deepgrind [PARAMETERS] <INPUT>

<INPUT> must be either a comma-separated list of MGF files or a dictionary containing MGF files. 
In the latter case, all files ending on '.mgf' will be used as input.

Required Parameters:

	-m <string>
	--msms-model <string>					Name of fragmentation prediction model to be used for the rescoring.
								Must be one of 'Prosit' or 'pDeep' (case insensitive).

	-r <string>
	--rt-model <string>					Name of retention time prediction model to be used for the rescoring.
								Right now, defaults to 'Prosit', but it is planned to also allow 'pDeep'
								once pDeep3 gets released, which will include a retention time model.

	-x <string>
	--msms-model-address <string>				URL of the server to request the MSMS prediction from.

	-y <string>
	--rt-model-address <string>				URL of the server to request the retention time prediction from.
								Required only if the model under the --msms-model-address does not return retention times 

	-n <boolean>
	--neutral-losses <boolean>				Include or Exclude Neutral Loss Predictions. Must be one of:
								'True' / 'Yes' / '1'
								'False' / 'No' / '0'
								Default: True

	-s <string>
	--scoring-function <string>				Scoring Function to be used for rescoring the experimental with the
								predicted fragmentation spectra. Must be one of:
								'dp' / 'dotproduct' --> Use Dot Product (Recommended for IT Spectra)
								'hs' / 'hyperscore' --> Use Hyperscore (Recommended for FT Spectra)

	-b <float>
	--binsize <float>					Binsize in m/z to be used in the rescoring. Recommended values are:
								0.005 [FT Spectra with Precursor Charges 2,3]
								0.01 [FT Spectra with Precursor Charge 1]
								0.1 [IT Spectra]

	-f <string>
	--fragmentation-mode <string>				Fragmentation mode used to acquire the MS2 spectra. Must be one of:
								CID [Default] 	[Collision Energy 35 will be used]
								HCD		[Collision Energy 25 will be used]

	-a <string>
	--analyzer <string>					Type of the Mass Analyzer used to acquire the MS2 spectra. Must be one of:
								'FT' --> Fourier Transform Mass Mass Spectrometer
								'IT' --> Ion Trap Mass Spectrometer		


Optional Parameters:

	-o <path>
	--output-directory <path> 				Path to the directory where the final percolator output will be sured.
								Please make sure you have read/write access to this directory.
								Default: Current Directory ('.')

	-e <path>
	--search-engine <path>					Path to an alternative search engine to be used instead of Comet. 
								Must return a tab-separated file of PSMs containing the following 
								mandatory columns:

	-p <path>
	--psms <path>						Path to a file (or a comma-separated list of files) of Search Engine 
								Results. If provided, these well be used for rescoring and the Comet 
								search will be skipped. The files must fulfil the same requirements as 
								the files returned by an alternative search engine (see above).

	-l <path>
	--library <path>					Path to an SQLite file containing libraries of predicted msms spectra 
								and retention times. If provided, only the predictions not present in
								the library will be predicted and added to the library afterwards.
								This can be used to speed up the pipeline. 

	-h 
	--help 							Show this help message and exit.
" 
}

parse_input()
{	
	#If input is directory, grep all '.mgf' files and return as array
	if [ -d $1 ]; then
		FILELIST=($1/*.mgf)
	else
		#If input is a string, split it by ',' and return as array
		IFS=',' read -r -a FILELIST <<< "$1"
	fi
	# #Check if all input elements are actually mgf files
	for FILE in ${FILELIST[@]}; do
		if [ ! -f ${FILE} ] || [[ ! ${FILE} == *.mgf ]]; then
			printf "Error in parsing input: '${FILE}' is not an MGF file.\nPlease provide either a comma-separated list of or a dictionary containing MGF files."
			exit 1
		fi	
	done

}

##### MAIN ROUTINE #####
printf "IMADeepGrind version ${VERSION}, released ${RELEASE_DATE}.\nWritten by JuMu.\n\n"

### 1. Parameter Initialization ###

#Set default values for some of the parameters
FRAGMENTATION_MODE="cid"
SEARCH_ENGINE="comet"
RT_MODEL="prosit"
USE_NEUTRAL_LOSSES=1
OUTPUT_DIRECTORY="."
#Parse Parameters
while [ "$1" != "" ]; do
	case $1 in
		-m | --msms-model )			shift
									MSMS_MODEL=$(printf $1 | tr '[:upper:]' '[:lower:]')
									;;
		-r | --rt-model )			shift
									RT_MODEL=$(printf $1 | tr '[:upper:]' '[:lower:]')
									;;
		-x | --msms-model-address )	shift
									MSMS_MODEL_ADDRESS=$(printf $1 | tr '[:upper:]' '[:lower:]')
									;;
		-y | --rt-model-address )	shift
									RT_MODEL_ADDRESS=$(printf $1 | tr '[:upper:]' '[:lower:]')
									;;
		-n | --neutral-losses )		shift
									USE_NEUTRAL_LOSSES=$(printf $1 | tr '[:upper:]' '[:lower:]')
									;;
		-s | --scoring-function )	shift
									SCORING_FUNCTION=$(printf $1 | tr '[:upper:]' '[:lower:]')
									;;
		-b | --binsize )			shift
									BINSIZE=$(printf $1 | tr '[:upper:]' '[:lower:]')
									;;
		-f | --fragmentation-mode )	shift
									FRAGMENTATION_MODE=$(printf $1 | tr '[:upper:]' '[:lower:]')
									;;
		-a | --analyzer )			shift
									ANALYZER=$(printf $1 | tr '[:upper:]' '[:lower:]')
									;;
		-o | --output-directory )	shift
									OUTPUT_DIRECTORY=$(printf $1 | tr '[:upper:]' '[:lower:]')
									;;
		-e | --search-engine )		shift
									SEARCH_ENGINE=$(printf $1 | tr '[:upper:]' '[:lower:]')
									;;
		-p | --psms )				shift
									PSMS=$(printf $1 | tr '[:upper:]' '[:lower:]')
									;;
		-l | --library ) 			shift
									LIBRARY=$(printf $1 | tr '[:upper:]' '[:lower:]')
									;;
		-h | --help )				show_help
									exit 1
									;;
		* )							parse_input $1
									
	esac
	shift
done

printf "MGF Input Files:\n"
for FILE in ${FILELIST[@]}; do
	printf "${FILE}\n"
done

printf "\nParameters:\n"
#If any mandatory parameter is missing/invalid, display a helpful string and exit
if [ "${MSMS_MODEL}" !=  'pdeep' ] && [ "${MSMS_MODEL}" !=  'prosit' ]; then
	if [  "${MSMS_MODEL}" == '' ]; then
		printf "Error: No MSMS model specified. Please specify a MSMS model using the '-m' option.\nThe accepted values are: 'pDeep' and 'Prosit'.\n"
	else
		printf "Error: '${MSMS_MODEL}' is not a valid MSMS model. The accepted values are: 'pDeep' and 'Prosit'.\n"
	fi
	exit 1
else
	printf "MSMS Model:\t\t${MSMS_MODEL}\n"
fi

if [ "${RT_MODEL}" !=  'prosit' ]; then
	printf "Error: '${RT_MODEL}' is not a valid retention time model. Right now, only 'Prosit' [Default] is accepted.\n"
	exit 1
else
	printf "RT Model:\t\t${RT_MODEL}\n"
fi

if [  ! "${MSMS_MODEL_ADDRESS}" ]; then
	printf "Error: MSMS model address not specified. Please specify the URL of the MSMS prediction model using the '-x' option.\nIf this model does NOT return RT predictions, you also need to specify the URL of an RT model using the '-y' option.\n"
	exit 1
fi

if [ ! "${RT_MODEL_ADDRESS}" ] && [ "${MSMS_MODEL}" != "${RT_MODEL}" ]; then
	printf "Error: Two different models were specified for MSMS and RT prediction, but only one address was given.\nPlease provide a separate address for the RT prediction using the '-y' option.\n"
	exit 1
fi

if [ "${USE_NEUTRAL_LOSSES}" ==  'true' ] || [ "${USE_NEUTRAL_LOSSES}" ==  'yes' ] || [ "${USE_NEUTRAL_LOSSES}" ==  1 ]; then
	USE_NEUTRAL_LOSSES=1
	printf "Predict Neutral Losses:\tYes\n"
elif [ "${USE_NEUTRAL_LOSSES}" ==  'false' ] || [ "${USE_NEUTRAL_LOSSES}" ==  'no' ] || [ "${USE_NEUTRAL_LOSSES}" ==  0 ]; then
	USE_NEUTRAL_LOSSES=0
	printf "Predict Neutral Losses:\tNo\n"
else
	printf "Error: '${USE_NEUTRAL_LOSSES}' is not a valid value for the '-n/--neutral-losses' option.\nThis is a boolean option. The accepted values are: 'true'/'yes'/'1' or 'false'/'no'/'0'.\n"
	exit 1
fi

if [ "${SCORING_FUNCTION}" == "dotproduct" ] || [ "${SCORING_FUNCTION}" == "dp" ]; then
	SCORING_FUNCTION='dp'
	printf "Scoring Function:\tDot Product\n"
elif [ "${SCORING_FUNCTION}" == "hyperscore" ] || [ "${SCORING_FUNCTION}" == "hs" ]; then
	SCORING_FUNCTION='hs'
	printf "Scoring Function:\tHyperscore\n"
else
	printf "Error: '${SCORING_FUNCTION}' is not a valid scoring function (option '-s/--scoring-function').\nThe accepted values are: 'dp'/'dotproduct' or 'hs'/'hyperscore'.\n"
	exit 1
fi

if  [[ ! "${BINSIZE}" =~ ^[+-]?[0-9]+\.?[0-9]*$ ]]; then
	printf "Error: Value of '-b/--binsize' option must be a floating point number.\n"
else
	printf "Bin Size:\t\t${BINSIZE}\n"
fi

if [ "${FRAGMENTATION_MODE}" !=  'cid' ] && [ "${FRAGMENTATION_MODE}" !=  'hcd' ]; then
	printf "Error: '${FRAGMENTATION_MODE}' is not a valid fragmentation mode (option -f/--fragmentation-mode).\nThe accepted values are: 'CID' [Default] or 'HCD'.\n"
	exit 1
else	
	printf "Fragmentation Mode:\t${FRAGMENTATION_MODE}\n"
fi


if [ "${ANALYZER}" == 'ft' ]; then
	COMET_PARAMS="comet_parameterfiles/FT.params"
	printf "Mass Analyzer:\t\t${ANALYZER}\n"
elif  [ "${ANALYZER}" != 'it' ]; then
	COMET_PARAMS="comet_parameterfiles/IT.params"
	printf "Mass Analyzer:\t\t${ANALYZER}\n"
else
	printf "Error: '${ANALYZER}' is not a valid mass analyzer type (option -a/--analyzer).\nThe accepted values are: 'FT' or 'IT'.\n"
fi

if [ "${SEARCH_ENGINE}" != 'comet' ]; then
	printf "Warning: 'Alternative search engine' option has not been implemented yet.\nFalling back to Comet.\n"
fi

if [ "${PSMS}" != '' ]; then
	printf "Warning: Processing of user-provided search engine results has not been implemented yet.\n"
fi

if [ "${LIBRARY}" != '' ]; then
	printf "Warning: Processing of user-provided prediction libraries has not been implemented yet.\n"
fi