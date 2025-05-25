
TEMPLATE_BATCH_FILE="/home/talmalu/thesis/projects/python/NetworkAnalysis/Scripts/template_batches/filter_only_of_taxonomy_network_batch_template.sbatch"
BATCH_RUNNER="/home/talmalu/thesis/projects/python/batch_runner.sh"

POSITIONAL_ARGS=()
TMP_TAXONOMY=()
TAXONOMY_LEVEL=""
XGMML_FILE=""
PROTEIN=""
OUTPUT=""
DEPENDENCY=""
MEMORY=""

while [[ $# -gt 0 ]]; do
	case "$1" in
	-p|--protein)
		export PROTEIN=$2
		shift
		shift
	;;
	-f|--xgmml_file)
		export XGMML_FILE=$2
		shift
		shift
	;;
	-t|--taxonomy)
		shift
		var="$1"
		while [[ ! -z "$var" && "$var" != -* ]] || [[ ! -z "$var" && "$var" != --* ]];
		do
			TMP_TAXONOMY+=("$var")
			shift
			var="$1"
		done
		export TAXONOMY=${TMP_TAXONOMY[@]}
	;;
	-db|--database)
		export DATABASE=$2
		shift
		shift
	;;
	-l|--taxonomy_level)
		export TAXONOMY_LEVEL=$2
		shift
		shift
	;;
	-o|--output)
		export OUTPUT=$2
		shift
		shift
	;;
	-m|--memory)
		MEMORY=$2
		shift
		shift
	;;
	-d|--dependency)
		DEPENDENCY=$2
		shift
		shift
	;;
	-*|--*)
		echo "Unknown option $1"
		exit 1
	;;
	*)
		var="$1"
		for v in ${var}
		do
			POSITIONAL_ARGS+=("$v") # save positional arg
		done
		shift # past argument
	;;
	esac
done

set -- "${POSITIONAL_ARGS[@]}" # restore positional parameters


#echo "PROTEIN               = $PROTEIN"
#echo "TAXONOMY_LEVEL        = $TAXONOMY_LEVEL"
#echo "TAXONOMY              = ${TAXONOMY[*]}"
#echo "XGMML_FILE            = $XGMML_FILE"
#echo "OUTPUT                = $OUTPUT"


SBATCH_COMMAND=("$BATCH_RUNNER")
TMP_FILE_PATH=$(mktemp)
cp $TEMPLATE_BATCH_FILE $TMP_FILE_PATH

#sed -i "s/###PROTEIN###/$PROTEIN/g" $TMP_FILE_PATH
#sed -i "s/###TAXONOMY_LEVEL###/$TAXONOMY_LEVEL/g" $TMP_FILE_PATH
#sed -i "s/###TAXONOMY###/${TAXONOMY[*]}/g" $TMP_FILE_PATH
#sed -i "s~###XGMML_FILE###~$XGMML_FILE~g" $TMP_FILE_PATH
#sed -i "s~###OUTPUT###~$OUTPUT~g" $TMP_FILE_PATH

if [[ ! -z "$DEPENDENCY" ]]
then
	SBATCH_COMMAND+=("--dependency $DEPENDENCY")
fi

if [[ -z "$MEMORY" ]]
then
	MEMORY=$(du -bcs $XGMML_FILE | grep total | awk '{printf "%dG",int($1/2^30 + 1) * 4 }')
fi

SBATCH_COMMAND+=("--memory $MEMORY")

SBATCH_COMMAND+=("$TMP_FILE_PATH")

JOB_ID=$(${SBATCH_COMMAND[@]}) #run command
echo $JOB_ID

