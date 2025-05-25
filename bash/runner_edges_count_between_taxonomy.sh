
TEMPLATE_BATCH_FILE="/home/talmalu/thesis/projects/python/NetworkAnalysis/Scripts/template_batches/edges_count_between_taxonomy_batch_template.sbatch"
BATCH_RUNNER="/home/talmalu/thesis/projects/python/batch_runner.sh"

POSITIONAL_ARGS=()
MIN_SIZE=""
RELEVANCE_TAXONOMY=""
TAXONOMY_LEVEL=""
XGMML_FILE=""
PROTEIN=""
OUTPUT=""
DEPENDENCY=""
DATABASE=""
MEMORY=""

while [[ $# -gt 0 ]]; do
	case "$1" in
	-p|--protein)
		PROTEIN=$2
		shift
		shift
	;;
	-f|--xgmml_file)
		XGMML_FILE=$2
		shift
		shift
	;;
	-m|--min_size_threshold)
		MIN_SIZE=$2
		shift
		shift
	;;
	-r|--relevance_taxonomy_xlsx)
		RELEVANCE_TAXONOMY=$2
		shift
		shift
        ;;
	-l|--taxonomy_level)
		TAXONOMY_LEVEL=$2
		shift
		shift
	;;
	-o|--output)
		OUTPUT=$2
		shift
		shift
	;;
	-db|--database)
		DATABASE=$2
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
#echo "MIN_SIZE              = $MIN_SIZE"
#echo "XGMML_FILE            = $XGMML_FILE"
#echo "OUTPUT                = $OUTPUT"
#echo "RELEVANCE_TAXONOMY    = $RELEVANCE_TAXONOMY"


SBATCH_COMMAND=("$BATCH_RUNNER")
TMP_FILE_PATH=$(mktemp)
cp $TEMPLATE_BATCH_FILE $TMP_FILE_PATH


sed -i "s/###PROTEIN###/$PROTEIN/g" $TMP_FILE_PATH
sed -i "s/###TAXONOMY_LEVEL###/$TAXONOMY_LEVEL/g" $TMP_FILE_PATH
sed -i "s/###MIN_SIZE###/$MIN_SIZE/g" $TMP_FILE_PATH
sed -i "s~###XGMML_FILE###~$XGMML_FILE~g" $TMP_FILE_PATH
sed -i "s~###OUTPUT###~$OUTPUT~g" $TMP_FILE_PATH
sed -i "s~###RELEVANCY_XLSX###~$RELEVANCE_TAXONOMY~g" $TMP_FILE_PATH

if [[ ! -z "$DATABASE" ]]
then
	sed -i "s~###DATABASE###~$DATABASE~g" $TMP_FILE_PATH
else
	sed -i "/###DATABASE###/d" $TMP_FILE_PATH
fi

if [[ ! -z "$DEPENDENCY" ]]
then
	SBATCH_COMMAND+=("--dependency $DEPENDENCY")
fi

if [[ -z "$MEMORY" ]]
then
	MEMORY=$(du -bcs $XGMML_FILE | grep total | awk '{printf "%d",int($1/2^30 + 1) * 6 }')
fi

SBATCH_COMMAND+=("--memory $MEMORY")

SBATCH_COMMAND+=("$TMP_FILE_PATH")

JOB_ID=$(${SBATCH_COMMAND[@]}) #run command
echo $JOB_ID
