
TEMPLATE_BATCH_FILE="/home/talmalu/thesis/projects/python/NetworkAnalysis/Scripts/template_batches/neighborhood_local_clustering_coefficient_batch_template.sbatch"
BATCH_RUNNER="/home/talmalu/thesis/projects/python/batch_runner.sh"

POSITIONAL_ARGS=()
PROTEIN=""
DATABASE=""
XGMML_FILE=""
OUTPUT=""
ALIGNMENT_WEIGHT=""
DEPENDENCY=""
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
	-db|--database)
		DATABASE=$2
		shift
		shift
	;;
	-o|--output)
		OUTPUT=$2
		shift
		shift
	;;
	-w|--alignment_weight)
		ALIGNMENT_WEIGHT=$2
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
#echo "DATABASE              = $DATABASE"
#echo "XGMML_FILE            = $XGMML_FILE"
#echo "IMAGE                 = $IMAGE"
#echo "IMAGE_TITLE           = $IMAGE_TITLE"
#echo "OUTPUT                = $OUTPUT"
#echo "POSITIONAL_ARGS       = ${POSITIONAL_ARGS[*]}"



SBATCH_COMMAND=("$BATCH_RUNNER")
TMP_FILE_PATH=$(mktemp)
cp $TEMPLATE_BATCH_FILE $TMP_FILE_PATH

sed -i "s/###PROTEIN###/$PROTEIN/g" $TMP_FILE_PATH
sed -i "s~###XGMML_FILE###~$XGMML_FILE~g" $TMP_FILE_PATH
sed -i "s~###OUTPUT###~$OUTPUT~g" $TMP_FILE_PATH


if [[ ! -z "$DATABASE" ]]
then
	sed -i "s~###DATABASE###~$DATABASE~g" $TMP_FILE_PATH
else
	sed -i "/###DATABASE###/d" $TMP_FILE_PATH
fi

if [[ ! -z "$ALIGNMENT_WEIGHT" ]]
then
	sed -i "s~###ALIGNMENT_WEIGHT###~$ALIGNMENT_WEIGHT~g" $TMP_FILE_PATH
else
	sed -i "/###ALIGNMENT_WEIGHT###/d" $TMP_FILE_PATH
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
