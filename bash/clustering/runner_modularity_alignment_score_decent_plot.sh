
TEMPLATE_BATCH_FILE="/home/talmalu/thesis/projects/python/NetworkAnalysis/Scripts/template_batches/clustering/modularity_alignment_score_decent_plot_batch_template.sbatch"
BATCH_RUNNER="/home/talmalu/thesis/projects/python/batch_runner.sh"

POSITIONAL_ARGS=()
PROTEIN=""
XGMML_FILE=""
DATABASE=""
OUT_FILE=""
TAXONOMY_DATA_FILE=""
GREEDY_DATA_FILE=""
LEIDEN_DATA_FILE=""
IMAGE_TITLE=""
REVERSE=""
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
	-t|--taxonomy_data_file)
		TAXONOMY_DATA_FILE=$2
		shift
		shift
	;;
	-g|--greedy_data_file)
		GREEDY_DATA_FILE=$2
		shift
		shift
	;;
	-l|--leiden_data_file)
		LEIDEN_DATA_FILE=$2
		shift
		shift
	;;
	-db|--database)
		DATABASE=$2
		shift
		shift
	;;
	-o|--out_file)
		OUT_FILE=$2
		shift
		shift
	;;
	-it|--image_title)
		IMAGE_TITLE=$2
		shift
		shift
	;;
	-r|--reverse)
		REVERSE=$2
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
#echo "DATA_FILE             = $DATA_FILE"
#echo "OUT_FILE              = $OUT_FILE"
#echo "XGMML_FILE            = $XGMML_FILE"
#echo "IMAGE_TITLE           = $IMAGE_TITLE"
#echo "REVERSE               = $REVERSE"
#echo "POSITIONAL_ARGS       = ${POSITIONAL_ARGS[*]}"



SBATCH_COMMAND=("$BATCH_RUNNER")
TMP_FILE_PATH=$(mktemp)
cp $TEMPLATE_BATCH_FILE $TMP_FILE_PATH

sed -i "s/###PROTEIN###/$PROTEIN/g" $TMP_FILE_PATH
sed -i "s~###XGMML_FILE###~$XGMML_FILE~g" $TMP_FILE_PATH
sed -i "s~###IMAGE_TITLE###~\"$IMAGE_TITLE\"~g" $TMP_FILE_PATH
sed -i "s~###OUT_FILE###~$OUT_FILE~g" $TMP_FILE_PATH
sed -i "s~###TAXONOMY_DATA_FILE###~$TAXONOMY_DATA_FILE~g" $TMP_FILE_PATH
sed -i "s~###GREEDY_DATA_FILE###~$GREEDY_DATA_FILE~g" $TMP_FILE_PATH
sed -i "s~###LEIDEN_DATA_FILE###~$LEIDEN_DATA_FILE~g" $TMP_FILE_PATH

if [[ "$REVERSE" == "TRUE" ]]; then
	sed -i "s~###REVERSE###~--reverse~g" $TMP_FILE_PATH
elif [[ "$REVERSE" == "FALSE" ]]; then
	sed -i "s~###REVERSE###~--no-reverse~g" $TMP_FILE_PATH
else
	sed -i "s~###REVERSE###~--no-reverse~g" $TMP_FILE_PATH
fi

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
