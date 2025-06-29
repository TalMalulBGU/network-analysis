
TEMPLATE_BATCH_FILE="/home/talmalu/thesis/projects/python/NetworkAnalysis/Scripts/template_batches/clustering/batch_modularity_leiden_algorithm_alignment_score_rule_batch_template.sbatch"
BATCH_RUNNER="/home/talmalu/thesis/projects/python/batch_runner.sh"

POSITIONAL_ARGS=()
ALIGNMENT_WEIGHT=""
ALIGNMENT_SCORE_RULE=""
OBJECTIVE_FUNCTION=""
N_ITERATIONS=-1
XGMML_FILE=""
PROTEIN=""
OUTPUT=""
DATABASE=""
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
	-w|--alignment_weight)
		ALIGNMENT_WEIGHT=$2
		shift
		shift
	;;
	-r|--alignment_score_rule)
		ALIGNMENT_SCORE_RULE=$2
		shift
		shift
	;;
	-i|--n_iterations)
		N_ITERATIONS=$2
		shift
		shift
	;;
	-h|--objective_function)
		OBJECTIVE_FUNCTION=$2
		shift
		shift
	;;
	-s|--start)
		START=$2
		shift
		shift
	;;
	-e|--end)
		END=$2
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
	-mem|--memory)
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
#echo "GRAPH_NAME            = $GRAPH_NAME"
#echo "ALGORITHM_FILE        = $ALGORITHM_FILE"
#echo "ALIGNMENT_WEIGHT      = $ALIGNMENT_WEIGHT"
#echo "MIN_IN_CLUSTER        = $MIN_IN_CLUSTER"
#echo "DATABASE              = $DATABASE"
#echo "XGMML_FILE            = $XGMML_FILE"
#echo "OUTPUT_FOLDER         = $OUTPUT_FOLDER"
#echo "POSITIONAL_ARGS       = ${POSITIONAL_ARGS[*]}"


SBATCH_COMMAND=("$BATCH_RUNNER")
TMP_FILE_PATH=$(mktemp)
cp $TEMPLATE_BATCH_FILE $TMP_FILE_PATH

sed -i "s/###PROTEIN###/$PROTEIN/g" $TMP_FILE_PATH
sed -i "s/###ALIGNMENT_WEIGHT###/$ALIGNMENT_WEIGHT/g" $TMP_FILE_PATH
sed -i "s/###MIN_ALIGNMENT_SCORE###/$START/g" $TMP_FILE_PATH
sed -i "s/###MAX_ALIGNMENT_SCORE###/$END/g" $TMP_FILE_PATH
sed -i "s~###OUTPUT###~$OUTPUT~g" $TMP_FILE_PATH
sed -i "s~###XGMML_FILE###~$XGMML_FILE~g" $TMP_FILE_PATH

if [[ ! -z "$DATABASE" ]]
then
	sed -i "s~###DATABASE###~$DATABASE~g" $TMP_FILE_PATH
else
	sed -i "/###DATABASE###/d" $TMP_FILE_PATH
fi

if [[ ! -z "$ALIGNMENT_SCORE_RULE" ]]
then
	sed -i "s~###ALIGNMENT_SCORE_RULE###~$ALIGNMENT_SCORE_RULE~g" $TMP_FILE_PATH
else
	sed -i "/###ALIGNMENT_SCORE_RULE###/d" $TMP_FILE_PATH
fi

if [[ ! -z "$N_ITERATIONS" ]]
then
	sed -i "s~###N_ITERATIONS###~$N_ITERATIONS~g" $TMP_FILE_PATH
else
	sed -i "/###N_ITERATIONS###/d" $TMP_FILE_PATH
fi

if [[ ! -z "$OBJECTIVE_FUNCTION" ]]
then
	sed -i "s~###OBJECTIVE_FUNCTION###~$OBJECTIVE_FUNCTION~g" $TMP_FILE_PATH
else
	sed -i "/###OBJECTIVE_FUNCTION###/d" $TMP_FILE_PATH
fi


if [[ ! -z "$DEPENDENCY" ]]
then
	SBATCH_COMMAND+=("--dependency $DEPENDENCY")
fi

if [[ -z "$MEMORY" ]]
then
	MEMORY=$(du -bcs $XGMML_FILE | grep total | awk '{printf "%dG",int($1/2^30 + 1) * 5 }')
fi

SBATCH_COMMAND+=("-a" "0-$((${END} - ${START}))")

SBATCH_COMMAND+=("--memory $MEMORY")

SBATCH_COMMAND+=("$TMP_FILE_PATH")

JOB_ID=$(${SBATCH_COMMAND[@]}) #run command
echo $JOB_ID
