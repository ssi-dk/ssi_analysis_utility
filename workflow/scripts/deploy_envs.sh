SCRPTDIR=$(dirname $0)
ENVSDIR="$SCRPTDIR"/../envs

find "$ENVSDIR" -type f -iname "*yaml" -exec conda env create -f {} --yes \;