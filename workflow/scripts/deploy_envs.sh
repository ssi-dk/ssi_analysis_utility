SCRPTDIR=$(dirname $0)
YAMLDIR="$SCRPTDIR"/../envs

# Find all yaml files in the envs/ dir and create environments from conda
find "$YAMLDIR" -type f -iname "*yaml" -exec conda env create -f {} \;
