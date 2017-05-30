bash_source="$BASH_SOURCE" 
printf "bash_source: %s\n" "$bash_source"

dir_name="$( dirname ${BASH_SOURCE[0]} )"
printf "dir_name: %s\n" "$dir_name"

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
printf "dir: %s\n" "$DIR"
