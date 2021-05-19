join_rec() {
    if [ $# -eq 1 ]; then
      join -1 1 -t $'\t' - "$1"
    else
      f=$1; shift
      join -1 1 -t $'\t' - "$f" | join_rec "$@"
    fi
}
if [ $# -eq 1 ]; then
    echo "$@"
elif [ $# -le 2 ]; then
    join -1 1 -t $'\t' "$@"
else
    f1=$1; f2=$2; shift 2
    join -1 1 -t $'\t' "$f1" "$f2" | join_rec "$@"
fi
