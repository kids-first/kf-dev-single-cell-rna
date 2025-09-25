tar xvf $1 \
--ignore-failed-read \
--wildcards "*_bc_matrix*" "*clustering*" \
--transform  's%.*/\([^/]*/count/.*\)%\1%' \
--transform 's%count%outs%' \
--transform 's%sample_\([filtered|raw]\)%\1%g' \
--exclude '*outs/multi*' \
--show-transformed-names
