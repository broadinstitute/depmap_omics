
for s in $(cut -f 6 -d, all.stats | sed 1d); do
if [ ! -f $(basename $s) ]; then
gsutil cp $s .
fi
done
