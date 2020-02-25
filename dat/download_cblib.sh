mkdir cblib
cat instances.txt | parallel -j1 "wget -P cblib http://cblib.zib.de/download/all/{}.cbf.gz"
gunzip cblib/*.cbf.gz