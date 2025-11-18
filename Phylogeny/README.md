
## Draft phylogeny for 153 induviduals of leatherback turtles
```
upgma_tree.R head.distance_matrix.tsv > dc_153_tree.upgma.nwk

 replace_tree_names_from_dict.py -t dc_153_tree.upgma.nwk -d <(awk -F"\t" '{print $1","$1"|"$14}' dc_153_samples_info.tsv) > renamed.dc_153_tree.upgma.nwk
```