# mapLoopLoci

## Apply conservation class labels to loop loci mapped across species and/or cell types.

Map loop loci, in bedpe format, across cells and/or species. Input loops are
labelled according to whether they are conserved (i.e., left and right anchors
in the query both map to left and right anchors in the same target loop,
allowing for possible inversions), partially conserved (one anchor is used in
both species/cells), or species/cell-specific (neither anchor is used in the
other species). Cross-species mappings are performed using the bnMapper
algorithm, to map features from the target species to the query species of a
chain alignment file.

## Usage
```mapLoopLoci.py [-h] [-s] [-o FILE] [-t FLOAT] [-g GAP] [-v {info,debug,silent}] [-k] [-m MIN_OVERLAP] [-w SLOP] query target alignment```


## Conservation classes
 | Class | Description |
 |-------|-------------|
 |C | Conserved. Both anchors map and are assigned to the same target loop.|
 |B2 | Both anchors map, but are assigned to different target loops.|
 |B1 | Both anchors map, but only one is assigned to a target loop.|
 |B0 | Both anchors map, but neither is assigned to a target loop.|
 |N1A | Semi-species-specific. Only one anchor maps to the target genome and the mapped anchor is used in target as a loop anchor).|
 |N1B | Species-specific, semi-mapping. Only one anchor maps to the target genome and the mapped anchor is not used in target as a loop anchor.|
 |N0 | Species-specific, non-mapping. Neither anchor maps to the target genome.|

## Positional Arguments
 | Argument | Description |
 |----------|-------------|
 |query | Input loops for the query species/cell.|
 |target | Input loops for the target species/cell.|
 |alignment | Alignment file (.chain or .pkl)|

## Optional Arguments
 | Short option | Long option | Argument Type |Description |
 |--------------|-------------|---------------|------------|
 | -h| --help | | Show help message and exit.|
 | -s| --same_species | | Query and target loops are from the same species. Cross-species mapping step will be skipped. (default: False)|
 | -o| --output | FILE | Output file. Mandatory if more than on file in input. (default: stdout)|
 | -t| --threshold | FLOAT | Mapping threshold i.e., (elem * threshold) <= (mapped_elem) (default: 0.0)|
 | -g | --gap | INT | Ignore elements with an insertion/deletion of this or bigger size. (default: -1 -- accept all gapped alignments)|
 | -v | --verbose | {info,debug,silent} | Verbosity level (default: info)|
 | -k | --drop_split | | If elements span multiple chains, silently drop instead of reporting the segment with the longest overlap. (This is the default behavior for bnMapper.) (default: False)|
 | -m | --min_overlap | INT | Minimum amount of overlap to consider a pair of query/target anchors as shared. (default: 1)|
 | -w | --slop | INT | Number of bases added up/downstream of query and target regions to enable flexible mapping. (default: 0)|

## Output
 | Column(s) | Description |
 |-----------|-------------|
 | 1-9 | Input BEDPE fields|
 | 10-13 | Orthologous coordinates of the LEFT anchor in the target genome.|
 | 14-17 | Orthologous coordinates of the RIGHT anchor in the target genome.|
 | 18-22 | Target loop anchor overlapping LEFT anchor of query loop. (name, anchor_in_target ("l"eft/"r"ight), chrom, chromStart, chromEnd)|
 | 23-27 | Target loop anchor overlapping RIGHT anchor of query loop. (name, anchor_in_target ("l"eft/"r"ight), chrom, chromStart, chromEnd)|
 | 28 | Conservation class assignment.|
