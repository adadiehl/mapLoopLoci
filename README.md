# mapLoopLoci

Apply conservation class labels to loop loci mapped across species and/or
cell types.

usage: mapLoopLoci.py [-h] [-s] [-o FILE] [-t FLOAT] [-g GAP]
                      [-v {info,debug,silent}] [-k] [-m MIN_OVERLAP] [-w SLOP]
                      query target alignment

Map loop loci, in bedpe format, across cells and/or species. Input loops are
labelled according to whether they are conserved (i.e., left and right anchors
in the query both map to left and right anchors in the same target loop,
allowing for possible inversions), partially conserved (one anchor is used in
both species/cells), or species/cell-specific (neither anchor is used in the
other species). Cross-species mappings are performed using the bnMapper
algorithm, to map features from the target species to the query species of a
chain alignment file.

Conservation classes:
C = Conserved. Both anchors map and are assigned to the same target loop.
B2 = Both anchors map, but are assigned to different target loops.
B1 = Both anchors map, but only one is assigned to a target loop
   (other present, but not used in target).
B0 = Both anchors map, but neither is assigned to a target loop
   (both present in target, but not used).
N0 = Species/cell-specific (Neither anchor maps to the target)
N1A = Species/cell-specific (One anchor does not map to the target,
    mapped anchor is used in target)
N1B = Species/cell-specific (One anchor does not map to the target,
    mapped anchor is not used in target)


positional arguments:
  query                 Input loops for the query species/cell.
  target                Input loops for the target species/cell.
  alignment             Alignment file (.chain or .pkl)

optional arguments:
  -h, --help            show this help message and exit
  -s, --same_species    Query and target loops are from the same species.
                        Cross-species mapping step will be skipped. (default:
                        False)
  -o FILE, --output FILE
                        Output file. Mandatory if more than on file in input.
                        (default: stdout)
  -t FLOAT, --threshold FLOAT
                        Mapping threshold i.e., |elem| * threshold <=
                        |mapped_elem| (default: 0.0)
  -g GAP, --gap GAP     Ignore elements with an insertion/deletion of this or
                        bigger size. (default: -1)
  -v {info,debug,silent}, --verbose {info,debug,silent}
                        Verbosity level (default: info)
  -k, --drop_split      If elements span multiple chains, silently drop
                        instead of reporting the segment with the longest
                        overlap. (This is the default behavior for bnMapper.)
                        (default: False)
  -m MIN_OVERLAP, --min_overlap MIN_OVERLAP
                        Minimum amount of overlap to consider a pair of
                        query/target anchors as shared. Default = 1. (default:
                        1)
  -w SLOP, --slop SLOP  Number of bases added up/downstream of query and
                        target regions to enable flexible mapping. Default =
                        0. (default: 0)


## Output:
Columns 1-9: Input BEDPE fields
Columns 10-13: Orthologous coordinates of the LEFT anchor in the target genome.
Columns 14-17: Orthologous coordinates of the RIGHT anchor in the target genome.
Columns 18-22: Target loop anchor overlapping LEFT anchor of query loop.
Columns 23-27: Target loop anchor overlapping RIGHT anchor of query loop.
Column 28: Conservation class assignment

For Columns 18-22 and 23-27, fields are arranged as follows:
name
anchor_in_target ("l"/"r", for left/right)
chrom
chromStart
chromEnd
