#!/usr/bin/env python

"""
Map loop loci, in bedpe format, across cells and/or species. Input loops are labelled
according to whether they are conserved (i.e., left and right anchors in the query 
both map to left and right anchors in the same target loop, allowing for possible
inversions), partially conserved (one anchor is used in both species/cells), or
species/cell-specific (neither anchor is used in the other species).

Cross-species mappings are performed using the bnMapper algorithm, to map
features from the target species to the query species of a chain alignment file.

Loop classifications are reported in the final column:
N0 = Species/cell-specific (Neither anchor maps to the target)
N1A = Species/cell-specific (One anchor does not map to the target, mapped anchor is used in target)
N1B = Species/cell-specific (One anchor does not map to the target, mapped anchor is not used in target)
C = Conserved. Both anchors map and are assigned to the same target loop.
B2 = Both anchors map, but are assigned to different target loops.
B1 = Both anchors map, but only one is assigned to a target loop (other present, but not used in target).
B0 = Both anchors map, but neither is assigned to a target loop (both present in target, but not used).

"""
import logging
import os
import sys
import re
from itertools import groupby
from operator import attrgetter, concat, itemgetter

import numpy as np
from six.moves import reduce

from bx.align import epo
from bx.align.epo import bed_union as elem_u
from bx.cookbook import argparse
from bx.intervals.intersection import IntervalTree, Interval

elem_t_bed4 = np.dtype([('chrom', np.str_, 30), ('start', np.int64), ('end', np.int64), ('id', np.str_, 500)])
elem_t = np.dtype([('chrom1', np.str_, 30), ('start1', np.int64), ('end1', np.int64),
                   ('chrom2', np.str_, 30), ('start2', np.int64), ('end2', np.int64),
                   ('id', np.str_, 500), ('IAB', np.int64), ('FDR', np.float),
                   ('mapped_chr_l', np.str_, 30), ('mapped_start_l', np.int64), ('mapped_end_l', np.int64), ('qStrand_l', np.str_, 1),
                   ('mapped_chr_r', np.str_, 30), ('mapped_start_r', np.int64), ('mapped_end_r', np.int64), ('qStrand_r', np.str_, 1),
                   ('id_l', np.str_, 500), ('anchor_l', np.str_, 1), ('chrom_l', np.str_, 30), ('start_l', np.int64), ('end_l', np.int64),
                   ('id_r', np.str_, 500), ('anchor_r', np.str_, 1), ('chrom_r', np.str_, 30), ('start_r', np.int64), ('end_r', np.int64),
                   ('class', np.str_, 2)])

LOG_LEVELS = {"info" : logging.INFO, "debug" : logging.DEBUG, "silent" : logging.ERROR}

logging.basicConfig()
log = logging.getLogger()

class GIntervalTree( IntervalTree ):
    """a set of IntervalTrees that is indexed by chromosomes"""

    def __init__(self, data=[]):
        self._trees = {}

    def add(self, chrom, element):
        """insert an element. use this method as the IntervalTree one.
        this will simply call the IntervalTree.add method on the right tree

        :param chrom: chromosome
        :param element: the argument of IntervalTree.insert_interval
        :return: None
        """

        self._trees.setdefault(chrom, IntervalTree()).insert_interval( element )

    def find(self, chrom, start, end):
        """find the intersecting elements

        :param chrom: chromosome
        :param start: start
        :param end: end
        :return: a list of intersecting elements"""

        tree = self._trees.get( chrom, None )
        if tree:
            return tree.find( start, end )
        #return always a list
        return []


def map_loops(ELEM_Q, ELEM_T, out_f, EPO, TREE, opt):
    """
    Map query loops to target loops, performing cross-species mapping
    if necessary.
    """
    FORMAT_STR = "%s\t%d\t%d\t%s\t%d\t%d\t%s\t%d\t%f\t%s\t%d\t%d\t%s\t%s\t%d\t%d\t%s\t%s\t%s\t%s\t%d\t%d\t%s\t%s\t%s\t%d\t%d\t%s\n"

    for from_elem in ELEM_Q:
        #sys.stderr.write("{}\n".format(from_elem))
        out_elem = list(from_elem)
        #sys.stderr.write("{}\n".format(out_elem))

        # If we are mapping across species, attempt to lift over the query coordinates to the target.
        if not opt.same_species:
            m = transform_elem(from_elem, EPO, TREE, opt)
            #sys.stderr.write("{}\n".format(m))
            # Check mappability and only proceed with loop comparison if at least one anchor mapped
            if m[0] == "." and m[6] == ".":
                out_f.write(FORMAT_STR % tuple(out_elem))
                continue
            else:
                elem_mapped = list(from_elem)
                elem_mapped[0] = m[0]
                elem_mapped[1] = m[1]
                elem_mapped[2] = m[2]
                elem_mapped[3] = m[6]
                elem_mapped[4] = m[7]
                elem_mapped[5] = m[8]

                # Store left-anchor mapping in fields 9-12
                if not m[0] == ".":
                    out_elem[9] = m[0]
                    out_elem[10] = m[1]
                    out_elem[11] = m[2]
                    out_elem[12] = m[3]
                # Store right-anchor mapping in fields 13-16
                if not m[6] == ".":
                    out_elem[13] = m[6]
                    out_elem[14] = m[7]
                    out_elem[15] = m[8]
                    out_elem[16] = m[9]

        else:
            elem_mapped = from_elem
            out_elem[9] = out_elem[10] = -1  # Use -1 as mappability indicator for same species

        # Compare (mapped) query loop to target loops to determine conservation/specificity
        for to_elem in ELEM_T:
            #sys.stderr.write("{}\n".format(to_elem))
            elem_q_l = [elem_mapped[0], elem_mapped[1], elem_mapped[2]]
            elem_q_r = [elem_mapped[3], elem_mapped[4], elem_mapped[5]]
            elem_t_l = [to_elem[0], to_elem[1], to_elem[2]]
            elem_t_r = [to_elem[3], to_elem[4], to_elem[5]]

            ll = overlaps(elem_q_l, elem_t_l, opt, opt.slop, 0)
            rr = overlaps(elem_q_r, elem_t_r, opt, opt.slop, 0)
            rl = overlaps(elem_q_r, elem_t_l, opt, opt.slop, 0)
            lr = overlaps(elem_q_l, elem_t_r, opt, opt.slop, 0)

            # If one/both loop anchor(s) maps to a target loop anchor, store the
            # target loop data in fields 17-21 (left mappings) and 22-26 (right mappings)
            if ll or rr or rl or lr:
                if ll or lr:
                    # Left anchor mapping
                    out_elem[17] = to_elem[6]  # Mapped loop ID
                    if ll:
                        out_elem[18] = "l"
                        out_elem[19] = to_elem[0]
                        out_elem[20] = to_elem[1]
                        out_elem[21] = to_elem[2]
                    else:
                        out_elem[18] = "r"
                        out_elem[19] = to_elem[3]
                        out_elem[20] = to_elem[4]
                        out_elem[21] = to_elem[5]
                        
                        
                if rr or rl:
                    # right anchro mapping
                    out_elem[22] = to_elem[6]  # Mapped loop ID                                                                                                                                         
                    if rl:
                        out_elem[23] = "l"
                        out_elem[24] = to_elem[0]
                        out_elem[25] = to_elem[1]
                        out_elem[26] = to_elem[2]
                    else:
                        out_elem[23] = "r"
                        out_elem[24] = to_elem[3]
                        out_elem[25] = to_elem[4]
                        out_elem[26] = to_elem[5]
            if out_elem[17] != "." and out_elem[22] != "." and out_elem[17] == out_elem[22]:
                """
                Break out of the loop if we've assigned both left and right anchors, and
                and the loop is conserved. Continue searching if the loop is not conserved.
                This gives priority to loop conservation.
                """
                out_elem[27] = "C"
                break

        # Figure out how to classify loops based on mapping and shared loop anchor(s)
        if (out_elem[9] != "." and out_elem[13] != ".") or opt.same_species:
            # Both ends map: both query anchor sequences are present in the target
            if out_elem[17] == "." and out_elem[22] == ".":
                # Neither query anchor overlaps a target anchor
                out_elem[27] = "B0"
            elif out_elem[17] == "." or out_elem[22] == ".":
                # One query anchor overlaps a target anchor
                out_elem[27] = "B1"
            elif out_elem[17] != out_elem[22]:
                # Both query anchors overlap target anchors, but the target anchors are from different loops
                out_elem[27] = "B2"
            else:
                # Conserved loop. This should be caught above, but just in case...
                out_elem[27] = "C"
        elif out_elem[9] == "." and out_elem[13] == ".":
            # Neither end maps: both query loop anchors are in query-specific sequence. Leave the default value.
            pass
        elif out_elem[9] == "." or out_elem[13] == ".":
            # Only one query anchor maps: one query loop anchor is present in the target, one is not.
            if out_elem[17] == "." and out_elem[22] == ".":
                # The mappable anchor does not map to a target loop anchor
                out_elem[27] = "N1B"
            else:
                # The mappable anchor does map to a target anchor
                out_elem[27] = "N1A"
        #sys.stderr.write("{}\n".format(out_elem))
        out_f.write(FORMAT_STR % tuple(out_elem))

def transform_elem(elem_qry, EPO, TREE, opt):

    elems = np.array([(elem_qry[0], elem_qry[1], elem_qry[2], "elem_1"),
                      (elem_qry[3], elem_qry[4], elem_qry[5], "elem_2")],
                     dtype=elem_t_bed4)

    elems_mapped = [".", -1, -1, ".", -1, -1, ".", -1, -1, ".", -1, -1]
    for chrom in set( elems['chrom'] ):
        els = (transform_by_chrom(EPO,
                                  elems[elems['chrom'] == chrom],
                                  TREE, chrom, opt))
        #elems_mapped = []
        if len(els) > 2:
            log.debug("%s: maps to more than two locations\n" % elem_qry)
        for el in els:
            #print el                                                                                                                                                         
            el_s = list(el[0])
            el_e = list(el[len(el)-1])
            #mapped_el = [el_s[0], el_s[1], el_e[2], el_s[3]]                                                                                                                 
            #print mapped_el                                                                                                                                                  
            #elems_mapped.append(mapped_el)                                                                                                                                   
            if el_s[3] == "elem_1":  # left anchor                                                                                                                            
                elems_mapped[0] = el_s[0]
                elems_mapped[1] = el_s[1]
                elems_mapped[2] = el_e[2]
                elems_mapped[3] = el_s[4]
                elems_mapped[4] = el_s[5]
                elems_mapped[5] = el_s[6]
            else:  # right anchor                                                                                                                                             
                elems_mapped[6] = el_s[0]
                elems_mapped[7] = el_s[1]
                elems_mapped[8] = el_e[2]
                elems_mapped[9] = el_s[4]
                elems_mapped[10] = el_s[5]
                elems_mapped[11] = el_s[6]
    #print elems_mapped                                                                                                                                                       

    return elems_mapped


def overlaps(elem_qry, elem_tgt, opt, slop = 0, debug_flag = 0):
    "See if two intervals overlap"

    "First check that chromosomes match"
    if elem_qry[0] == elem_tgt[0]:
        if debug_flag:
            log.debug("\toverlaps: %s\t%s\t", elem_qry, elem_tgt)
        """ Next see if intervals overlap. slop allows for some                                                                                                        
        flexibility in exact positioning of overlaps. """
        if (elem_qry[1] + slop >= elem_tgt[1] - slop and \
            elem_qry[1] - slop <= elem_tgt[2] + slop) or \
            \
           (elem_qry[2] + slop >= elem_tgt[1] - slop and \
            elem_qry[2] - slop <= elem_tgt[2] + slop) or \
            \
           (elem_qry[1] <= elem_tgt[1] + slop and \
            elem_qry[2] >= elem_tgt[2] - slop) or \
            \
           (elem_qry[1] + slop >= elem_tgt[1] and \
            elem_qry[2] - slop <= elem_tgt[2]):
           # check overlap length                                                                                                                                             
            ol = min(elem_qry[2] + slop, elem_tgt[2] + slop) - max(elem_qry[1] - slop, elem_tgt[1] - slop)
            if debug_flag:
                log.debug("\t%s\t%s\t%s\n", elem_qry, elem_tgt, ol)
            "Finally, see if the overlap passes the minimum threshold"
            if ol >= opt.min_overlap:
                return 1
    "Return 0 if any conditions are not met"
    return 0

    
    
def transform(elem, chain_CT_CQ, max_gap):
    """transform the coordinates of this elem into the other species.

    elem intersects this chain's ginterval.
    :return: a list of the type [(to_chr, start, end, elem[id]) ... ]"""
    (chain, CT, CQ) = chain_CT_CQ
    start, end = max(elem['start'], chain.tStart) - chain.tStart, min(elem['end'], chain.tEnd) - chain.tStart

    assert np.all( (CT[:,1] - CT[:,0]) == (CQ[:,1] - CQ[:,0]) )
    to_chrom = chain.qName
    to_gab_start = chain.qStart

    start_idx = np.where( CT[:,1] > start )[0][0]
    end_idx = np.where( CT[:,0] < end )[0][-1]

    if start_idx > end_idx: #maps to a gap region on the other species
        return []

    ## apply the gap threshold
    if max_gap >= 0 and start_idx < end_idx - 1:
        if np.max(CT[(start_idx+1):end_idx,0] - CT[start_idx:(end_idx-1),1]) > max_gap or np.max(CQ[(start_idx+1):end_idx,0] - CQ[start_idx:(end_idx-1),1]) > max_gap:
            return []

    assert start < CT[start_idx, 1]
    assert  CT[end_idx, 0] < end
    to_start = CQ[start_idx, 0] + max(0, start - CT[start_idx,0]) # correct if on middle of interval
    to_end = CQ[end_idx, 1] - max(0, CT[end_idx, 1] - end)        # idem

    if start_idx == end_idx: #elem falls in a single run of matches
        slices = [(to_start, to_end)]
    else:
        slices = [(to_start, CQ[start_idx,1])]
        slices += [(CQ[i,0], CQ[i,1]) for i in range(start_idx+1, end_idx)]
        slices.append( (CQ[end_idx,0], to_end) )
    if chain.qStrand == '-':
        Sz = chain.qEnd - chain.qStart
        slices =  [(Sz-t[1], Sz-t[0]) for t in slices]
    return [(to_chrom, to_gab_start + t[0], to_gab_start + t[1], elem['id'], chain.qStrand, chain.tStart, chain.tEnd) for t in slices]

def union_elements(elements):
    """elements = [(chr, s, e, id), ...], this is to join elements that have a
    deletion in the 'to' species
    """

    if len(elements) < 2: return elements
    assert set( [e[3] for e in elements] ) == set( [elements[0][3]] ), "more than one id"
    el_id = elements[0][3]
    el_strand = elements[0][4]
    el_tStart = elements[0][5]
    el_tEnd = elements[0][6]

    unioned_elements = []
    for ch, chgrp in groupby(elements, key=itemgetter(0)):
        for (s, e) in elem_u( np.array([itemgetter(1, 2)(_) for _ in chgrp], dtype=np.uint) ):
            if (s < e):
                unioned_elements.append( (ch, s, e, el_id, el_strand, el_tStart, el_tEnd) )
    assert len(unioned_elements) <= len(elements)
    return unioned_elements

def transform_by_chrom(all_epo, from_elem_list, tree, chrom, opt):
    assert len( set(from_elem_list['chrom']) ) <= 1

    elems_mapped = []

    mapped_elem_count = 0
    for from_elem in from_elem_list:
        
        matching_block_ids = [attrgetter("value")(_) for _ in tree.find(chrom, from_elem['start'], from_elem['end'])]

        # do the actual mapping
        to_elem_slices = [_ for _ in (transform(from_elem, all_epo[i], opt.gap) for i in matching_block_ids) if _]
        """ # Original version: silently discard split alignments
        if len(to_elem_slices) > 1 or len(to_elem_slices) == 0:
            log.debug("%s no match or in different chain/chromosomes" % (str(from_elem)))
            continue
        to_elem_slices = to_elem_slices[0]
        """
        """ Modified version below allows liftOver-like behavior of
        keeping the longest alignment when alignments are split across
        multiple chains. Added by Adam Diehl (adadiehl@umich.edu)
        """
        max_elem_idx = 0
        if len(to_elem_slices) == 0:
            log.debug("%s: no match in target: discarding." % (str(from_elem)))
            continue
        elif len(to_elem_slices) > 1 and not opt.drop_split:
            log.debug("%s spans multiple chains/chromosomes. Using longest alignment." % (str(from_elem)))
            max_elem_len = 0
            for i in xrange(len(to_elem_slices)):
                elem_len = to_elem_slices[i][-1][2] - to_elem_slices[i][0][2]
                if elem_len > max_elem_len:
                    max_elem_len = elem_len
                    max_elem_idx = i
        elif len(to_elem_slices) > 1:
            log.debug("%s spans multiple chains/chromosomes: discarding." % (str(from_elem)))
            continue
        to_elem_slices = to_elem_slices[max_elem_idx]
        """ End AGD modifications """

        # apply threshold
        if (from_elem[2] - from_elem[1]) * opt.threshold > reduce(lambda b,a: a[2]-a[1] + b, to_elem_slices, 0):
            log.debug("%s did not pass threshold" % (str(from_elem)))
            continue

        # if to_species had insertions you can join elements
        to_elem_list = sorted(union_elements(to_elem_slices), key=lambda a: a[1])
        if to_elem_list:
            mapped_elem_count += 1
            log.debug("\tjoined to %d elements" % (len(to_elem_list)))
            start = to_elem_list[0][1]
            end = to_elem_list[-1][2]
            elems_mapped.append(to_elem_list)
            
    log.debug("%s: %d of %d elements mapped" % (chrom, mapped_elem_count, from_elem_list.shape[0]))
    return elems_mapped

    
def loadChains(path):
    "name says it."

    EPO = epo.Chain._parse_file(path, True)
    ## convert coordinates w.r.t the forward strand (into slices)
    ## compute cummulative intervals
    for i in range( len(EPO) ):
        ch, S, T, Q = EPO[i]
        if ch.tStrand == '-':
            ch = ch._replace(tEnd = ch.tSize - ch.tStart,
                    tStart = ch.tSize - ch.tEnd)
        if ch.qStrand == '-':
            ch = ch._replace(qEnd = ch.qSize - ch.qStart,
                    qStart = ch.qSize - ch.qEnd)
        EPO[i] = (ch,
                epo.cummulative_intervals(S, T),
                epo.cummulative_intervals(S, Q)
                )
    ##now each element of epo is (chain_header, target_intervals, query_intervals)
    assert all( t[0].tStrand == '+' for t in EPO ), "all target strands should be +"
    return EPO

def loadFeatures(path, opt):
    """
    Load features. For BED, only BED4 columns are loaded.
    For narrowPeak, all columns are loaded.
    """

    log.info("loading from %s ..." % path)
    data = []
    with open(path) as fd:
        for line in fd:
            cols = line.split()
            # Correct any chromosome records that don't have the "chr" prefix
            if re.match('(?!chr)', cols[0]):
                cols[0] = "%s%s" % ("chr", cols[0])
            if re.match('(?!chr)', cols[3]):
                cols[3] = "%s%s" % ("chr", cols[3])
            data.append( (cols[0], cols[1], cols[2],
                          cols[3], cols[4], cols[5],
                          cols[6], cols[7], cols[8],
                          ".", -1, -1, ".",
                          ".", -1, -1, ".",
                          ".", ".", ".", -1, -1,
                          ".", ".", ".", -1, -1,
                          "N0")
            )
    data = np.array(data, dtype=elem_t)
    return data

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__, epilog="Adam Diehl (Boyle Lab)",
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("query",
            help="Input loops for the query species/cell.")
    parser.add_argument("target",
            help="Input loops for the target species/cell.")
    parser.add_argument("alignment", help="Alignment file (.chain or .pkl)")

    parser.add_argument("-s", '--same_species', default=False, action='store_true',
            help="Query and target loops are from the same species. Cross-species mapping step will be skipped.")
    parser.add_argument("-o", '--output', metavar="FILE", default='stdout',
            type=lambda s: ((s in ('stdout', '-') and "/dev/stdout") or s),
            help="Output file. Mandatory if more than on file in input.")
    parser.add_argument("-t", '--threshold', metavar="FLOAT", default=0., type=float,
            help="Mapping threshold i.e., |elem| * threshold <= |mapped_elem|")
    parser.add_argument('-g', '--gap', type=int, default=-1,
            help="Ignore elements with an insertion/deletion of this or bigger size.")
    parser.add_argument('-v', '--verbose', type=str, choices=list(LOG_LEVELS.keys()), default='info',
            help='Verbosity level')
    parser.add_argument("-k", '--drop_split', default=False, action='store_true',
            help="If elements span multiple chains, silently drop instead of reporting the segment with the longest overlap. (This is the default behavior for bnMapper.)")
    parser.add_argument("-m", '--min_overlap', type=int, default=1,
            help="Minimum amount of overlap to consider a pair of query/target anchors as shared. Default = 1.")
    parser.add_argument("-w", '--slop', type=int, default=0,
            help="Number of bases added up/downstream of query and target regions to enable flexible mapping. Default = 0.")

    opt = parser.parse_args()
    log.setLevel(LOG_LEVELS[opt.verbose])

    EPO = []
    TREE = []
    if not opt.same_species:
        #loading alignments from opt.alignment
        EPO = dict( (ch[0].id, ch) for ch in loadChains(opt.alignment) )

        ## create an interval tree based on chain headers (from_species side)
        ## for fast feature-to-chain_header searching
        log.info("indexing %d chains ..." % (len(EPO),))
        TREE = GIntervalTree()
        for gabid in EPO:
            chain, t, q = EPO[gabid]
            TREE.add(chain.tName, Interval(chain.tStart, chain.tEnd, chain.id))

    # Load up query and target files
    ELEM_Q = loadFeatures( opt.query, opt )
    ELEM_T = loadFeatures( opt.target, opt )

    # Map loops from query file to target file.
    with open(opt.output, 'w') as out_f:
        map_loops(ELEM_Q, ELEM_T, out_f, EPO, TREE, opt)
