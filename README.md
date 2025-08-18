# nVennPy
This package adds a Python interface to the `nVenn2` algorithm to create generalized, quasi-proportional Venn diagrams. 

## The problem
We have several `sets` composed of `elements`, like gene symbols. Each gene symbol can belong to one or more sets, which define `regions`. A region is defined by the sets it belongs to and the sets it does not belong to. 

A proportional Venn diagrams shows sets inside lines that may intersect and define regions. The area of each region is approximately proportional to the number of gene symbols that belong to that region. The nVenn2 algorithm is also generalized, as it can be used on any number of sets. In practice, a Venn diagram with more that six sets is not practical. However, more sets can be used if most regions are empty (see example 2).

## Input


# Example 1

    import nvenn2
    n = nvenn2.diagram("Set1 TP53 SF3B1 POT1 \nSet2 TP53 KRAS NRAS\nSet3 SF3B1 POT1 LMNA\nSet4 TP53 KRAS SF3B1", 2)
    n.simulate
    print(n.tosvg())
