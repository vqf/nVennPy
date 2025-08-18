# nVennPy
This package adds a Python interface to the `nVenn2` algorithm to create generalized, quasi-proportional Venn diagrams.

# Example

    import nvenn2
    n = nvenn2.diagram("a b c\n1 2 3\n2 4", 0)
    n.simulate
    print(n.tosvg())
