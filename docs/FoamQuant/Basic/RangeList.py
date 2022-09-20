def RangeList(i1, i2, verbose=False):
    import numpy as np
    List = np.uint8(np.linspace(i1,i2,i2-i1+1))
    if verbose:
        print(List)
    return List