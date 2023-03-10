def strindex(i, n0):
    """
    Return str index written on n0 digit
    
    :param i: index 
    :type i: int
    :param n0: number of 0 digit
    :type n0: int
    :return: str index
    """    
    
    istr = str(i)
    istrlen = len(istr)
    fullistr = (n0-istrlen)*'0'+istr
    return fullistr


def RangeList(i1, i2, verbose=False):
    """
    Return an index range
    
    :param i1: first index 
    :type i1: int
    :param i2: last index 
    :type i2: int
    :param verbose: If True, print the range
    :type verbose: Bool
    :return: int numpy array
    """
    
    import numpy as np
    List = np.uint8(np.linspace(i1,i2,i2-i1+1))
    if verbose:
        print(List)
    return List
