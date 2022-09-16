""
Library of SPAM functions for manipulating images
Copyright (C) 2020 SPAM Contributors

This program is free software: you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by the Free
Software Foundation, either version 3 of the License, or (at your option)
any later version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
more details.

You should have received a copy of the GNU General Public License along with
this program.  If not, see <http://www.gnu.org/licenses/>.
"""

import numpy
import tifffile
import spam.label
import progressbar

[docs]
def stackToArray(prefix, sufix='.tif', stack=range(10), digits='05d', erosion=False, verbose=False):
    """
    Convert of stack of 2D sequential tif images into a 3D array.

    Parameters
    ----------
        prefix : string
            The common name of the 2D images files before the sequential number

        sufix : string, default='.tif'
            The common name and extension of the 2D images after the sequential number

        stack : sequence, default=range(10)
            The numbers of the slices with no formating (with no leading zeros)

        digits : string, default='05d'
            The format (number of digits) of the numbers (add leading zeros).

        erosion : bool, default=None
            Apply an erosion of 1px to the mask in order to avoid border noise.

        verbose: bool, default=False
            Verbose mode if True.

    Returns
    -------
        im : array
            The 3D image

        mask : array
            The 3D mask (full of 1 if ``createMask=None``)
    """

    if verbose:
        print("spam.helpers.imageManipulation.stackToArray:")
        print('\tfrom: {p}{first:{d}}{s}'.format(p=prefix, s=sufix, first=stack[0], d=digits))
        print('\tto:   {p}{last:{d}}{s}'.format(p=prefix, s=sufix, last=stack[-1], d=digits))

    # Step 1 If nBytes is not defined: we open the first slice just for the dimensions
    slice_name = '{h}{s:{d}}{t}'.format(h=prefix, t=sufix, s=stack[0], d=digits)
    slice_im = tifffile.imread(slice_name)

    # Step 2 compute the dimension and create the 3D array
    ny = slice_im.shape[0]
    nx = slice_im.shape[1]
    nz = len(stack)
    im = numpy.zeros([nz, ny, nx], dtype=slice_im.dtype)

    # Step 2.1 create empty mask (8 bits)
    mask = numpy.ones([nz, ny, nx], dtype='<u1')

    # Step 3 stack all the slices
    for i, s in enumerate(stack):
        slice_name = '{h}{s:{d}}{t}'.format(h=prefix, t=sufix, s=s, d=digits)
        if verbose:
            print('\tStack slice number {}/{} ({s:{d}})'.format(i + 1, len(stack), s=s, d=digits))
        im[i, :, :] = tifffile.imread(slice_name)
        # we fill the mask
        #if createMask:
            #mask[i, :, :] = _mask2D(im[i, :, :], erosion=erosion)

    return im, mask
