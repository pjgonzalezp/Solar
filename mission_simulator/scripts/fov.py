"""
FOV class
"""
import numpy as np

__all__ = ['FOV']


class FOV(object):
    """
    Field of View (FOV) class

    The shape of the FOV dictates the allowed size parameter(s), which is any number type that can
    be compared (`~astropy.units.Quantity` is recommended).

    The `in` operator will check whether a provided 2-element offset from the center of the FOV is
    contained within the FOV.

    Parameters
    ----------
    shape : str
        Must be 'circle', 'rectangle', or 'square'.  Defaults to 'circle'.
    diameter :
        Valid only if ``shape`` is 'circle'. Must be a named argument.
    width :
        Valid only if ``shape`` is 'rectangle' or 'square'. Must be a named argument.
    height :
        Valid only if ``shape`` is 'rectangle' or 'square'. Must be a named argument.

    Examples
    --------
    >>> circle_fov = FOV('circle', diameter=6*u.arcmin)
    >>> square_fov = FOV('square', width=9.8*u.arcmin)
    >>> (4*u.arcmin, -3*u.arcmin) in circle_fov
    False
    >>> (4*u.arcmin, -3*u.arcmin) in square_fov
    True
    """
    def __init__(self, shape='circle', *, diameter=None, width=None, height=None):
        self.shape = shape
        if self.shape == 'circle':
            if diameter is None:
                raise ValueError("For a circular FOV, `diameter` must be provided")
            self.diameter = diameter
        elif self.shape == 'rectangle':
            if width is None or height is None:
                raise ValueError("For a rectangular FOV, both `width` and `height` must be provided")
            self.width, self.height = width, height
        elif self.shape == 'square':
            if width is None and height is None:
                raise ValueError("For a rectangular FOV, either `width` or `height` must be provided")
            if width is not None and height is not None and width != height:
                raise ValueError("For a square FOV, if both `width` and `height` are provided, they must be equal")
            self.width, self.height = (width, width) if width is not None else (height, height)
        else:
            raise ValueError("Unknown FOV type")

    def __contains__(self, offset_xy):
        if self.shape == 'circle':
            return offset_xy[0]**2 + offset_xy[1]**2 < (self.diameter/2)**2
        elif self.shape in ['rectangle', 'square']:
            return np.abs(offset_xy[0]) < self.width/2 and np.abs(offset_xy[1]) < self.height/2

    def __repr__(self):
        if self.shape == 'circle':
            description = f"Circular FOV of diameter {self.diameter}"
        elif self.shape == 'rectangle':
            description = f"Rectangular FOV with width {self.width} and height {self.height}"
        elif self.shape == 'square':
            description = f"Square FOV with width/height {self.width}"

        return description
