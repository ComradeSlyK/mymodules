# Author: Silvio Gregorini (silviogregorini@openforce.it)

import sympy as sym


class Point:

    def __init__(self, *coordinates, **kwargs):
        self._coordinates = self.convert_coordinates(coordinates)
        constrains = kwargs.get('constrains')
        if callable(constrains):
            constrains = (constrains, )
        self._constrains = tuple(constrains or [])
        self.check_constrains()

    def get_constrains(self):
        return self._constrains

    def get_coordinates(self):
        return self._coordinates

    def get_dimension(self):
        return len(self._coordinates)

    def set_constrains(self, constrains):
        self._constrains += tuple([c for c in constrains])
        self.check_constrains()

    def set_coordinates(self, coordinates):
        self._coordinates = self.convert_coordinates(coordinates)
        self.check_constrains()

    constrains = property(get_constrains, set_constrains)
    coordinates = property(get_coordinates, set_coordinates)
    dimension = property(get_dimension)

    def __abs__(self):
        return self.get_norm()

    def __add__(self, other):
        if isinstance(other, Point):
            raise ValueError("Can't sum two points.")
        # This might raise an error, let it be
        return other + self

    def __bool__(self):
        return any(self._coordinates)

    def __eq__(self, other):
        if not isinstance(other, Point):
            return False
        if self.dimension == other.dimension:
            return self._coordinates == other._coordinates
        return all(
            self[n] == other[n]
            for n in range(max(self.dimension, other.dimension))
        )

    def __getitem__(self, key):
        if isinstance(key, int):
            try:
                return self._coordinates[key]
            except IndexError:
                return 0
        elif isinstance(key, slice):
            return self._coordinates[key]
        return getattr(self, key)

    def __hash__(self):
        return hash(frozenset(self._coordinates))

    def __iter__(self):
        for c in self._coordinates:
            yield c

    def __len__(self):
        return self.dimension

    def __mul__(self, other):
        raise NotImplementedError("`__mul__` not implemented on Points")

    def __repr__(self):
        return "Point{}".format(self._coordinates)

    def __sub__(self, other):
        raise NotImplementedError("`__sub__` not implemented on Points")

    __radd__ = __add__
    __rmul__ = __mul__
    __rsub__ = __sub__
    __str__ = __repr__

    def check_constrains(self):
        for constrain in getattr(self, '_constrains', ()):
            constrain(self)

    def convert_coordinates(self, coordinates):
        """
        Converts `False`, `None`, etc to zeros, checks whether every coordinate
        fits its relative constrain (if there's any).
        :return: list of coordinates
        """
        return tuple([v if v else 0 for v in coordinates])

    def distance(self, other, p=2):
        """
        Minkowski distance.
        For p=2, we get the Euclidean distance (default behaviour).
        For p=1, we get the Manhattan distance (or discrete).
        For p=+oo or p=-oo, we get Chebyshev distance.
        """
        if not isinstance(other, Point):
            raise ValueError(
                "Distances can be computed only between two points"
            )

        if p == 0:
            raise ValueError(
                "`p` parameter can't be 0 when calculating distance"
            )
        elif p == sym.oo:
            return max(
                abs(self[n] - other[n])
                for n in range(max(self.dimension, other.dimension))
            )
        elif p == -sym.oo:
            return min(
                abs(self[n] - other[n])
                for n in range(max(self.dimension, other.dimension))
            )
        return pow(
            sum(abs(self[n] - other[n]) ** p
                for n in range(max(self.dimension, other.dimension))),
            (1 / p)
        )

    def get_norm(self, p=2):
        return self.distance(Point(), p)

    def to_vector(self):
        from mymodules.mygeometry.vectors import Vector
        return Vector(*self._coordinates, constrains=self._constrains)
