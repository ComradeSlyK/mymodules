# Author: Silvio Gregorini (silviogregorini@openforce.it)

from sympy import Expr, Symbol, solve

from mymodules.mygeometry.points import Point


def linear_independent(*vectors):
    """
    Checks whether list of vectors are linearly independent from each other
    """
    xs = [Symbol('x{}'.format(n)) for n in range(len(vectors))]
    v = Vector.vector_sum(*[xs[n] * vectors[n] for n in range(len(vectors))])
    solutions = solve(v.coordinates, xs)
    return solutions and all(solutions.get(x) == 0 for x in xs)


def scalar_product(v1, v2):
    if all(isinstance(v, Vector) for v in (v1, v2)):
        if v1.dimension != v2.dimension:
            raise ValueError(
                "Cannot compute scalar product for vectors of different"
                " dimensions"
            )
        return sum(
            v1._coordinates[n] * v2._coordinates[n]
            for n in range(v1.dimension)
        )

    raise NotImplementedError(
        "Cannot compute scalar product between {} and {}"
        .format(v1, v2)
    )


# Compatibility for nomenclature's sake
dot_product = scalar_product


def vector_sum(*vectors):
    """ Returns a new vector as a result of the sum of given vectors """
    if not vectors:
        return Vector()
    max_dim = max(v.dimension for v in vectors)
    coordinates = [
        sum(v[n] for v in vectors)
        for n in range(max_dim)
    ]
    constrains = ()
    for v in vectors:
        constrains += v._constrains
    return Vector(*coordinates, constrains=tuple(set(constrains)))


class Vector(Point):

    def __init__(self, *coordinates, **kwargs):
        super().__init__(*coordinates, **kwargs)

    def __add__(self, other):
        if not isinstance(other, (Vector, Point)):
            raise NotImplementedError(
                "Can't sum vectors to anything but other vectors or points."
            )
        if not (self and other):
            return Vector(
                *[c for c in (self or other)._coordinates],
                constrains=(self or other)._constrains
            )
        elif self == other:
            return self * 2
        return vector_sum(self, other)

    def __mul__(self, other):
        # a * V = (aV1, aV2, ..., aVn)
        if isinstance(other, (int, float, Expr)):
            return Vector(
                *[other * v for v in self._coordinates],
                constrains=self._constrains
            )
        # V * U = V1*U1 + V2*U2 + ... + Vn*Un
        return scalar_product(self, other)

    def __repr__(self):
        return "Vector{}".format(self._coordinates)

    def __sub__(self, other):
        if not isinstance(other, (Vector, Point)):
            raise NotImplementedError(
                "Can't subtract anything but vectors from other vectors"
            )
        return self + (other * -1)

    __radd__ = __add__
    __rmul__ = __mul__
    __rsub__ = __sub__
    __str__ = __repr__

    def is_parallel_to(self, other):
        """ Returns whether `self` is parallel to `other` """
        if not (self and other):
            raise ValueError(
                "Null vector can't be compared for parallelism"
            )
        if not isinstance(other, Vector):
            raise ValueError(
                "Can compute parallelism only between vectors"
            )
        return not linear_independent(self, other)

    def is_perpendicular_to(self, other):
        """ Returns whether `self` is perpendicular to `other` """
        if not (self and other):
            raise ValueError(
                "Null vector can't be compared for perpendicularity"
            )
        if not isinstance(other, Vector):
            raise ValueError(
                "Can compute perpendicularity only between vectors"
            )
        return scalar_product(self, other) == 0

    @staticmethod
    def linear_independent(*vectors):
        return linear_independent(*vectors)

    def scalar_product(self, other):
        return scalar_product(self, other)

    def to_point(self):
        return Point(*self.coordinates, constrains=self.constrains)

    @staticmethod
    def vector_sum(*vectors):
        return vector_sum(*vectors)
