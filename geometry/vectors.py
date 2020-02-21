# Author: Silvio Gregorini (silviogregorini@openforce.it)

from sympy import Expr, Symbol, solve
from mymodules.geometry.points import Point


class Vector(Point):

    def __init__(self, *coordinates, **kwargs):
        super().__init__(*coordinates, **kwargs)

    def __add__(self, other):
        if not isinstance(other, (Vector, Point)):
            raise NotImplementedError(
                "Can't sum vectors to anything but other vectors."
            )
        if not (self and other):
            return Vector(
                *[c for c in (self or other)._coordinates],
                constrains=(self or other)._constrains
            )
        elif self == other:
            return 2 * self
        return Vector(
            *[self.get_nth_coordinate(n) + other.get_nth_coordinate(n)
              for n in range(max(self.dimension, other.dimension))],
            constrains=tuple(
                set(tuple(self._constrains) + tuple(other._constrains))
            )
        )

    def __mul__(self, other):
        return self.scalar_product(other)

    def __repr__(self):
        return "Vector{}".format(self._coordinates)

    def __sub__(self, other):
        if not isinstance(other, (Vector, Point)):
            raise NotImplementedError(
                "Can't subtract anything but vectors from other vectors"
            )
        return self + -1 * other

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
        return not Vector.linear_independent(self, other)

    def is_perpendicular_to(self, other):
        """ Returns whether `self` is perpendicular to `other` """
        if not isinstance(other, Vector):
            raise ValueError(
                "Can compute perpendicularity only between vectors"
            )
        return self.scalar_product(other) == 0

    @staticmethod
    def linear_independent(*vectors):
        """
        Checks whether list of vectors are linearly independent from each other
        """
        coefficients = [Symbol('x{}'.format(n)) for n in range(len(vectors))]
        vectors = [coefficients[n] * vectors[n] for n in range(len(vectors))]
        v = Vector.sum(*vectors)
        equations = [v.get_nth_coordinate(n) for n in range(v.dimension)]
        solutions = solve(equations, coefficients)
        return solutions and all(solutions.get(c) == 0 for c in coefficients)

    def scalar_product(self, other):
        # a * V = (aV1, aV2, ..., aVn)
        if isinstance(other, (int, float, Expr)):
            return Vector(
                *[other * v for v in self._coordinates],
                constrains=self._constrains
            )

        # V * U = V1*U1 + V2*U2 + ... + Vn*Un
        elif isinstance(other, Vector):
            if self.dimension != other.dimension:
                raise ValueError(
                    "Cannot compute scalar product for vectors of different"
                    " dimensions"
                )
            return sum(
                self._coordinates[n] * other._coordinates[n]
                for n in range(self.dimension)
            )

        raise NotImplementedError(
            "Cannot compute scalar product between {} and {}"
            .format(self, other)
        )

    @staticmethod
    def sum(*vectors):
        """ Returns a new vector as a result of the sum of given vectors """
        if not vectors:
            return Vector()
        max_dim = max(v.dimension for v in vectors)
        coordinates = [
            sum(v.get_nth_coordinate(n) for v in vectors)
            for n in range(max_dim)
        ]
        constrains = ()
        for v in vectors:
            constrains += v._constrains
        return Vector(*coordinates, constrains=tuple(set(constrains)))
