#!/usr/bin/env python3

# Library for polynomial based functions


from collections import defaultdict
import math
import tqdm


# Helper for ** method
def phelp(poly, n):
    if n == 0:
        return Polynomial(nvars=poly.nvars, coeffs={tuple([0] * poly.nvars): 1})
    if n == 1:
        return poly
    if n % 2 == 0:
        return phelp(poly, n // 2) * phelp(poly, n // 2)
    return poly * phelp(poly, n // 2) * phelp(poly, n // 2)


LEX_CACHE = dict()


# lex_index python rewrite
def lex(arr, n, d):
    if arr in LEX_CACHE:
        return LEX_CACHE[arr]
    op = 0
    s = 0
    for j in range(1, n):
        if j > 1:
            s += int(arr[j - 2])
        for i in range(1, int(arr[j - 1]) + 1):
            op += math.factorial(n - j + d + 1 - s - i - 1) / (
                math.factorial(d + 1 - s - i) * math.factorial(n - j - 1)
            )
    op += 1
    LEX_CACHE[arr] = op
    return op


class Polynomial:
    def __init__(self, nvars, coeffs=None) -> None:
        self.nvars = nvars  # eg. num variables: 3
        self.coeffs = defaultdict(int)
        if coeffs is not None:
            for k, v in coeffs.items():  # eg. dict: {(0,0,0):4,(0,1,8):2}
                self.coeffs[k] = v

    def get_dim(self):
        """Returns the degree of the polynomial"""
        if len(self.coeffs.keys()) == 0:
            self.dim = 0
        else:
            self.dim = max(sum(k) for k in self.coeffs.keys())
        return self.dim

    def get_active_vars(self):
        """Returns list of variables with non-zero coefficients"""
        s = set()
        for i in self.coeffs.keys():
            if self.coeffs[i] == 0:
                continue
            for j, v in enumerate(i):
                if v == 0:
                    continue
                s.add(j)
        return s

    def numberadd(self, num):
        """Add a number to the polynomial and return a new polynomial"""
        a = self.coeffs
        c = Polynomial(nvars=self.nvars, coeffs=a)
        zero = tuple([0 for _ in range(self.nvars)])
        c.coeffs[zero] += num
        return c

    def inumberadd(self, num):
        """Add a number to a polynomial in place"""
        a = self.coeffs
        zero = tuple([0 for _ in range(self.nvars)])
        self.coeff[zero] += num
        return self

    def __add__(self, other):
        """Add two polynomials and return a new one"""
        if not isinstance(other, Polynomial):
            return self.numberadd(other)
        a = self.coeffs
        b = other.coeffs
        c = Polynomial(nvars=self.nvars)
        for k in a.keys():
            c.coeffs[k] += a[k]
        for k in b.keys():
            c.coeffs[k] += b[k]
        return c

    def __iadd__(self, other):
        """Add a polynomial in place"""
        if not isinstance(other, Polynomial):
            return self.inumberadd(other)
        for k in other.coeffs:
            self.coeffs[k] += other.coeffs[k]
        return self

    def __sub__(self, other):
        """Subtract two polynomials"""
        if not isinstance(other, Polynomial):
            return self.numberadd(-other)
        a = self.coeffs
        b = other.coeffs
        c = Polynomial(nvars=self.nvars)
        for k in a.keys():
            c.coeffs[k] += a[k]
        for k in b.keys():
            c.coeffs[k] -= b[k]
        return c

    def numbermul(self, num):
        """Multiply polynomial by a number"""
        a = self.coeffs
        c = Polynomial(nvars=self.nvars)
        if num == 0:
            return c
        for k in a:
            c.coeffs[k] = num * a[k]
        return c

    def __mul__(self, other):
        """Multiply two polynomials"""
        if not isinstance(other, Polynomial):
            return self.numbermul(other)
        a = self.coeffs
        b = other.coeffs
        c = Polynomial(nvars=self.nvars)
        for ai in a:
            for bj in b:
                tupsum = [0] * c.nvars
                for i in range(c.nvars):
                    tupsum[i] = ai[i] + bj[i]
                c.coeffs[tuple(tupsum)] += a[ai] * b[bj]
        return c

    def __repr__(self) -> str:
        def sign(x):
            if x >= 0:
                return "+"
            return ""

        s = ""
        for k, coeff in self.coeffs.items():
            if coeff == 0:
                continue
            a = ""
            if sum(k) != 0:
                for i, ex in enumerate(k):
                    if ex == 0:
                        continue
                    a += f"x_{{{i+1}}}^{ex}"
            s += sign(coeff) + str(coeff) + a
        return s if s != "" else "0"

    def __pow__(self, n):
        return phelp(self, n)

    def to_file(self, filename, verbose=False, i=None, j=None):
        """Write full polynomial to a file"""
        with open(filename, "w") as f:
            runner = self.coeffs.keys()
            d = 4  # self.get_dim()
            if verbose:
                i = 0
                s = ""
                for k in runner:
                    if self.coeffs[k] != 0:
                        s += f"{lex(k,int(self.nvars),d)},{self.coeffs[k]},{k}" + "\n"
                    i += 1
                    if i == 1_000:  # Chunked data write
                        f.write(s)
                        s = ""
                        i = 0
                f.write(s)
            else:
                inds = [
                    (a, b, c)
                    for c in range(d + 1)
                    for b in range(d + 1)
                    for a in range(d + 1)
                    if a + b + c <= d
                ]
                i = 0
                s = ""
                for i in inds:
                    g = False
                    for k in runner:
                        if k == i:
                            s += f"{self.coeffs[k]}" + "\n"
                            g = True
                            break
                    if not g:
                        s += f"0.0\n"
                f.write(s)


def test():
    p = Polynomial(
        nvars=2,
        coeffs={
            (0, 0): 3,
            (1, 0): -3,
            (2, 0): 5,
        },
    )
    q = Polynomial(
        nvars=2,
        coeffs={
            (0, 0): -3,
            (1, 0): 3,
            (2, 2): -5,
        },
    )
    print(p, q, p + q, p * q, p**2, p - q, p - p, q - q, sep="\n")
    p += p
    print(p)
    # +3-3x_{1}^1+5x_{1}^2
    # -3+3x_{1}^1-5x_{1}^2x_{2}^2
    # +5x_{1}^2-5x_{1}^2x_{2}^2
    # -9+18x_{1}^1-15x_{1}^2x_{2}^2-24x_{1}^2+15x_{1}^3x_{2}^2+15x_{1}^3-25x_{1}^4x_{2}^2
    # +9-18x_{1}^1+39x_{1}^2-30x_{1}^3+25x_{1}^4
    # +6-6x_{1}^1+10x_{1}^2


if __name__ == "__main__":
    test()
    # main()
