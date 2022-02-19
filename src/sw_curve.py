import sympy as sym
import numpy as np


# We start with just one example for a Seiberg-Witten curve/Mirror curve for (C^{*})^{3}
# More examples may be added later.
# TODO: Move this to an extra file (?)
# FAQ: Do we need to consider this multiplied by 1/x ?
def H_c3(x, y):
    return -x + y**2 + y


# In version one there is a parameter "comp" that tells
# the Program 
class sw_curve():
    # The SW curve (in our case) only depends on the defining function H
    def __init__(self, H):
        # Initialize formal variables
        self.x_sym, self.y_sym  = sym.symbols('x y')
        # Symbolic defining equation
        self.H_sym = sym.simplify(H(self.x_sym, self.y_sym))
        # Decendents of H_sym
        self.dHx_sym = sym.diff(self.H_sym, self.x_sym)
        self.dHy_sym = sym.diff(self.H_sym, self.y_sym)
        self.d2Hy2_sym = sym.diff(self.dHy_sym, self.y_sym)

        self.dHx = sym.lambdify([self.x_sym, self.y_sym],
                                self.dHx_sym, "numpy")
        self.dHy = sym.lambdify([self.x_sym, self.y_sym],
                                self.dHy_sym, "numpy")
        self.d2Hy2 = sym.lambdify([self.x_sym, self.y_sym],
                                  self.d2Hy2_sym, "numpy")
        self.n_y, self.d_y = sym.fraction(sym.together(self.H_sym), self.y_sym)
        # Lambdify symbolic functions for later use
        self.H = sym.lambdify([self.x_sym, self.y_sym],
                              self.H_sym, "numpy")
        # Call method to find branch points and singularities
        self.branch_points, self.sing_points = self.__branch_singular_points()

    # FAQ: Is it smart to use one method for
    # spectral and exponential differential?
    def sw_differential(self, pt, theta, expo=False):
        x = pt[0]
        y1 = pt[1]
        y2 = pt[2]
        if not(expo):
            dx = x * np.exp(1j * theta) / (y2 - y1)
            dy1 = -self.dHx(x, y1) / self.dHy(x, y1) * dx
            dy2 = -self.dHx(x, y2) / self.dHy(x, y2) * dx
        else:
            if self.dHy(x, np.exp(y2)) == 0 or self.dHy(x, np.exp(y1)) == 0 or np.isnan(self.dHy(x, np.exp(y1))) or np.isnan(self.dHy(x, np.exp(y2))) or np.exp(y1) == 0 or np.exp(y2) == 0 or (y1 - y2) == 0:
                return np.array([0,0,0])
            dx = x * np.exp(1j * theta) / (y2 - y1)
            dy1 = -self.dHx(x, np.exp(y1)) / self.dHy(x, np.exp(y1)) * dx / np.exp(y1)
            dy2 = -self.dHx(x, np.exp(y2)) / self.dHy(x, np.exp(y2)) * dx / np.exp(y2)
            if np.isinf(dy2):
                print(-self.dHx(x, np.exp(y2)), self.dHy(x, np.exp(y2)), dx, np.exp(y2))
        return np.array([dx, dy1, dy2])
    

     # Find branch points and singularities (Private)
    def __branch_singular_points(self):
        disc_poly_sym = sym.Poly(
            self.H_sym*self.d_y, self.y_sym).discriminant()
        n, d = sym.fraction(sym.together(disc_poly_sym))
        sing_sym = sym.solve(sym.Poly(d, self.x_sym))
        sing = np.array([sing_sym[i].evalf() for i in range(len(sing_sym))])
        branch_coeffs = sym.Poly(
            sym.simplify(d*disc_poly_sym), self.x_sym).all_coeffs()
        branch = np.roots(branch_coeffs)
        return branch, sing

    def get_fiber(self, x):
        rts_dict = sym.roots(
                (self.d_y*self.H_sym).subs(
                    self.x_sym, x))
        rts = []
        for rt in rts_dict:
                for i in range(rts_dict[rt]):
                    rts.append(complex(rt))

        return rts

    def match_fiber_point(self, x, Y):
        rts = self.get_fiber(x)
        fiber_Y = []
        for y in Y:
            j = np.argmin(np.array([abs(rts[j] - y) for j in range(len(rts))]))
            fiber_Y.append(complex(rts[j]))
        return fiber_Y

