# coding=utf-8
# Numerical methods for computer engineering
# Author: Miguel Blanco Godón, Computer Engineering, 2021

import matplotlib.pyplot as pyplot
import numpy


class ConvDiff:

    def __init__(self, steps: int):
        self.__steps = steps

    def compute(self):
        dx = (1 - (-1)) / self.__steps
        dy = (1 - (-1)) / self.__steps
        x = []
        for i in range(self.__steps + 1):
            x.append(-1 + i * dx)
        y = []
        for j in range(self.__steps + 1):
            y.append(-1 + j * dy)
        # domain
        u = []
        for i in range(self.__steps + 1):
            u.append([])
            for j in range(self.__steps + 1):
                u[i].append(0)
        for i in range(1, self.__steps):
            for j in range(1, self.__steps):
                f = (x[i] + 1)*(x[i]+1) + y[j]*y[j]
                if f < (1 / 20):
                    u[j][i] = 60
        # speed vectorial field
        b = [[], []]
        for i in range(self.__steps + 1):
            b[0].append(-3 * y[i])
            b[1].append(-3 * x[i])
        a = []
        for i in range(self.__steps + 1):
            a.append([])
            for j in range(self.__steps + 1):
                # diffusion coefficient
                a[i].append(0.06 - 0.05 * x[i] * x[i] - 0.002 * y[j]*y[j])

        c = []
        for i in range(self.__steps + 1):
            c.append([])
            for j in range(self.__steps + 1):
                # reaction coefficient
                c[i].append(10 * x[i] * x[i] * y[j] * y[j])

        pyplot.subplot(321)
        pyplot.imshow(u, extent=[-1, 1, -1, 1])
        pyplot.title("Initial state Pixel")
        pyplot.axhline(y=0, color="black")
        pyplot.axvline(x=0, color="black")
        pyplot.grid(True, which="both")
        pyplot.colorbar()

        pyplot.subplot(322)
        pyplot.contour(x, y, u)
        pyplot.title("Initial state Contour")
        pyplot.axhline(y=0, color="black")
        pyplot.axvline(x=0, color="black")
        pyplot.grid(True, which="both")
        pyplot.colorbar()

        p = numpy.zeros((self.__steps**2, self.__steps**2))
        res = numpy.zeros(self.__steps**2)

        # inner nodes
        for j in range(self.__steps+1):
            for i in range(1, self.__steps-1):
                k = i*self.__steps+j
                p[k][k] = ((2*a[i][j])/(dx*dx)+(2*a[i][j])/(dy*dy)+c[i][j])
                p[k][k-1] = ((0.1*x[i])/(2*dx)+(3*y[j])/(2*dx)-a[i][j]/(dx*dx))
                p[k][k+1] = ((-0.1*x[i])/(2*dx)+(-3*y[j])/(2*dx)-a[i][j]/(dx*dx))
                p[k][k-(self.__steps-2)] = ((0.004*y[j]-3*x[i])/(2*dy)-a[i][j]/(dy*dy))
                p[k][k+(self.__steps-2)] = ((-0.004*y[j]+3*x[i])/(2*dy)-a[i][j]/(dy*dy))
                # second member of system
                res[k] = u[i][j]
        # side conditions
        for i in range(self.__steps):
            p[i][i] = 1
            p[-i][-i] = 1
            res[i] = 0
            res[-i] = 0
            p[i][0] = 1
            p[i][-1] = 1

        new_p = numpy.linalg.solve(p, res)

        new_u = []
        for i in range(self.__steps+1):
            new_u.append([])
            for j in range(self.__steps+1):
                new_u[i].append(0)
        for j in range(1, self.__steps-1):
            for i in range(1, self.__steps):
                k = i*self.__steps+j
                new_u[i][j] = new_p[k]
        # plotting
        pyplot.subplot(323)
        pyplot.imshow(new_u, extent=[-1, 1, -1, 1])
        pyplot.title("Solution (a) Pixel")
        pyplot.axhline(y=0, color="black")
        pyplot.axvline(x=0, color="black")
        pyplot.grid(True, which="both")
        pyplot.get_current_fig_manager().full_screen_toggle()
        pyplot.colorbar()

        pyplot.subplot(324)
        pyplot.contour(x, y, new_u)
        pyplot.title("Solution (a) Contour")
        pyplot.axhline(y=0, color="black")
        pyplot.axvline(x=0, color="black")
        pyplot.grid(True, which="both")
        pyplot.get_current_fig_manager().full_screen_toggle()
        pyplot.colorbar()

        p_neumann = numpy.zeros((self.__steps**2, self.__steps**2))
        res_neumann = numpy.zeros(self.__steps**2)

        # inner nodes
        for j in range(self.__steps+1):
            for i in range(1, self.__steps-1):
                k = i*self.__steps+j
                p_neumann[k][k] = ((2*a[i][j])/(dx*dx)+(2*a[i][j])/(dy*dy)+c[i][j])
                p_neumann[k][k-1] = ((0.1*x[i])/(2*dx)+(3*y[j])/(2*dx)-a[i][j]/(dx*dx))
                p_neumann[k][k+1] = ((-0.1*x[i])/(2*dx)+(-3*y[j])/(2*dx)-a[i][j]/(dx*dx))
                p_neumann[k][k-(self.__steps-2)] = ((0.004*y[j]-3*x[i])/(2*dy)-a[i][j]/(dy*dy))
                p_neumann[k][k+(self.__steps-2)] = ((-0.004*y[j]+3*x[i])/(2*dy)-a[i][j]/(dy*dy))
                # second member of system
                res_neumann[k] = u[i][j]
        # side conditions
        for i in range(self.__steps):
            p_neumann[i][i] = 1
            p_neumann[-i][-i] = 1
            res_neumann[i] = 0
            res_neumann[-i] = 0
            p_neumann[i][0] = 1
            p_neumann[i][-1] = 1

        # Neumann conditions
        for i in range(int(self.__steps/2)):
            p_neumann[-i][self.__steps] = 1/dy
            p_neumann[-i][-i] = -1/dy
            res_neumann[-i] = 0.08

        neumann_p = numpy.linalg.solve(p_neumann, res_neumann)

        neumann_u = []
        for i in range(self.__steps+1):
            neumann_u.append([])
            for j in range(self.__steps+1):
                neumann_u[i].append(0)
        for j in range(1, self.__steps-1):
            for i in range(1, self.__steps):
                k = i*self.__steps+j
                neumann_u[i][j] = neumann_p[k]
        # plotting
        pyplot.subplot(325)
        pyplot.imshow(neumann_u, extent=[-1, 1, -1, 1])
        pyplot.title("Solution (b) Pixel")
        pyplot.axhline(y=0, color="black")
        pyplot.axvline(x=0, color="black")
        pyplot.grid(True, which="both")
        pyplot.get_current_fig_manager().full_screen_toggle()
        pyplot.colorbar()

        pyplot.subplot(326)
        pyplot.contour(x, y, neumann_u)
        pyplot.title("Solution (b) Contour")
        pyplot.axhline(y=0, color="black")
        pyplot.axvline(x=0, color="black")
        pyplot.grid(True, which="both")
        pyplot.get_current_fig_manager().full_screen_toggle()
        pyplot.colorbar()

        pyplot.show()
