#!/usr/bin/python
# coding=utf-8
# Numerical methods for computer engineering
# Author: Miguel Blanco Godón, Computer Engineering, 2021

from matplotlib import pyplot


class Oscillator:

    def __init__(self, w: float, psi: float, y0: float, v0: float, t: float, n: int):
        self.__w = w
        self.__psi = psi
        self.__y0 = y0
        self.__v0 = v0
        if t <= 0:
            raise ValueError("Simulation time must be positive number!")
        self.__t = t
        if n <= 0:
            raise ValueError("Step number must be a positive integer!")
        self.__n = n

    def get_w(self):
        return self.__w

    def set_w(self, w: float):
        self.__w = w

    def get_psi(self):
        return self.__psi

    def set_psi(self, psi: float):
        self.__psi = psi

    def get_y0(self):
        return self.__y0

    def set_y0(self, y0: float):
        self.__y0 = y0

    def get_v0(self):
        return self.__v0

    def set_v0(self, v0: float):
        self.__v0 = v0

    def get_t(self):
        return self.__t

    def set_t(self, t: float):
        if t <= 0:
            raise ValueError("Simulation time must be a positive number")
        self.__t = t

    def oscillate_explicit(self):
        delta = self.__t /self.__n
        t = []
        for n in range(self.__n):
            t.append(n * delta)
        y = [self.__y0]
        v = [self.__v0]
        for n in range(1, self.__n):
            v.append(v[n-1]+delta*(-(self.__w*self.__w)*y[n-1] - 2*self.__psi*v[n-1]))
            y.append(y[n-1]+delta*v[n-1])

        pyplot.subplot(221)
        pyplot.plot(t, y, 'r-')
        pyplot.axhline(y=0, color="black")
        pyplot.axvline(x=0, color="black")
        pyplot.title("Position with PSI="+str(self.__psi)+" (Euler explicit)")
        pyplot.xlabel("t")
        pyplot.ylabel("y")
        pyplot.legend(["y(t)"], loc="lower right")
        pyplot.grid(True, which="both")

        pyplot.subplot(222)
        pyplot.plot(v, y)
        pyplot.axhline(y=0, color="black")
        pyplot.axvline(x=0, color="black")
        pyplot.xlabel("dy(t)/dt")
        pyplot.ylabel("y(t)")
        pyplot.title("Phase diagram with PSI="+str(self.__psi)+" (Euler explicit)")
        pyplot.grid(True, which="both")
        pyplot.get_current_fig_manager().full_screen_toggle()

        pyplot.subplot(223)
        pyplot.plot(t, v, 'b-')
        pyplot.axhline(y=0, color="black")
        pyplot.axvline(x=0, color="black")
        pyplot.xlabel("t")
        pyplot.ylabel("dy(t)/dt")
        pyplot.title("Velocity with PSI="+str(self.__psi)+" (Euler explicit)")
        pyplot.legend(["dy/dt"], loc="lower right")
        pyplot.grid(True, which="both")
        pyplot.show()

    def oscillate_implicit(self):
        delta = self.__t / self.__n
        t = []
        for n in range(self.__n):
            t.append(n * delta)
        y = [self.__y0]
        v = [self.__v0]
        for n in range(1, self.__n):
            v.append((v[n-1] - self.__w*self.__w*y[n-1]*delta)/(1+2*self.__psi*delta))
            y.append(y[n-1]+delta*v[n])

        pyplot.subplot(221)
        pyplot.plot(t, y, 'r-')
        pyplot.axhline(y=0, color="black")
        pyplot.axvline(x=0, color="black")
        pyplot.title("Position with PSI="+str(self.__psi)+" (Euler implicit)")
        pyplot.xlabel("t")
        pyplot.ylabel("y")
        pyplot.legend(["y(t)"], loc="lower right")
        pyplot.grid(True, which="both")

        pyplot.subplot(222)
        pyplot.plot(v, y)
        pyplot.axhline(y=0, color="black")
        pyplot.axvline(x=0, color="black")
        pyplot.xlabel("dy(t)/dt")
        pyplot.ylabel("y(t)")
        pyplot.title("Phase diagram with PSI="+str(self.__psi)+" (Euler implicit)")
        pyplot.grid(True, which="both")
        pyplot.get_current_fig_manager().full_screen_toggle()

        pyplot.subplot(223)
        pyplot.plot(t, v, 'b-')
        pyplot.axhline(y=0, color="black")
        pyplot.axvline(x=0, color="black")
        pyplot.xlabel("t")
        pyplot.ylabel("dy(t)/dt")
        pyplot.title("Velocity with PSI="+str(self.__psi)+" (Euler implicit)")
        pyplot.legend(["dy/dt"], loc="lower right")
        pyplot.grid(True, which="both")
        pyplot.show()
