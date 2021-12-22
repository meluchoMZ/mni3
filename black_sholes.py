#!/usr/bin/python
# coding=utf-8
# Numerical methods for computer engineering
# Author: Miguel Blanco Godón, Computer Engineering, 2021

import thomas
from matplotlib import pyplot


class BlackSholes:

    def __init__(self, sigma: float, r: float, e: float, theta: float, tf: float, a: float, b: float):
        self.__sigma = sigma
        self.__r = r
        self.__e = e
        self.__theta = theta
        self.__tf = tf
        self.__a = a
        self.__b = b

    def get_sigma(self):
        return self.__sigma

    def set_sigma(self, sigma: float):
        self.__sigma = sigma

    def get_r(self):
        return self.__r

    def set_t(self, r: float):
        self.__r = r

    def get_e(self):
        return self.__e

    def set_e(self, e: float):
        self.__e = e

    def get_theta(self):
        return self.__theta

    def set_theta(self, theta: float):
        self.__theta = theta

    def compute(self):
        nn = 500
        delta_x = (4*self.__e) / nn
        tt = 120
        delta_t = self.__tf / tt
        x = []
        for n in range(nn):
            x.append(self.__a + n*delta_x)
        t = []
        for _t in range(tt):
            t.append(_t*delta_t)
        u = [max(self.__e - xi, 0) for xi in x]
        lt = len(t)
        lu = len(u) - 1
        # du/dx(x = 0, t) = -1
        left = [0]
        middle = [-1/delta_x]
        right = [1/delta_x]
        for i in range(1, lu):
            # u(i-1)
            left.append(-(self.__sigma * self.__sigma * x[i] * x[i] * self.__theta) / (2 * delta_x * delta_x))
            # u(i)
            middle.append((1 / delta_t) + (self.__sigma * self.__sigma * x[i] * x[i] * self.__theta) / (
                    delta_x * delta_x) + (self.__r*x[i]*self.__theta)/delta_x + self.__r*self.__theta)
            # u(i+1)
            right.append(-(self.__sigma * self.__sigma * x[i] * x[i] * self.__theta) / (2 * delta_x * delta_x) - (
                self.__r*x[i]*self.__theta)/delta_x)
        # u(x = 4E, t) = 0
        left.append(0)
        middle.append(1)
        right.append(0)
        middle, left, right = thomas.thomas01(middle, left, right)
        for t in range(lt):
            # δu/δx (x = 0, t) = -1
            fu = [-1+((1-self.__theta)*u[0])/delta_x - ((1-self.__theta)*u[1])/delta_x]
            for i in range(1, lu):
                fu.append\
                    (
                        # u[i-1]
                        (u[i-1]*((self.__sigma*self.__sigma*x[i]*x[i]*(1-self.__theta))/(2*delta_x*delta_x))) +
                        # u[i]
                        (u[i]*((1/delta_t)-((1-self.__theta)*self.__r)-(self.__sigma*self.__sigma*x[i]*x[i]*(
                                1-self.__theta))/(delta_x*delta_x)-((self.__r*x[i]*(1-self.__theta))/delta_x))) +
                        # u[i+1]
                        (u[i+1]*((self.__sigma*self.__sigma*x[i]*x[i]*(1-self.__theta))/(2*delta_x*delta_x)+((
                                self.__r*x[i]*(1-self.__theta))/delta_x)))
                    )
            # u(x = 4E, t) 0, 1, = 0
            fu.append(0)
            sol = thomas.thomas02(middle, left, right, fu)
            if t == 0:
                pyplot.plot(x, sol)
            for i in range(lu):
                u[i] = sol[i+1]
        pyplot.plot(x, sol)
        pyplot.show()
