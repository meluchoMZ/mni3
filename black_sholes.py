#!/usr/bin/python
# coding=utf-8
# Numerical methods for computer engineering
# Author: Miguel Blanco Godón, Computer Engineering, 2021

import thomas
from matplotlib import pyplot


class BlackSholes:

    def __init__(self, sigma: float, r: float, e: float, theta: float, tf: float, x_steps: int, t_steps: int):
        self.__sigma = sigma
        self.__r = r
        self.__e = e
        self.__theta = theta
        self.__tf = tf
        if x_steps <= 0 or t_steps <= 0 :
            raise ValueError("Step number must be a positive integer")
        self.__x_steps = x_steps
        self.__t_steps = t_steps

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

    def get_x_steps(self):
        return self.__x_steps

    def set_x_steps(self, x_steps):
        self.__x_steps = x_steps

    def get_t_steps(self):
        return self.__t_steps

    def set_t_steps(self, t_steps):
        self.__t_steps = t_steps

    def compute(self):
        delta_x = (4*self.__e) / self.__x_steps
        delta_t = self.__tf / self.__t_steps
        x = []
        for n in range(self.__x_steps):
            x.append(n*delta_x)
        t = []
        for _t in range(self.__t_steps):
            t.append(_t*delta_t)
        u = [max(self.__e - xi, 0) for xi in x]
        lt = len(t)
        lu = len(u) - 1
        # du/dx(x = 0, t) = -1
        left = [0]
        middle = [-(self.__theta/delta_x)]
        right = [self.__theta/delta_x]
        for i in range(1, lu):
            # u(i-1)
            left.append(-((self.__sigma * self.__sigma * x[i] * x[i] * self.__theta) / (2 * delta_x * delta_x)))
            # u(i)
            middle.append((1 / delta_t) + (self.__sigma * self.__sigma * x[i] * x[i] * self.__theta) / (
                    delta_x * delta_x) + (self.__r*x[i]*self.__theta)/delta_x + self.__r*self.__theta)
            # u(i+1)
            right.append(-((self.__sigma * self.__sigma * x[i] * x[i] * self.__theta) / (2 * delta_x * delta_x)) - ((
                self.__r*x[i]*self.__theta)/delta_x))
        # u(x = 4E, t) = 0
        left.append(0)
        middle.append(1)
        right.append(0)
        middle, left, right = thomas.thomas01(middle, left, right)
        sol = []
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
            for i in range(lu):
                u[i] = sol[i+1]
        pyplot.plot(x, sol)
        pyplot.title("Theta = " + str(self.__theta) + "; Sigma = " + str(self.__sigma) +  "; r = " + str(self.__r))
        pyplot.xlabel("t")
        pyplot.ylabel("Price")
        pyplot.axhline(y=0, color="black")
        pyplot.axvline(x=0, color="black")
        pyplot.grid(True, which="both")
        pyplot.get_current_fig_manager().full_screen_toggle()
        pyplot.show()
