#!/usr/bin/python
# coding=utf-8
# Numerical methods for computer engineering
# Author: Miguel Blanco Godón, Computer Engineering, 2021


class BlackSholes:

    def __init__(self, sigma: float, r: float, e: float, theta: float):
        self.__sigma = sigma
        self.__r = r
        self.__e = e
        self.__theta = theta

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
        print(self.__r)
