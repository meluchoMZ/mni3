# coding=utf-8
# Numerical methods for computer engineering
# Author: Miguel Blanco Godón, Computer Engineering, 2021

from colours import Colours
from oscillator import Oscillator
from black_sholes import BlackSholes


def handle():
    valid_inputs = ['0', '1', '2', '3']
    while True:
        print("    (1): Damped oscillator")
        print("    (2): Sell option (Black-Sholes model)")
        print("    (3): Convection-Diffusion equation")
        print("    (0):" + Colours.red(" Quit"))
        x = input("Please select an operation: ")
        while x not in valid_inputs:
            x = input(Colours.red("Error. Please select a valid operation: "))
        if x == '0':
            print(Colours.green("Bye!"))
            break
        if x == '1':
            m = 2
            w = 1
            psi = [0, 0.1, 0.9, 1.0, 1.2, 2.0]
            y0 = -1
            v0 = 0
            t = float(input("Please insert simulation time: "))
            for i in psi:
                o = Oscillator(m, w, i, y0, v0, t)
                o.oscillate_explicit()
                o.oscillate_implicit()
        elif x == '2':
            tf = float(input("Please insert simulation time: "))
            sigma = 0.85
            r = 0.02
            e = 10
            for theta in [1, 0.5, 0]:
                BlackSholes(sigma, r, e, theta, tf, 0, 0).compute()
        else:
            print("Third problem:")
