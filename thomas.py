# encoding: utf-8

import numpy as np

# ----------------------------------------------------------------------------------------

def thomas01 (a,b,c):

    '''
    THOMAS01: Factorización LU de una matriz tridiagonal.
              La matriz está almacenada en tres vectores.

    Entradas:
    - a: diagonal principal de la matriz.
    - b: diagonal inferior.
    - c: diagonal superior.

    Salidas:
    - alpha: diagonal principal de U.
    - beta:  diaongal inferior de L.
    - gamma: diagonal superior de U (no se calcula, es igual a C).
    '''

    n = len (a)
    alpha = np.zeros (n, 'f')
    beta = np.zeros (n, 'f')

    alpha [0] = a [0]
    for i in range (1, n-1):
        beta [i] = b [i] / alpha [i-1]
        alpha [i] = a [i] - beta [i] * c [i-1]
    beta [n-1] = b [n-1] / alpha [n-2]
    alpha [n-1] = a [n-1] - beta [n-1] * c [n-2]

    return [alpha, beta, c]

# ----------------------------------------------------------------------------------------

def thomas02 (alpha,beta,gamma,f):

    '''
    THOMAS02: Resolución (descenso y remonte) de un sistema con matriz tridiagonal.
              La matriz ha sido previamente factorizada en forma LU.
              La matriz factorizada está almacenada en tres vectores.

    Entradas:
    - alpha: diagonal principal de U.
    - beta:  diaongal inferior de L.
    - gamma: diagonal superior de U.
    - f:     segundo miembro.

    Salida:
    - x: solución del sistema.
    '''

    n = len (alpha)
    x = np.zeros (n)
    y = np.zeros (n)

    # Descenso.

    y [0] = f [0]
    for i in range (1,n):
        y [i] = f [i] - beta [i] * y [i-1]

    # Remonte.

    x [n-1] = y [n-1] / alpha [n-1]
    for i in range (n-2,-1,-1):
        x [i] = (y [i] - gamma [i] * x [i+1]) / alpha [i]

    return x

# ----------------------------------------------------------------------------------------

if (__name__ == '__main__'):

    '''
    # Datos, construyendo los vectores.

    m = np.dot (np.array ([[1, 0, 0, 0, 0], [-1, 1, 0, 0, 0], [0, -2, 1, 0, 0], [0, 0, 1, 1, 0], [0, 0, 0, 1, 1]]) , \
                np.array ([[4, 1, 0, 0, 0], [0, 6, -1, 0, 0], [0, 0, 4, 2, 0], [0, 0, 0, 5, 1], [0, 0, 0, 0, 4]]))
    s = np.array ([1, 2, 3, 4, 5])

    n = len (s)
    a = np.zeros (n, 'f')
    b = np.zeros (n, 'f')
    c = np.zeros (n, 'f')
    a [0] = m [0,0]
    c [0] = m [0,1]
    for i in range (1,n-1):
        a [i] = m [i,i]                    # Diagonal principal.
        b [i] = m [i,i-1]                  # Diagonal inferior.
        c [i] = m [i,i+1]                  # Diagonal superior.
    a [n-1] = m [n-1,n-1]
    b [n-1] = m [n-1,n-2]

    d = np.dot (m,s)

    # Datos.

    a = np.array ([4., 5., 6., 7., 5.])             # Diagonal principal.
    b = np.array ([0., -4., -12., 4., 5.])          # Diagonal inferior.
    c = np.array ([1., -1., 2., 1., 0.])            # Diagonal superior.
    d = np.array ([6., 3., 2., 45., 45.])

    #####################################################
    n = len (d)
    m = np.zeros ((n,n), 'f')
    for i in range (n):
        m [i,i] = a [i]
    for i in range (n-1):
        m [i, i+1] = c [i]
    for i in range (1,n):
        m [i, i-1] = b [i]

    print ' Matriz: '
    print m
    x = np.array ([1., 2., 3., 4., 5.])
    print ' M x:'
    print np.dot (m,x)
    #####################################################

    # Resolución.

    [alpha, beta, gamma] = thomas01 (a,b,c)
    x = thomas02 (alpha, beta, gamma, d)

    print
    print '     Solución: ', x
    print

    '''
# ----------------------------------------------------------------------------------------

