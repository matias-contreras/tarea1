
# Tarea 1 fisica de la tierra
# encontrar la masa de la tierra y la gravedad que produce
# dados los intervalos de densidades.

import numpy as np
from matplotlib import pyplot as plt
from scipy import integrate


# constantes
Rt = 6371  # en kilometros
Mt = 5.973 * 10**24  # en kilogramos
G = 6.6732 * 10**-11 * 10**-6  # en N * Km**2 / kg**2

# profundidades y densidades
rydens = [[0, 2.84],
          [15, 2.84],
          [15, 3.313],
          [60, 3.332],
          [100, 3.348],
          [200, 3.387],
          [300, 3.424],
          [350, 3.441],
          [350, 3.7],
          [400, 3.775],
          [413, 3.795],
          [500, 3.925],
          [600, 4.075],
          [650, 4.15],
          [650, 4.2],
          [800, 4.38],
          [984, 4.529],
          [1000, 4.538],
          [1200, 4.655],
          [1400, 4.768],
          [1600, 4.877],
          [1800, 4.983],
          [2000, 5.087],
          [2200, 5.188],
          [2400, 5.288],
          [2600, 5.387],
          [2800, 5.487],
          [2878, 5.527],
          [2878, 9.927],
          [3000, 10.121],
          [3200, 10.421],
          [3400, 10.697],
          [3600, 10.948],
          [3800, 11.176],
          [4000, 11.383],
          [4200, 11.57],
          [4400, 11.737],
          [4600, 11.887],
          [4800, 12.017],
          [4982, 12.121],
          [5000, 12.13],
          [5121, 12.197],
          [5200, 12.229],
          [5400, 12.301],
          [5600, 12.36],
          [5800, 12.405],
          [6000, 12.437],
          [6200, 12.455],
          [6371, 12.46]]

rydens = np.array(rydens)

# se cambia de profundidad desde la
# superficie a radio desde el centro terrestre

radios = (6371 - rydens[0:, 0])
densidades = rydens[0:, 1]

rydens[0:, 0] = radios  # en kilometros
rydens[0:, 1] = densidades * 10**12  # en kg/km**3

rydens = rydens[::-1]  # para dejar el cero al inicio del vector


def coef_recta(r1, r2):
    """calcula [m, n] de una recta  y = mx + n, dados 2 puntos
     [r1, r2] = [(x1, y1), (x2, y2)] """
    m = (r2[1] - r1[1])/(r2[0] - r1[0])
    n = r1[1] - m*r1[0]
    return [m, n]


def func_recta(x, m, n):
    """dados los coeficientes y un x, encuentra el valor de la funcion"""
    y = m*x + n
    return y


def densidad(radio):
    """encuentra la densidad de la tierra en cierto radio"""
    if radio == 0:
        densidad = rydens[0, 1]
    else:
        for i in range(len(rydens)):
            if rydens[i][0] < radio <= rydens[i + 1][0]:
                coefs = coef_recta(rydens[i], rydens[i + 1])
                densidad = func_recta(radio,coefs[0], coefs[1])

                return densidad
    return densidad


def func_a_integrar_masa(radio):
    """funcion a integrar para encontrar la masa"""
    func = 4 * np.pi * densidad(radio) * radio**2
    return func


def masa_total(radio):
    """masa en kilogramos resultante al integrar, usando el metodo scipy.quad"""

    if radio == 0:
        integral = 0
    else:
        integral, err = integrate.quad(func_a_integrar_masa, 0, radio)
    return np.float(integral)


def gravedad(radio):
    """encuentra la gravedad usando la ley de gauss"""
    if radio == 0:
        g = 0
    else:
        g = - G * masa_total(radio) / radio**2
    return g


def masa_cascaron(radio, dens, delta):
    """calcula la masa de una cascara de esfera dadas la densidad, el radio menor y mayor"""
    masa_cascaron = dens * 4/3 * np.pi * ((radio)**3 - (radio - delta)**3)
    return masa_cascaron


def masa_total_cascaron(radio, delta=1):
    """calcula la masa en un radio, usando una discretización de tamaño delta, 1 km por defecto"""
    x = np.arange(0, radio + delta, delta)
    masas = []
    for r in x:
        dens = densidad(r)

        masa = masa_cascaron(r, dens, delta)

        masas.append(masa)

    return np.sum(masas)


radios = rydens[0:, 0]
densidades_plot = []
masas_plot = []
gravedades_plot = []

for r in radios:
    dens = densidad(r)
    masa = masa_total(r)
    grav = gravedad(r)
    densidades_plot.append(dens)
    masas_plot.append(masa)
    gravedades_plot.append(grav)


fig1 = plt.figure(1)
plt.plot(radios, densidades_plot)
plt.title('Gráfico densidad v/s radio')
plt.ylabel('Densidad $\left[\dfrac{Kg}{Km^{3}}\\right]$')
plt.xlabel('Radio [Km]')
plt.savefig('grafico_densidad')
fig1.show()

fig2 = plt.figure(2)
plt.plot(radios, masas_plot)
plt.title('Gráfico masa v/s radio')
plt.ylabel('Masa [Kg]')
plt.xlabel('Radio [Km]')
plt.savefig('grafico_masa')
fig2.show()

fig3 = plt.figure(3)
plt.plot(radios, np.abs(gravedades_plot))
plt.title('Gráfico modulo de la Ac. de gravedad v/s radio')
plt.ylabel('Módulo de la gravedad $\left[\dfrac{m}{s^{2}}\\right]$')
plt.xlabel('Radio [Km]')
plt.savefig('grafico_gravedad')
fig3.show()
