import numpy as np
import matplotlib.pyplot as plt



beta = 1.0
sigma = 1.0
gamma = 0.1

s_0 = 0.99
e_0 = 0.01
i_0 = 0
r_0 = 0


def sus(s, e, i):
    return -beta * i * s


def exp(s, e, i):
    return beta * i * s - sigma * e


def inf(e, i):
    return sigma * e - gamma * i


def rem(i):
    return gamma * i


def metoda_rungego_kutty(s, e, i, r, h):
    k1_s = sus(s, e, i)
    k1_e = exp(s, e, i)
    k1_i = inf(e, i)
    k1_r = rem(i)

    k2_s = sus(s + 0.5 * h * k1_s, e + 0.5 * h * k1_e, i + 0.5 * h * k1_i)
    k2_e = exp(s + 0.5 * h * k1_s, e + 0.5 * h * k1_e, i + 0.5 * h * k1_i)
    k2_i = inf(e + 0.5 * h * k1_e, i + 0.5 * h * k1_i)
    k2_r = rem(i + 0.5 * h * k1_i)

    k3_s = sus(s + 0.5 * h * k2_s, e + 0.5 * h * k2_e, i + 0.5 * h * k2_i)
    k3_e = exp(s + 0.5 * h * k2_s, e + 0.5 * h * k2_e, i + 0.5 * h * k2_i)
    k3_i = inf(e + 0.5 * h * k2_e, i + 0.5 * h * k2_i)
    k3_r = rem(i + 0.5 * h * k2_i)

    k4_s = sus(s + h * k3_s, e + h * k3_e, i + h * k3_i)
    k4_e = exp(s + h * k3_s, e + h * k3_e, i + h * k3_i)
    k4_i = inf(e + h * k3_e, i + h * k3_i)
    k4_r = rem(i + h * k3_i)

    s_new = s + (h / 6) * (k1_s + 2 * k2_s + 2 * k3_s + k4_s)
    e_new = e + (h / 6) * (k1_e + 2 * k2_e + 2 * k3_e + k4_e)
    i_new = i + (h / 6) * (k1_i + 2 * k2_i + 2 * k3_i + k4_i)
    r_new = r + (h / 6) * (k1_r + 2 * k2_r + 2 * k3_r + k4_r)
    

    return s_new, e_new, i_new, r_new


t_max = 50
h = 0.1
steps = int(t_max / h) + 1


time = np.linspace(0, t_max, steps)
susceptible = np.zeros(steps)
exposed = np.zeros(steps)
infected = np.zeros(steps)
removed = np.zeros(steps)


susceptible[0] = s_0
exposed[0] = e_0
infected[0] = i_0
removed[0] = r_0


for i in range(1, steps):
    susceptible[i], exposed[i], infected[i], removed[i] = metoda_rungego_kutty(susceptible[i - 1],
                                               exposed[i - 1], infected[i - 1], removed[i - 1], h)


plt.plot(time, susceptible, label="Susceptible")
plt.plot(time, exposed, label="Exposed")
plt.plot(time, infected, label="Infected")
plt.plot(time, removed, label="Removed")
plt.xlabel("Czas")
plt.ylabel("Populacja")
plt.legend()
plt.show()













