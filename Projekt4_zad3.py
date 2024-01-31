import numpy as np
import matplotlib.pyplot as plt


s_0 = 0.99
e_0 = 0.01
i_0 = 0
r_0 = 0
sigma = 1.0
gamma = 0.1


def susceptible(s, e, i, beta):
    return -beta * i * s

def exposed(s, e, i, beta, sigma):
    return beta * i * s - sigma * e

def infected(e, i, sigma, gamma):
    return sigma * e - gamma * i

def removed(i, gamma):
    return gamma * i

def metoda_rungego_kutty(s, e, i, r, beta, sigma, gamma, h):
    k1_s = susceptible(s, e, i, beta)
    k1_e = exposed(s, e, i, beta, sigma)
    k1_i = infected(e, i, sigma, gamma)
    k1_r = removed(i, gamma)

    k2_s = susceptible(s + 0.5 * h * k1_s, e + 0.5 * h * k1_e, i + 0.5 * h * k1_i, beta)
    k2_e = exposed(s + 0.5 * h * k1_s, e + 0.5 * h * k1_e, i + 0.5 * h * k1_i, beta, sigma)
    k2_i = infected(e + 0.5 * h * k1_e, i + 0.5 * h * k1_i, sigma, gamma)
    k2_r = removed(i + 0.5 * h * k1_i, gamma)

    k3_s = susceptible(s + 0.5 * h * k2_s, e + 0.5 * h * k2_e, i + 0.5 * h * k2_i, beta)
    k3_e = exposed(s + 0.5 * h * k2_s, e + 0.5 * h * k2_e, i + 0.5 * h * k2_i, beta, sigma)
    k3_i = infected(e + 0.5 * h * k2_e, i + 0.5 * h * k2_i, sigma, gamma)
    k3_r = removed(i + 0.5 * h * k2_i, gamma)

    k4_s = susceptible(s + h * k3_s, e + h * k3_e, i + h * k3_i, beta)
    k4_e = exposed(s + h * k3_s, e + h * k3_e, i + h * k3_i, beta, sigma)
    k4_i = infected(e + h * k3_e, i + h * k3_i, sigma, gamma)
    k4_r = removed(i + h * k3_i, gamma)

    s_new = s + (h / 6) * (k1_s + 2 * k2_s + 2 * k3_s + k4_s)
    e_new = e + (h / 6) * (k1_e + 2 * k2_e + 2 * k3_e + k4_e)
    i_new = i + (h / 6) * (k1_i + 2 * k2_i + 2 * k3_i + k4_i)
    r_new = r + (h / 6) * (k1_r + 2 * k2_r + 2 * k3_r + k4_r)


    return s_new, e_new, i_new, r_new


def model_R0(beta, sigma, gamma, title):
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
                        exposed[i - 1], infected[i - 1], removed[i - 1], beta, sigma, gamma, h)
        

    plt.plot(time, susceptible, label="Susceptible")
    plt.plot(time, exposed, label="Exposed")
    plt.plot(time, infected, label="Infected")
    plt.plot(time, removed, label="Removed")
    plt.xlabel("Czas")
    plt.ylabel("Populacja")
    plt.legend()
    plt.title(title)
    plt.show()
    
    

R0 = 0.5
beta = R0 * gamma
model_R0(beta, sigma, gamma, f'Model SEIR dla R0 = {R0}')

#R0 ustawiono na 0.5, co oznacza, że średnio każda zainfekowana osoba zaraża 0.5 innych osób.
#R0<1, nie dochodzi do rozwoju epidemii, a epidemia wygasa.


R1 = 50
beta1 = R1 * gamma
model_R0(beta1, sigma, gamma, f'Model SEIR dla R0 = {R1}')

#R0 ustawiono na 50, co wskazuje na znacznie gwałtowniejsze rozprzestrzeniania się epidemii.
#Choroba ma istotny wpływ na populację, na co wskazuje liczba osób zainfekowanych oraz wystawionych na działanie wirusa.
#Wysoka wartość R0 sugeruje wyższą zdolność do przenoszenia się choroby.
#Prowadzi to do bardziej rozległej epidemii.





















