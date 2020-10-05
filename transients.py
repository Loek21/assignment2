import math
import random
import numpy as np
import matplotlib.pyplot as plt
import bisect
from scipy.integrate import odeint
from scipy.signal import lfilter
from scipy.fft import fft

N = 1000
T = 1000
t = 0.0
beta = 1.5
gamma = 0.5
mu = 0.02
nu = 0.02
alpha = 0.1
I = 100
S = N - I
R = 0

def prob_sum(prob):
    """Calculates the sum of the probabilities"""
    total = sum(prob)
    result = []
    val = 0
    for p in prob:
        val += p
        result.append(val / total)
    return result

def disc_event(t, S, I, R):
    data = []
    data.append((t, S, I, R))

    while t < T:
        if I == 0:
            print("Extinction")
            break

        event1 = beta * S * I / N
        event2 = gamma * I
        event3 = mu * N
        event4 = nu * S
        event5 = nu * I
        event6 = nu * R
        event7 = alpha * math.sqrt(N)
        E_tot = event1 + event2 + event3 + event4 + event5 + event6 +event7

        dt = -math.log(1-random.uniform(0.0, 1.0)) / E_tot
        t += dt

        event_list = prob_sum([event1, event2, event3, event4, event5, event6, event7])

        selected_event = bisect.bisect(event_list, random.uniform(0.0, 1.0))

        if selected_event == 0:
            S -= 1
            I += 1
        elif selected_event == 1:
            R += 1
            I -= 1

        elif selected_event == 2:
            S += 1

        elif selected_event == 3:
            S -= 1

        elif selected_event == 4:
            I -= 1

        elif selected_event == 5:
            R -= 1

        else:
            I += 1

        data.append((t, S, I, R))

    t_list = list(map(lambda x: data[x][0], np.arange(0, len(data), 1)))
    S_list = list(map(lambda x: data[x][1], np.arange(0, len(data), 1)))
    I_list = list(map(lambda x: data[x][2], np.arange(0, len(data), 1)))
    R_list = list(map(lambda x: data[x][3], np.arange(0, len(data), 1)))
    return t_list, S_list, I_list, R_list

def calc_ode3(b_val, g_val, death_rate, birth_rate, S, I, R, N):

    # initial conditions for the problem
    x0 = [S, I, R]

    # up to 50 days after incident
    t = np.linspace(0, 100, 10000)

    # pop size
    # N = 1000

    def SIRmodel3(x,t):
        s = x[0]    # susceptible
        i = x[1]    # infected
        r = x[2]    # recovered

        # define the parameters for the model
        b = b_val    # rate of infection
        g = g_val    # rate of recovery
        d_r = death_rate
        b_r = birth_rate

        # get the ODE's, this time with birth and death rates
        dsdt = (-b * s * i) / N + b_r * N - d_r * s
        didt = (b * s * i) / N - g * i - d_r * i + alpha * math.sqrt(N)
        drdt = g * i - d_r * r

        return [dsdt, didt, drdt]

    # calculate
    x = odeint(SIRmodel3, x0, t)

    # get the values for s, i and r from the array x
    s = x[:, 0]
    i = x[:, 1]
    r = x[:, 2]


    return s, i, r

S_D, I_D, R_D = calc_ode3(beta, gamma, mu, nu, S, I, R, N)

t_list, S_list, I_list, R_list = disc_event(t, S, I, R)

# Phase plots for transients
S_dict = {}
I_dict = {}
for i in range(5):
    t_list, S_list, I_list, R_list = disc_event(t, S, I, R)
    S_dict[i] = S_list[:20000]
    I_dict[i] = I_list[:20000]

av_S_list = []
av_I_list = []
for i in range(20000):
    tempS = []
    tempI = []
    for j in range(5):
        tempS.append(S_dict[j][i])
        tempI.append(I_dict[j][i])
    av_S_list.append(np.mean(tempS))
    av_I_list.append(np.mean(tempI))

s_mean = np.mean(S_list)
i_mean = np.mean(I_list)

S_D, I_D, R_D = calc_ode3(beta, gamma, mu, nu, S, I, R, N)

fig = plt.figure()
ax1 = fig.add_subplot(1,1,1)
# ax1.plot(S_dict[0], I_dict[0])
# ax1.plot(S_dict[1], I_dict[1])
# ax1.plot(S_dict[2], I_dict[2])
# ax1.plot(S_dict[3], I_dict[3])
# ax1.plot(S_dict[4], I_dict[4])
ax1.plot(av_S_list, av_I_list, linewidth=1, label="Stochastic")
ax1.plot(S_D, I_D, label="Deterministic")
ax1.legend(loc="best")
plt.show()
