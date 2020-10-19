import math
import random
import numpy as np
import matplotlib.pyplot as plt
import bisect
from scipy.integrate import odeint
from scipy.signal import lfilter
from scipy.fft import fft
#import seaborn as sns
#import pandas as pd

N = 1000
T = 500
t = 0.0
beta = [1.5, 2.5]
gamma = [0.5, 0.5]
mu = 0.02
nu = 0.02
alpha = 0#0.01
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

def disc_event(t, S, I, R, pops, beta, gamma):
    data = []
    data.append((t, S, I, R))

    while t < T:
        if I == 0:
            break

        event1 = beta * S * I / pops
        event2 = gamma * I
        event3 = mu * pops
        event4 = nu * S
        event5 = nu * I
        event6 = nu * R
        # event7 = alpha * math.sqrt(N)
        E_tot = event1 + event2 + event3 + event4 + event5 + event6

        dt = -math.log(1-random.uniform(0.0, 1.0)) / E_tot
        t += dt

        event_list = prob_sum([event1, event2, event3, event4, event5, event6])

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

        else:
            R -= 1

        data.append((t, S, I, R))

    t_list = list(map(lambda x: data[x][0], np.arange(0, len(data), 1)))
    S_list = list(map(lambda x: data[x][1], np.arange(0, len(data), 1)))
    I_list = list(map(lambda x: data[x][2], np.arange(0, len(data), 1)))
    R_list = list(map(lambda x: data[x][3], np.arange(0, len(data), 1)))

    # plt.plot(t_list, I_list, label=f'beta={beta}, gamma={gamma}, pop={pops}')
    # plt.plot(t_list, S_list)
    # plt.legend()
    # plt.show()
    return t_list, S_list, I_list, R_list

#plt.show()
fig = plt.figure(figsize=(9,4))
ax1 = fig.add_subplot(1,1,1)

for k in range(len(beta)):
    # Phase plots for transients
    covariance = {}
    pop_sizes = [2000, 2500, 3000, 3500, 4000, 4500, 5000, 5500, 6000]
    for j in range(10):
        for i in pop_sizes:
            
            t_list, S_list, I_list, R_list = disc_event(t, S, I, R, i, beta[k], gamma[k])

            while t_list[-1] < 400:
                t_list, S_list, I_list, R_list = disc_event(t, S, I, R, i, beta[k], gamma[k])


            starting_index = 0
            # for t in range(len(t_list)):
            #     if t_list[t] > 30:
            #         starting_index = t 
            #         print(t)
            #         break

            #stop = int(i*(15))
            s_mean = np.mean(S_list[starting_index:])
            i_mean = np.mean(I_list[starting_index:])

            covariance_si = 0

            

            for i in range(len(S_list[starting_index:])):
                covariance_si += (S_list[i+starting_index] - s_mean) * (I_list[i+starting_index] - i_mean)

            covariance_si /= len(S_list[starting_index:])
            if j in covariance.keys():
                covariance[j].append(covariance_si)
            else:
                covariance[j] = [covariance_si]

    av_Covar = []
    std_Covar = []
    for j in range(len(pop_sizes)):
        temp = []
        for i in range(10):
            temp.append(covariance[i][j])
        av_Covar.append(np.mean(temp))
        std_Covar.append(np.std(temp))

    ax1.errorbar(pop_sizes, av_Covar, std_Covar, capsize=5, label=f"Beta: {beta[k]}, Gamma: {gamma[k]}")

time_list = np.arange(300,5500,1)
null_list = np.zeros((len(time_list), 1))
#ax1.plot(time_list, null_list, "--", color="black")
ax1.legend(loc="best")
ax1.set_xlabel("Population Size")
ax1.set_ylabel("Covariance")
ax1.set_title("Covariance between Infected and Susceptibles against Population Size")
ax1.set_xlim(1900,6100)
plt.show()
