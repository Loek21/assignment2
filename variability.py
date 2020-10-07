import math
import random
import numpy as np
import matplotlib.pyplot as plt
import bisect
from scipy.integrate import odeint
from scipy.signal import lfilter
from scipy.fft import fft
import seaborn as sns
import pandas as pd

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

t_list, S_list, I_list, R_list = disc_event(t, S, I, R)

# Phase plots for transients
S_dict = {}
I_dict = {}
R_dict = {}
t_dict = {}
for i in range(5):
    t_list, S_list, I_list, R_list = disc_event(t, S, I, R)
    S_dict[i] = S_list[:20000]
    I_dict[i] = I_list[:20000]
    R_dict[i] = R_list[:20000]
    t_dict[i] = t_list[:20000]

av_S_list = []
av_I_list = []
av_R_list = []
av_t_list = []
std_S_list = []
std_I_list = []
std_R_list = []

for i in range(20000):
    tempS = []
    tempI = []
    tempR = []
    tempt = []
    for j in range(5):
        tempS.append(S_dict[j][i])
        tempI.append(I_dict[j][i])
        tempR.append(R_dict[j][i])
        tempt.append(t_dict[j][i])
    av_S_list.append(np.mean(tempS))
    av_I_list.append(np.mean(tempI))
    av_R_list.append(np.mean(tempR))
    av_t_list.append(np.mean(tempt))
    std_S_list.append(np.std(tempS))
    std_I_list.append(np.std(tempI))
    std_R_list.append(np.std(tempR))


s_mean = np.mean(S_list)
i_mean = np.mean(I_list)

plt.show()
fig = plt.figure()
ax1 = fig.add_subplot(1,1,1)
markers, caps, bars = ax1.errorbar(av_t_list, av_S_list, std_S_list, color="blue", linewidth=0.4, ecolor="orange", label="Susceptible")
[bar.set_alpha(0.8) for bar in bars]
markers, caps, bars = ax1.errorbar(av_t_list, av_I_list, std_I_list, color="red", linewidth=0.4, ecolor="orange", label="Infected")
[bar.set_alpha(0.8) for bar in bars]
markers, caps, bars = ax1.errorbar(av_t_list, av_R_list, std_R_list, color="green", linewidth=0.4, ecolor="orange", label="Recovered")
[bar.set_alpha(0.8) for bar in bars]
ax1.legend(loc="best")
ax1.set_xlabel("Time")
ax1.set_ylabel("Population")
ax1.set_title("Variability in Stochastic models")
plt.show()
