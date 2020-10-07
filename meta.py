import math
import random
import numpy as np
import matplotlib.pyplot as plt
import bisect

N = 1000
T = 1000
t = 0.0
beta = 1.5
gamma = 0.5
mu = 0.02
nu = 0.02
alpha = 0.1
I1 = 10
S1 = N - I1
R1 = 0
I2 = 0
S2 = N - I2
R2 = 0
rho12 = 0.1
rho21 = 0.1

def prob_sum(prob):
    """Calculates the sum of the probabilities"""
    total = sum(prob)
    result = []
    val = 0
    for p in prob:
        val += p
        result.append(val / total)
    return result

def disc_event(t, S1, I1, R1, rho12, S2, I2, R2, rho21):
    data = []
    data.append((t, S1, I1, R1, S2, I2, R2))

    while t < T:
        if I1 == 0:
            print("Extinction")
            break

        event1 = beta * S1 * I1 / N
        event2 = beta * S1 * I2 * rho12 / N
        event3 = beta * S2 * I2 / N
        event4 = beta * S2 * I1 * rho21 / N
        event5 = gamma * I1
        event6 = gamma * I2
        E_tot = event1 + event2 + event3 + event4 + event5 + event6

        dt = -math.log(1-random.uniform(0.0, 1.0)) / E_tot
        t += dt

        event_list = prob_sum([event1, event2, event3, event4, event5, event6])

        selected_event = bisect.bisect(event_list, random.uniform(0.0, 1.0))

        if selected_event == 0:
            S1 -= 1
            I1 += 1
        elif selected_event == 1:
            S1 -= 1
            I1 += 1

        elif selected_event == 2:
            S2 -= 1
            I2 += 1

        elif selected_event == 3:
            S2 -= 1
            I2 += 1

        elif selected_event == 4:
            I1 -= 1
            R1 += 1

        elif selected_event == 5:
            I2 -= 1
            R2 += 1

        data.append((t, S1, I1, R1, S2, I2, R2))

    t_list = list(map(lambda x: data[x][0], np.arange(0, len(data), 1)))
    S1_list = list(map(lambda x: data[x][1], np.arange(0, len(data), 1)))
    I1_list = list(map(lambda x: data[x][2], np.arange(0, len(data), 1)))
    R1_list = list(map(lambda x: data[x][3], np.arange(0, len(data), 1)))
    S2_list = list(map(lambda x: data[x][4], np.arange(0, len(data), 1)))
    I2_list = list(map(lambda x: data[x][5], np.arange(0, len(data), 1)))
    R2_list = list(map(lambda x: data[x][6], np.arange(0, len(data), 1)))

    return t_list, S1_list, I1_list, R1_list, S2_list, I2_list, R2_list

t_list, S1_list, I1_list, R1_list, S2_list, I2_list, R2_list = disc_event(t, S1, I1, R1, rho12, S2, I2, R2, rho21)

max1 = np.argmax(I1_list)
max2 = np.argmax(I2_list)
delay_time = t_list[max2] - t_list[max1]

print(delay_time)

fig = plt.figure()
ax1 = fig.add_subplot(1,1,1)
ax1.plot(t_list, S1_list, label="Susceptible 1")
ax1.plot(t_list, I1_list, label="Infected 1")
ax1.plot(t_list, R1_list, label="Recovered 1")
ax1.plot(t_list, S2_list, label="Susceptible 2")
ax1.plot(t_list, I2_list, label="Infected 2")
ax1.plot(t_list, R2_list, label="Recovered 2")
ax1.legend(loc="best")
plt.show()
