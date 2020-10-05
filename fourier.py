import math
import random
import matplotlib.pyplot as plt
import numpy as np
import bisect
from scipy.integrate import odeint
from numpy.fft import fft, fftfreq
import scipy.signal as signal

def prob_sum(prob):
    """Calculates the sum of the probabilities"""
    total = sum(prob)
    result = []
    val = 0
    for p in prob:
        val += p
        result.append(val / total)
    return result

N = 2100
T = 1000
t = 0.0

beta = 1.5
gamma = 0.5
mu = 0.02
nu = 0.02
alpha = 0.002

I = 100
S = N - I
R = 0

data = []
data.append((t, S, I, R))
extinctions = 0

switch = 0
while t < T:
  if I == 0 and switch == 0:
    print("Extinction")
    extinctions += 1
    switch = 1

  if I > 0:
      switch = 0
    
  event1 = beta * S * I / N
  event2 = gamma * I
  event3 = mu * N
  event4 = nu * S
  event5 = nu * I
  event6 = nu * R
  event7 = alpha * math.sqrt(N)
  E_tot = event1 + event2 + event3 + event4 + event5 + event6 + event7
  
  dt = -math.log(random.uniform(0.0, 1.0)) / E_tot
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
s_list = list(map(lambda x: data[x][1], np.arange(0, len(data), 1)))
i_list = list(map(lambda x: data[x][2], np.arange(0, len(data), 1)))
r_list = list(map(lambda x: data[x][3], np.arange(0, len(data), 1)))
tot_list = list(map(lambda x: data[x][1] + data[x][2] + data[x][3], np.arange(0, len(data), 1)))

##### LOWPASS FILTER FOR INFECTEDS#####
raw_data = i_list

N_order = 3
Wn = 0.003
B, A = signal.butter(N_order, Wn, output='ba')
smooth_data = signal.filtfilt(B,A, raw_data)

########################################
  
# print(t_list)
# print(len(data))
plt.plot(t_list, smooth_data)
#plt.plot(t_list, i_list)
plt.plot(t_list, r_list)
plt.plot(t_list, s_list)
plt.plot(t_list, tot_list)
#plt.show()

s_mean = np.mean(s_list)
i_mean = np.mean(i_list)

covariance_si = 0

for i in range(len(s_list)):
    covariance_si += (s_list[i] - s_mean) * (i_list[i] - i_mean)

covariance_si /= len(s_list)

print(covariance_si)

def calc_ode3(b_val, g_val, death_rate, birth_rate):
    
    # initial conditions for the problem
    x0 = [N-100, 100, 0]
    
    # up to 50 days after incident
    t = np.linspace(0, 1000, 1000)
    
    # pop size
    #N = 1000
    
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



data_ode = calc_ode3(beta,gamma,mu,nu)
total_pop = []
for i in range(len(data_ode[0])):
    total_pop.append(data_ode[0][i] + data_ode[1][i] + data_ode[2][i])
time_list = np.linspace(0,1000,1000)
plt.plot(time_list, data_ode[1], '--b', label='I')
plt.plot(time_list, data_ode[2], linestyle='dashed', color='orange', label='R')
plt.plot(time_list, data_ode[0], '--g', label='S')
plt.plot(time_list, total_pop, '--r', label='T')
plt.legend(loc='upper right')
plt.show()

print(extinctions)

cutoff = 200

# add fourier analysis
intervals = cutoff
t = np.linspace(0, 1, intervals)
fourier = fft(data_ode[1][:cutoff])
fourier_abs = (fourier.real**2 + fourier.imag**2)**(1/2)
fourier2 = fft(smooth_data[:cutoff])
fourier_abs2 = (fourier2.real**2 + fourier2.imag**2)**(1/2)

plt.subplot(1,2,1)
plt.bar(t[1:intervals // 2], np.abs(fourier)[1:intervals // 2]*(1/intervals), width=0.003)

plt.subplot(1,2,2)
plt.bar(t[1:intervals // 2], np.abs(fourier2)[1:intervals // 2]*(1/intervals), width=0.003)
plt.show()

