from multimetaclasses import Subpopulation
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
import math
import bisect
import time

GRIDLEN = 6

def prob_sum(prob):
    """Calculates the sum of the probabilities"""
    total = sum(prob)
    result = []
    val = 0
    for p in prob:
        val += p
        result.append(val / total)
    return result

# we save the subpopulations in a list
subpopulations = [] 

random_rho = [0.3, 0.4, 0.5, 0.6, 0.7]

# initiate the subpopulations
for i in range(GRIDLEN):
    for j in range(GRIDLEN):
        subpopulation = Subpopulation(1000, i, j, 1000, 0, 0, 0.4, 0.1, np.random.choice(random_rho))
        subpopulations.append(subpopulation)

# time span
t = 0.0
#T = 11.1
T = 10.5

time_list = []

# make dict to save the population data
data_dict = {}
for q in range(GRIDLEN**2):
    data_dict[f'pop{q}'] = []

# change the first population to have some infecteds
subpopulations[0].S -= 5
subpopulations[0].I += 5

# use this if we want to change rho values for specific subpopulations
# subpopulations[10].rho = 0.7
# subpopulations[25].rho = 0.7

while t < T:
    print(t)
    eventlist = []
    eventnames = []

    # create all possible events for every subpopulation
    for i in range(GRIDLEN**2):
        
        pop_a = subpopulations[i]
        x1, y1 = pop_a.x, pop_a.y

        # from I to R for chosen subpopulation
        eventname_IR = f'IR_{i}'
        eventweight_IR = pop_a.gamma * pop_a.I

        eventlist.append(eventweight_IR)
        eventnames.append(eventname_IR)

        for j in range(GRIDLEN**2):

            pop_b = subpopulations[j]
            x2, y2 = pop_b.x, pop_b.y

            # calculate coupling between pops based on distance
            #rho = abs(1 - np.random.normal(0,1)*np.sqrt((abs(x1-x2))**2+(abs(y1-y2))**2) / (2*np.sqrt(2)))
            if np.sqrt((abs(x1-x2))**2+(abs(y1-y2))**2) < 0.1:
                rho = 1
            elif np.sqrt((abs(x1-x2))**2+(abs(y1-y2))**2) <= 1.5:
                rho = pop_a.rho
            else:
                rho = 0
            #rho = abs(1 - np.sqrt((abs(x1-x2))**2+(abs(y1-y2))**2) / (2*np.sqrt(2)))

            eventname_SI = f'SI_{i}_{j}'
            eventweight_SI = pop_a.beta * pop_a.S * pop_b.I * rho / pop_a.size
            eventlist.append(eventweight_SI)
            eventnames.append(eventname_SI)

    # decide which event will take place
    event_list = prob_sum(eventlist)

    selected_event = bisect.bisect(event_list, np.random.uniform(0.0, 1.0))
    
    # total event rates for the time iterator
    event_total = sum(eventlist)
    dt = -math.log(np.random.uniform(0.0, 1.0)) / event_total
    t += dt

    name_event = eventnames[selected_event]

    
    pop_nr = int(name_event[3:5].strip('_'))

    if name_event[0:2] == 'SI':
        subpopulations[pop_nr].S -= 1
        subpopulations[pop_nr].I += 1

    elif name_event[0:2] == 'IR':
        subpopulations[pop_nr].I -= 1
        subpopulations[pop_nr].R += 1

    time_list.append(t)
    #print(t)

    for k in range(GRIDLEN**2):

        SIRdata = (subpopulations[k].S, subpopulations[k].I, subpopulations[k].R)
        data_dict[f'pop{k}'].append(SIRdata)

max_infections = 0
for w in range(GRIDLEN**2):
    test_dataS = []
    test_dataI = []
    test_dataR = []

    for i in range(len(data_dict[f'pop{w}'])):
        test_dataS.append(data_dict[f'pop{w}'][i][0])
        test_dataI.append(data_dict[f'pop{w}'][i][1])
        test_dataR.append(data_dict[f'pop{w}'][i][2])

        # get max infections for later
        if data_dict[f'pop{w}'][i][1] > max_infections:
            max_infections = data_dict[f'pop{w}'][i][1]

    #plt.plot(time_list, test_dataS, label=f'S{w}')
    plt.plot(time_list, test_dataI, label=f'I{w}')
    #plt.plot(time_list, test_dataR, label=f'R{w}')
plt.legend()
plt.show()


# heat map
datamatrix = np.random.random((GRIDLEN, GRIDLEN))

fig = plt.figure(figsize=(7.5,10))

for i in range(1,13):
    ax = fig.add_subplot(4,3,i)
    

    count = 0
    for w in range(GRIDLEN):
        for v in range(GRIDLEN):
            datamatrix[w][v] = data_dict[f'pop{count}'][math.floor((i-1) * len(time_list)/12)][1] / 1000

            coupling = subpopulations[count].rho
            ax.annotate(f'{coupling}', (w-0.25,v+0.1), fontsize=8)
            ax.axis('off')
            count += 1

            if count > GRIDLEN**2:
                count = 0
    im = ax.imshow(datamatrix, cmap=plt.get_cmap('OrRd'), vmin=0, vmax=max_infections/1000)
    #fig.colorbar(im)
fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
cbar = fig.colorbar(im, cax=cbar_ax)
cbar.set_label('Fraction of population infected', rotation=90)
plt.subplots_adjust(wspace=0.05, hspace=0.05)
plt.show()

            