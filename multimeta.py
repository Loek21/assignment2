from multimetaclasses import Subpopulation
import matplotlib.pyplot as plt
import numpy as np
import math
import bisect


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


# initiate the subpopulations
for i in range(3):
    for j in range(3):
        subpopulation = Subpopulation(1000, i, j, 1000, 0, 0, 0.1, 0.5)
        subpopulations.append(subpopulation)

# time span
t = 0.0
T = 10

time_list = []

# make dict to save the population data
data_dict = {}
for q in range(9):
    data_dict[f'pop{q}'] = []

# change the first population to have some infecteds
subpopulations[0].S -= 5
subpopulations[0].I += 5

while t < T:
    eventlist = []
    eventnames = []

    # create all possible events for every subpopulation
    for i in range(9):
        
        pop_a = subpopulations[i]
        x1, y1 = pop_a.x, pop_a.y

        # from I to R for chosen subpopulation
        eventname_IR = f'IR_{i}'
        eventweight_IR = pop_a.gamma * pop_a.I

        eventlist.append(eventweight_IR)
        eventnames.append(eventname_IR)

        for j in range(9):

            pop_b = subpopulations[j]
            x2, y2 = pop_b.x, pop_b.y

            # calculate coupling between pops based on distance
            #rho = abs(1 - np.random.normal(0,1)*np.sqrt((abs(x1-x2))**2+(abs(y1-y2))**2) / (2*np.sqrt(2)))
            rho = abs(1 - np.sqrt((abs(x1-x2))**2+(abs(y1-y2))**2) / (2*np.sqrt(2)))

            eventname_SI = f'SI_{i}_{j}'
            eventweight_SI = pop_a.beta * pop_a.S * pop_b.I * rho * pop_a.size
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

    pop_nr = int(name_event[3])
    if name_event[0:2] == 'SI':
        subpopulations[pop_nr].S -= 1
        subpopulations[pop_nr].I += 1

    elif name_event[0:2] == 'IR':
        subpopulations[pop_nr].I -= 1
        subpopulations[pop_nr].R += 1

    time_list.append(t)
    #print(t)

    for k in range(9):

        SIRdata = (subpopulations[k].S, subpopulations[k].I, subpopulations[k].R)
        data_dict[f'pop{k}'].append(SIRdata)




for w in range(1):
    test_dataS = []
    test_dataI = []
    test_dataR = []

    for i in range(len(data_dict[f'pop{w}'])):
        test_dataS.append(data_dict[f'pop{w}'][i][0])
        test_dataI.append(data_dict[f'pop{w}'][i][1])
        test_dataR.append(data_dict[f'pop{w}'][i][2])

    plt.plot(time_list, test_dataS, label=f'S{w}')
    plt.plot(time_list, test_dataI, label=f'I{w}')
    plt.plot(time_list, test_dataR, label=f'R{w}')
plt.legend()
#plt.show()
print(data_dict['pop0'][0:50])