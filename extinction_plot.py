import matplotlib.pyplot as plt 
import csv
import numpy as np

# population sizes used in the extinction experiments
pop_sizes = [1000, 1250, 1500, 1750, 2000, 2250, 2500, 2750, 3000]


extinctions_mean = []
extinctions_std = []
with open('extinctions.csv') as csv_file:
    csv_reader = csv.reader(csv_file, delimiter=';')

    data = []
    for row in csv_reader:
        for number in row:
            data.append(float(number))

        npdata = np.array(data)
        extinctions_mean.append(npdata.mean())
        extinctions_std.append(npdata.std())
        data = []

print(extinctions_mean)
print(extinctions_std)

fig, ax = plt.subplots()

ax.errorbar(pop_sizes, extinctions_mean,yerr=extinctions_std, capsize=5)
ax.set_ylabel('Extinctions per 1000 days')
ax.set_xlabel('Population size')
plt.show()


ftps_mean = []
ftps_std = []
with open('first_passage.csv') as csv_file:
    csv_reader = csv.reader(csv_file, delimiter=';')

    data = []
    for row in csv_reader:
        for number in row:
            if float(number) > 0:
                data.append(float(number))

        npdata = np.array(data)
        ftps_mean.append(npdata.mean())
        ftps_std.append(npdata.std())
        data = []

print(ftps_mean)
print(ftps_std)

fig, ax = plt.subplots()

ax.errorbar(pop_sizes, ftps_mean,yerr=ftps_std, capsize=5)
ax.set_ylabel('First passage time')
ax.set_xlabel('Population size')
plt.show()


        