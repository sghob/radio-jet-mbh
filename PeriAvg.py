import matplotlib.pyplot as plt

# read the data from the txt file
with open('possible.txt', 'r') as f:
    data = f.readlines()

# extract x and y values from data
x = []
y = []
for line in data:
    values = line.split()
    x.append(float(values[0])/(3600 * 24 * 365 * 1000000000)) #/(3600 * 24 * 365 * 1000000000)
    y.append(float(values[1]))


#construction of average a0 data using integration
with open('avgset.txt', 'w') as f:
    pstart = 0
    isum = 0

    for i in range(1, len(y) - 1):

        if y[i] > y[i - 1] and y[i + 1] < y[i]:
            tau = x[i] - x[pstart]
            a_avg = isum / tau
            time = x[pstart] + tau / 2
            f.write("%E\t%E\n" % (time, a_avg))
            pstart = i
            isum = 0

        dt = x[i] - x[i - 1]
        isum += dt * y[i]

with open('avgset.txt', 'r') as f:
    avg_data = f.readlines()  

x_avg = []
y_avg = []

for line in avg_data:
    vals = line.split()
    x_avg.append(float(vals[0]))
    y_avg.append(float(vals[1]))


# create a scatter plot
plt.plot(x, y, '-')

# set the axis labels
plt.xlabel('Time (Gyr)')
plt.ylabel('Semi-major Axis (kpc)')

# set the plot title
plt.title('Orbital Evolution')

# show the plot
plt.show()