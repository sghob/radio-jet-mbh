import matplotlib.pyplot as plt

# Open the input file in read mode and read the contents into a list
with open("ithink.txt", "r") as f:
    data = [line.strip() for line in f]

# Create a set to store the distinct pairs
unique_pairs = set()

# Loop through each line in the data list and add distinct pairs to the set
for line in data:
    pair = tuple(line.split())
    if pair not in unique_pairs:
        unique_pairs.add(pair)

# Open the output file in write mode and write the unique pairs to it
with open("output.txt", "w") as f:
    for pair in unique_pairs:
        f.write(f"{pair[0]}\t{pair[1]}\n")

# read the data from the txt file
with open('output.txt', 'r') as f:
    data = f.readlines()

# extract x and y values from data
x = []
y = []
for line in data:
    values = line.split()
    x.append(float(values[0]))
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
plt.scatter(x, y)

# set the axis labels
plt.xlabel('X Axis Label')
plt.ylabel('Y Axis Label')

# set the plot title
plt.title('Scatter Plot of Two Columns of Data')

# show the plot
plt.show()