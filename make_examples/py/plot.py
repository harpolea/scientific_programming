import numpy as np
import matplotlib.pyplot as plt
plt.switch_backend('agg')
import sys

def plot(input_filename='../data/data.csv', output_name='../figures/fig.pdf'):

    # read input file
    dat = np.loadtxt(input_filename, dtype=None, delimiter=',', skiprows=1)

    print(dat)
    plt.plot(dat[:,0], dat[:,1])
    plt.savefig(output_name)


if __name__ == '__main__':
    if len(sys.argv) > 1:
        input_filename = sys.argv[1]
        filename = sys.argv[2]

        plot(input_filename=input_filename, filename=filename)

    else:
        plot()
