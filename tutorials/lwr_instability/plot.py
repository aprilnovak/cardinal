import matplotlib.pyplot as plt
import csv
import os
import numpy as np
import common as specs

plt.rcParams.update({'font.size': 14})

def generate_colors(n: int):
    cmap = plt.get_cmap("rainbow")
    return cmap(np.linspace(1, 0, max(n, 1)))

files = [filename for filename in os.listdir('.') if filename.startswith('openmc_out_axial_power')]
files = sorted(files)

x = np.linspace(0, specs.heated_length, specs.n_layers)
colors = generate_colors(len(files))

for i, f in enumerate(files):
  if (i == 0):
    continue

  power = []
  with open(f) as file:
    reader = csv.reader(file)
    next(reader) # skip header
    for row in reader:
      power.append(float(row[0]))

  plt.plot(x, power, label='Iteration {}'.format(i), color=colors[i])

plt.legend()
plt.grid()
plt.xlabel('Axial Position [m]')
plt.ylabel('Power Density [W/m$^3$]')
plt.savefig('power.png')
plt.close()
