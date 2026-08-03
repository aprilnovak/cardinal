import csv
import os
import matplotlib.pyplot as plt
import numpy as np

plt.rcParams.update({
    "font.size": 9,
    "axes.titlesize": 10,
})

# mdot_values = [1.0, 3.2, 5.4]
# fuel_values = [0.5, 1.0, 5.0]
mdot_values = [1.2, 10.0]
fuel_values = [2.0, 20.0]


nrows = len(fuel_values)
ncols = len(mdot_values)

fig, axs = plt.subplots(
    nrows=nrows,
    ncols=ncols,
    figsize=(12, 10),
    sharex=True,
    sharey=True
)

if nrows == 1 and ncols == 1:
    axs = np.array([[axs]])
else:
    axs = np.atleast_2d(axs)
    if nrows == 1:
        axs = axs.reshape(1, -1)
    elif ncols == 1:
        axs = axs.reshape(-1, 1)


im = None

for i, fuel in enumerate(fuel_values):
    for j, mdot in enumerate(mdot_values):
        ax = axs[i, j]
        folder = f"results_m{mdot}_f{fuel}"

        files = sorted([f for f in os.listdir(folder) if f.startswith("openmc_out_axial_power")])

        profiles = []

        for file in files:
            if file.endswith("_0000.csv") or "_0000" in file:
                continue

            power = []
            with open(os.path.join(folder, file)) as f:
                reader = csv.reader(f)
                next(reader)

                for row in reader:
                    power.append(float(row[0]))

            power = np.array(power)

            if np.mean(power) > 0:
                power /= np.mean(power)

            profiles.append(power)

        heatmap = np.array(profiles).T

        im = ax.imshow(heatmap,origin="lower",aspect="auto",cmap="jet",interpolation="nearest")

        if i == 0:
            ax.set_title(rf"$\dot{{m}} = {mdot}$")

        if j == ncols - 1:
            ax.text(1.05, 0.5, f"$k_f = {fuel}$", transform=ax.transAxes, va="center", ha="left", fontsize=11)

        if i == nrows - 1:
            ax.set_xlabel("Iteration")
            ax.set_xticks(np.arange(heatmap.shape[1]))
            ax.set_xticklabels(np.arange(1, heatmap.shape[1] + 1))

        if j == 0:
            ax.set_ylabel("Axial Cell")
            yticks = np.arange(0, heatmap.shape[0], 5)
            ax.set_yticks(yticks)
            ax.set_yticklabels(yticks + 1)

if im is not None:
    cbar_ax = fig.add_axes([0.15, 0.05, 0.7, 0.02])
    fig.colorbar(im, cax=cbar_ax, orientation="horizontal", label="Normalized Power")

plt.subplots_adjust(bottom=0.15, hspace=0.2, wspace=0.1)
plt.savefig("grid_power_heatmaps.png", dpi=300)
plt.close()
