import matplotlib.pyplot as plt
import numpy as np


def load_data(filename):
    with open(filename) as f:
        lines = [line for line in f if not line.startswith("#")]
        # skip first 2 non-comment lines
        data = np.loadtxt(lines[2:])
    return data


c1_data = load_data(
    "/Users/kmaitreys/Projects/ProsiaLAB/scorpia/examples/default/kappa_opacity_comp_1.dat"
)  # Diana
c2_data = load_data(
    "/Users/kmaitreys/Projects/ProsiaLAB/scorpia/examples/default/kappa_opacity_comp_2.dat"
)  # Dsharp
c3_data = load_data(
    "/Users/kmaitreys/Projects/ProsiaLAB/scorpia/examples/default/kappa_opacity_comp_3.dat"
)  # Custom


# Plot the data
fig, ax = plt.subplots(1, 3, figsize=(12, 4))

ax[0].plot(c1_data[:, 0], c1_data[:, 1], label="Diana")
ax[0].plot(c2_data[:, 0], c2_data[:, 1], label="Dsharp")
ax[0].plot(c3_data[:, 0], c3_data[:, 1], label="Custom")
ax[0].set_xlabel("Wavelength [micron]")
ax[0].set_ylabel("Absorption coefficient")
ax[0].set_xscale("log")
ax[0].set_yscale("log")
ax[0].legend()

ax[1].plot(c1_data[:, 0], c1_data[:, 2], label="Diana")
ax[1].plot(c2_data[:, 0], c2_data[:, 2], label="Dsharp")
ax[1].plot(c3_data[:, 0], c3_data[:, 2], label="Custom")
ax[1].set_xlabel("Wavelength [micron]")
ax[1].set_ylabel("Scattering coefficient")
ax[1].set_xscale("log")
ax[1].set_yscale("log")
ax[1].legend()

ax[2].plot(c1_data[:, 0], c1_data[:, 3], label="Diana")
ax[2].plot(c2_data[:, 0], c2_data[:, 3], label="Dsharp")
ax[2].plot(c3_data[:, 0], c3_data[:, 3], label="Custom")
ax[2].set_xlabel("Wavelength [micron]")
ax[2].set_ylabel("Asymmetry factor")
ax[2].set_xscale("log")
ax[2].set_yscale("log")
ax[2].legend()


plt.show()
