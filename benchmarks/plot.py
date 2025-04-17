import matplotlib.pyplot as plt
import numpy as np

# Read the data (optool)
optool_data = np.loadtxt(
    "/Users/kmaitreys/Documents/software/optool/dustkappa.dat", skiprows=24
)
optool_lam = optool_data[:, 0]
optool_kabs = optool_data[:, 1]
optool_ksca = optool_data[:, 2]
optool_g = optool_data[:, 3]

# Read the data (kappa)
kappa_data = np.loadtxt(
    "/Users/kmaitreys/Documents/ProsiaLAB/kappa/output/kappa_opacity.dat", skiprows=23
)
kappa_lam = kappa_data[:, 0]
kappa_kabs = kappa_data[:, 1]
kappa_ksca = kappa_data[:, 2]
kappa_g = kappa_data[:, 3]

# Plot the data
fig, ax = plt.subplots(1, 3, figsize=(12, 4))

ax[0].plot(optool_lam, optool_kabs, label="optool")
ax[0].plot(kappa_lam, kappa_kabs, label="kappa")
ax[0].set_xlabel("Wavelength [micron]")
ax[0].set_ylabel("Absorption coefficient")
ax[0].set_xscale("log")
ax[0].set_yscale("log")
ax[0].legend()

ax[1].plot(optool_lam, optool_ksca, label="optool")
ax[1].plot(kappa_lam, kappa_ksca, label="kappa")
ax[1].set_xlabel("Wavelength [micron]")
ax[1].set_ylabel("Scattering coefficient")
ax[1].set_xscale("log")
ax[1].set_yscale("log")
ax[1].legend()

ax[2].plot(optool_lam, optool_g, label="optool")
ax[2].plot(kappa_lam, kappa_g, label="kappa")
ax[2].set_xlabel("Wavelength [micron]")
ax[2].set_ylabel("Asymmetry factor")
ax[2].set_xscale("log")
ax[2].set_yscale("log")
ax[2].legend()


plt.show()
