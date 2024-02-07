#!/usr/bin/python

import os
import sys
import math
import numpy as np
from scipy import signal 
import popFOAM as pop
import matplotlib.pyplot as plt


forces_file = "postProcessing/forces/0/forces.dat"

if not os.path.isfile(forces_file):
	print("Forces file not found at ", forces_file)
	print("Be sure that the case has been run and you have the right directory!")
	print("Exiting.")
	sys.exit()

time = np.array([])
forceX = np.array([])
forceY = np.array([])

time, forceX, forceY, _ = pop.createForceFile(forces_file)

average_forceX = np.mean(forceX)
average_forceY = np.mean(forceY)

rho = 998
draft = 5
motionAmp = 0.5
wprime = 0.3
w = np.sqrt(wprime*9.81/draft)
truncMax = 75
truncMin = 0.2

min_truncate_index = np.argmax(time >= truncMin)
max_truncate_index = np.argmax(time >= truncMax)

time_truncated = time[min_truncate_index:max_truncate_index]
forceX_truncated = forceX[min_truncate_index:max_truncate_index]
forceY_truncated = forceY[min_truncate_index:max_truncate_index]

R = draft
alpha = np.arccos(motionAmp/R)
Awp = 2 * R * np.sin(alpha)


forceY_truncated = forceY_truncated - np.average(forceY_truncated)

dt = time_truncated[1] - time_truncated[0] 
n = len(forceY_truncated)
forceY_truncated_hat = np.fft.fft(forceY_truncated,n) # inverse FFT
PSD = forceY_truncated_hat * np.conj(forceY_truncated_hat) / n # Power spectrum density
freq = 1 / (dt * n) * np.arange(n) # create x axis on frequencies

indices = PSD > 0.5 * 10**10
PSD_clean = PSD * indices
forceY_truncated_hat = indices * forceY_truncated_hat
forceY_truncated_filtered = np.fft.ifft(forceY_truncated_hat)

fig, axs = plt.subplots(2, 1)

plt.sca(axs[0])
plt.plot(time_truncated, forceY_truncated, color = 'c', linewidth=1.5, label = 'Noisy')
plt.plot(time_truncated, forceY_truncated_filtered, color = 'r', linewidth=1.5, label = 'Clean')
plt.legend()

plt.sca(axs[1])
plt.semilogy(freq, PSD, color = 'c', linewidth=2, label = 'Noisy')
plt.semilogy(freq, PSD_clean, color = 'r', linewidth=2, label = 'Clean')
plt.legend()
plt.show()

print(np.max(PSD))
'''

# Filter the data
window_length = 500
poly_order = 4
hydrostaticforceAVG = np.average(forceY_truncated)
#hydrostaticforceEmp = 9.81 * rho * math.pi * R**2 / 2
#print(hydrostaticforceAVG/hydrostaticforceEmp*100)

forceY_truncated = forceY_truncated - hydrostaticforceAVG

forceY_filtered = signal.savgol_filter(forceY_truncated, window_length, poly_order)
# Example: Finding peaks in the forceY_truncated array
peaks, _ = signal.find_peaks(np.absolute(forceY_filtered), prominence=50, threshold=1)

# Getting the time and force values at these peaks
peak_times = time_truncated[peaks]
peak_forceY = forceY_filtered[peaks]

print(np.average(np.absolute(peak_forceY)), '\n')

# Printing or processing the peak values
for t, f in zip(peak_times, peak_forceY):
    print(f"Peak at time {t}: Force Y = {f}")

plt.plot(time_truncated, forceY_truncated, alpha = 0.5, marker = '', color = 'grey')
plt.plot(time_truncated, forceY_filtered, alpha = 0.5, marker = '', color = 'yellow')
plt.plot(peak_times, peak_forceY, alpha = 0.5, marker = 'o', color = 'red', linestyle='')
plt .show()

'''

#motion_signal = motionAmp * np.sin(w * trunc_time)
'''
# Calculate the cross-correlation
cross_correlation = np.correlate(trunc_forceY - np.mean(trunc_forceY), motion_signal - np.mean(motion_signal), 'full')

# Calculate the delays for each cross-correlation value
delays = np.arange(1-len(motion_signal), len(motion_signal))

# Find the delay corresponding to the maximum cross-correlation value
delay = delays[np.argmax(cross_correlation)]

print(f"The phase lag (delay) between the force and motion is: {delay} time units")

						
# Identify the indices of the peaks in the truncated force data
peak_indices = argrelextrema(np.array(trunc_forceY), np.greater)[0]
peak_forces = [trunc_forceY[i] for i in peak_indices]

for peak_force in peak_forces:
    # Instantiate the HydroCoeff class object
    hydro_coeff = pop.HydroCoeff(delay, peak_force, motionAmp, w)
    damping = hydro_coeff.damping/(np.pi() / 2 * 998.2 * draft**2)
    added_mass = hydro_coeff.addedmass/(np.pi() / 2 * 998.2 * draft**2)

    # Print the calculated damping for each peak
    print(f"For peak force {peak_force}, the calculated damping is: {damping} units")			

'''