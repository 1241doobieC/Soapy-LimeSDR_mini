#Written by Levi Powell, July 2023
#This code plots the data saved by the code in measureDelay.c


import time
import struct
import array
import matplotlib.pyplot as plt
import numpy as np
from numpy.fft import fft
from numpy.fft import fftshift

f_config_name = 'config.txt'
f_data_name = 'data.bin'
f_nulldata_name = 'nulldata.bin'
plotPhase = False
window = True
interpolate = False
fftLength = 100
intData = True # true => reading integer 


print("Reading configuration file...")

f_config = open(f_config_name, 'r')
line = f_config.read()
row = line.split(',')
sampleRate = int(row[0])
bandwidth = int(row[1])
numSamples = int(row[2])
numChirps = int(row[3])

print(f"\tSample rate was set to {sampleRate / 1e6} MS/s")
print(f"\tBandwidth was set to {bandwidth / 1e6} MHz")
print(f"\tNumber of samples was set to {numSamples}")
print(f"\tNumber of chirps was set to {numChirps}")

dataMatrix = np.empty([numChirps,numSamples], np.csingle)
realMatrix = np.empty([numChirps,numSamples], np.float32)
FFTMagMatrix = np.empty([numChirps,fftLength], np.float32)

dataMatrix2 = np.empty([numChirps,numSamples], np.csingle)
realMatrix2 = np.empty([numChirps,numSamples], np.float32)
FFTMagMatrix2 = np.empty([numChirps,fftLength], np.float32)

# dataMatrixDiff = np.empty([numChirps,numSamples], complex)
realMatrixDiff = np.empty([numChirps,numSamples], np.float32)
FFTMagMatrixDiff = np.empty([numChirps,fftLength], np.float32)

print("Opening data files...")
f_nulldata = open(f_nulldata_name, 'rb')
f_data = open(f_data_name, 'rb')

windowArray = np.float32(np.hanning(numSamples))



print(f"Reading {f_nulldata_name}...")
startTime = time.time()

for row in range(numChirps):
    # print(f"\tProcessing chirp #{row}...")
    
    if intData:
        reals = array.array('i')
        imags = array.array('i')
    else:
        reals = array.array('f')
        imags = array.array('f')
    reals.fromfile(f_nulldata, numSamples)
    imags.fromfile(f_nulldata, numSamples)
    reals = np.float32(reals)
    imags = np.float32(imags)

    
    # dataMatrix[row] = [complex(r,i) for r,i in zip(reals, imags)]
    dataMatrix[row] = reals + 1j*imags
    realMatrix[row] = reals

    del reals, imags
if window:
    dataMatrix = dataMatrix * windowArray
    realMatrix = realMatrix * windowArray

endTime = time.time()
print(f"\tTotal time was {endTime - startTime} s")



print(f"Calculating FFT for {f_nulldata_name}...")
startTime = time.time()

#calculate FFT
# for i in range(numChirps):
#     dataMatrix[i] = fftshift(fft(dataMatrix[i]))
# largeFFTMagMatrix = np.empty([numChirps,numSamples], np.csingle)
# print("created large matrix)")
# largeFFTMagMatrix = fftshift(fft(dataMatrix))
# print("took fft")
# FFTMagMatrix = largeFFTMagMatrix[:, numSamples//2 - fftLength//2:numSamples//2 + fftLength//2]
# print("copied part of fft")
# FFTMagMatrix = 10*np.log10(np.abs(FFTMagMatrix))
# print("took log")
# del largeFFTMagMatrix

# FFTMagMatrix = 10*np.log10(np.abs(fftshift(fft(dataMatrix))))[:, numSamples//2 - fftLength//2:numSamples//2 + fftLength//2]

for i in range(len(dataMatrix)):
    sr = sampleRate
    
    X = fftshift(fft(dataMatrix[i]))
    # N = len(X)
    # k = np.arange(N)
    # T = N/sr
    # freq = k/T

    X = X[numSamples//2 - fftLength//2:numSamples//2 + fftLength//2]

    FFTMagMatrix[i] = np.abs(X)
    FFTMagMatrix[i] = 10*np.log10(FFTMagMatrix[i])
    if plotPhase:
        FFTPhaseMatrix[i] = np.angle(X)
        # FFTPhaseMatrix[i] = np.angle(dataMatrix[i])
        for j in range(len(FFTPhaseMatrix[i])):
            if FFTPhaseMatrix[i][j] < 0:
                FFTPhaseMatrix[i][j] = 2*np.pi + FFTPhaseMatrix[i][j]

del dataMatrix

endTime = time.time()
print(f"\tTotal time was {endTime - startTime} s")



print(f"Reading {f_data_name}...")
startTime = time.time()

for row in range(numChirps):
    # print(f"\tProcessing chirp #{row}...")
    
    if intData:
        reals = array.array('i')
        imags = array.array('i')
    else:
        reals = array.array('f')
        imags = array.array('f')
    reals.fromfile(f_data, numSamples)
    imags.fromfile(f_data, numSamples)
    reals = np.float32(reals)
    imags = np.float32(imags)

    
    # dataMatrix2[row] = [complex(r,i) for r,i in zip(reals, imags)]
    dataMatrix2[row] = reals + 1j*imags
    realMatrix2[row] = reals

    del reals, imags
if window:
    dataMatrix2 = dataMatrix2 * windowArray
    realMatrix2 = realMatrix2 * windowArray

endTime = time.time()
print(f"\tTotal time was {endTime - startTime} s")



print(f"Calculating FFT for {f_data_name}...")
startTime = time.time()

for i in range(len(dataMatrix2)):
    sr = sampleRate
    
    X = fftshift(fft(dataMatrix2[i]))
    # N = len(X)
    # k = np.arange(N)
    # T = N/sr
    # freq = k/T

    X = X[numSamples//2 - fftLength//2:numSamples//2 + fftLength//2]

    FFTMagMatrix2[i] = np.abs(X)
    FFTMagMatrix2[i] = 10*np.log10(FFTMagMatrix2[i])
    if plotPhase:
        FFTPhaseMatrix2[i] = np.angle(X)
        # FFTPhaseMatrix2[i] = np.angle(dataMatrix2[i])
        for j in range(len(FFTPhaseMatrix2[i])):
            if FFTPhaseMatrix2[i][j] < 0:
                FFTPhaseMatrix2[i][j] = 2*np.pi + FFTPhaseMatrix2[i][j]

del dataMatrix2

endTime = time.time()
print(f"\tTotal time was {endTime - startTime} s")



del windowArray

# Interpolation - multiply length of matricies by 10x
if interpolate:
    extDataMatrix = np.zeros([numChirps, 100000], complex)
    extDataMatrix[:, 40000:60000] = dataMatrix
    dataMatrix = extDataMatrix
    del extDataMatrix

    extDataMatrix2 = np.zeros([numChirps, 100000], complex)
    extDataMatrix2[:, 40000:60000] = dataMatrix2
    dataMatrix2 = extDataMatrix2
    del extDataMatrix2

    FFTMagMatrix = np.empty([numChirps, 100000], float)
    FFTMagMatrix2 = np.empty([numChirps, 100000], float)
    FFTMagMatrixDiff = np.empty([numChirps, 100000], float)

if plotPhase:
    FFTPhaseMatrix = np.empty([numChirps, len(FFTMagMatrix[0])], float)
    FFTPhaseMatrix2 = np.empty([numChirps, len(FFTMagMatrix[0])], float)
    FFTPhaseMatrixDiff = np.empty([numChirps, len(FFTMagMatrix[0])], float)



print(f"Calculating average FFT of {f_nulldata_name}...")
startTime = time.time()

# Average the targetless chirps in FFT
for i in range(fftLength):
    FFTMagMatrix[:,i] = np.average(FFTMagMatrix[:,i])

endTime = time.time()
print(f"\tTotal time was {endTime - startTime} s")



print(f"Calculating difference matrix...")
startTime = time.time()

realMatrixDiff = realMatrix2 - realMatrix
FFTMagMatrixDiff = FFTMagMatrix2 - FFTMagMatrix
if plotPhase:
    FFTPhaseMatrixDiff = FFTPhaseMatrix2 - FFTPhaseMatrix
    for i in range(len(FFTPhaseMatrixDiff)):
        for j in range(len(FFTPhaseMatrixDiff[i])):
            if FFTPhaseMatrixDiff[i][j] < 0:
                FFTPhaseMatrixDiff[i][j] = 2*np.pi + FFTPhaseMatrixDiff[i][j]

endTime = time.time()
print(f"\tTotal time was {endTime - startTime} s")



print(f"Calculating distance...")
startTime = time.time()

# Calculate distance
c = 2*10**8
T = numSamples/sampleRate
k = np.arange(numSamples)
freq = k/T
# freq = 6600
chirpsPerSec = sampleRate/numSamples
slope = bandwidth*chirpsPerSec
t = freq/slope
distance = c*t
distance = distance[numSamples//2 - fftLength//2:numSamples//2 + fftLength//2]
distance = distance-distance[fftLength//2]
distance = np.round(distance)
distance = distance.astype(int)

endTime = time.time()
print(f"\tTotal time was {endTime - startTime} s")



print("Plotting...")

vmaxT = '40'
vminT = '-40'
vmaxF = '50'
vminF = '-10'
# xlim = [49975,50025]
# xlim = [499950,500050]
# xlim = [0,numSamples]
# xlim = [numSamples//2 - fftLength//2, numSamples//2 + fftLength//2]

# if interpolate:
#     xlim = [9750,10250]

subplot1 = 220
subplot2 = 210
if plotPhase:
    subplot1 = 320
    subplot2 = 310

plt.figure(1)

# Plot the two received signals (before subtraction)
plt.subplot(subplot1+1)
# plt.subplot(221)
plt.title('Dataset 1')
plot = plt.matshow(realMatrix, cmap='hot', aspect="auto", fignum=False)
# plot = plt.matshow(realMatrix, cmap='hot', vmax=vmaxT, vmin=vminT, aspect="auto", fignum=False)
plt.colorbar(plot)
plt.xlabel('Sample Number')
plt.ylabel('Chirp Number')

plt.subplot(subplot1+2)
# plt.subplot(222)
plt.title('Dataset 2')
plot = plt.matshow(realMatrix2, cmap='hot', aspect="auto", fignum=False)
# plot = plt.matshow(realMatrix2, cmap='hot', vmax=vmaxT, vmin=vminT, aspect="auto", fignum=False)
plt.colorbar(plot)
plt.xlabel('Sample Number')
plt.ylabel('Chirp Number')

plt.subplot(subplot1+3)
# plt.subplot(223)
plot = plt.matshow(FFTMagMatrix, cmap='hot', aspect="auto", fignum=False)
# plot = plt.matshow(FFTMagMatrix, cmap='hot', vmax=vmaxF, vmin=vminF, aspect="auto", fignum=False)
plt.colorbar(plot)
plt.xlabel('Frequency Bin')
plt.ylabel('Chirp Number')
# plt.xlim(xlim)

plt.subplot(subplot1+4)
# plt.subplot(224)
plot = plt.matshow(FFTMagMatrix2, cmap='hot', aspect="auto", fignum=False)
# plot = plt.matshow(FFTMagMatrix2, cmap='hot', vmax=vmaxF, vmin=vminF, aspect="auto", fignum=False)
plt.colorbar(plot)
plt.xlabel('Frequency Bin')
plt.ylabel('Chirp Number')
# plt.xlim(xlim)

if plotPhase:
    plt.subplot(subplot1+5)
    # plot = plt.matshow(FFTMagMatrix, vmax='4', vmin='-4', aspect="auto", fignum=False)
    plot = plt.matshow(FFTPhaseMatrix, cmap='twilight', aspect="auto", fignum=False)
    plt.colorbar(plot)
    plt.xlabel('Frequency Bin')
    plt.ylabel('Chirp Number')
    # plt.xlim(xlim)

    plt.subplot(subplot1+6)
    # plot = plt.matshow(FFTMagMatrix, vmax='4', vmin='-4', aspect="auto", fignum=False)
    plot = plt.matshow(FFTPhaseMatrix2, cmap='twilight', aspect="auto", fignum=False)
    plt.colorbar(plot)
    plt.xlabel('Frequency Bin')
    plt.ylabel('Chirp Number')
    # plt.xlim(xlim)

# Free up memory
del FFTMagMatrix, FFTMagMatrix2
del realMatrix, realMatrix2

# Plot the signals after subtraction
plt.figure(2)
plt.subplot(subplot2+1)
# plt.subplot(211)
plt.title('Difference')
plot = plt.matshow(realMatrixDiff, cmap='hot', aspect="auto", fignum=False)
# plot = plt.matshow(realMatrixDiff, cmap='hot', vmax=vmaxT, vmin=vminT, aspect="auto", fignum=False)
plt.colorbar(plot)
plt.xlabel('Sample Number')
plt.ylabel('Chirp Number')

plt.subplot(subplot2+2)
# plt.subplot(212)
plot = plt.matshow(FFTMagMatrixDiff, cmap='hot', aspect="auto", fignum=False)
# plot = plt.matshow(FFTMagMatrixDiff, cmap='hot', vmax=vmaxF, vmin=vminF, aspect="auto", fignum=False)
plt.colorbar(plot)
# plt.xlabel('Frequency Bin')
plt.xlabel('Distance (m)')
plt.ylabel('Chirp Number')
# plt.xlim(xlim)
# plt.xticks(np.arange(0, fftLength), distance)
# h2 = plt.twiny()
# h2.set_xticks(np.arange(0, fftLength), distance)
# h2.set_xlabel('Distance (m)')

if plotPhase:
    plt.subplot(subplot2+3)
    # plot = plt.matshow(FFTMagMatrix, vmax='4', vmin='-4', aspect="auto", fignum=False)
    plot = plt.matshow(FFTPhaseMatrixDiff, cmap='twilight', aspect="auto", fignum=False)
    plt.colorbar(plot)
    plt.xlabel('Frequency Bin')
    plt.ylabel('Chirp Number')
    # plt.xlim(xlim)


# Free up memory
del FFTMagMatrixDiff
del realMatrixDiff




plt.show()

f_config.close()
f_nulldata.close()
f_data.close()
print("Done")