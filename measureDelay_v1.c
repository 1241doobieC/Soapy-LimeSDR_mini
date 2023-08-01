/*
Written by Levi Powell, September 2022
This code implements a basic radar using the LimeSDR Mini SDR.

Compile with: gcc -std=c99 measureDelay.c -lSoapySDR -lm -o measureDelay.out
Compile with memory checking: gcc -std=c99 measureDelay.c -fsanitize=address -fno-omit-frame-pointer -lSoapySDR -lm -o measureDelay.out
Compile with speed optimization: gcc -std=c99 measureDelay.c -lSoapySDR -lm -Ofast -o measureDelay.out
*/


#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <stdlib.h>
#include <stdint.h>
#include <unistd.h>
#include <stdint.h>

#include <SoapySDR/Device.h>
#include <SoapySDR/Formats.h>
#include <SoapySDR/Constants.h>

#define SAVE_TO_FILE true
#define ZERO_DELAY_SAMPLES 96

#define BUFF_TX_LENGTH (size_t) 1e3
#define BUFF_RX_LENGTH (size_t) 1e3

// #define FREQUENCY 222.500e6 //Hz
#define FREQUENCY 0.5e9 //Hz
#define CLOCK_RATE 100e6 //Hz
#define SAMPLE_RATE_TX 30e6 //Hz
#define SAMPLE_RATE_RX 30e6 //Hz
#define BANDWIDTH_TX 50e6
#define BANDWIDTH_RX 50e6
#define CHANNEL_TX 0
#define CHANNEL_RX 0

#define STREAM_DELAY_TX (long long) 0.25e9 //nanoseconds
#define STREAM_DELAY_RX (long long) 0.25e9 //nanoseconds
#define STREAM_TIMEOUT_TX 3e6 //microseconds
#define STREAM_TIMEOUT_RX 3e6 //microseconds
#define FLAGS_TX SOAPY_SDR_HAS_TIME | SOAPY_SDR_END_BURST
#define FLAGS_RX SOAPY_SDR_HAS_TIME | SOAPY_SDR_END_BURST

#define CONTIGUOUS_BUFF_TX_LENGTH (int) 1e5
#define CONTIGUOUS_BUFF_RX_LENGTH (int) 1e5

#define CHIRP_BANDWIDTH SAMPLE_RATE_RX
#define NUM_CHIRPS 10
#define CHIRP_DELAY (long long) 0.3e9 //nanoseconds

#define PI 3.1415926535

//Define complex integer data type
//(C includes complex floats, but not 16-bit integers)
struct cs16_struct
{
    int16_t real;
    int16_t imag;
};
typedef struct cs16_struct cs16;

//Function declarations
struct SoapySDRDevice* Setup(void);
void DeviceInfo(struct SoapySDRDevice* sdr);
void SetParameters(SoapySDRDevice* sdr);
SoapySDRStream* MakeStream(SoapySDRDevice* sdr, const int direction);
// void SaveData(FILE* fp, int* sampleNumber, const complex float* buffer, const int length, const bool saveToFile);
void SaveData(FILE* fp, int* sampleNumber, const cs16* buffer, const int length, const bool saveToFile);
// void FillBuffer(complex float* buff, int length);
void FillBuffer(cs16* buff, int length);
// void TransmitReceive(SoapySDRDevice* sdr, SoapySDRStream* txStream, SoapySDRStream* rxStream, complex float* bufferTx, complex float* bufferRx, complex float* contBufferTx, complex float* contBufferRx, long long transmitTime, long long receiveTime, int* firstSampleIndex);
void TransmitReceive(SoapySDRDevice* sdr, SoapySDRStream* txStream, SoapySDRStream* rxStream, cs16* bufferTx, cs16* bufferRx, cs16* contBufferTx, cs16* contBufferRx, long long transmitTime, long long receiveTime, int* firstSampleIndex);
// void TrimBuffer(complex float* contBufferRx, int firstSampleIndex);
void TrimBuffer(cs16* contBufferRx, int firstSampleIndex);
// void MixSignals(complex float* contBufferTx, complex float* contBufferRx);
void MixSignals(cs16* contBufferTx, cs16* contBufferRx);


int main()
{
    //Setup data and configuration files
    FILE* fp_config;
    FILE* fp_data;
    FILE* fp_nulldata;
    fp_config = fopen("config.txt", "w");
    fp_data = fopen("data.bin", "wb");
    fp_nulldata = fopen("nulldata.bin", "wb");
    int sampleNumber = 0;

    //Load constants into configuration file
    if(fp_config == NULL)
    {
        printf("The configuration file cannot be opened.\n");
    }
    else
    {
        fprintf(fp_config, "%d,%d,%d,%d\n", (int)SAMPLE_RATE_RX, (int)CHIRP_BANDWIDTH, (int)CONTIGUOUS_BUFF_RX_LENGTH, (int)NUM_CHIRPS);
    }

    //Setup SDR
    struct SoapySDRDevice* sdr = Setup();

    //Display SDR information
    // DeviceInfo(sdr);

    //Set SDR parameters
    SetParameters(sdr);

    //A buffer to hold the data to be transmitted
    // complex float contBufferTx[CONTIGUOUS_BUFF_TX_LENGTH];
    // complex float* contBufferTx = malloc(8*CONTIGUOUS_BUFF_TX_LENGTH);
    // complex float* bufferTx = contBufferTx;
    cs16* contBufferTx = malloc(4*CONTIGUOUS_BUFF_TX_LENGTH);
    cs16* bufferTx = contBufferTx;
    printf("Length of Tx buffer: %d\n", CONTIGUOUS_BUFF_TX_LENGTH);

    //A buffer to hold the data to be received
    // complex float contBufferRx[CONTIGUOUS_BUFF_RX_LENGTH];
    // complex float* contBufferRx = malloc(8*CONTIGUOUS_BUFF_RX_LENGTH);
    // complex float* bufferRx = contBufferRx;
    cs16* contBufferRx = malloc(4*CONTIGUOUS_BUFF_RX_LENGTH);
    cs16* bufferRx = contBufferRx;
    printf("Length of Rx buffer: %d\n", CONTIGUOUS_BUFF_RX_LENGTH);
    //Initialize to -1 (easy to detect errors)
    for(int i = 0; i < CONTIGUOUS_BUFF_RX_LENGTH; i++)
    {
        // contBufferRx[i] = -1 - 1*I;
        contBufferRx[i].real = -1;
        contBufferRx[i].imag = -1;
    }

    //Create rx and tx streams
    SoapySDRStream* txStream = MakeStream(sdr, SOAPY_SDR_TX);
    SoapySDRStream* rxStream = MakeStream(sdr, SOAPY_SDR_RX);
    printf("Created RX and TX streams.\n");

    //Fill Tx buffer with data to transmit
    FillBuffer(contBufferTx, CONTIGUOUS_BUFF_TX_LENGTH);
    
    //Prapare timed streaming
    const char* timeSource = SoapySDRDevice_getTimeSource(sdr);
    long long currentHardwareTime;
    long long transmitTime;
    long long receiveTime;

    //Print start time
    currentHardwareTime = SoapySDRDevice_getHardwareTime(sdr, timeSource);
    printf("Current hardware time is %lf s.\n", currentHardwareTime / 1.0e9);

    //Store the index of the first valid sample received
    int firstSampleIndex = 0;

    // sleep(1);

    //Transmit and receive
    for(int k = 0; k < NUM_CHIRPS; k++)
    {
        printf("\nReady to transmit/receive chirp %d. Press any key to continue...\n", k);
        // getchar();

        currentHardwareTime = SoapySDRDevice_getHardwareTime(sdr, timeSource);
        printf("Started chirp code at time %lf s\n", currentHardwareTime / 1.0e9);

        transmitTime = CHIRP_DELAY * (k+1);
        receiveTime = CHIRP_DELAY * (k+1);
        
        TransmitReceive(sdr, txStream, rxStream, bufferTx, bufferRx, contBufferTx, contBufferRx, transmitTime, receiveTime, &firstSampleIndex);

        //remove garbage values from the beginning of the contiguous buffer
        long long startCleanTime = SoapySDRDevice_getHardwareTime(sdr, timeSource);
        TrimBuffer(contBufferRx, firstSampleIndex);
        long long endCleanTime = SoapySDRDevice_getHardwareTime(sdr, timeSource);

        //Mix Tx/Rx signals
        long long startMixTime = SoapySDRDevice_getHardwareTime(sdr, timeSource);
        // MixSignals(contBufferTx, contBufferRx);
        long long endMixTime = SoapySDRDevice_getHardwareTime(sdr, timeSource);

        //Save data to file
        long long startSaveTime = SoapySDRDevice_getHardwareTime(sdr, timeSource);
        SaveData(fp_nulldata, &sampleNumber, contBufferRx, CONTIGUOUS_BUFF_RX_LENGTH, SAVE_TO_FILE);
        long long endSaveTime = SoapySDRDevice_getHardwareTime(sdr, timeSource);

        bufferTx = contBufferTx;
        bufferRx = contBufferRx;

        // printf("Chirp Rx timestamp was %lf s\n", firstRxTimestamp / 1.0e9);
        printf("\tTime to clean signal: %lf s\n", (endCleanTime - startCleanTime) / 1.0e9);
        printf("\tTime to mix signals: %lf s\n", (endMixTime - startMixTime) / 1.0e9);
        printf("\tTime to save data to file: %lf s\n", (endSaveTime - startSaveTime) / 1.0e9);
        
        currentHardwareTime = SoapySDRDevice_getHardwareTime(sdr, timeSource);
        printf("Finished chirp code at time %lf s\n", currentHardwareTime / 1.0e9);
    }


    printf("\nCollected null data. Ready to collect target data. Press any key to continue...\n");
    getchar();

    long long currentHardwareTime2 = SoapySDRDevice_getHardwareTime(sdr, timeSource);

    //Transmit and receive
    for(int k = 0; k < NUM_CHIRPS; k++)
    {
        printf("\nReady to transmit/receive chirp %d. Press any key to continue...\n", k);
        // getchar();

        currentHardwareTime = SoapySDRDevice_getHardwareTime(sdr, timeSource);
        printf("Started chirp code at time %lf s\n", currentHardwareTime / 1.0e9);

        transmitTime = CHIRP_DELAY * (k+1) + currentHardwareTime2;
        receiveTime = CHIRP_DELAY * (k+1) + currentHardwareTime2;
        
        TransmitReceive(sdr, txStream, rxStream, bufferTx, bufferRx, contBufferTx, contBufferRx, transmitTime, receiveTime, &firstSampleIndex);

        //remove garbage values from the beginning of the contiguous buffer
        long long startCleanTime = SoapySDRDevice_getHardwareTime(sdr, timeSource);
        TrimBuffer(contBufferRx, firstSampleIndex);
        long long endCleanTime = SoapySDRDevice_getHardwareTime(sdr, timeSource);

        //Mix Tx/Rx signals
        long long startMixTime = SoapySDRDevice_getHardwareTime(sdr, timeSource);
        // MixSignals(contBufferTx, contBufferRx);
        long long endMixTime = SoapySDRDevice_getHardwareTime(sdr, timeSource);

        //Save data to file
        long long startSaveTime = SoapySDRDevice_getHardwareTime(sdr, timeSource);
        SaveData(fp_data, &sampleNumber, contBufferRx, CONTIGUOUS_BUFF_RX_LENGTH, SAVE_TO_FILE);
        long long endSaveTime = SoapySDRDevice_getHardwareTime(sdr, timeSource);

        bufferTx = contBufferTx;
        bufferRx = contBufferRx;

        // printf("Chirp Rx timestamp was %lf s\n", firstRxTimestamp / 1.0e9);
        printf("\tTime to clean signal: %lf s\n", (endCleanTime - startCleanTime) / 1.0e9);
        printf("\tTime to mix signals: %lf s\n", (endMixTime - startMixTime) / 1.0e9);
        printf("\tTime to save data to file: %lf s\n", (endSaveTime - startSaveTime) / 1.0e9);
        
        currentHardwareTime = SoapySDRDevice_getHardwareTime(sdr, timeSource);
        printf("Finished chirp code at time %lf s\n", currentHardwareTime / 1.0e9);
    }

    //Clean up memory
    printf("\nFreeing memory...\n");
    free(contBufferTx);
    free(contBufferRx);
    
    printf("Closing streams...\n");
    SoapySDRDevice_deactivateStream(sdr, txStream, 0, 0);
    SoapySDRDevice_deactivateStream(sdr, rxStream, 0, 0);
    SoapySDRDevice_closeStream(sdr, txStream);
    SoapySDRDevice_closeStream(sdr, rxStream);
    SoapySDRDevice_unmake(sdr);

    printf("Closing files...\n");
    fclose(fp_config);
    fclose(fp_data);
    
    printf("DONE\n");
    return 0;
}

struct SoapySDRDevice* Setup(void)
{
    //enumerate devices
    size_t length;
    SoapySDRKwargs* results = SoapySDRDevice_enumerate(NULL, &length);
    for (size_t i = 0; i < length; i++)
    {
        printf("[Setup] Found device #%d: ", (int)i);
        for (size_t j = 0; j < results[i].size; j++)
        {
            printf("%s=%s, ", results[i].keys[j], results[i].vals[j]);
        }
        printf("\n");
    }
    SoapySDRKwargsList_clear(results, length);
    //create device instance
    //args can be user defined or from the enumeration result
    SoapySDRKwargs args = {};
    SoapySDRKwargs_set(&args, "driver", "lime");
    SoapySDRDevice* sdr = SoapySDRDevice_make(&args);
    SoapySDRKwargs_clear(&args);
    if (sdr == NULL)
    {
        printf("[Setup] SoapySDRDevice_make fail: %s\n", SoapySDRDevice_lastError());
        //return EXIT_FAILURE;
        return NULL;
    }

    return sdr;
}

void DeviceInfo(struct SoapySDRDevice* sdr)
{
    //query device info
    size_t length;
    char** names = SoapySDRDevice_listAntennas(sdr, SOAPY_SDR_RX, 0, &length);
    printf("[DeviceInfo] Rx antennas: ");
    for (size_t i = 0; i < length; i++) printf("%s, ", names[i]);
    printf("\n");
    SoapySDRStrings_clear(&names, length);
    names = SoapySDRDevice_listGains(sdr, SOAPY_SDR_RX, 0, &length);
    printf("[DeviceInfo] Rx gains: ");
    for (size_t i = 0; i < length; i++) printf("%s, ", names[i]);
    printf("\n");
    SoapySDRStrings_clear(&names, length);
    SoapySDRRange* ranges = SoapySDRDevice_getFrequencyRange(sdr, SOAPY_SDR_RX, 0, &length);
    printf("[DeviceInfo] Rx freq ranges: ");
    for (size_t i = 0; i < length; i++) printf("[%g Hz -> %g Hz], ", ranges[i].minimum, ranges[i].maximum);
    printf("\n");
    free(ranges);
}

void SetParameters(SoapySDRDevice* sdr)
{
    //Set clock rate
    if(SoapySDRDevice_hasHardwareTime(sdr, NULL))
    {
        SoapySDRDevice_setMasterClockRate(sdr, CLOCK_RATE);
        printf("[SetParameters] Master Clock rate was set to: %e Hz\n", CLOCK_RATE);
        printf("[SetParameters] Master Clock rate was read back as: %e Hz\n", SoapySDRDevice_getMasterClockRate(sdr));
    }
    else
    {
        printf("[SetParameters] This device does not support timed streaming.\n"); 
    }

    //Setup Tx parameters
    if (SoapySDRDevice_setSampleRate(sdr, SOAPY_SDR_TX, 0, SAMPLE_RATE_TX) != 0)
    {
        printf("[SetParameters] Tx setSampleRate fail: %s\n", SoapySDRDevice_lastError());
    }
    if (SoapySDRDevice_setAntenna(sdr, SOAPY_SDR_TX, 0, "Auto") != 0)
    {
        printf("[SetParameters] Tx setAntenna fail: %s\n", SoapySDRDevice_lastError());
    }
    if (SoapySDRDevice_setGain(sdr, SOAPY_SDR_TX, 0, 40) != 0)
    // if(SoapySDRDevice_setGainElement(sdr, SOAPY_SDR_TX, 0, "PGA", 40) != 0)
    {
        printf("[SetParameters] Tx setGain fail: %s\n", SoapySDRDevice_lastError());
    }
    if (SoapySDRDevice_setFrequency(sdr, SOAPY_SDR_TX, 0, FREQUENCY, NULL) != 0)
    {
        printf("[SetParameters] Tx setFrequency fail: %s\n", SoapySDRDevice_lastError());
    }
    if (SoapySDRDevice_setBandwidth(sdr, SOAPY_SDR_TX, CHANNEL_TX, BANDWIDTH_TX) != 0)
    {
        printf("[SetParameters] Tx setBandwidth fail: %s\n", SoapySDRDevice_lastError());
    }

    //Setup Rx parameters
    if (SoapySDRDevice_setSampleRate(sdr, SOAPY_SDR_RX, 0, SAMPLE_RATE_RX) != 0)
    {
        printf("[SetParameters] Rx setSampleRate fail: %s\n", SoapySDRDevice_lastError());
    }
    if (SoapySDRDevice_setAntenna(sdr, SOAPY_SDR_RX, 0, "Auto") != 0)
    {
        printf("[SetParameters] Rx setAntenna fail: %s\n", SoapySDRDevice_lastError());
    }
    if (SoapySDRDevice_setGain(sdr, SOAPY_SDR_RX, 0, 40) != 0)
    // if(SoapySDRDevice_setGainElement(sdr, SOAPY_SDR_RX, 0, "PGA", 40) != 0)
    {
        printf("[SetParameters] Rx setGain fail: %s\n", SoapySDRDevice_lastError());
    }
    if (SoapySDRDevice_setFrequency(sdr, SOAPY_SDR_RX, 0, FREQUENCY, NULL) != 0)
    {
        printf("[SetParameters] Rx setFrequency fail: %s\n", SoapySDRDevice_lastError());
    }
    if (SoapySDRDevice_setBandwidth(sdr, SOAPY_SDR_RX, CHANNEL_RX, BANDWIDTH_RX) != 0)
    {
        printf("[SetParameters] Rx setBandwidth fail: %s\n", SoapySDRDevice_lastError());
    }
}

SoapySDRStream* MakeStream(SoapySDRDevice* sdr, const int direction)
{
    SoapySDRKwargs args = SoapySDRKwargs_fromString("WIRE=CS16"); // 設定設備與主機之間的Sample Format為複數int16

    //Create and activate a stream
    // SoapySDRStream* stream = SoapySDRDevice_setupStream(sdr, direction, SOAPY_SDR_CF32, NULL, 0, NULL);
    SoapySDRStream* stream = SoapySDRDevice_setupStream(sdr, direction, SOAPY_SDR_CS16, NULL, 0, &args);
    if (stream == NULL)
    {
        printf("[MakeStream] Setup stream fail: %s\n", SoapySDRDevice_lastError());
    }
    // if(SoapySDRDevice_activateStream(sdr, stream, 0, 0, 0) != 0)
    // {
    //     printf("[MakeStream] Activate stream fail: %s\n", SoapySDRDevice_lastError());
    // }
    SoapySDRKwargs_clear(&args); // 釋放args的資源
    return stream;
}

// void SaveData(FILE* fp, int* sampleNumber, const complex float* buffer, const int length, const bool saveToFile)
void SaveData(FILE* fp, int* sampleNumber, const cs16* buffer, const int length, const bool saveToFile)
{
    if(saveToFile)
    {
        printf("[SaveData] Saving data to file...\n");
        
        if(fp == NULL)
        {
            printf("[SaveData] The file cannot be opened.\n");
        }
        // else
        // {
        //     float real;
        //     float imag;
        //     for(int i = 0; i < length; i++)
        //     {
        //         // fprintf(fp, "%d,%f,%f\n", *sampleNumber, crealf(buffer[i]), cimagf(buffer[i]));
        //         // (*sampleNumber)++;
                
        //         // fprintf(fp, "%+f%+fj,", creal(buffer[i]), cimag(buffer[i]));
        //         real = creal(buffer[i]);
        //         imag = cimag(buffer[i]);
        //         fwrite(&real, sizeof(real), 1, fp);
        //         fwrite(&imag, sizeof(imag), 1, fp);
        //     }
        //     // fprintf(fp, "\n");
        // }
        else
        {
            // float real[length];
            // float imag[length];
            // float* real = malloc(length * sizeof(float));
            // float* imag = malloc(length * sizeof(float));
            int16_t* real = malloc(length * sizeof(int16_t));
            int16_t* imag = malloc(length * sizeof(int16_t));
            if(real == NULL || imag == NULL) printf("[SaveData] Failed to allocate memory!\n");
            for(int i = 0; i < length; i++)
            {
                // real[i] = creal(buffer[i]);
                // imag[i] = cimag(buffer[i]);
                real[i] = buffer[i].real;
                // real[i] = 0;
                /**
                 * @test  原本長如下:
                 * @brief real[i] = buffer[i].real;
                 * @brief imag[i] = buffer[i].imag;
                 * @date 07.31.2023 16:10
                */
                imag[i] = buffer[i].imag;
                // imag[i] = 0;
            }
            // fwrite(&real, sizeof(float)*length, 1, fp);
            // fwrite(&imag, sizeof(float)*length, 1, fp);
            // size_t writtenReal = fwrite(real, sizeof(float), length, fp);
            // size_t writtenImag = fwrite(imag, sizeof(float), length, fp);
            size_t writtenReal = fwrite(real, sizeof(int16_t), length, fp);
            size_t writtenImag = fwrite(imag, sizeof(int16_t), length, fp);

            if(writtenReal != length || writtenImag != length)
            {
                printf("[SaveData] Failed to write data to file!\n");
            }
            free(real);
            free(imag);
        }
    }
    else
    {
        printf("[SaveData] Printing data to screen...\n");
        for(int i = 0; i < length; i++)
        {
            // printf("%d: [%f, %fi]\n", *sampleNumber, crealf(buffer[i]), cimagf(buffer[i]));
            printf("%d: [%d, %di]\n", *sampleNumber, buffer[i].real, buffer[i].imag);
            (*sampleNumber)++;
        }
    }
}

// void FillBuffer(complex float* buff, int length)
void FillBuffer(cs16* buff, int length)
{
    float amplitude = 2046; //Center to peak
    int period = 100; //Number of samples
    for(int i = 0; i < length; i++)
    {
        //Square wave
        // buff[i] = (i % period < period/2) ? amplitude*(1 + 1*I) : amplitude*(-1 + -1*I);

        //Continuous wave
        // buff[i] = amplitude*(1 + 1*I);

        //Sawtooth wave
        // buff[i] = amplitude*(2*(i % period) / (float)period - 1)*(1 + 1*I);

        //Sine wave
        // buff[i] = amplitude*cos(2*PI*i/period) + amplitude*sin(2*PI*i/period)*I;
        // buff[i].real = amplitude*cos(2*PI*i/period);
        // buff[i].imag = amplitude*sin(2*PI*i/period);

        //Chirp
        // int chirpsPerSec = 10;
        float chirpsPerSec = ((float) SAMPLE_RATE_TX) / CONTIGUOUS_BUFF_TX_LENGTH; //Enough to have one chirp in the contiguous buffer.
        int chirpStartFreq = 0;
        int chirpEndFreq = CHIRP_BANDWIDTH;
        double time = i / SAMPLE_RATE_TX;
        double chirpSlope = (chirpEndFreq - chirpStartFreq) * chirpsPerSec;
        const double angle = (2*PI*time)*(chirpStartFreq + time*chirpSlope/2);
        // buff[i] = cos(angle) + sin(angle)*I;
        buff[i].real = cos(angle);
        buff[i].imag = sin(angle);

        //Impulse
        // if(i == length / 2)
        // {
        //     buff[i] = 1+1*I;
        // }
        // else
        // {
        //     buff[i] = 0+0*I;
        // }
    }
}

// void TransmitReceive(SoapySDRDevice* sdr, SoapySDRStream* txStream, SoapySDRStream* rxStream, complex float* bufferTx, complex float* bufferRx, complex float* contBufferTx, complex float* contBufferRx, long long transmitTime, long long receiveTime, int* firstSampleIndex)
void TransmitReceive(SoapySDRDevice* sdr, SoapySDRStream* txStream, SoapySDRStream* rxStream, cs16* bufferTx, cs16* bufferRx, cs16* contBufferTx, cs16* contBufferRx, long long transmitTime, long long receiveTime, int* firstSampleIndex)
{
    const void* buffsTx[] = {bufferTx};
    void* buffsRx[] = {bufferRx};

    long long timestampRx;

    printf("[TransmitReceive] Scheduled to transmit/receive a chirp at %lf s\n", transmitTime / 1.0e9);

    // double firstRxTimestamp = -1;

    //Activate streams
    int flags = FLAGS_TX;
    if(SoapySDRDevice_activateStream(sdr, txStream, 0, 0, 0) != 0)
    {
        printf("[TransmitReceive] Activate Tx stream fail: %s\n", SoapySDRDevice_lastError());
    }
    if(SoapySDRDevice_activateStream(sdr, rxStream, flags, receiveTime, 0) != 0)
    {
        printf("[TransmitReceive] Activate Rx stream fail: %s\n", SoapySDRDevice_lastError());
    }

    //Begin transmitting
    int txStreamStatus;
    for(int i = 0; i < CONTIGUOUS_BUFF_TX_LENGTH / BUFF_TX_LENGTH; i++)
    {
        buffsTx[0] = bufferTx;
        txStreamStatus = SoapySDRDevice_writeStream(sdr, txStream, buffsTx, BUFF_TX_LENGTH, &flags, transmitTime, STREAM_TIMEOUT_TX);
        transmitTime += BUFF_TX_LENGTH / SAMPLE_RATE_TX * 1e9; //Begin next transmission immediately after this one finishes.
        bufferTx += BUFF_TX_LENGTH; //Move buffer pointer to the next section of the contiguous buffer to be transmitted.

        if(txStreamStatus != BUFF_TX_LENGTH)
        {
            printf("[TransmitReceive] Write stream failed: %d\n", txStreamStatus);
        }
    }
    
    //Begin receiving
    int rxStreamStatus;
    for(int i = 0; i < CONTIGUOUS_BUFF_RX_LENGTH / BUFF_RX_LENGTH; i++)
    {
        buffsRx[0] = bufferRx;
        rxStreamStatus = SoapySDRDevice_readStream(sdr, rxStream, buffsRx, BUFF_RX_LENGTH, &flags, &timestampRx, STREAM_TIMEOUT_RX);
        bufferRx += BUFF_RX_LENGTH; //Move buffer pointer to the next section of the contiguous buffer to be received.
        
        if(rxStreamStatus != BUFF_RX_LENGTH)
        {
            printf("[TransmitReceive] Read stream failed. Samples captured: %d\n", rxStreamStatus);
            *firstSampleIndex = BUFF_RX_LENGTH - rxStreamStatus + ZERO_DELAY_SAMPLES;
        }

        // if(i == 0)
        // {
        //     firstRxTimestamp = timestampRx;
        // }

        if(rxStreamStatus < 0)  // Detect a failed chirp
        {
            printf("[TransmitReceive] Detected failed chirp! Filling buffer with zersos.\n");
            for(int j = 0; j < CONTIGUOUS_BUFF_RX_LENGTH; j++)
            {
                // contBufferRx[j] = 0*(1+I);
                contBufferRx[j].real = 0;
                contBufferRx[j].imag = 0;
            }
            break;  // Skip this chirp, move on to the next one.
        }
    }
}

// void TrimBuffer(complex float* contBufferRx, int firstSampleIndex)
void TrimBuffer(cs16* contBufferRx, int firstSampleIndex)
{
    // complex float* contBufferRxCleaned = malloc(8*CONTIGUOUS_BUFF_RX_LENGTH);
    for(int i = 0; i < CONTIGUOUS_BUFF_RX_LENGTH; i++)
    {
        if(i < CONTIGUOUS_BUFF_RX_LENGTH-firstSampleIndex)
        {
            contBufferRx[i] = contBufferRx[i+firstSampleIndex];
        }
        else
        {
            // contBufferRx[i] = 0;
            contBufferRx[i].real = 0;
            contBufferRx[i].imag = 0;
        }
    }
}

// void MixSignals(complex float* contBufferTx, complex float* contBufferRx)
void MixSignals(cs16* contBufferTx, cs16* contBufferRx)
{
    // complex float* mixedSignal = malloc(8*CONTIGUOUS_BUFF_RX_LENGTH);
    for(int i = 0; i < CONTIGUOUS_BUFF_RX_LENGTH; i++)
    {
        // contBufferRx[i] = (creal(contBufferTx[i]) + cimag(contBufferTx[i])*I) * (creal(contBufferRx[i]) - cimag(contBufferRx[i])*I);

    }
}
