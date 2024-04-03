# skrf_tools
tools based on skrf to help with RF-design and s-parameter files.

## filter_identification.py
Will help you identify your s-parameter files.
```
$ ./filter_identification.py bandpassfilterbank/filter*.s2p
filter_0
	 Isolation -79.69[dB]
filter_1
	 Bandpass filter. Center frequency 350.00[MHz] bandwidth 200.00[MHz]
filter_2
	 Bandpass filter. Center frequency 848.50[MHz] bandwidth 503.00[MHz]
filter_3
	 Bandpass filter. Center frequency 1295.50[MHz] bandwidth 551.00[MHz]
filter_4
	 Bandpass filter. Center frequency 2044.50[MHz] bandwidth 1151.00[MHz]
filter_5
	 Bandpass filter. Center frequency 3312.00[MHz] bandwidth 1674.00[MHz]
filter_6
	 Lowpass filter. 3dB frequency 2108.00[MHz]
filter_7
	 Bandpass filter. Center frequency 4726.00[MHz] bandwidth 3212.00[MHz]
filter_8
	 Bandpass filter. Center frequency 12581.50[MHz] bandwidth 13003.00[MHz]
filter_9
	 Through, loss -3.85[dB]
```
Also contains useful functions that can identify pass bands in S-parameter files.
Help:
```
usage: filter_identification.py [-h] [--passband-threshold PASSBAND_THRESHOLD]
                                [--isolation-threshold ISOLATION_THRESHOLD]
                                [--merge-distance MERGE_DISTANCE] [--plot]
                                [touchstone [touchstone ...]]

Automatic identification of filters in s-parameters

positional arguments:
  touchstone            S2P files to identify

optional arguments:
  -h, --help            show this help message and exit
  --passband-threshold PASSBAND_THRESHOLD, -pt PASSBAND_THRESHOLD
                        Threshold for determining passband relative to max
                        value of S21. Default -3dB
  --isolation-threshold ISOLATION_THRESHOLD, -it ISOLATION_THRESHOLD
                        Threshold for determining isolation relative to max
                        value of s21. Default -10dB
  --merge-distance MERGE_DISTANCE, -md MERGE_DISTANCE
                        For denoising. If passbands are closer (1+md)*bw they
                        are merged in to one. Default 0.05
  --plot, -p            Plot filter data

```


## mixer_sim.py
Simulate spectral components from a mixer. Adds and keeps track of harmonics. Can be used with single frequency simulations or s21 can be read from a s2p-file to simulate spectral response.
Multiple s2p-files can be given for if and rf filter, which will simply cascade to form combined s-parameters. Note that there will be a lot of spectral components som the legend will be very comperhensive, use -nl or --no-legend to hide.

Example usage:
```
$ ./mixer_sim.py --rf-filter bandpassfilterbank/denoised/filter_1.s2p filters/denoised/filter_circtel_lm9_5400.s2p --if-filter if_filters/denoised/filter_ewt_14_2_to_14_55.s2p -lau -nl
```

Will show three figures, one with input parameters, one with all spectral components on IF and one with all spectral components on IF filtered by IF-filter.
Help:
```
usage: mixer_sim.py [-h] [--rf-filter [RF_FILTER ...]] [--rf-frequency RF_FREQUENCY] [--if-filter [IF_FILTER ...]] [--lo-frequency LO_FREQUENCY] [--lo-auto-lsb] [--lo-auto-usb]
                    [--lo-harmonics {1,2,3,4,5}] [--rf-harmonics {1,2,3,4,5}] [--no-legend]

Calculate harmonic transfer functions for mixers

options:
  -h, --help            show this help message and exit
  --rf-filter [RF_FILTER ...], -rffi [RF_FILTER ...]
                        S2P file containing rf (preselector) filter. S21 will be used
  --rf-frequency RF_FREQUENCY, -rffr RF_FREQUENCY
                        Single tone representing RF signal
  --if-filter [IF_FILTER ...], -iff [IF_FILTER ...]
                        S2P file containing if (roofing) filter. S21 will be used
  --lo-frequency LO_FREQUENCY, -lof LO_FREQUENCY
                        Frequency of LO
  --lo-auto-lsb, -lal   Automatically determine LO frequency so rf lsb will end up in middle of IF-filter
  --lo-auto-usb, -lau   Automatically determine LO frequency so rf usb will end up in middle of IF-filter
  --lo-harmonics {1,2,3,4,5}, -lh {1,2,3,4,5}
                        Number of lo harmonics to consider, default = 5
  --rf-harmonics {1,2,3,4,5}, -rh {1,2,3,4,5}
                        Number of rf harmonics to consider, default = 3
  --no-legend, -nl      Do not display legend in plots containing spectral components

```

## remove_noise.py
Remove noise from s-parameter files. 
Example usage:
```
$ ./remove_noise.py noisy.s2p --write denoised.s2p --plot
```
Help:
```
usage: remove_noise.py [-h] [--write WRITE] [--epsilon [EPSILON]]
                       [--filter-length [FILTER_LENGTH]]
                       [--threshold [THRESHOLD]]
                       [--dilation [{[0.000000-1.000000]}]] [--plot]
                       input

Condition s2p-files by removing trace noise du to low signal level

positional arguments:
  input                 S2P file to condition

optional arguments:
  -h, --help            show this help message and exit
  --write WRITE, -w WRITE
                        Write conditioned data to this s2p-file
  --epsilon [EPSILON], -e [EPSILON]
                        Replace noise with this value in dB. Default -130 dB
  --filter-length [FILTER_LENGTH], -f [FILTER_LENGTH]
                        Number of taps in averaging filter. Default 50 taps
  --threshold [THRESHOLD], -t [THRESHOLD]
                        Threshold in dB for noise. Applied to second
                        derivative of s-parameter filtered using filter above.
                        Default 4 dB
  --dilation [{[0.000000-1.000000]}], -d [{[0.000000-1.000000]}]
                        Dilation threshold to make contiguous noise regions.
                        Default 0.1
  --plot, -p            Plot result

```
