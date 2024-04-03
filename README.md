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

## mixer_sim.py
Simulate spectral components from a mixer. Adds and keeps track of harmonics. Can be used with single frequency simulations or s21 can be read from a s2p-file to simulate spectral response.
Multiple s2p-files can be given for if and rf filter, which will simply cascade to form combined s-parameters. Note that there will be a lot of spectral components som the legend will be very comperhensive, use -nl or --no-legend to hide.

Example usage:
```
$ ./mixer_sim.py --rf-filter bandpassfilterbank/denoised/filter_1.s2p filters/denoised/filter_circtel_lm9_5400.s2p --if-filter if_filters/denoised/filter_ewt_14_2_to_14_55.s2p -lau -nl
```

Will show three figures, one with input parameters, one with all spectral components on IF and one with all spectral components on IF filtered by IF-filter.

## remove_noise.py
Remove noise from s-parameter files. 
