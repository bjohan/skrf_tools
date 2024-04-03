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
