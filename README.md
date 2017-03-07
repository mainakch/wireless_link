# wireless_link

Wireless_link aims to be a free, open source ray tracing software
for modeling wireless propagation in radio frequencies.  It
supports point transmitters and receivers, rectangular reflecting
surfaces, and perfect reflections off of the scattering
environment.  This project grew out of my research on
understanding the physical layer of multiple input multiple output
(MIMO) systems better.

## Installation

Tested with a modern gcc compiler.  Please make sure you have
cmake build system and GNU Scientific Library installed.
Then change the binary directory to a location in your PATH and type in

```
$bash INSTALL
```

## Running 
There are some example input files in `examples/`.  Pick a file and then type in 
```
wireless_link -i <example_input_file> -o <outputfile>
```

The code will output a time series of IQ measurements at each of the receiver antennas.  Note that due to finite speed of light, the receiver does not see the signal immediately when the transmit signal changes.

## Input format
