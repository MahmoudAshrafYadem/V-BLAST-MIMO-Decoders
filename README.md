# V-BLAST-MIMO-Decoders
MATLAB implementation and performance analysis of V-BLAST MIMO detection algorithms (ZF, ZF-SIC, MD, MMSE, MMSE-SIC)

## Overview
This project presents a comprehensive comparative study of five V-BLAST (Vertical Bell Labs Layered Space-Time) MIMO detection algorithms, analyzing their performance in terms of Bit Error Rate (BER) versus number of receive antennas and SNR.

## Detection Algorithms
This repository implements and compares the following detection algorithms:

1. **Zero Forcing (ZF)** - Linear detector with simple implementation
2. **Zero Forcing with Successive Interference Cancellation (ZF-SIC)** - Improved ZF with interference cancellation
3. **Maximum Likelihood / Minimum Distance (MD/ML)** - Optimal detector with exponential complexity
4. **Minimum Mean Square Error (MMSE)** - Linear detector with better noise handling
5. **MMSE with Successive Interference Cancellation (MMSE-SIC)** - Best performance-complexity trade-off

## System Configuration
- 3×3 MIMO system (3 transmit antennas, 3 receive antennas)
- BPSK modulation
- Rayleigh fading channel
- Performance metrics: BER vs SNR and BER vs number of receive antennas

## Key Findings
- MD detector provides the best BER performance but with exponential computational complexity
- MMSE-SIC achieves near-optimal performance with reasonable complexity
- MMSE-SIC is recommended for practical wireless communication systems
- Performance improves with increasing number of receive antennas

## Repository Contents
- `Task1.m` - MATLAB implementation of detection algorithms
- `Task2.m` - Performance analysis and simulation
- `Report.pdf` - Detailed project report with theoretical background, simulation results, and analysis

## How to Run
1. Open MATLAB
2. Run `Task1.m` for basic algorithm implementation
3. Run `Task2.m` for performance analysis and plotting

## Documentation
For detailed theoretical background, methodology, simulation parameters, and comprehensive analysis, please refer to the included `Report.pdf`.

## Author
Mahmoud Ashraf Mahmoud

## References
See the project report for complete references and citations.
