# TimedHN
This is the Matlab implementation of TimedHN, a time-aware oncogenetic modeling method.<br />
## Description of demo scripts
1. The script `TimedHN_demo.m` is for the inferrence of the hazard network.<br />
2. The script `pseudo_order_demo.m` computes the conditional expectation of progression times using a learned hazard network.<br />
3. The script `simulation_demo.m` is for the generation of synthetic data.<br />

## A visualization of the accumulation process
TimedHN progress on 7-cube in layered layout. The animation shows the accumulation process of 7 events. The observation probabilities of all $2^7$ possible states are shown in redish gradient color.<br />
![layered_7_cube](https://user-images.githubusercontent.com/45474252/202872997-e7966248-d583-4981-926b-b3f2c5779af6.gif)


## Reference
1. [Jian Chen. "Inferring time-aware models of cancer progression
using Timed Hazard Networks." bioRxiv (2022): 2022-10.](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0283004)<br />


