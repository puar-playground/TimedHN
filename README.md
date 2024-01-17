# TimedHN
This is the Matlab implementation of [TimedHN](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0283004), a time-aware oncogenetic modeling method.<br />
## Description of demo scripts
1. The script `TimedHN_demo.m` is for the inferrence of the hazard network.<br />
2. The script `pseudo_order_demo.m` computes the conditional expectation of progression times using a learned hazard network.<br />
3. The script `simulation_demo.m` is for the generation of synthetic data.<br />

## A visualization of the accumulation process
TimedHN progress on 7-cube in layered layout. The animation shows the accumulation process of 7 events. The observation probabilities of all $2^7$ possible states are shown in redish gradient color.<br />
![layered_7_cube](https://user-images.githubusercontent.com/45474252/202872997-e7966248-d583-4981-926b-b3f2c5779af6.gif)


## Reference
```
@article{chen2023timed,
  title={Timed hazard networks: Incorporating temporal difference for oncogenetic analysis},
  author={Chen, Jian},
  journal={Plos one},
  volume={18},
  number={3},
  pages={e0283004},
  year={2023},
  publisher={Public Library of Science San Francisco, CA USA}
}
```
