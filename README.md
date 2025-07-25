# Analytical Reconstruction of Dark-Field CT

In this work, we proposed the analytical reconstruction of dark-field CT.
The dark-field CT is modeled as a weighted Radon transform due to its positional dependency, 
and can be reconstructed based on the corresponding inversion formula.

This repository provides the implementation in the following paper:

**Analytical Reconstruction of Human-Scale Dark-Field CT** <br/>
IEEE Transactions on Medical Imaging ([TMI](https://ieeexplore.ieee.org/)) <br/>
[paper link](links-to-be-updated)


## Documentation & Demonstration

The principles of the inversion formula can be found in the [references attached here](articles/).

See the code samples in `demo/` directory for implementation of the reconstruction functions.
Unfortunately github doesn't support uploading large files, so some of the data used in demos are provided in [releases](https://github.com/GPY1999/DFCT-Reconstruction/releases/tag/data).

The definition of geometry parameters can refer to Fig. 4 of [our paper](#analytical-reconstruction-of-dark-field-ct).

### Principles of the inversion formula

I would recommend these references:

1. Boman, J., Strömberg, JO. Novikov’s inversion formula for the attenuated Radon transform—A new approach. J Geom Anal 14, 185–198 (2004). https://doi.org/10.1007/BF02922067
2. Huang Q, Zeng G L, You J, et al. An FDK-like cone-beam SPECT reconstruction algorithm for non-uniform attenuated projections acquired using a circular trajectory. Phys. Med. Biol. 50, 2329 (2005). https://doi.org/10.1088/0031-9155/50/10/010

## Dependencies

- MATLAB, in which the codes are written

- [ASTRA toolbox](https://astra-toolbox.com/): Optional, only used for comparison with the proposed analytical algorithm.

## Citation

If you are inspired by our work or use these codes for your research, we would appreciate it if you would cite our paper:

```
to be updated
```

## Acknowledgements

This study is supported in part by National Natural Science Foundation of China (No. 62227804), National key research and development program of China (No. 2023YFC2605802) and Beijing Natural Science Foundation (No. L232115).

## Contact

GUO Peiyuan: [guopy21@mails.tsinghua.edu.cn](guopy21@mails.tsinghua.edu.cn)

WANG Zhentian: [wangzhentian@tsinghua.edu.cn](wangzhentian@tsinghua.edu.cn)
