# SPaCR
---
## Spatial Peak-Aware Collaborative Representation for Hyperspectral Imagery Classification

***Chengle Zhou, Bing Tu, Qi Ren, Siyuan Chen***

*IEEE Geoscience and Remote Sensing, Early Access, 11 June 2021.*

[DOI:10.1109/LGRS.2021.3083416](https://ieeexplore.ieee.org/document/9452044)

---

![FIg.1](https://github.com/chengle-zhou/MY-IMAGE/raw/main/SPaCR/Fig.1.png)

Fig.1 Metrics of flat-topped hexagonal grids.



## Abstract

In this letter, a novel spatial peak-aware collaborative representation (SPaCR) method is proposed for hyperspectral imagery (HSI) classification, which introduces spectral-spatial information among superpixel clusters into egularization terms to construct a new CR-based closed-form solution. The proposed method is composed of the following key steps. First, the raw HSI is clustered into many superpixels according to an oversegmentation strategy. Then, cluster pixels  are determined based on spectral-spatial correlation between pixels within each superpixel. Next, spectral distance and  spatial coherence of superpixel clusters corresponding to training samples and testing pixels are fused to define  differences between pixels. Finally, the difference information between clusters as a spectral-spatial feature induced regularization term is incorporated into the objective function. Experimental results on the Indian Pines and University of Pavia HSIs indicated that the proposed SPaCR method, without any pre-processing and post-processing, outperforms well-known and state-of-the-art classifiers on the limited labeled samples.

---

The code is for the work, if you find it is useful, please cite our work:
>   @article{SPaCR2021,  
	title={Spatial Peak-Aware Collaborative Representation for Hyperspectral Imagery Classification},  
	author={Zhou, Chengle and Tu, Bing and Ren, Qi and Chen, Siyuan},  
	journal={IEEE Geoscience and Remote Sensing},  
	year={2021},  
    volume={},  
	number={},  
	pages={1-5},  
	doi={10.1109/LGRS.2021.3083416}}

If you need another two datasets (i.e., University of Pavia and Salinas Valley), please feel free to contact me. Or you can download them from http://www.ehu.eus/ccwintco/index.php/Hyperspectral_Remote_Sensing_Scenes

University of Pavia: http://www.ehu.eus/ccwintco/uploads/e/ee/PaviaU.mat, http://www.ehu.eus/ccwintco/uploads/5/50/PaviaU_gt.mat

Salinas Valley: http://www.ehu.eus/ccwintco/uploads/a/a3/Salinas_corrected.mat, http://www.ehu.eus/ccwintco/uploads/f/fa/Salinas_gt.mat
