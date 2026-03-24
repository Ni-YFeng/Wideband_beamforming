
This simulation code package is mainly used to reproduce the results of the following paper [1]:

[1] L. Dai, J. Tan, Z. Chen, and H. Vincent Poor, “Delay-phase precoding for wideband THz massive MIMO,” IEEE Trans. Wireless Commun., vol. 21, no. 9, pp. 7271-7286, Sep. 2022.

*********************************************************************************************************************************
If you use this simulation code package in any way, please cite the original paper [1] above. 
 
The author in charge of this simulation code pacakge is: Jingbo Tan (email: tanjb17@mails.tsinghua.edu.cn).

Reference: We highly respect reproducible research, so we try to provide the simulation codes for our published papers (more information can be found at: 
http://oa.ee.tsinghua.edu.cn/dailinglong/publications/publications.html). 

Please note that the MATLAB R2018b is used for this simulation code package,  and there may be some imcompatibility problems among different MATLAB versions. 

Copyright reserved by the Broadband Communications and Signal Processing Laboratory (led by Dr. Linglong Dai), Beijing National Research Center for Information Science and Technology (BNRist), Department of Electronic Engineering, Tsinghua University, Beijing 100084, China. 

*********************************************************************************************************************************
Abstract of the paper: 

Benefiting from tens of GHz bandwidth, terahertz (THz) communication has been considered as a promising technology to provide ultra-high speed data rates for future 6G wireless systems. To compensate for the serious propagation attenuation of THz signals, massive multiple-input multipleoutput (MIMO) with hybrid precoding can be utilized to generate directional beams with high array gains. However, the classical hybrid precoding architecture based on frequency-independent phase-shifters cannot cope with the beam split effect caused by the very large bandwidth in THz massive MIMO systems, where the directional beams will split into different physical directions at different subcarrier frequencies. The beam split effect will result in a serious array gain loss over the whole bandwidth, which has not been well investigated in THz massive MIMO systems to the best of our knowledge. In this paper, we first reveal and quantify the seriousness of the beam split effect in THz massive MIMO systems by analyzing the array gain loss caused by it. Then, we propose a new precoding architecture called delay-phase precoding (DPP) to mitigate this effect. Specifically, the proposed DPP introduces a time delay network as a new precoding layer between radio-frequency chains and phaseshifters in the classical hybrid precoding architecture. In this way, the conventional phase-controlled analog beamforming can be converted into the delay-phase controlled analog beamforming. Unlike the frequency-independent phase shifts, the time delay network introduced in the DPP can realize frequencydependent phase shifts, which can be elaborately designed to generate frequency-dependent beams towards the target physical direction over the whole THz bandwidth. Due to the joint control of delay and phase, the proposed DPP can significantly relieve the array gain loss caused by the beam split effect. Furthermore, we propose a hardware structure by using true-time-delayers to realize frequencydependent phase shifts for realizing the concept of DPP. Theoretical analysis and simulation results show that the proposed DPP can significantly mitigate the beam split effect over the whole THz bandwidth, and it can achieve the near-optimal achievable rate performance with higher energy efficiency than theclassical hybrid precoding architecture.
*********************************************************************************************************************************
How to use this simulation code package?

Fig. 8, Fig. 9, and Fig. 11 can be derived by running the corresponding m file.

The simulation results in Fig. 8  can be directly obtained by running Fig8.m.

The simulation results in Fig. 9  can be directly obtained by running Fig9.m.

The simulation results in Fig. 11  can be directly obtained by running Fig11.m. 
*********************************************************************************************************************************
Enjoy the reproducible research!