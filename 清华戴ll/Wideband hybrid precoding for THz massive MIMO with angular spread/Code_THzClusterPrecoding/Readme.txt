This simulation code package is mainly used to reproduce the results of the following paper [1]:

[1] Mingyao Cui, Jingbo Tan, and Linglong Dai, "Wideband hybrid precoding for THz massive MIMO with angular spread (in Chinese)," Sci Sin Inform, Aug. 2022. 

*********************************************************************************************************************************
If you use this simulation code package in any way, please cite the original paper [1] above. 
 
The author in charge of this simulation code pacakge is: Mingyao Cui (email: cmy20@mails.tsinghua.edu.cn).

Reference: We highly respect reproducible research, so we try to provide the simulation codes for our published papers 
( more information can be found at: 
http://oa.ee.tsinghua.edu.cn/dailinglong/publications/publications.html )

Please note that the MATLAB R2020b is used for this simulation code package,  
and there may be some imcompatibility problems among different MATLAB versions. 

Copyright reserved by the Broadband Communications and Signal Processing Laboratory (led by Dr. Linglong Dai), 
Beijing National Research Center for Information Science and Technology (BNRist), 
Department of Electronic Engineering, Tsinghua University, Beijing 100084, China. 

*********************************************************************************************************************************
Abstract of the paper: 

Terahertz (THz) wideband hybrid precoding is a promising technology for future wireless communications. 
Since the bandwidth in THz is very large, a severe beam split effect will be induced to limit the system capacity. 
Existing THz wideband precoding schemes mainly focus on overcoming the beam split effect in line-of-sight (LOS) scenarios, 
while these schemes will suffer from unacceptable performance loss in non-line-of-sight (NLOS) scenarios. 
To address this challenging problem, we first analyze and compare the difference of beam split effect in LOS scenarios and NLOS scenarios. 
Then, a delay-time joint optimization algorithm is proposed to overcome the beam split effect in NLOS scenarios. 
Finally, simulations are provided to verify the effectiveness of the proposed method.
*********************************************************************************************************************************
How to use this code:

1. By setting parameters 'sigma_mt' and 'sigma_bs' in the file "AverageRate_vs_SNR.m" to 0 and running this file, we can obtain Fig. 5.
2. By setting parameters 'sigma_mt' and 'sigma_bs' in the file "AverageRate_vs_SNR.m" to 5*pi/180 and running this file, we can obtain Fig. 6.
3. Run "AverageRate_vs_AS.m" and obtain Fig. 7.
4. Run "EE_vs_K.m" and obtain Fig. 8.
5. Run "AverageRate_vs_K.m" and obtain Fig. 9.

*********************************************************************************************************************************
Enjoy the reproducible research!