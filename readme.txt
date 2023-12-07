This code is for receiver function velocity analysis technique (RFVAT).
For more details, please refer to:
Shi, J., Wang, T., & Chen, L. (2020).Receiver function velocity analysis technique and its application to remove multiples. Journal of GeophysicalResearch: Solid Earth, 125,e2020JB019420. https://doi.org/10.1029/2020JB019420
##########################
Here are the instructions for environmental configuration and the parameters in the paraRFVAT.cfg.
##########################
To successfully run this code, please ensure the python version is above 3.0, simply because the author dose not test it in python2.0. Besides, additional libraries are also needed, including numpy, obspy, tqdm, and matplotlib. For reference, the author currently runs in python 3.7.3 and the versions of such libraries are numpy1.16.2, obspy1.1.1, tqdm4.31.1, and matplotlib3.1.1.
##########################
The author prepares an example file which contains an input data file and an output image file, so the reader can directly run this code through this example.
[python RFVAT.py -r paraRFVAT.cfg]
##########################
In particular:
1. As the example file, the reader needs a data and an image files, then assign their relative paths to the data_path and out_path variables.In data file, the type of receiver function (rf) data is SAC, whose header variable of user0 is ray parameter and its unite is s/km.
2. nsample: the sampling number of the rf trace.
3. shift:  the time shift of the rf trace (the time before the direct P wave). (unit: second)
4. dt: the sampling interval of the rf trace. (unit: second)
5. t_max: the maximum t0 of scanning. (unit: second)
6. e_max: the maximum E of scanning. (E = Vp*Vs)
7. t1: all of the upper time boundaries for the extrema in the Ps stacking amplitude spectrum.
8. t2: all of the lower time boundaries for the extrema in the Ps stacking amplitude spectrum.
9. k1: the ratio fo TPpPs(p=0)/TPs(p=0) in the crust.
10. k2: the ratio fo TPpPs(p=0)/TPs(p=0) in the mantle.
11. c: the arrival time bound of TPs(p=0) between the crust and mantle. (unit: second)
12. Thre: the threshold of stacking amplitude for the extrema in the PpPs stacking amplitude spectrum.  
Note: t1 and t2 are confirmed by the reader, between which the reader considers there may exist extrema. The author recommends reading the automatic scanning strategy in the paper mentioned above.
##########################
For the calculation of receiver function, the author suggest browsing seispy in https://xumijian.me.
##########################
Contact me, if you have any questions.
Copyright, please cite the above paper, if you use any part of this code.
Jing Shi
shij@smail.nju.edu.cn
 



 
