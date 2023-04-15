// sin: 原始麦克风当前帧输入
// ref: 参考信号当前帧输入

1. 将当前帧ref更新到历史XmRef缓存中；

2. 根据XmRef和sin的相干性计算时延大小（因果系统不需要缓存sin的历史帧）；

3. 根据2中计算的时延大小，获取对应的历史参考信号；

4. 根据3中获取的历史信号做滤波，获取回声估计；
   4.1 双滤波结构+后滤波误差信号,还可以设置较小的滤波器作为快速滤波；
   4.2 将双滤波器估计回声大的值作为最终的回声估计值，作为双滤波器的回声估计输出值；

5. 根据4中计算的回声估计，结合sin计算残余信号；

6. 滤波器发散估计；
   6.1 kalman误差信号大于100倍的sin能量，则判断为滤波器发散；
   6.2 NLMS
   6.3 需要防止误判断导致滤波器高频率被重置，静音或者小语音场景很容易导致误判发散；

7. 进入NLP非线性处理

8. 根据kamlan的滤波器权重系数来重新选择历史参考信号
   iDex = 0;
   iMax = 0;
   for (i = 0; i < numFilterBlock; ++i) {
    sum = 0;
    for (j = 0; j < 257; ++j) {
        sum += kalman->weight[i][j];
    }
    if (sum > iMax) {
        iMax = sum;
        iDex = i;
    }
   }

9. 扬声器模型，根据相干性做NLP回声估计，跟之前滤波器回声估计取大值

10. 计算各种场景近端单讲，远端单讲，双讲等参数，计算底噪，相干性等； 

11. 对9估计获取到的回声估计信号，做谐波增强。并根据相干性对回声估计进一步调整,echoPow；

12. 根据8中重新得到的参考信号，结合ERLsm获取回声估计，echoPow3；

13. 根据自适应对ERL, SER进行平滑ERLsm, SERsm；

14. 根据cohxdMean, coheecMean等参数，进行单双讲判断。
    并在双讲的时候计算afAlpha对回声echoPow调整；（双讲效果很重要）
    afAlpha参数根据alph1,2,3,4等进行计算;

15. 在远端单讲和回声起始段；
    for (i = 0; i < 257; ++i) {
        echoPow[i] = max(echoPow[i], echoPow3[i]);
    }

16. 使用谐波增强进行回声降噪gain；
     errEstPow *= gain

17. 根据16获取的errEstPow使用倒谱降噪进行回声降噪;

18. 做低频指导高频回声消除；

20. 可以根据kalman的ABDf进行回声gain进一步处理；
    for (i = 0; i < 257; ++i) {
        gain[i] *= kalman->ABDf[i];
    }
