// 1. 低频指导高频
void HighProc(float *gain)
{
    int bins = 257;
    const int FREQ_400HZ = 300;
    const int FREQ_1200HZ = 1200;
    const int FREQ_8000HZ = 8000;

    int lowFreq = FREQ_400HZ / FREQ_8000HZ * bins;
    int highFreq = FREQ_1200HZ / FREQ_8000HZ * bins;
    int band = (highFreq - lowFreq) / 4 + 0.5;

    int maxId = bins;
    int lenght = 0;
    int lastLength = 0;
    for (int i = lowFreq; i <= highFreq; ++i) {
        if (gain[i] < 0.05) {
            lenght += 1;
            lastLength += 1;
        } else {
            if (lastLength > band) {
               maxId = i;
            } else {
                lastLength = 0;
            }

        }

    }

    if ((lenght > (highFreq - lowFreq) / 2) && (lastLength > band)) {
        for (int i = maxId; i < bins; ++i) {
            gain[i] = 0;
        }
    }
}

// 2. 计算ERL， SER
float xPowSmooth[257]; // 远端信号平滑能量谱 [0-257]
float dPowSmooth[257]; // 近端信号平滑能量谱 [0-257]
float echoEstPow[257]; // 回声估计
float errEstPow[257]; // 残差估计

float curERL[257];
float curSER[257];
float smoothERL[257];
float smoothSER[257];

for (int i = 0; i < 257; ++i) {
    curERL[i] = dPowSmooth[i] / (xPowSmooth[i] + EPS);
    curSER[i] = errEstPow[i] / (echoEstPow[i] + EPS);
}

SMOOTH(smoothERL, curERL, 0.98);
SMOOTH(smoothSER, curSER, 0.98);


// 3. 谐波增强消除回声
const int FFT_SIZE = 514;
float echoEstPow[257]; // 回声估计能量谱
float errEst[514];   // 滤波器+后滤波误差信号估计
float tmpErrEst[514]; // 误差估计增强信号
float snrEchoEst[257]; // 信回比

copy(errEst, tmpErrEst, FFT_SIZE); // errEst->tmpErrEst

IFFT(tmpErrEst);
for (int i = 0; i < 257; ++i) {
    tmpErrEst[2 * i] = max(tmpErrEst[2 * i], 0);
}
FFT(tmpErrEst);

for (int i = 0; i < 257; ++i) {
    float tmpR1 = tmpErrEst[2 * i];
    float tmpI1 = tmpErrEst[2 * i + 1];

    float tmpR2 = errEst[2 * i] - tmpR1;
    float tmpI2 = errEst[2 * i] - tmpI1;

    float tmp = ((tmpR1 * tmpR1 + tmpI1 * tmpI1) + (tmpR2 * tmpR2 + tmpI2 * tmpI2)) / 2;

    snrEchoEst[i] = 0.2 * snrEchoEst[i] + 0.8 * (tmp / (echoEstPow[i] + EPS));
}

// 3.1 根据snr计算回声最终增益值
for (i = 0; i < 257; ++i) {
    tmp1 = snrEchoEst[i];
    tmp2 = tmp1 + 1;
    gain[i] = tmp1 / tmp2;
    gain[i] = min(gain[i], aec->gainMinLimit);
}


// 4. 倒谱消回声
for () {


}

// 5. 相干性计算
   X: 远端频域数据
   Y: 近端频域数据
   errEst: 误差估计/残差信号
   echoEst: 回声估计
   
   X_power： 远端信号能量 X_power_smooth 能量平滑信号
   Y_power: 近端信号能量
   err_power: 误差信号能量
   echo_power: 回声信号能量
   
   corr_xy: 远近端互相关
   corr_xerr: 远端和残差互相关
   
  // cohxd
  // 1. 计算分母能量
  for (i = 0; i < 257; ++i) {
    X_power[i] = X[2 * i] * X[2 * i] + X[2 * i + 1] * X[2 * i + 1];
    Y_power[i] = Y[2 * i] * Y[2 * i] + Y[2 * i + 1] * Y[2 * i + 1];
    
    X_power_smooth[i] = 0.9 * X_power_smooth[i] + 0.1 * X_power[i];
  }
  
  // 2. 互相关
  for (i = 0; i < 257; ++i) {
    corr_xy[2 * i] = X[2 * i] * Y[2 * i] + X[2 * i + 1] * Y[2 * i + 1];
    corr_xy[2 * i + 1] = X[2 * i + 1] * Y[2 * i] - X[2 * i] * Y[2 * i + 1]; 
  }
   // 3. 计算相干性
  for (i = 0; i < 257; ++i) {
    tmp = corr_xy[2 * i] * corr_xy[2 * i] + corr_xy[2 * i + 1] * corr_xy[2 * i + 1];
    cohxd[i] = tmp / (X_power_smooth[i] * Y_power_smooth[i] + EPS);
  }

// 6. 根据cohxd, coheec计算NLP增益
   for (i = 0; i < 400Hz; ++i) {
       cohxd[i] = min(cohxd[i], 1);
       coheec[i] = min(coheec[i], 1);
       gain = ((cohxd[i] - coheec[i]) * 1.5 + 1.0) * coheec[i];
       gain = min(gain, 1.0);
   }
   


