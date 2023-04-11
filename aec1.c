// 1. 低频指导高频
void HighProc(float *gain)
{
    int bins = 257;
    const int FREQ_300HZ = 300;
    const int FREQ_1200HZ = 1200;
    const int FREQ_8000HZ = 8000;

    int lowFreq = FREQ_300HZ / FREQ_8000HZ * bins;
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
    curERL[i] = (dPowSmooth[i] - noiseEst[i]) / (xPowSmooth[i] + EPS);
    curSER[i] = (errEstPow[i]  - noiseEst[i]) / (echoEstPow[i] + EPS);
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
for (int i = 0; i < 514; ++i) {
    tmpErrEst[i] = max(tmpErrEst[i], 0);
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
  // 分成0~4K低频做和4~8K
  // cemp8k: 0~4K倒谱
  // cemp16k: 4~8K倒谱
  // echoEst: 回声估计
  // 先求倒谱
    for (i = 0; i < 123; ++i) {
        tmpPow[i] = echoEst[2 * i] * echoEst[2 * i] + 
                    echoEst[2 * i + 1] * echoEst[2 * i + 1];
    }

    tmpPow = Log(tmpPow);

    tmpPow[123] = 0;
    for (i = 123; i >= 0; --i) {
        tmpPow[2 * i ] = tmp[i];
        tmpPow[2 * i + 1] = 0;
    }

    IFFT(tmpPow);

    for (i = 0; i < 257; ++i) {
        tmpPow[i] = max(tmpPow[i], 0);
    }

    cemp8k = FFT(tmpPow);

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
  // a. 计算分母能量
  for (i = 0; i < 257; ++i) {
    X_power[i] = X[2 * i] * X[2 * i] + X[2 * i + 1] * X[2 * i + 1];
    Y_power[i] = Y[2 * i] * Y[2 * i] + Y[2 * i + 1] * Y[2 * i + 1];
    
    X_power_smooth[i] = 0.9 * X_power_smooth[i] + 0.1 * X_power[i];
  }
  
  // b. 互相关
  for (i = 0; i < 257; ++i) {
    corr_xy[2 * i] = X[2 * i] * Y[2 * i] + X[2 * i + 1] * Y[2 * i + 1];
    corr_xy[2 * i + 1] = X[2 * i + 1] * Y[2 * i] - X[2 * i] * Y[2 * i + 1]; 
    smooth_corr_xy[2 * i] = 0.1 * corr_xy[2 * i] + 0.9 * smooth_corr_xy[2 * i];
    smooth_corr_xy[2 * i + 1] = 0.1 * corr_xy[2 * i + 1] + 0.9 * smooth_corr_xy[2 * i + 1]
  }
   // c. 计算相干性
  for (i = 0; i < 257; ++i) {
    tmp = smooth_corr_xy[2 * i] * smooth_corr_xy[2 * i] + smooth_corr_xy[2 * i + 1] * smooth_corr_xy[2 * i + 1];
    cohxd[i] = tmp / (X_power_smooth[i] * Y_power_smooth[i] + EPS);
  }

  // d. 计算相干性ERL
  for(i = 0; i < 257; ++i) {
    coheec[i] = min(1, coheec[i]);
    cohed[i] = min(1, cohed[i]);
    cohERL[i] = coheec * (1 - cohed[i]); // 不做平滑
  }

// 7. 根据cohxd, cohxe计算NLP增益
   for (i = 0; i < 400Hz; ++i) {
       gain = ((cohxd[i] - cohxe[i]) * 1.5 + 1.0) * cohxe[i];
       gain = max(gain, cohxe[i]);
       gain = min(gain, 0.95); // 0.95 不同频域设置不同值，频带越大设置越大
       gain = 1 - gain;
   }
   

// 8. NLMS步长更新
  // errEst: 误差信号估计
  // curErrPow：当前帧误差能量平滑值
  // constErrPow: 0.98平滑系数的误差信号能量平滑值
  // smoothErrPow: 0.85平滑系数的误差能量平滑值
  // refPow: 参考信号能量
  // FilterBlockNum： 滤波器块数

  float mu = 0.1;
  for (i = 0; i < 257; ++i) {
    tmp1 = min(sqrt(constErrPow[i] / max(EPS, curErrPow[i])), 1);
    tmp2 = refPow[i] * FilterBlockNum + smoothErrPow[i];
    R = tmp1 * errEst[2 * i] / (tmp2 + EPS);
    I = tmp1 * errEs[2 * i + 1] / (tmp2 + EPS);

    alpha = 0.002f / max(sqrt(R * R + I * I), 0.002f);
    R *= alpha;
    I *= alpha;
    
    step[2 * i] = mu * R;
    step[2 * i + 1] = mu * I;
  }

// 9. 基频检测
  // mag: 幅度谱

  for (i = 0; i < 257; ++i) {
    cemp[i] = mag[2 * i] * mag[2 * i] + mag[2 * i + 1] * mag[2 * i + 1];
  }
  cemp = log(cemp);

  for (i = 257; i >= 0; --i) {
    cemp[2 * i] = cemp[i];
    cemp[2 * i + 1] = 0;
  }

  cemp = IFFT(cemp)

  L = Fs/400 - 1;
  H = Fs / 70 - 1;

  pitchId = -1;
  for (L:H) {
    if (cemp[i] >= 0.36) {
        pitch = i;
        break;
    }
  }

// 10. 能量求dB
  // energy: 能量谱
  // [L, H]: 求dB的能量谱范围
  
  for (L: H) {
    sumEnergy += energy[i];
  }
  coff = (H - L) * 0.5 * stft_size;
  tmp = sumEnergy / coff;
  dB = LOG(tmp) * 0.434294481903252;
  dB = max(10 * dB, -20);

// 11. 根据cohxd, cohxe计算NLP增益
   for (i = 0; i < 400Hz; ++i) {
       gain = ((cohxd[i] - cohxe[i]) * 1.5 + 1.0) * cohxe[i];
       gain = max(gain, cohxe[i]);
       gain = min(gain, 0.95); // 0.95 不同频域设置不同值，频带越大设置越大
       gain = 1 - gain;
   }

 
