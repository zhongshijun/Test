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
float echoEstPow[257]; // 回声估计平滑能量[0-257]
float errEstPow[257]; // 残差估计平滑能量[0-257]

float curERL[257];
float curSER[257];
float smoothERL[257];
float smoothSER[257];

for (int i = 0; i < 257; ++i) {
    curERL[i] = (dPowSmooth[i] - noiseEst[i]) / max(xPowSmooth[i], 1);
    curERL[i] = max(EPS, curERL[i]);
    curSER[i] = (errEstPow[i]  - noiseEst[i]) / max(echoEstPow[i], 1);
    curSER[i] = max(EPS, curSER[i]);
}

// alpha需要根据单双讲，当前帧回声大小动态调整
SMOOTH(smoothERL, curERL, alpha);
SMOOTH(smoothSER, curSER, alpha);

// 可以根据在[400Hz, 2000Hz]范围内，就散curSER, curERL, smoothERL, smoothSER来作为
// 单双讲的辅助判断依据

// 2.1 ERL回声估计
  for (i = 0; i < 257; ++i) {
    X_pow[i] = echoEst[2 * i] * echoEst[2 * i] + echoEst[2 * i + 1] + echoEst[2 * i + 1];
  }

  for (i = 0; i < 257; ++i) {
    tmp3 = tmp * smoothERL;
    Y_pow = echoPow[2 * i ] * echoPow[2 * i ] + echoPow[2 * i + 1] * echoPow[2 * i + 1]; // 谐波增强后的回声估计

    echoPow3[i] = tmp3 + Y_pow - 2 * (X_pow[i] * Y_pow);
    echoPow3[i] = min(echoPow3[i], tmp3);
    echoPow3[i] = max(echoPow3[i], 0);
  }

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
  // errEst误差谱
  // 先求倒谱
    for (i = 0; i < 128; ++i) {
        tmpPow[i] = errEst[2 * i] * errEst[2 * i] + 
                    errEst[2 * i + 1] * errEst[2 * i + 1];
    }

    tmpPow = Log(tmpPow);
    for (i = 128; i >= 0; --i) {
        tmpPow[2 * i ] = tmpPow[i];
        tmpPow[2 * i + 1] = 0;

        cemp16k[2 * i] = tmpPow[i];
        cemp16k[2 * i + 1] = tmpPow[i];
    }
    tmpPow[258] = tmpPow[256] - 4;
    tmpPow[257] = 0;

    cemp8k = IFFT(tmpPow);

    // 16k
    for (i = 128; i < 257; ++i) {
      tmpPow[i - 128] = errEst[2 * i] * errEst[2 * i] + 
            errEst[2 * i + 1] * errEst[2 * i + 1];
    }

    tmpPow = Log(tmpPow);

    for (i = 128; i >= 0; --i) {
      cemp16k[2 * (i + 128)] = tmpPow[i];
      cemp16k[2 * (i + 128) + 1] = 0;
    }
    tmpPow[514] = tmpPow[512] - 4;
    tmpPow[513] = 0;

    cemp16k = IFFT(cemp16k);

    
    //  倒谱平滑值 
    alpha[16] = [0, 0, 0, 0, 0, 0, 0.1, 0.4, 0.6, 0.8, 0.8, 0.85, 0.9, 0.95, 0.95, 0.997]
    for (i = 0; i < 16; ++i) {
        smooth_cemp8k = SMOOTH(cemp8k, alpha[i]);
    }
    alpha1 = alpha[14];
        for (i = 16; i < 126; ++i) {
        smooth_cemp8k = SMOOTH(cemp8k, alpha1);
    }

    cempMod: // 原始倒谱值 [2k-4k]
    cempModAmp: //平滑倒谱值 [0-2k]
    for (i = 0; i < 126; ++i) {
        tmp1 = cemp8k[i];
        tmp2 = smooth_cemp8k[i];
        cempMod = tmp2;
        if (tmp1 > 0.36 && i >= 15) {
            smooth_cemp8k[i] = 0;
            cempMod = tmp1;
            tmp2 = cempMod * 5;
        }
        cempMod[i] = cempMod; // 2~4k
        cempModAmp[i] = tmp2; // 0~2k

        cempMod[i + 126] = cempMod;
        cempModAmp[i + 126] = cempMod;
    }

    cempMod = FFT(cempMod); // 0 ~ 256
    cempModAmp = FFT(cempModAmp); // 0 ~ 256
    
    idx = 126 / 2;
    for (i = 0; i < idx; ++i) {
      epxCemp[i] = cempModAmp[2 * i];
      epxCemp[i + idx] = cempMod[2 * (i + idx)];
    }

    epxCemp = Exp(epxCemp);
    for (i = 0; i < 126; ++i) {
        snr = epxCemp[i];
        gain[i] = snr / (echoEst + EPS);
    }

    // 16K
    for (i = 0; i < 16; ++i) {
        smooth_cemp16k = SMOOTH(cemp8k, alpha[i]);
    }
    alpha1 = alpha[14];
        for (i = 16; i < 126; ++i) {
        smooth_cemp16k = SMOOTH(cemp8k, alpha1);
    }

    for (i = 0; i < 126; ++i) {
      epxCemp[i] = smooth_cemp16k[i];
    }
    epxCemp = EXP(epxCemp);
    for (i = 0; i < 126; ++i) {
        snr = epxCemp[i];
        gain[i + 126] = snr / (echoEst + EPS);
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

// 5.1 计算相干性ERL
  // coheec: 回声平滑估计和误差平滑估计相干性
  // cohed: 误差平滑估计和近端信号平滑估计相干性
  for(i = 0; i < 257; ++i) {
    coheec[i] = min(1, coheec[i]);
    cohed[i] = min(1, cohed[i]);
    cohERL[i] = coheec * (1 - cohed[i]); // 不做平滑
  }

// 7. 根据cohxd, cohxe计算NLP增益
   for (i = 0; i < 96; ++i) {
       gain = ((cohxd[i] - cohxe[i]) * 1.5 + 1.0) * cohxe[i];
       gain = max(gain, cohxe[i]);
       gain = min(gain, 0.95); // 0.95(0~3K)/0.98(3~8K) 不同频域设置不同值，频带越大设置越大
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

// 11. 回声谐波增强
  // enhance
  for (i = 0; i < 257; ++i) {
    echoPow[i] = echoEst[2 * i] * echoEst[2 * i] + echoEst[2 * i + 1] * echoEst[2 * i + 1];
  }

  for (i = 0; i < 514; ++i) {
    tmpFft[2 * i] = echoPow[i];
    tmpFft[2 * i + 1] = 0;
  }
  tmpFft = IFFT(tmpFft);

  for (i = 0; i < 514; ++i) {
    tmpFft[i] = max(0, tmpFft[i]);
  }

  tmpFft = FFT(tmpFft);

  for (i = 0; i < 257; ++i) {
      enhance_echoPow[i] = fabs(tmpFft[2 * i]);
  }

  // 回声估计
  for (i = 0; i < 257; ++i) {
    echoEstPow[i] = max(enhance_echoPow[i] * 0.2, echoEstPow[i]); // 0.2可以根据回声情况调整
  }

  for (i = 0; i < 257; ++i) {
    echoEstPow[i] =  echoEstPow[i] * cohERL[i]; // 相干性ERL
  }

  if (frame == 1) {
    for (i = 0; i < 257; ++i) {
      smooth_echoEstPow[i] = echoEstPow[i];
    }
  }

  alpha = 0.6;
  if (近端单讲) {
    alpha = 0.8;
  }

  for (i = 0; i < 257; ++i) {
    tmp = echoEstPow[i];
    smooth_echoEstPow[i] = alpha * smooth_echoEstPow[i] + (1.0 - alpha) * tmp;
    echoEstPow[i] = max(tmp, smooth_echoEstPow[i]);
  }


// 12. kalman filter
  // a. 获取当前回声估计值
  echoEst, errEst = filter(kalman.weight, X, Y);
  
  // b. 获取S
  for (i = 0; i < numFilterBlock; ++i) {
    P = kalman->P + i * 257;
    x_pow = X_pow + i * 257; // 当前帧参考能量
    for (j = 0; j < 257; ++j) {
      S[i * 257 + j] = P[j] * x_pow[j]; 
    }
  }

  for (i = 0; i < 257; ++i) {
    S[i] = 1.0 / (S[i] + smoothErrEstPow[i] + EPS); // smoothErrEstPow: alhap=0.85, frame=1,alpha=0
    M2[i] = 1.0;
  }

  for (0~3K) {
    Pmin[i] = 0.01;
  }
  
  for (3K~8K) {
    Pmin[i] = 0.01;
  }

  // kalman参数更新
  for (i = 0; i < FilterBlockNum; ++i) {
    P = P + i * 257;
    S = S + i * 257;
    W = W + i * 257;
    Q = Q + i * 257;
    M2 = M2 + i * 257;
    // X = X + i * 257; ??
    // errEst = errEst + i * 257;
    for (j = 0; j < 257; ++j) {
      mu = P[j] * S[j];

      R = errEst[2 * j] * X[2 * j] + errEst[2 * j + 1] * X[2 * j + 1]; // ??
      I = errEst[2 * j + 1] * X[2 * j] - errEst[2 * j] * X[2 * j + 1]; // ??

      // 更新W
      W[2 * j] = W[2 * j] + R * mu;
      W[2 * j + 1] = W[2 * j + 1] + I * mu;

      // 更新P
      P[j] = P[j] + (Q[j] - P[j] * mu * X2); // X2??
      P[j] = max(P[j], Pmin);
      P[j] = min(P[j], 0.99);

      // 更新M2
      M2[j] = M2[j] * (1 - mu * X2); // X2?

      // 更新Q
      R = R * mu;
      I = I * mu;
      Q[j] = 0.98 * Q[j] + 0.02 * (R * R + I * I);
    }
  }

  for (i = 0; i < FilterBlockNum; ++i) {
    M2 = M2 + i * 257;
    for (j = 0; j < 257; ++j) {
      kalman->ABDf[i].p[j] = min(1.0, max(0.2, M2[j]));
    }
  }
