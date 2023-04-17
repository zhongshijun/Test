// 能量求dB
  // energy: 能量谱
  // [L, H]: 求dB的能量谱范围
  int PowdB(float *energy, int H, int L, int stft_size)
  {
    for (L: H) {
      sumEnergy += energy[i];
    }
    coff = (H - L) * 0.5 * stft_size;
    tmp = sumEnergy / coff;
    dB = LOG(tmp) * 0.434294481903252;
    dB = max(10 * dB, -20);

    return dB;
  }

// 【1】先做双线性滤波
// 1. NLMS + Kalman 滤波处理;
// 2. 根据[1]的滤波参数，再次做后滤波;
// 3. 选择两者回声能量大的作为最终回声估计输出;

// [1.1]. kalman filter
  // a. 获取当前回声估计值
  echoEst = filter(kalman.weight, X);
  errEst = Errestimation(Y, echoEst);

  float *x_pow[numFilterBlock]; // 缓存numFilterBlock历史参考信号
  for (i = 0; i < numFilterBlock; ++i) {
    x_pow_tmp = x_pow[i];
    X_tmp = X[i];
    for(j = 0; j < 257; ++j) {
        x_pow_tmp[j] = (X_tmp[2 * j] * X_tmp[2 * j] + X_tmp[2 * i + 1] * X_tmp[2 * i + 1]); 
    }

  }

  // b. 获取S
  for (i = 0; i < numFilterBlock; ++i) {
    P = kalman->P + i * 257;
    x_pow = X_pow + i * 257; // 当前帧参考能量
    for (j = 0; j < 257; ++j) {
      S[j] = S[j] + P[j] * x_pow[j]; 
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
    X = X[i];  // 缓存的历史参考信号
    X2 = x_pow[i]; // 第i帧历史参考能量
    for (j = 0; j < 257; ++j) {
      mu = P[j] * S[j];

      R = errEst[2 * j] * X[2 * j] + errEst[2 * j + 1] * X[2 * j + 1]; // errEst当前计算出来的误差信号，X对应时延帧的信号
      I = errEst[2 * j + 1] * X[2 * j] - errEst[2 * j] * X[2 * j + 1];

      // 更新W
      W[2 * j] = W[2 * j] + R * mu;
      W[2 * j + 1] = W[2 * j + 1] + I * mu;

      // 更新P
      P[j] = P[j] + (Q[j] - P[j] * mu * X2);
      P[j] = max(P[j], Pmin);
      P[j] = min(P[j], 0.99);

      // 更新M2
      M2[j] = M2[j] * (1 - mu * X2);

      // 更新Q
      R = R * mu;
      I = I * mu;
      Q[j] = 0.98 * Q[j] + 0.02 * (R * R + I * I);
    }
  }

  for (j = 0; j < 257; ++j) {
    kalman->ABDf.p[j] = min(1.0, max(0.2, M2[j]));
  }

// [1.2] NLMS步长更新
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

// [1.3] kalman录波器发散判断
    errPow = 0.
    sinPow = 0.
    for (j = 0; j < bins; ++j) {
      sinPow += micInPow[j];
      errPow += errEstPow[j];
    }
    if (errPow > 100 * sinPow) {
      // 发散，滤波器系数重置
    }
    
    // 注意为了防止遇到假发散，应该判断连续的发散帧后才重设滤波器系数；
    // 注意当远近端语音能量很小的时候很容易引起能量差异判断错误，导致假发散；所以要特别判断该场景；

// 【2】非线性处理
   // 1. 取kalman的权重能量最大的滤波块作为时延帧，获取参考信号X；
   // 2. 根据相关性对【1】中估计的回声进一步处理；（弱NLP）
   // 3. 计算基本信息；相干性，噪声估计，单双讲判断；
   // 4. 对回声信号进行谐波增强处理；并根据相干性对估计的回声进一步处理；
   // 5. 根据ERL对回声进行处理；
   // 6. 做后滤波处理；谐波增强降噪+倒谱降噪；
   // 7. 增益输出；

// [2.1] 根据kalman系数获取当前参考信号X
   iDex = 0;
   iMax = 0;
   for (i = 0; i < numFilterBlock; ++i) {
    sum = 0;
    for (j = 257Hz; j < 3600Hz; ++j) {
        sum += kalman->weight[i][2 * j] * kalman->weight[i][2 * j] + 
               kalman->weight[i][2 * j + 1] * kalman->weight[i][2 * j + 1];
    }
    if (sum > iMax) {
        iMax = sum;
        iDex = i;
    }
   }

// [2.2] 基础信息计算
    //    X: 远端频域数据(是根据kalman滤波系数获取到的吗？？？)
    //    Y: 近端频域数据
    //    errEst: 误差估计/残差信号
    //    echoEst: 回声估计
    
    //    X_power： 远端信号能量 X_power_smooth 能量平滑信号
    //    Y_power: 近端信号能量
    //    err_power: 误差信号能量
    //    echo_power: 回声信号能量
    
    //    corr_xy: 远近端互相关
    //    corr_xerr: 远端和残差互相关
   
  // [2.2.1] cohxd
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

  // 掐头去尾计算cohedMean,coheecMean, cohxdMean, cohxd600HzMean; bandsize = 800Hz, [300Hz, 3800Hz]
    // 可以用来辅助单双讲判断条件,从而控制回声估计对双讲的保护和单讲的抑制
    // 800Hz bandsize
    bandsize = 800Hz
   for (300Hz: 300Hz + bandsize) {
    sum += cohxd[i];
   }

   tmpSum = 0;
   for (i = 300Hz; i < 1500Hz; i++) {
    sum -= cohxd[i]; // 掐头
    sum += cohxd[i + bandsize]; // 增加新尾
    tmpSum += tmp;
   }
   cohxdMean = tmpSum / (800Hz - 1500Hz + 1);

    // 600Hz bandsize

   bandsize = 600Hz
  for (300Hz: 300Hz + bandsize) {
    sum += cohxd[i];
   }

   tmpCohxdMax600Hz = -100;
   tmpCohxdMin600Hz = 100;
   for (i = 300Hz; i < 3800Hz; i++) {
    tmp = sum - cohxd[i]; // 掐头
    tmp += cohxd[i + bandsize]; // 增加新尾
    tmpCohxdMax600Hz = max(tmpCohxdMax600Hz, tmp);
    tmpCohxdMin600Hz = min(tmpCohxdMin600Hz, tmp);
   }
   aec->cohxdMax600Hz = tmpCohxdMax600Hz / bandsize;
   aec->cohxdMin600Hz = tmpCohxdMin600Hz / bandsize;

    // [2.2.2]基频检测
    // errEst: 误差估计
    for (i = 0; i < 257; ++i) {
        cemp[i] = errEst[2 * i] * errEst[2 * i] + 
                  errEst[2 * i + 1] * errEst[2 * i + 1] + EPS;
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

    // [2.2.3] 近端，远端，误差信号能量平滑值均是 0.9
    // 求XdB（参考信号）, EdB(误差信号),能量范围250Hz~3600Hz

// [2.3] 弱NLP计算
   for (i = 0; i < 96; ++i) {
       gain = ((cohxd[i] - cohxe[i]) * 1.5 + 1.0) * cohxe[i];
       gain = max(gain, cohxe[i]);
       gain = min(gain, 0.95); // 0.95(0~3K)/0.98(3~8K) 不同频域设置不同值，频带越大设置越大
       gain = 1 - gain;
   }

// 【3. 相关性回声估计】
  // [3.1 计算相干性ERL]
  // coheec: 回声平滑估计和误差平滑估计相干性
  // cohed: 误差平滑估计和近端信号平滑估计相干性
  for(i = 0; i < 257; ++i) {
    aec->coheec[i] = min(1, aec->coheec[i]);
    aec->cohed[i] = min(1, aec->cohed[i]);
    aec->cohERL[i] = aec->coheec[i] * (1 - aec->cohed[i]);
  }

  // [3.2 回声谐波增强]
  // enhance
  for (i = 0; i < 257; ++i) {
    echoPow[i] = echoEst[2 * i] * echoEst[2 * i] + 
                  echoEst[2 * i + 1] * echoEst[2 * i + 1];
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
  for (i = 450Hz; i < 257; ++i) {
    echoEstPow[i] = max(enhance_echoPow[i] * 0.2, echoEstPow[i]); // 0.2可以根据回声情况调整;比如远端单讲就设置为1
  }

  for (i = 0; i < 257; ++i) {
    echoEstPow[i] =  echoEstPow[i] * aec->cohERL[i]; // 相干性ERL
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


// 【4. 计算ERL, SER 回声估计】 
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

    // ERL回声估计
    for (i = 0; i < 257; ++i) {
        X_pow[i] = X[2 * i] * X[2 * i] + X[2 * i + 1] + X[2 * i + 1]; // 对应时延帧出的参考信号
    }

    for (i = 0; i < 257; ++i) {
        tmp3 = X_pow[i] * smoothERL[i];
        Y_pow = echoPow[2 * i ] * echoPow[2 * i ] + echoPow[2 * i + 1] * echoPow[2 * i + 1]; // 谐波增强后的回声估计

        echoPow3[i] = tmp3 + Y_pow - 2 * sqrt(tmp3[i] * Y_pow);
        echoPow3[i] = min(echoPow3[i], tmp3);
        echoPow3[i] = max(echoPow3[i], 0);
    }

    if (aec->xbegin > 0) {
        for (i = 0; i < 257; ++i) {
            echoPow3[i] = max(echoPow3[i], 0.1 * X_pow[i]);
        }
    }

    // ERL, SER平滑系数更新; 需要根据单双讲，当前帧回声大小动态调整 
    sdr = (远端单讲 & xdB > 10 & fabs(xdbSmooth - xdB) < 15);
    alpha = sdr ? 0.9 : 1;
    alpha = (alpha == 0.9 && aec->SERcnt == 0) ? 0 : alpha; // 回声刚开始，滤波器还没收敛
    if (aec->noiseEstFlag == 0) { // 非回声
        alpha = 1
    } else {
        if (非双讲 & xdB > 10 & xbegin == 0) {
            alpha = 0.6；
        }
        else if(aec->SERcnt > 10) {
            tmp1 = ERLdB - ERLSmoothdB;
            tmp2 = SERdB - SERSmoothdB;
            if (tmp1 < 0 & alpha != 1) { // 远端回声
                alpha = 0.6;
            }
            else if (tmp2 > 5) {
                alpha = 1;
            }
            else if (tmp2 > 3) {
                alpha = max(0.98, alhap);
            }
        }
    }

    aec->SERCnt = alpha != 0? aec->SERCnt + 1 : aec->SERCnt;

    if (严格远端单讲) {
        for (0:400Hz) {
            ERLSmooth = SMOOTH(ERL, alpha);
            SERSmooth = SMOOTH(SER, alhap);
        }
    }
    for (400Hz:8K) {
        ERLSmooth = SMOOTH(ERL, alpha);
        SERSmooth = SMOOTH(SER, alhap);
    }

    // 求ERLdBsm, SERdBsm, ERLdB, SERdB
    // 可以根据ERLdBsm, SERdBsm, ERLdB, SERdB来作为单双讲的辅助判断依据
    ERLdBsm =PowerdB(ERLsm, 400, 2000);
    SERdBsm =PowerdB(SERsm, 400, 2000);

    ERLdB =PowerdB(ERL, 400, 2000);
    SERdB =PowerdB(SER, 400, 2000);

    // 严格远端单讲判断
    if (aec->SERcnt < 10) {
      ERLth = -30;
    } else {
      ERLth = ERLdBsm;
    }

    if (ERLdB > ERLdBsm & aec->F0 == false & 近端vad == 0 & XdB > 10 & ) {
       aec->STDR = true;
    }

 /////////////////////////[双讲场景回声估计系数调整]/////////////////////////////////////////////
    // 计算如果在双讲场景下，回声估计的调整系数
    // SER_2sm计算
    for (600Hz:2000Hz) {
        E1 += SERsmooth[i] * echoEstPow[i]; // 相关性回声估计
        E2 += echoPow3[i]; // ERL 回声估计
    }
    aec->DTSERCnt = (E1 > 2 * E2) ? aec->DTSERCnt + 1 : 0;
    SER_2sm = aec->DTSERCnt > 10 ? min(1.0, 2 * E2 / E1) : 1; 
    if (非双讲场景) { // 根据实际情况再考虑吧！！
        SER_2sm = 1;
    }

    // alpha1,2,3,4值计算
    if (frameNum == 1) {
        alpha1 = 1.0;
        alpha2 = 1.0;
        alpha3 = 1.0;
        alpha4 = 1.0;
    }

    if (严格远端单讲 & 非双讲 & xdB > max(10, xdBsm- xdB) & 非回声起始段) {
        F0 = 1200Hz;
        tmpPowSm = 0;
        for (1200Hz:6800Hz) {
            tmpPowSm[i]  = aec->se[i] - noiseEstimation[i]; // se:近端误差信号平滑能量
        }

        // alph1 1200~2400Hz
        ErrMeandB =  PowerdB(tmpPowSm, 1200Hz, 2400Hz);
        fXdB = PowerdB(XPow, 1200, 2400)
        aec->alph1 = if (fXdB > 10 ){
            if (ErrMeandB > aec->smoothdB) {
            aec->smoothdB = 0.9 * aec->smoothdB + 0.1 * fXdB;
            } else {
            aec->smoothdB = 0.99 * aec->smoothdB + 0.01 * fXdB;
            }
            alpha = 0.2 * Max(5, aec->smoothdB - 35);
        }

        // alpha2 2400~3000Hz
        ErrMeandB =  PowerdB(tmpPowSm, 2400Hz, 3000Hz);
        fXdB = PowerdB(XPow, 2400Hz, 3000Hz)
        aec->alph2 = if (fXdB > 10 ){
            if (ErrMeandB > aec->smoothdB) {
            aec->smoothdB = 0.9 * aec->smoothdB + 0.1 * fXdB;
            } else {
            aec->smoothdB = 0.99 * aec->smoothdB + 0.01 * fXdB;
            }
            alpha = 0.2 * Max(5, aec->smoothdB - 35);
        }

        // alpha3 3000~3600Hz
        ErrMeandB =  PowerdB(tmpPowSm, 3000Hz, 3000Hz);
        fXdB = PowerdB(XPow, 3000Hz, 3000Hz)
        aec->alph3 = if (fXdB > 10 ){
            if (ErrMeandB > aec->smoothdB) {
            aec->smoothdB = 0.9 * aec->smoothdB + 0.1 * fXdB;
            } else {
            aec->smoothdB = 0.99 * aec->smoothdB + 0.01 * fXdB;
            }
            alpha = 0.2 * Max(5, aec->smoothdB - 35);
        }

        // alpha4 3600~6800Hz
        ErrMeandB =  PowerdB(tmpPowSm, 3600Hz, 6800Hz);
        fXdB = PowerdB(XPow, 3600Hz, 6800Hz)
        aec->alph4 = if (fXdB > 10 ){
            if (ErrMeandB > aec->smoothdB) {
            aec->smoothdB = 0.9 * aec->smoothdB + 0.1 * fXdB;
            } else {
            aec->smoothdB = 0.99 * aec->smoothdB + 0.01 * fXdB;
            }
            alpha = 0.2 * Max(5, aec->smoothdB - 35);
        }
    }

    // 13.2 调整双讲场景回声估计系数afalph的系数alpha
    // a. 300Hz以下根据F0设置afAlpah的值
    // b. 300, 600, 1200, 2400, 3000, 3600, 6800
    alpha = 0.5 * SER_2sm
    for (300Hz: 600Hz) {
    aec->afAlpha[i] = min(1, alpha * SERSmooth[i]);
    }

    alpha = 1 * SER_2sm
    for (600Hz: 1200Hz) {
    aec->afAlpha[i] = min(1, alpha * SERSmooth[i]);
    }

    alpha = aec->alpha1 * SER_2sm
    for (1200Hz: 2400Hz) {
    aec->afAlpha[i] = min(1, alpha * SERSmooth[i]);
    }

    alpha = aec->alpha2 * SER_2sm
    for (2400Hz: 3000Hz) {
    aec->afAlpha[i] = min(1, alpha * SERSmooth[i]);
    }

    alpha = aec->alpha3 * SER_2sm
    for (3000Hz: 3600Hz) {
    aec->afAlpha[i] = min(1, alpha * SERSmooth[i]);
    }

    alpha = aec->alpha4 * SER_2sm
    for (3600Hz: 6800Hz) {
    aec->afAlpha[i] = min(1, alpha * SERSmooth[i]);
    }
///////////////////////////////////////////////////////////////////////////////////////////////

// 5.3 根据计算的coheec,计算conheecMeanCnt，作为双讲回声估计系数调整的重要依据
   cnt = 0
   for (100Hz: 3800Hz) {
    cnt = (aec->coheec[i] > 0.4) ? cnt+1:cnt;
   }    
   aec->coheecMeanCnt = Max(0, aec->coheecMeanCnt - 1);

   if (aec->coheecMean > 0.2 || cnt > 40) {
      aec->coheecMeanCnt = Max(aec->coheecMeanCnt, 30);
   } else if (aec->coheecMean > 0.15 || cnt > 30) {
      aec->coheecMeanCnt = Max(aec->coheecMeanCnt, 20);
   } else if (aec->coheecMean > 0.06 || cnt > 20) {
      aec->coheecMeanCnt = Max(aec->coheecMeanCnt, 10);
   } else if (aec->coheecMean > 0.015 || cnt > 10) {
      aec->coheecMeanCnt = Max(aec->coheecMeanCnt, 5);
   }

    // 相干性回声估计和ERL回声估计结果融合
    if (远端单讲 & aec->SERCnt > 20) {
        for (600Hz:8KHz) {
            echoPow[i] = max(echoPow[i], echoPow3[i]);
        }
    }
    if (aec->xbegin > 0) {
        for (0: 8KHz) {
            echoPow[i] = max(echoPow[i], echoPow3[i]);
        }
    }

// 【3. 回声增益计算】
// [3.1] 谐波增强消除回声
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

// [3.2 求Gain]
   Gain(snrEchoEst, gain);
   for (i = 0; i < 257; ++i) {
        errEst[2 * i] *= gain[i];
        errEst[2 * i + 1] *= gain[i];
   }

// [3.3 倒谱消回声]
  // 分成0~4K低频做和4~8K
  // cemp8k: 0~4K倒谱
  // cemp16k: 4~8K倒谱
  // errEst误差谱

    for (i = 0; i < 128; ++i) {
        tmpPow[i] = errEst[2 * i] * errEst[2 * i] + 
                    errEst[2 * i + 1] * errEst[2 * i + 1] + EPS;
    }

    tmpPow = Log(tmpPow);
    for (i = 127; i >= 0; --i) {
        tmpPow[2 * i ] = tmpPow[i];
        tmpPow[2 * i + 1] = 0;

        cemp16k[2 * i] = tmpPow[i];
        cemp16k[2 * i + 1] = 0;
    }
    tmpPow[256] = tmpPow[254] - 4;
    tmpPow[255] = 0;

    cemp8k = IFFT(tmpPow);

    // 16k
    for (i = 128; i < 257; ++i) {
      tmpPow[i - 128] = errEst[2 * i] * errEst[2 * i] + 
            errEst[2 * i + 1] * errEst[2 * i + 1];
    }

    tmpPow = Log(tmpPow + 128);

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
    for (i = 0; i < 128; ++i) {
        tmp1 = cemp8k[i];
        cempMod = smooth_cemp8k[i];
        cempModAmp = cempMod;
        if (tmp1 > 0.36 && i >= 15) {
            smooth_cemp8k[i] = 0;
            cempMod = tmp1;
            cempModAmp = cempMod * 5;
        }
        cempMod[i] = cempMod; // 2~4k
        cempModAmp[i] = cempModAmp; // 0~2k

        cempMod[257 - i] = cempMod;
        cempModAmp[257 - i] = cempModAmp;
    }

    cempMod = FFT(cempMod); // 0 ~ 256
    cempModAmp = FFT(cempModAmp); // 0 ~ 256
    
    idx = 128 / 2;
    for (i = 0; i < idx; ++i) {
      epxCemp[i] = cempModAmp[2 * i];
      epxCemp[i + idx] = cempMod[2 * (i + idx)];
    }

    epxCemp = Exp(epxCemp);
    for (i = 0; i < 128; ++i) {
        snr = epxCemp[i];
        snrGain[i] = snr / (echoEst + EPS);
    }

    // 16K
    for (i = 0; i < 16; ++i) {
      tmpR = smooth_cemp16k[2 * i] * alpha[i] + cemp16k[2 * i] * (1.0 - alpha[i]);
      tmpI = smooth_cemp16k[2 * i + 1] * alpha[i] + cemp16k[2 * i + 1] * (1.0 - alpha[i]);

      smooth_cemp16k[2 * i] = tmpR;
      smooth_cemp16k[2 * i + 1] = tmpI;

      // 下面FFT需要继续使用
      cemp16k[2 * i] = tmpR;
      cemp16k[2 * i + 1] = tmpI;
      cemp16k[514 - 2 * i] = tmpR;
      cemp16k[514 - (2 * i + 1)] = tmpI;
    }

    alpha1 = alpha[15];
    for (i = 32; i < 257; ++i) {
        smooth_cemp16k = SMOOTH(cemp16k, alpha1);
        cemp16k[i] = smooth_cemp16k[i];
        cemp16k[514 - 2 * i] = cemp16k[i];
    }

    cemp16k = FFT(cemp16k);

    for (i = 128; i < 256; ++i) {
      epxCemp[i] = smooth_cemp16k[2 * i];
    }
    epxCemp = EXP(epxCemp + 128);
    for (i = 128; i < 256; ++i) {
        snr = epxCemp[i];
        snrGain[i + 128] = snr / (echoEst + EPS);
    }

// [3.4 求倒谱降噪输出增益]
    Gain(snrGain, gain);

// [3.3 根据snr计算回声最终增益值] 
void Gain(float *snrEchoEst, float *gain)
{
    for (i = 0; i < 257; ++i) {
        tmp1 = snrEchoEst[i];
        tmp2 = tmp1 + 1;
        gain[i] = tmp1 / tmp2;
        gain[i] = min(gain[i], aec->gainMinLimit);
    }
}


// 【5. 低频指导高频】
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
