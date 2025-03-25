// interface/QCDReweightInfo.h
#ifndef L1SCOUTINGANALYZER_QCDREWEIGHTINFO_H
#define L1SCOUTINGANALYZER_QCDREWEIGHTINFO_H

struct QCDReweightInfo {
  float weight;      // legacy weight (if needed)
  float weight_v2;   // computed weight from QCDWeightCalc
};

#endif // L1SCOUTINGANALYZER_QCDREWEIGHTINFO_H