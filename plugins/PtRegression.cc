#include "NtuplizerGenPart/L1ScoutingAnalyzer/interface/PtRegression.h"

// Sigmoid function
inline float sigmoid(float x){
    return (1.0/(1.0 + std::exp(-1.*x)));
}


// Constructor
PtRegression::PtRegression(std::string model_path){
    XGBoosterCreate(NULL, 0, &booster_);
    XGBoosterLoadModel(booster_, model_path.c_str());
}

// Destructor
PtRegression::~PtRegression(){}

// Get regressed pt
float PtRegression::get_regressed_pt(std::vector<float> features){
    float result;
    float values[1][features.size()];
    int ivar=0;

    for(auto& var: features){
        values[0][ivar] = var;
        ivar++;
    }
    DMatrixHandle dvalues;
    XGDMatrixCreateFromMat(reinterpret_cast<float*>(values), 1, features.size(), -9999., &dvalues);

    // Output prediction
    bst_ulong out_dim;

    float const *out_result = NULL;
    auto ret = XGBoosterPredict(booster_, dvalues, 0, 0, 0, &out_dim, &out_result);

    XGDMatrixFree(dvalues);

    if(ret == 0){
        result = out_result[0];
    }

    return result;
}