#ifndef ptregression_h
#define ptregression_h

// Standard libraries
#include <vector>
#include <string>
#include <cmath>
#include <iostream>

// XGBoost
#include <xgboost/c_api.h>


class PtRegression{
    public:
        PtRegression();
        PtRegression(std::string model_path);
        ~PtRegression();
        float get_regressed_pt(std::vector<float> features);

    private:
        BoosterHandle booster_;
};

#endif