#include <cmath>
#include <limits.h>
#include "gtest/gtest.h"
#include "include/binning.h"

class BinningTest : public testing::Test { 
protected:
    Binning p_bin, d_bin, he_bin, he3_bin;

    BinningTest():
              p_bin(0.938, 1, 1), 
              d_bin(1.875, 2, 1), 
             he_bin(3.727, 2, 4), 
            he3_bin(2.814, 2, 3) {}

    virtual void SetUp() { }
    virtual void TearDown() { }
};

float R_from_beta(float beta, float mq) { return mq * beta / sqrt(1 - beta * beta);  }
float beta_from_R(float R,    float mq) { return R / sqrt(R * R + mq * mq);  }
float Ek_from_beta(float beta, float m, float  n) { return (m / sqrt(1 - beta * beta) - m) / n; }

TEST_F(BinningTest, TestBinsLogarithmic) { 
    p_bin.setBinsFromEk(10,1,10);

    EXPECT_EQ(10, p_bin.size());

    float dlog = log(p_bin.EkBins()[1]) - log(p_bin.EkBins()[0]);
    for(int i = 2; i < 10; i++){
        float current = log(p_bin.EkBins()[i]) - log(p_bin.EkBins()[i-1]);
        float tolerance = fabs( dlog / current - 1); 
        EXPECT_LE(tolerance, 0.001 );
    }
}

TEST_F(BinningTest, TestBinsParticle) { 
    Binning   p_bin(0.938, 1, 1); 
    Binning   d_bin(1.875, 1, 2); 
    Binning  he_bin(3.727, 2, 4); 
    Binning he3_bin(2.814, 2, 3);

    EXPECT_FLOAT_EQ( 0.938,   p_bin.getParticle().getMass() );
    EXPECT_FLOAT_EQ( 1.875,   d_bin.getParticle().getMass() );
    EXPECT_FLOAT_EQ( 3.727,  he_bin.getParticle().getMass() );
    EXPECT_FLOAT_EQ( 2.814, he3_bin.getParticle().getMass() );
    
    EXPECT_EQ( 1,   p_bin.getParticle().getZ() );
    EXPECT_EQ( 1,   d_bin.getParticle().getZ() );
    EXPECT_EQ( 2,  he_bin.getParticle().getZ() );
    EXPECT_EQ( 2, he3_bin.getParticle().getZ() );
    
    EXPECT_EQ( 1,   p_bin.getParticle().getA() );
    EXPECT_EQ( 2,   d_bin.getParticle().getA() );
    EXPECT_EQ( 4,  he_bin.getParticle().getA() );
    EXPECT_EQ( 3, he3_bin.getParticle().getA() );
}
