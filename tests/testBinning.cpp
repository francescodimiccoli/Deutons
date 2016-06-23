#include <cmath>
#include <limits.h>
#include "gtest/gtest.h"
#include "include/binning.h"

class BinningTest : public testing::Test { 
protected:
    Binning p_bin, d_bin, he_bin, he3_bin;

    BinningTest():
              p_bin(0.938, 1, 1), 
              d_bin(1.875, 1, 2), 
             he_bin(3.727, 2, 4), 
            he3_bin(2.814, 2, 3) {}

    virtual void SetUp() { }
    virtual void TearDown() { }
};

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

TEST_F(BinningTest, TestBinsCopyConstr) { 
    p_bin.setBinsFromEk(10,1,10);

    Binning copy(p_bin);

    EXPECT_FLOAT_EQ( 0.938, copy.getParticle().getMass() );
    EXPECT_FLOAT_EQ(     1, copy.getParticle().getZ()    );
    EXPECT_FLOAT_EQ(     1, copy.getParticle().getA()    );

    EXPECT_EQ(10, copy.size());
    for(int i = 0; i < 10; i++){
        EXPECT_EQ(p_bin.EkBins  ()[i], copy.EkBins  ()[i]);
        EXPECT_EQ(p_bin.EtotBins()[i], copy.EtotBins()[i]);
        EXPECT_EQ(p_bin.MomBins ()[i], copy.MomBins ()[i]);
        EXPECT_EQ(p_bin.RigBins ()[i], copy.RigBins ()[i]);
        EXPECT_EQ(p_bin.BetaBins()[i], copy.BetaBins()[i]);
    }
}



float R_from_beta(float beta, float mq) { return mq * beta / sqrt(1 - beta * beta);  }
float beta_from_R(float R,    float mq) { return R / sqrt(R * R + mq * mq);  }
float Ek_from_beta(float beta, float m, float  n) { return (m / sqrt(1 - beta * beta) - m) / n; }


//float Ek_from_R(float R, int z, float m) { return sqrt( z*z*R*R + m*m ) - m; }
//TEST_F(BinningTest, TestEkToRigidity) { 
//      p_bin.setBinsFromEk(2,1,10);
//      d_bin.setBinsFromEk(2,1,10);
//     he_bin.setBinsFromEk(2,1,10);
//    he3_bin.setBinsFromEk(2,1,10);
//    
//    EXPECT_FLOAT_EQ( 1, Ek_from_R(  p_bin.RigBins()[0], 1, 0.938) );
//    EXPECT_FLOAT_EQ( 1, Ek_from_R(  d_bin.RigBins()[0], 1, 1.875) );
//    EXPECT_FLOAT_EQ( 1, Ek_from_R( he_bin.RigBins()[0], 2, 3.727) );
//    EXPECT_FLOAT_EQ( 1, Ek_from_R(he3_bin.RigBins()[0], 2, 2.814) );
//    
//    EXPECT_FLOAT_EQ( 10, Ek_from_R(  p_bin.RigBins()[1], 1, 0.938) );
//    EXPECT_FLOAT_EQ( 10, Ek_from_R(  d_bin.RigBins()[1], 1, 1.875) );
//    EXPECT_FLOAT_EQ( 10, Ek_from_R( he_bin.RigBins()[1], 2, 3.727) );
//    EXPECT_FLOAT_EQ( 10, Ek_from_R(he3_bin.RigBins()[1], 2, 2.814) );
//}
