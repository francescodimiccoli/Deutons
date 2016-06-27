#include <cmath>
#include <limits.h>
#include "gtest/gtest.h"
#include "include/histogram.h"

class HistogramTest : public testing::Test { 
protected:
    Histogram hist;
    HistogramTest() { }
    virtual void SetUp() { }
    virtual void TearDown() { }
};

TEST_F(HistogramTest, TestLogNormGeneration) { 

    // I want three bins: 1--2----4------------8
    hist.genMCLogNormFlux(3, 1, 8);

    std::vector<float> edges   = hist.getEdges();
    std::vector<float> content = hist.getContent();

    EXPECT_EQ( 4, edges.size() );
    EXPECT_FLOAT_EQ( 1, edges[0] );
    EXPECT_FLOAT_EQ( 2, edges[1] );
    EXPECT_FLOAT_EQ( 4, edges[2] );
    EXPECT_FLOAT_EQ( 8, edges[3] );

    EXPECT_FLOAT_EQ( 1, content[0] );
    EXPECT_FLOAT_EQ( 1, content[1] );
    EXPECT_FLOAT_EQ( 1, content[2] );
}


TEST_F(HistogramTest, TestNormalize1) { 

    // I want three bins: 1--2----4------------8
    hist.genMCLogNormFlux(3, 1, 8);
    hist.normalize();

    std::vector<float> content = hist.getContent();
    // The integral should be 1 + 2 + 4 = 7
    EXPECT_FLOAT_EQ( 1.0 / 7.0, content[0] );
    EXPECT_FLOAT_EQ( 1.0 / 7.0, content[1] );
    EXPECT_FLOAT_EQ( 1.0 / 7.0, content[2] );
}
