#include <cmath>
#include <limits.h>
#include "gtest/gtest.h"
#include "include/histogram.h"

class HistogramTest : public testing::Test { 
protected:
    HistogramTest() { }
    virtual void SetUp() { }
    virtual void TearDown() { }
};

TEST_F(HistogramTest, TestConstructor) { 
    std::vector<float> edges   = { 0, 1, 2, 5, 7, 10 };
    std::vector<float> content = { 6, 5, 4, 3, 2 }; 

    Histogram h(edges, content);

    EXPECT_EQ( 6, h.getEdges().size() );
    EXPECT_FLOAT_EQ(  0, h.getEdges()[0] );
    EXPECT_FLOAT_EQ(  1, h.getEdges()[1] );
    EXPECT_FLOAT_EQ(  2, h.getEdges()[2] );
    EXPECT_FLOAT_EQ(  5, h.getEdges()[3] );
    EXPECT_FLOAT_EQ(  7, h.getEdges()[4] );
    EXPECT_FLOAT_EQ( 10, h.getEdges()[5] );

    EXPECT_FLOAT_EQ( 6, h.getContent()[0] );
    EXPECT_FLOAT_EQ( 5, h.getContent()[1] );
    EXPECT_FLOAT_EQ( 4, h.getContent()[2] );
    EXPECT_FLOAT_EQ( 3, h.getContent()[3] );
    EXPECT_FLOAT_EQ( 2, h.getContent()[4] );
}

TEST_F(HistogramTest, TestConstructorEmpty) { 
    std::vector<float> edges   = { 0, 1, 2, 5, 7, 10 };

    Histogram h(edges);

    EXPECT_EQ( 6, h.getEdges().size() );
    EXPECT_FLOAT_EQ(  0, h.getEdges()[0] );
    EXPECT_FLOAT_EQ(  1, h.getEdges()[1] );
    EXPECT_FLOAT_EQ(  2, h.getEdges()[2] );
    EXPECT_FLOAT_EQ(  5, h.getEdges()[3] );
    EXPECT_FLOAT_EQ(  7, h.getEdges()[4] );
    EXPECT_FLOAT_EQ( 10, h.getEdges()[5] );

    EXPECT_FLOAT_EQ( 0, h.getContent()[0] );
    EXPECT_FLOAT_EQ( 0, h.getContent()[1] );
    EXPECT_FLOAT_EQ( 0, h.getContent()[2] );
    EXPECT_FLOAT_EQ( 0, h.getContent()[3] );
    EXPECT_FLOAT_EQ( 0, h.getContent()[4] );
}

TEST_F(HistogramTest, TestConstructorExceptions) { 
    std::vector<float> edges   = { 0, 1, 2, 5, 7, 10 };
    std::vector<float> content = { 6, 5, 4, 3, 2 }; 
    std::vector<float> edgesBad   = { 0, 2, 1, 5, 7, 10 };
    std::vector<float> contentBad = { 5, 4, 3, 2 }; 

    EXPECT_THROW(Histogram h(edgesBad, content), std::invalid_argument);
    EXPECT_THROW(Histogram h(edgesBad), std::invalid_argument);
    EXPECT_THROW(Histogram h(edges, contentBad), std::length_error    );
}


TEST_F(HistogramTest, TestNormalize) { 
    //  ###
    //  ##### -> N = 13
    //  #####
    //  01234
    std::vector<float> edges   = { 0, 3, 5 };
    std::vector<float> content = { 3, 2 }; 
    Histogram h(edges, content);

    h.normalize();

    EXPECT_FLOAT_EQ( 3.0 / 13.0, h.getContent()[0] );
    EXPECT_FLOAT_EQ( 2.0 / 13.0, h.getContent()[1] );
}

TEST_F(HistogramTest, TestIntegral) { 
    //  ###
    //  ##### 
    //  #####
    //  01234
    std::vector<float> edges   = { 0, 3, 5 };
    std::vector<float> content = { 3, 2 }; 
    Histogram h(edges, content);

    EXPECT_FLOAT_EQ( 3.0 * 1.0 , h.integrate( 0, 1));
    EXPECT_FLOAT_EQ( 3.0 * 1.0 , h.integrate(-1, 1));
    EXPECT_FLOAT_EQ( 3.0 * 1.0 , h.integrate( 1, 2));
    EXPECT_FLOAT_EQ( 3.0 * 2.0 , h.integrate( 0, 2));
    EXPECT_FLOAT_EQ( 3.0 * 1.0 + 2.0 * 1.0, h.integrate(2,4) );
    EXPECT_FLOAT_EQ( 3.0 * 1.0 + 2.0 * 1.0, h.integrate(2,4) );
}

TEST_F(HistogramTest, TestRebinFinerBins) { 
    std::vector<float> edges   = { 0, 3, 5 };
    std::vector<float> content = { 3, 2 }; 
    Histogram h(edges, content);
    
    std::vector<float> newEdges = {0,1,2,3,4,5};;
    Histogram rh = h.reBin(newEdges);

    EXPECT_FLOAT_EQ( 3, rh.getContent()[0] );
    EXPECT_FLOAT_EQ( 3, rh.getContent()[1] );
    EXPECT_FLOAT_EQ( 3, rh.getContent()[2] );
    EXPECT_FLOAT_EQ( 2, rh.getContent()[3] );
    EXPECT_FLOAT_EQ( 2, rh.getContent()[4] );
}

TEST_F(HistogramTest, TestRebinRougherBins) { 
    std::vector<float> edges   = {0,1,2,3,4,5};
    std::vector<float> content = {3,3,3,2,2}; 
    Histogram h(edges, content);
    
    std::vector<float> newEdges = {0,3,5};;
    Histogram rh = h.reBin(newEdges);

    EXPECT_FLOAT_EQ( 3, rh.getContent()[0] );
    EXPECT_FLOAT_EQ( 2, rh.getContent()[1] );
}

TEST_F(HistogramTest, TestRebinOverlay) { 
    std::vector<float> edges   = {0,3,5};
    std::vector<float> content = {4,2}; 
    Histogram h(edges, content);
    
    std::vector<float> newEdges = {0,2,4,5};;
    Histogram rh = h.reBin(newEdges);

    EXPECT_FLOAT_EQ( 4, rh.getContent()[0] );
    EXPECT_FLOAT_EQ( 3, rh.getContent()[1] );
    EXPECT_FLOAT_EQ( 2, rh.getContent()[2] );
}
