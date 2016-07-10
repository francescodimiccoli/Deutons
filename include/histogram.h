#ifndef HISTOGRAM_H
#define HISTOGRAM_H

#include "printmatrix.h"
#include <string>
#include <stdexcept>
#include <algorithm>

/** A histogram is a combination of two vectors, edges and contents
 * Edges is sorted, contents has one element less (both are checked by CheckIntegrity)
 *
 * */

class Histogram {
public:
    std::vector<float> getEdges()   {return edges  ;}
    std::vector<float> getContent() {return content;}

    inline Histogram(std::vector<float> edges);
    inline Histogram(std::vector<float> edges, std::vector<float> content);
    inline float at(float x) const;
    inline void fill(float x, float w=1.0);
    inline void fillWithCounts (std::vector<float> counts);
    inline void normalize();

    inline float integrate(float min, float max) const;
    inline Histogram reBin(std::vector<float> targetEdges);
    inline void printContent() const;


private:
    std::vector<float> edges ;
    std::vector<float> content ;
};


Histogram::Histogram(std::vector<float> e): edges(e), content(e.size()-1, 0) 
{
    if( !std::is_sorted(edges.begin(),edges.end()) ) 
        throw std::invalid_argument("Histogram::Histogram sorted bin edges are expected.");
}

Histogram::Histogram(std::vector<float> e, std::vector<float> c): edges(e), content(c) 
{
    if( !std::is_sorted(edges.begin(),edges.end()) ) 
        throw std::invalid_argument("Histogram::Histogram sorted bin edges are expected.");
    if( edges.size() != content.size() + 1)            
        throw std::length_error("Histogram::Histogram contents vector length should correspond to the number of bins.");
}

float Histogram::at(float x) const {
    if( x < edges.front() ) return 0;
    if( edges.back() > x  ) return 0;
    for( int i = 0; i < content.size(); i++ ) {
        if(edges[i] < x && x < edges[i+1]) return content[i];
    }
}

void Histogram::fill(float x, float w) {
    if( x < edges.front() ) return;
    if( edges.back() > x  ) return;
    for( int i = 0; i < content.size(); i++ ) {
        if(edges[i] < x && x < edges[i+1]) {
            content[i] += w / (edges[i+1] - edges[i]);
            return;
        } 
    }
}

void Histogram::fillWithCounts(std::vector<float> counts) {
    if( edges.size() != counts.size() + 1)            
        throw std::length_error("Histogram::fillWithCounts counts vector length should correspond to the number of bins.");
    for( int i = 0; i < content.size(); i++ ) 
        content[i] = counts[i] / (edges[i+1] - edges[i]);
}

void Histogram::normalize() {
    float integral = 0;
    for(int i = 0; i < content.size(); i++)  
        integral += content[i] * (edges[i+1] - edges[i]);
    for (int i = 0; i < content.size(); i++)  content[i] /= integral;
}

float Histogram::integrate(float min, float max) const {
    // This is for integrals with "backward" limits
    int sign = 1;
    if( min > max ) { std::swap(min,max); sign = -1; }

    float integral = 0;
    for( int i = 0; i < content.size(); i++ ) {
        if( edges[i+1] < min || max < edges[i] ) // bin is not inside the range, skip
            continue;
        else if( edges[i] <= min && max <= edges[i+1] ) // range is fully inside the bin, break 
            return sign * content[i] * (max - min);
        else if( min <= edges[i] && edges[i+1] <= max ) // bin is fully inside the range, add
            integral += content[i] * (edges[i+1] - edges[i]);
        else if( edges[i] <= min && edges[i+1] <= max ) // lower bound is inside the bin, add
            integral += content[i] * (edges[i+1] - min);
        else if( min <= edges[i] && max <= edges[i+1] ) // upper bound is inside the bin, add 
            integral += content[i] * (max - edges[i]);
    }
    return sign * integral; 
}


inline Histogram Histogram::reBin(std::vector<float> targetEdges){
    std::vector<float> newContent(targetEdges.size() - 1, 0);
    for( int i = 0; i < targetEdges.size() - 1; i++ ) {
        float dx = targetEdges[i+1] - targetEdges[i];
        newContent[i] = integrate(targetEdges[i], targetEdges[i+1]) / dx; 
    }
    return Histogram(targetEdges, newContent);
}


void Histogram::printContent() const { printMatrix::print( {edges, content} ); }

#endif
