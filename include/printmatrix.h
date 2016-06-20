#ifndef PRINTMAT_H
#define PRINTMAT_H

#include <vector>
#include <string>
#include <iostream>
#include <iomanip>


class printMatrix {

   public:

      static void print (std::vector<std::vector <float>> matrix)
      {
         uint ncolumns = matrix.size();
         uint nlines=0;
         for (auto column : matrix)
            if (column.size() > nlines)  nlines = column.size();

         printPadding(ncolumns);

         for (uint iline=0; iline<nlines; iline++) {
            std::cout << std::right << std::setw (4)  << std::setfill (' ') << iline;
            for (uint icol=0; icol<ncolumns; icol++) {
               std::cout << std::right << std::setw (11) << std::setfill (' ');
               if (matrix[icol].size() >iline)
                  std::cout << std::setprecision (4) <<   matrix[icol][iline] ;
               else
                  std::cout << " " ;
            }
            std::cout << std::endl;
         }
         printPadding(ncolumns);
         return;
      }

      static void printTransposed(std::vector<std::vector <float>> matrix)
      {
         uint ncolumns = matrix.size();
         uint nlines=0;
         for (auto column : matrix)
            if (column.size() > nlines)  nlines = column.size();

         printPadding(nlines);

         for (uint icol=0; icol<ncolumns; icol++) {
            std::cout << std::right << std::setw (4)  << std::setfill (' ') << icol;
                for (uint iline=0; iline<nlines; iline++) {     
               std::cout << std::right << std::setw (11) << std::setfill (' ');
               if (matrix[icol].size() >iline)
                  std::cout << std::setprecision (4) <<   matrix[icol][iline] ;
               else
                  std::cout << " " ;
            }
            std::cout << std::endl;
         }
         printPadding(nlines);
         return;
      }

      static void print ( std::vector< std::vector<float> > matrix, std::vector<std::string> titles)
      {
         printPadding(titles.size());
         std::cout << std::right << std::setw (4)  << std::setfill (' ') << "Bin";
         for (auto title : titles)
            std::cout << std::right << std::setw (11) << std::setfill (' ') << title;
         std::cout << std::endl;
         print (matrix);
      }

   private:
      static void printPadding (int ncolumns) {
         int nchar=4 + 11*ncolumns + 1;
         std::cout << std::string (nchar, '-') << std::endl;
         return;
      }
};


#endif
