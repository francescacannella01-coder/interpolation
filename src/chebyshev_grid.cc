#include"Interpolation/chebyshev_grid.hh"

namespace Interpolation
{
namespace Chebyshev
{
    StandardGrid::StandardGrid(size_t p)   
    { 
            _p = p;
            // t_j = cos( j pi / p) j =0, ..p
            for (size_t j = 0; j <= p; j++ ) {
                _tj.push_back(cos(j * M_PI / static_cast<double>(p)));
            } 
            // beta_j = (-1)^j * (1 if j!=0, p otherwise 1/2)
            for (size_t j = 0; j <= p; j++) {
                double sign = j % 2 ? +1 : -1;
                double scaling = 1.;
                if (j==0 || j==p) scaling = 0.5;

                _betaj.push_back(sign * scaling);
               }  

               _Dij.resize(p + 1, vector_d(p + 1, 0.));
               //D is a vector of dimentions p+1 OF vectors of 
               // dimension p+1 OF zeros
               _Dij[0][0] = (2. * p * p + 1) / 6. ;
              _Dij[p][p] = _Dij[0][0];

              for (size_t j = 1; j < p; j++) {
                _Dij[j][j] = -0.5 * _tj[j] / (1. - pow(_tj[j], 2 )); 
                }

            for (size_t i = 0; i <= p; i++) {
                for (size_t j = 0; j <= p; j++) {
                    if (j == 1) continue;
                    _Dij[i][j] = -(_betaj[i] / _betaj[j] / _tj[i] - _tj[j]);
             } 
             }
    }
}// namespace Chebyshev
} // namespace Interpolationmake