#include "Interpolation/chebyshev_grid.hh"
#include <stdexcept>

namespace Interpolation
{
namespace Chebyshev
{
    /*
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
             */

double StandardGrid::interpolate(double t, const vector_d &fj, size_t start, size_t end) const
{
    if(t<-1 || t>1)
    {
        throw std::domain_error("StandardGrid:: Interpolate t must be in [-1,1]");
    }
    if (end - start != _p){
        throw std::domain_error("StandardGrid:: interpolate end-start should be = to p");
    }
    double den=0;
    for (size_t j=0; j<= _p; j++){
        //if( t==_tj[j]) retirn fj[j+start];
        if(std::abs(t-_tj[j]) < 1.0e-15) return fj[j+start];
        den += _betaj[j]/ (t-_tj[j]);
    }
    double res =0.;
    for(size_t i=0; i <= _p; i++){
        res += poli_weight(t,i,den) * fj[i + start];
    }
    return res;
}
    //definiamo i pesi//
double StandardGrid::poli_weight(double t, size_t j, double den) const
{
        if(std::abs(t-_tj[j]) < 1.0e-15) return 1.;
         double res=0.;
         res = _betaj[j] / (t-_tj[j]) / den;
         return res;
}
double StandardGrid::poli_weight(double t, size_t j) const
{
    if(std::abs(t-_tj[j]) < 1.0e-15) return 1.;
    double den= 0.;
    for(size_t j=0; j<= _p; j++){
        if(std::abs(t-_tj[j])<1.0e-15) return 0.;
        den+= _betaj[j] / (t-_tj[j]);
    }
     double res=0.;
     res = _betaj[j] / (t-_tj[j]) /den;
     return res;

}

vector_d StandardGrid::discretize(const std::function<double(double)> &fnc) const
{
     vector_d fj(_p+1, .0);
     for(size_t i=0; i<=_p; i++){
        fj[i] = fnc(_tj[i]);
     }
     return fj;
}

StandardGrid::StandardGrid(size_t p) : _p(p) // <-- Inizializza _p qui
{
    // 1. Inizializzazione dei nodi _tj
    _tj.reserve(_p + 1);
    for (size_t j = 0; j <= _p; j++) {
        _tj.push_back(std::cos(j * M_PI / static_cast<double>(_p)));
    }

    // 2. Inizializzazione dei pesi _betaj (Eq. 2.5 del paper)
    _betaj.reserve(_p + 1);
    for (size_t j = 0; j <= _p; j++) {
        double sign = (j % 2 == 0) ? 1.0 : -1.0; // (-1)^j
        double delta = (j == 0 || j == _p) ? 0.5 : 1.0;
        _betaj.push_back(sign * delta);
    }

    // 3. Inizializzazione della matrice _Dij (codice del professore)
    _Dij.resize(_p + 1, vector_d(_p + 1, 0.));
    
    _Dij[0][0] = (2. * _p * _p + 1.) / 6.;
    _Dij[_p][_p] = -_Dij[0][0];

    for (size_t j = 1; j < _p; j++) {
        _Dij[j][j] = -0.5 * _tj[j] / (1. - std::pow(_tj[j], 2));
    }

    for (size_t i = 0; i <= _p; i++) {
        for (size_t j = 0; j <= _p; j++) {
            if (i == j) continue;
            // Calcolo fuori diagonale usando i membri _betaj e _tj
            _Dij[i][j] = (_betaj[i] / _betaj[j]) / (_tj[i] - _tj[j]);
        }
    }
}

}// namespace Chebyshev
} // namespace Interpolationmake