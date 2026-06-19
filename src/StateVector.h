#pragma once

struct StateVector {
    double rho      ;
    double rho_u    ;
    double rho_e    ;

    //   A + B  entre deux StateVector
    StateVector operator+(const StateVector& other) const {
        return { rho + other.rho, rho_u + other.rho_u, rho_e + other.rho_e };
    }
    //  A - B
    StateVector operator-(const StateVector& other) const {
        return { rho - other.rho, rho_u - other.rho_u, rho_e - other.rho_e };
    }
    //  U * s   (vecteur fois scalaire)
    StateVector operator*(double s) const {
        return { rho * s, rho_u * s, rho_e * s };
    }
};
// autoriser aussi  s * U  (scalaire à gauche) :
inline StateVector operator*(double s, const StateVector& U) { return U * s; }
