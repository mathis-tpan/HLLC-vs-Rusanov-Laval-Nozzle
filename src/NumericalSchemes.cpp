#include "NumericalSchemes.h"
#include <algorithm>    // std::max
#include <cmath>        // std::abs
#include <array>        // std::array

namespace NumericalSchemes {

    StateVector  getRusanov(Cellule Gauche, Cellule Droite, double gamma) {
        StateVector F_interface ;

        //Récupération des vitesses du son (gauche et droite) "Ci"
        double  c_g = Gauche.getSoundSpeed(gamma);
        double  c_d = Droite.getSoundSpeed(gamma);

        //Récupération des vitesses  (gauche et droite) "Ui"
        double u_g = Gauche.getFluidVelocity();
        double u_d = Droite.getFluidVelocity();

        // Récupération des flux des cellules (gauche et droite) "Fi"
        StateVector F_g = Gauche.getFlux(gamma);
        StateVector F_d = Droite.getFlux(gamma);

        // Calcul de Smax
        double S_max = std::max(std::abs(u_g) + c_g , std::abs(u_d) + c_d);


        // Calcul du flux à l'interface
        F_interface.rho = 0.5 * (F_d.rho + F_g.rho) - 0.5 * S_max * (Droite.U.rho - Gauche.U.rho);
        F_interface.rho_u = 0.5 * (F_d.rho_u + F_g.rho_u) - 0.5 * S_max * (Droite.U.rho_u - Gauche.U.rho_u);
        F_interface.rho_e = 0.5 * (F_d.rho_e + F_g.rho_e) - 0.5 * S_max * (Droite.U.rho_e - Gauche.U.rho_e);

        return F_interface;
    }

    StateVector getHLLC(Cellule Gauche, Cellule Droite, double gamma) {

        // Primitives gauche / droite
        double c_g = Gauche.getSoundSpeed(gamma), c_d = Droite.getSoundSpeed(gamma);
        double u_g = Gauche.getFluidVelocity(),   u_d = Droite.getFluidVelocity();
        double p_g = Gauche.getPressure(gamma),   p_d = Droite.getPressure(gamma);

        // États conservatifs et flux physiques : ce sont déjà des StateVector, zéro conversion
        StateVector U_g = Gauche.U,              U_d = Droite.U;
        StateVector F_g = Gauche.getFlux(gamma), F_d = Droite.getFlux(gamma);

        // Vitesses d'ondes
        double S_g = std::min(u_g - c_g, u_d - c_d);
        double S_d = std::max(u_g + c_g, u_d + c_d);
        double S_c = (p_d - p_g + U_g.rho*u_g*(S_g - u_g) - U_d.rho*u_d*(S_d - u_d))
                   / (U_g.rho*(S_g - u_g) - U_d.rho*(S_d - u_d));

        // État étoile U*_k puis saut de Rankine-Hugoniot, en UNE ligne lisible
        auto getHLLCFlux = [&](const StateVector& U_k, const StateVector& F_k,
                               double u_k, double p_k, double S_k)
        {
            double coeff = U_k.rho * (S_k - u_k) / (S_k - S_c);

            StateVector U_HLLC = {
                coeff,
                coeff * S_c,
                coeff * ( U_k.rho_e / U_k.rho
                        + (S_c - u_k) * (S_c + p_k / (U_k.rho * (S_k - u_k))) )
            };

            return F_k + S_k * (U_HLLC - U_k);   // F*_k = F_k + S_k (U*_k - U_k)
        };


        StateVector F_interface;
        if      (S_g >= 0.0) F_interface = F_g;
        else if (S_c >= 0.0) F_interface = getHLLCFlux(U_g, F_g, u_g, p_g, S_g);
        else if (S_d >= 0.0) F_interface = getHLLCFlux(U_d, F_d, u_d, p_d, S_d);
        else                 F_interface = F_d;
        return F_interface;
    }
}
