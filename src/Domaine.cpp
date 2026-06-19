#include "Domaine.h"
#include <iostream>
#include <fstream>
#include <limits>    // std::numeric_limits
#include <algorithm> // std::min
#include <cmath>     // std::abs

double Domaine::getDeltat(const Cellule &cell, double gamma, double CFL) const {
    // initialisation de c (la vietsse du son)
    const double c = cell.getSoundSpeed(gamma) ;

    // récupération de la vitesse
    const double u = cell.getFluidVelocity() ;

    // Calcul du delta t
    double delta_t = CFL * cell.dx / (std::abs(u) + c); ;

    return delta_t ;
}

double Domaine::getGlobalDeltat(std::vector<Cellule> tube, double CFL, double gamma) const {
    double dt_min = 1000000.0;
    int N = tube.size();

    for (int i = 0; i < N; i++) {
        double dt_local = getDeltat(tube[i], gamma, CFL);
        if (dt_local < dt_min) {
            dt_min = dt_local;
        }
    }
    return dt_min;
}

void Domaine::exporterInstantane(int numero_photo) {
    // Chemin complet vers results/ depuis cmake-build-debug/ en créant un fichier nommé resultat_0.csv, resultat_1.csv, etc.
    std::string chemin = "../results/resultat_"
                         + std::to_string(numero_photo)
                         + ".csv";

    std::ofstream fichier(chemin);

    // En-tête du fichier CSV
    fichier << "x,rho,vitesse,pression\n";

    int N = tuyere.size() - 2;
    for (int i = 1; i <= N; i++) {
        double u = tuyere[i].getFluidVelocity();
        // Pression = (rho_E - 0.5 * rho * u^2) * (gamma - 1)
        double p = tuyere[i].getPressure(gamma);

        // On écrit les colonnes séparées par des virgules
        fichier << tuyere[i].x << ","
                << tuyere[i].U.rho << ","
                << u << ","
                << p << "\n";
    }
    fichier.close();
}

void Domaine::appliquerConditionsLimites(double r) {
    int N = tuyere.size() - 2;
    tuyere[0] = tuyere[1];
    tuyere[N+1] = tuyere[N];

    // Entrée subsonique : conditions imposées

    double p_inlet = 101325;  // condition initiale p0 qui remplie la tuyere
    double T_inlet = 300;     // condition initiale T0 qui remplie la tuyere
    double u_inlet = 34.7;    // condition initiale u0 qui remplie la tuyere (initialisation à Mach = 0.1 pour eviter probleme d'ecoulement nul au départ)
    double rho_inlet = p_inlet/(r*T_inlet);
    tuyere[0].U.rho = rho_inlet;
    tuyere[0].U.rho_u = rho_inlet * u_inlet;
    tuyere[0].U.rho_e = p_inlet / (gamma - 1) + 0.5 * rho_inlet * u_inlet * u_inlet;
}

Domaine::Domaine(int L, int N, double gamma, double r, double t_max) {
    this -> gamma = gamma ;
    this -> t_max = t_max ;


    tuyere.resize(N+2); // prendre en compte les cellules fantomes pour le flux
    double dx = static_cast<double>(L)/N;
    for (int i =1; i <= N; i++) {
        tuyere[i].dx = dx;
        tuyere[i].x = (i -0.5)*dx;
        tuyere[i].A = 1 + 2.2*(tuyere[i].x -1.5)*(tuyere[i].x -1.5) ;
        tuyere[i].dAdx = 4.4*(tuyere[i].x -1.5) ;

        double p = 101325;  // condition initiale p0 qui remplie la tuyere
        double T = 300;     // condition initiale T0 qui remplie la tuyere
        double u = 34.7;    // condition initiale u0 qui remplie la tuyere (initialisation à Mach = 0.1 pour eviter probleme d'ecoulement nul au départ)
        tuyere[i].U.rho = p / (r*T) ;
        tuyere[i].U.rho_u = tuyere[i].U.rho * u ;
        tuyere[i].U.rho_e = p/ (gamma -1 ) ;
    }
    appliquerConditionsLimites(r);
}

void Domaine::run(double r) {
    int N = tuyere.size() - 2; // Le nombre de vraies cellules
    std::vector<StateVector> flux_interfaces(N + 1);
    double t = 0.0;

    // Pour le fichier
    int compteur_iterations = 0;
    int compteur_photos = 0;

    // Photo initiale à t = 0
    exporterInstantane(compteur_photos);
    compteur_photos++;

    while (t < t_max) {
        appliquerConditionsLimites(r); // creer les cellules fantomes de valeurs egale à aux autres pour le flux

        double dt = getGlobalDeltat(tuyere, 0.5,gamma);

        // Étape 1 : La photographie des flux aux interfaces
        for (int i = 0; i <= N; i++) {
            flux_interfaces[i] = NumericalSchemes::getHLLC(tuyere[i], tuyere[i+1], gamma);
        }

        for (int i = 1 ; i <= N; i++) {
            StateVector S = tuyere[i].getSource(gamma);
            tuyere[i].U.rho   = tuyere[i].U.rho - (dt/tuyere[i].dx)*(1/tuyere[i].A)*(flux_interfaces[i].rho * tuyere[i].A - flux_interfaces[i-1].rho * tuyere[i-1].A) +dt*S.rho ;
            tuyere[i].U.rho_u = tuyere[i].U.rho_u - (dt/tuyere[i].dx)*(1/tuyere[i].A)*(flux_interfaces[i].rho_u * tuyere[i].A - flux_interfaces[i-1].rho_u * tuyere[i-1].A) +dt*S.rho_u ;
            tuyere[i].U.rho_e = tuyere[i].U.rho_e - (dt/tuyere[i].dx)*(1/tuyere[i].A)*(flux_interfaces[i].rho_e * tuyere[i].A - flux_interfaces[i-1].rho_e * tuyere[i-1].A) +dt*S.rho_e ;
        }

        // Le temps avance !
        t = t + dt;

        compteur_iterations++;

        // Toutes les 15 itérations, on prend une photo !
        if (compteur_iterations % 15 == 0) {
            exporterInstantane(compteur_photos);
            compteur_photos++;
        }
    }
    std::cout << "Simulation terminee avec succes a t = " << t << " s !" << std::endl;

    // Une dernière photo à la fin
    exporterInstantane(compteur_photos);
    std::cout << "Simulation terminee ! " << compteur_photos << " fichiers CSV generes !" << std::endl;

}

