#pragma once
#include <cmath>

class RegulatorPID
{
public:
    enum LiczCalk { PROSTOKATNY, TRAPEZOWY, Wew, Zew, ZERO };

private:
    double Kp;
    double Ti;
    double Td;
    double T;

    double e_prev;
    double sum_internal;
    double sum_external;

    LiczCalk typCalki;
    static constexpr double EPS = 1e-12;

public:
    RegulatorPID(double Kp_ = 1.0, double Ti_ = 0.0, double Td_ = 0.0, double T_ = 1.0)
        : Kp(Kp_), Ti(Ti_), Td(Td_), T(T_),
          e_prev(0.0), sum_internal(0.0), sum_external(0.0),
          typCalki(PROSTOKATNY)
    {}

    void reset()
    {
        e_prev = 0.0;
        sum_internal = 0.0;
        sum_external = 0.0;
    }

    void setKp(double kp) { Kp = kp; }
    void setTi(double ti) { Ti = ti; }
    void setTd(double td) { Td = td; }
    void setT(double t)  { T = t; }
    void setStalaCalk(double ti) { Ti = ti; }

    void setLiczCalk(LiczCalk metoda)
    {
        if (typCalki != metoda) {
            // Synchronizacja pamięci całki
            if (typCalki == Zew && metoda != Zew) {
                // Zew -> Wew: przenieś stan
                sum_internal = sum_external / (Ti > EPS ? Ti : 1.0);
            }
            else if (typCalki != Zew && metoda == Zew) {
                // Wew -> Zew: przenieś stan
                sum_external = sum_internal * (Ti > EPS ? Ti : 1.0);
            }
        }
        typCalki = metoda;
    }

    LiczCalk getLiczCalk() const { return typCalki; }

    double symuluj(double e)
    {
        // --- P ---
        double P = Kp * e;

        // --- I ---
        double I = 0.0;
        bool integral_enabled = (typCalki != ZERO) && (Ti > EPS);

        if (integral_enabled)
        {
            switch (typCalki)
            {
                case PROSTOKATNY:
                case Wew:
                    sum_internal += (T / Ti) * e;
                    I = sum_internal;
                    break;

                case TRAPEZOWY:
                    sum_internal += (T / (2.0 * Ti)) * (e + e_prev);
                    I = sum_internal;
                    break;

                case Zew:
                    sum_external += T * e;
                    I = sum_external / Ti;
                    break;

                default:
                    I = 0.0;
                    break;
            }
        }

        // --- D ---
        double D = 0.0;
        if (Td > EPS)
            D = Td * (e - e_prev) / T;

        double u = P + I + D;
        e_prev = e;

        return u;
    }
};
