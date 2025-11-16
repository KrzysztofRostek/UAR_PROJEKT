#pragma once
#include <cmath>
#include <cstdint>

class GeneratorSygnalu
{
public:
    enum Tryb { SINUS, PROSTOKAT };

private:
    Tryb tryb;

    double A;        // amplituda
    double S;        // sk³adowa sta³a
    double p;        // wype³nienie (0..1)
    double TRZ;      // okres rzeczywisty [s]
    int TT;          // taktowanie symulacji [ms]

    int T;           // okres dyskretny [próbki]

public:
    GeneratorSygnalu()
        : tryb(SINUS), A(1.0), S(0.0), p(0.5), TRZ(1.0), TT(100), T(10)
    {}

    void ustawTryb(Tryb t) { tryb = t; }
    void ustawA(double a) { A = a; przeliczT(); }
    void ustawS(double s) { S = s; }
    void ustawP(double pp) { p = pp; }
    void ustawTRZ(double trz) { TRZ = trz; przeliczT(); }
    void ustawTT(int tt) { TT = tt; przeliczT(); }

    int getT() const { return T; }
    Tryb getTryb() const { return tryb; }

    // OBOWI¥ZKOWE przeliczenie okresu dyskretnego
    void przeliczT()
    {
        if (TT <= 0) TT = 1;
        double Td = (1000.0 * TRZ) / TT;
        T = static_cast<int>(std::round(Td));
        if (T < 1) T = 1;
    }

    // Zwraca w(i)
    double generuj(int i) const
    {
        if (T <= 0) return S;

        int k = i % T;

        switch (tryb)
        {
        case SINUS:
            return A * std::sin(k * (2.0 * M_PI / T)) + S;

        case PROSTOKAT:
            return (k < p * T) ? (A + S) : S;

        default:
            return S;
        }
    }
};
