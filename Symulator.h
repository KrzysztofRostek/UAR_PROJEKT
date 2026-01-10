#pragma once

#include "GeneratorSygnalu.h"
#include "RegulatorPID.h"
#include "ModelARX.h"
#include "ProstyUAR.h"

class SymulatorUAR
{
private:
    // ===== WARSTWA DANYCH =====
    GeneratorSygnalu generator;
    RegulatorPID pid;
    ModelARX arx;
    ProstyUAR uar;

    // ===== STAN SYMULACJI =====
    int k;          // numer kroku
    double w;       // wartość zadana
    double e;       // uchyb
    double u;       // sterowanie (PID)
    double y;       // wyjście obiektu

public:
    // ===== KONSTRUKTOR =====
    SymulatorUAR(const GeneratorSygnalu& gen,
                 const RegulatorPID& pid_,
                 const ModelARX& arx_)
        : generator(gen),
          pid(pid_),
          arx(arx_),
          uar(arx, pid),
          k(0),
          w(0.0),
          e(0.0),
          u(0.0),
          y(0.0)
    {}

    // ===== RESET =====
    void reset()
    {
        k = 0;
        w = e = u = y = 0.0;
        uar.reset();
    }

    // Generator
    void setGeneratorTryb(GeneratorSygnalu::Tryb t) { generator.ustawTryb(t); }
    void setGeneratorA(double a) { generator.ustawA(a); }
    void setGeneratorS(double s) { generator.ustawS(s); }
    void setGeneratorP(double p) { generator.ustawP(p); }
    void setGeneratorTRZ(double trz) { generator.ustawTRZ(trz); }
    void setGeneratorTT(int tt) { generator.ustawTT(tt); }

    //  Regulator PID
    void setPID_Kp(double kp) { pid.setKp(kp); }
    void setPID_Ti(double ti) { pid.setStalaCalk(ti); }
    void setPID_Td(double td) { pid.setTd(td); }
    void setPID_T(double t) { pid.setT(t); }
    void setPID_TypCalki(RegulatorPID::LiczCalk typ) { pid.setLiczCalk(typ); }
    // ARX
    void setARX(const std::vector<double>& a,
            const std::vector<double>& b,
            int opoznienie,
            double szum)
{
    arx = ModelARX(a, b, opoznienie, szum);
    reset();
}

    void krokSymulacji()
    {
        uar.krok(w, e, u, y, generator, k);
        k++;
    }

    int getKrok() const { return k; }

    double getWartoscZadana() const { return w; }   // w(k)
    double getUchyb() const { return e; }           // e(k)
    double getSterowanie() const { return u; }      // u(k)
    double getWyjscie() const { return y; }         // y(k)
};
