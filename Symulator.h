#pragma once

#include "GeneratorSygnalu.h"
#include "RegulatorPID.h"
#include "ModelARX.h"
#include "ProstyUAR.h"

class SymulatorUAR
{
private:
    // WARSTWA DANYCH
    GeneratorSygnalu generator;
    RegulatorPID pid;
    ModelARX arx;
    ProstyUAR uar;

    // STAN SYMULACJI
    int k;
    double w, e, u, y;

    // ZEGAR
    bool symuluj;
    int interwalMs;

public:
    //KONSTRUKTOR
    SymulatorUAR(const GeneratorSygnalu& gen,
                 const RegulatorPID& pid,
                 const ModelARX& arx);

    //STEROWANIE SYMULACJĄ
    void start();
    void stop();
    void reset();
    void krokSymulacji();

    //USTAWIENIA GENERATORA
    void setGeneratorTryb(GeneratorSygnalu::Tryb t);
    void setGeneratorA(double a);
    void setGeneratorS(double s);
    void setGeneratorP(double p);
    void setGeneratorCzestotliwosc(double f);
    void setGeneratorTT(int tt);

    //PID
    void setPID_Kp(double kp);
    void setPID_Ti(double ti);
    void setPID_Td(double td);
    void setPID_T(double t);
    void setPID_TypCalki(RegulatorPID::LiczCalk typ);

    //ARX
    void setARX(const std::vector<double>& a,
                const std::vector<double>& b,
                int opoznienie,
                double szum);

    //GETTERY DLA GUI
    int getKrok() const;
    double getWartoscZadana() const;
    double getUchyb() const;
    double getSterowanie() const;
    double getWyjscie() const;

    int getInterwalMs() const;
    void setInterwalMs(int ms);
    bool czysymuluj() const;
};
