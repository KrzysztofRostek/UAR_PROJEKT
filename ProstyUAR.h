#pragma once

#include "GeneratorSygnalu.h"
#include "ModelARX.h"
#include "RegulatorPID.h"

class ProstyUAR
{
private:
    ModelARX &obiekt;  //referemcka dp ARX
    RegulatorPID &pid; //referencja do PID

    double wartwyjsc_poprzedni; //poprzednie wyjście y(k-1)
    double k_generator;

public:
    ProstyUAR(ModelARX &arx, RegulatorPID &pid)
        : obiekt(arx), pid(pid), wartwyjsc_poprzedni(0.0), k_generator(0)
    {}

    void reset()
    {
        wartwyjsc_poprzedni = 0.0;
        k_generator = 0;   // reset licznika generatora
        pid.reset();
        obiekt.reset();
    }

    void krok(double &generatorVal, double &uchyb, double &PID, double &WartWyjsc, GeneratorSygnalu &gen, int)
    {
        // Używamy własnego licznika generatora
        generatorVal = gen.generuj(k_generator);
        k_generator++;

        uchyb = generatorVal - wartwyjsc_poprzedni;
        PID = pid.symuluj(uchyb);
        WartWyjsc = obiekt.symuluj(PID);
        wartwyjsc_poprzedni = WartWyjsc;
    }
};
