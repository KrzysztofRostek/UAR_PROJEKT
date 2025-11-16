#pragma once
#include "GeneratorSygnalu.h"
#include "RegulatorPID.h"
#include "ModelARX.h"

class Symulator
{
private:
    GeneratorSygnalu gen;
    RegulatorPID pid;
    ModelARX arx;

    double y_prev;     // ŷ(i) = y(i-1)
    int i_krok;

public:
    Symulator(const GeneratorSygnalu& g,
              const RegulatorPID& r,
              const ModelARX& m)
        : gen(g), pid(r), arx(m), y_prev(0.0), i_krok(0)
    {}

    void reset()
    {
        pid.reset();
        arx.reset();
        y_prev = 0.0;
        i_krok = 0;
    }

    // wykonuje jeden krok symulacji
    void krok(double& w_out, double& e_out,
              double& u_out, double& y_out)
    {
        // 1. wartość zadana
        double w = gen.generuj(i_krok);

        // 2. uchyb
        double e = w - y_prev;

        // 3. regulator PID
        double u = pid.symuluj(e);

        // 4. obiekt ARX
        double y = arx.symuluj(u);

        // 5. zapisz aktualną wartość jako poprzednią
        y_prev = y;
        i_krok++;

        // przekazujemy wyniki na zewnątrz
        w_out = w;
        e_out = e;
        u_out = u;
        y_out = y;
    }
};
