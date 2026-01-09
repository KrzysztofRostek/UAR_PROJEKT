#pragma once

#include "GeneratorSygnalu.h"
#include "RegulatorPID.h"
#include "ModelARX.h"
#include "ProstyUAR.h"

class SymulatorUAR
{
private:
    GeneratorSygnalu generator;
    RegulatorPID pid;
    ModelARX arx;
    ProstyUAR uar;

    int k; // numer kroku

public:
    SymulatorUAR(const GeneratorSygnalu& gen,
                 const RegulatorPID& pid_,
                 const ModelARX& arx_)
        : generator(gen),
          pid(pid_),
          arx(arx_),
          uar(arx, pid),
          k(0)
    {}

    void reset()
    {
        k = 0;
        uar.reset();
    }

    void krok(double& w, double& e, double& u, double& y)
    {
        w = generator.generuj(k);
        y = uar.symuluj(w);
        e = w - y;

        // u jest pośrednie – możemy je odtworzyć logicznie
        u = pid.symuluj(e); // (opcjonalne – jeśli GUI chce u)

        k++;
    }

    // ===== SETTERY dla GUI =====
    void setKp(double kp) { pid.setKp(kp); }
    void setTi(double ti) { pid.setStalaCalk(ti); }
    void setTd(double td) { pid.setTd(td); }

    void setGeneratorTryb(GeneratorSygnalu::Tryb t) { generator.ustawTryb(t); }
    void setGeneratorA(double a) { generator.ustawA(a); }
    void setGeneratorTRZ(double trz) { generator.ustawTRZ(trz); }
};
