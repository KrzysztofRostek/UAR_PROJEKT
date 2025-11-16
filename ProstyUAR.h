#pragma once
// ProstyUAR.h  -- podmieniæ wszystkie stare pliki ProstyUAR

#include "ModelARX.h"
#include "RegulatorPID.h"

class ProstyUAR
{
private:
    ModelARX& obiekt;
    RegulatorPID* pid;
    double y_prev;   // pomiar opóŸniony o 1 próbkê

public:
    // Konstruktor - przyjmujemy referencje do istniej¹cych obiektów
    ProstyUAR(ModelARX& ob, RegulatorPID& r)
        : obiekt(ob), pid(&r), y_prev(0.0)
    {}

    // Reset: resetujemy pamiêæ regulatora i obiektu
    void reset()
    {
        y_prev = 0.0;
        if (pid) pid->reset();
        obiekt.reset();
    }

    // ----------------------------
    // Signtura zgodna z Twoim main:
    //   uar.krok(w,e,u,y);
    // gdzie w,e,u,y s¹ zmiennymi w main (przez referencjê)
    // ----------------------------
    void krok(double& w, double& e, double& u, double& y)
    {
        // y_prev to pomiar opóŸniony (y(k-1))
        double y_hat = y_prev;
        e = w - y_hat;

        // wylicz sygna³ steruj¹cy przez regulator
        if (pid)
            u = pid->symuluj(e);
        else
            u = 0.0;

        // podaj sterowanie do modelu
        y = obiekt.symuluj(u);

        // zapisz aktualne y jako poprzednie dla nastêpnego kroku
        y_prev = y;
    }

    // ----------------------------
    // Alternatywna metoda u¿ywana w niektórych testach:
    //   y = uar.symuluj(w);
    // ----------------------------
    double symuluj(double w)
    {
        double e, u, y;
        // wywo³ujemy krok z lokalnymi zmiennymi
        krok(w, e, u, y);
        return y;
    }
};
