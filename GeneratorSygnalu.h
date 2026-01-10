#pragma once
#include <cmath>
#include <cstdint>
using namespace std;
class GeneratorSygnalu
{
public:
    enum Tryb { SINUS, PROSTOKAT };//określenie typu pracy generatora

private:
    Tryb tryb;

    double A;        // amplituda
    double S;        // skŁadowa staŁa
    double p;        // wypeŁnienie
    double TRZ;      // okres rzeczywisty [s]
    int TT;          // taktowanie symulacji [ms]
    int T;           // ile próbek na okres sygbału

public:
    GeneratorSygnalu()
        : tryb(SINUS), A(1.0), S(0.0), p(0.5), TRZ(1.0), TT(100), T(10)
    {}//Konstruktor

    void ustawTryb(Tryb t) {//seter ustawienie trybu pracy generatora
         tryb = t;
    }
    void ustawA(double a) {//setter ustawienie Amplitudy
        A = a;
        przeliczT();
    }
    void ustawS(double s) {//seter ustawienie składniowej stałej
        S = s;
    }
    void ustawP(double pp) {//seter ustawienie wypełnienia
        p = pp;
    }
    void ustawTRZ(double trz) {//seter ustawienie okresu
        TRZ = trz;
        przeliczT();
    }
    void ustawTT(int tt) {//seter ustawienie taktowania
        TT = tt;
        przeliczT();
    }

    int getT() const {//geter okres w próbkach
        return T;
    }
    Tryb getTryb() const {//geter tryb sygnału
        return tryb;
    }

    // OBOWI¥ZKOWE przeliczenie okresu dyskretnego
    void przeliczT()
    {
        if (TT <= 0) {//zapobieganie dzieleniu przez 0
            TT = 1;
        }
        double Td = (1000.0 * TRZ) / TT;//zamiana z s na ms
        T = static_cast<int>(round(Td));//zaokrąglanie liczby do całkowitej
        if (T < 1) T = 1;//ustawienie okresu minimalnego
    }

    double generuj(int i) const
    {
        if (T <= 0){//zapobieganie błędnemu T
            return S;
        }
        int k = i % T;//przekształcenie numeru próbki na pozycję wewnątrz okresu sygnału

        switch (tryb)
        {
        case SINUS:{
            return A * sin(k * (2.0 * M_PI / T)) + S;
        }
        case PROSTOKAT:{
            if(k<p*T){
                return A + S;//stan wysoki
            }else{
                return S;//stan niski
            }
        }
        default:{//wypadek nieznanego trybu
            return S;
            }
        }
    }
};
