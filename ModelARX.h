using namespace std;
#include <vector>
#include <cstdlib>
#include <ctime>
#include <cmath>

class ModelARX
{
private:
    double A[10];//zależność od poprzednich wartości wyjścia
    double B[10];//zależność od wcześniejszych wartości sterowania
    int IleA;//liczba wartości wyjścia
    int IleB;//liczba wartości sterujących
    int opoznienie;//jakie opóźnienie między sygnałem u a reakcją y
    double odchylenie;//amplituda możliwego szumu
    double HistWyjsc[10];//tablica histori wyjść y[k-1]
    double HistSter[20];//tablica histori watości sterujących u[k-1]

public:
    ModelARX(const vector<double>& a,//wektor współczynników A
             const vector<double>& b,//weksto współczynników B
             int opoznienie = 1,
             double odchylenie = 0.0)
        : odchylenie(odchylenie), opoznienie(opoznienie)//inicializacja pól
    {
        IleA = a.size();//zapisanie długości wektora a i b
        IleB = b.size();
        for (int i = 0; i < IleA; ++i){
            A[i] = a[i];//kopiowanie wspczynników wekoru a do tablicy A
        }
        for (int i = 0; i < IleB; ++i){
            B[i] = b[i];//kopiowanie wspczynników wekoru b do tablicy B
        }
        srand(time(0));
        reset();
    }
    void reset()//resetowanie histori wyjść i histori wartości sterujących
    {
        for (int i = 0; i < IleA; ++i) {
            HistWyjsc[i] = 0.0;}
        for (int i = 0; i < IleB+opoznienie+1; ++i) {
            HistSter[i] = 0.0;}
    }
    double losowySzum()//generowanie losowego szumu w zakresie [-odchylenie, odchylenie]
    {
        if (odchylenie == 0.0) return 0.0;
        double r = (double)rand() / RAND_MAX; // [0,1]
        r = 2.0 * r - 1.0; // [-1,1]
        return r * odchylenie;
    }
    double symuluj(double WartSter)//symulacja jednego przejścia
    {
        for (int i = IleB+opoznienie-1; i > 0; --i){//przesunięcie tablicy z historią wartości sterujących o jeden
            HistSter[i] = HistSter[i - 1];
            }
        HistSter[0] = WartSter;//zapisanie ostatniej wartości sterującej na początku tablicy

        double WartWyjsc = 0.0;

        for (int i = 0; i < IleB; ++i)
        {
            int OpozSter = i + opoznienie;//obliczenie którą komórkę histori wartości sterujących użyć do każdego współczynnika B
            if (OpozSter < 20)//zabezpieczenie przed wyjściem poza zakres
                WartWyjsc += B[i] * HistSter[OpozSter];//sumowanie wgodnie ze wzorem  (b0ui-k + b1ui-k-1 + … + bdBui-k-dB)
        }

        for (int i = 0; i < IleA; ++i){
            WartWyjsc -= A[i] * HistWyjsc[i];//odejmowanie od wartości wyjścia zgodznie ze wzorem sumy (a1yi-1 + a2yi-2 + … + adAyi-dA)
        }

        WartWyjsc += losowySzum(); // szum dodany do wyjścia składając wszystko otrzymujemy cały wzór yi = (b0ui-k + b1ui-k-1 + … + bdBui-k-dB) - (a1yi-1 + a2yi-2 + … + adAyi-dA) + zi,

        for (int i = IleA - 1; i > 0; --i){//przesunięcie tablicy z historią wartości wyjść o jeden
            HistWyjsc[i] = HistWyjsc[i - 1];
            }
        HistWyjsc[0] = WartWyjsc;//zapisanie ostatniej wartości wyjścia na początku tablicy

        return WartWyjsc;
    }
};

