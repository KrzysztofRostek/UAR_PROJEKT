#include "mainwindow.h"
#include "Symulator.h"
#include "arxwindow.h"
#include "qvalueaxis.h"
#include "ui_mainwindow.h"
#include <QFileDialog>
#include <QMessageBox>
#include <QtCharts/QChart>
#include <QtCharts/QChartView>
#include <QtCharts/QLineSeries>

#include <QDebug>


MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent)
    , ui(new Ui::MainWindow)
    , symulator(GeneratorSygnalu(), RegulatorPID(), ModelARX({0}, {0}))
    , aktualnyWektorA({0})
    , aktualnyWektorB({0})
    , aktualneOpoznienie(1)
    , aktualnySzum(0.0)
{
    ui->setupUi(this);
    this->showMaximized();

    //----Serie----//
    seriaP = new QLineSeries();
    seriaI = new QLineSeries();
    seriaD = new QLineSeries();
    seriaUchyb = new QLineSeries();
    seriaRegulator = new QLineSeries();
    seriaZad = new QLineSeries();
    seriaRegulowana = new QLineSeries();

    //----Tytuły----//
    seriaZad->setName("Wartość zadana");
    seriaRegulowana->setName("Wartość regulowana");
    seriaP->setName("P");
    seriaI->setName("I");
    seriaD->setName("D");
    seriaUchyb->setName("Uchyb");
    seriaRegulator->setName("Sterowanie u");

    //----Połączenie z symulatorem----//
    connect(&symulator, &SymulatorUAR::krokWykonany, this, &MainWindow::onKrokWykonany);

    //----Konfiguracja Wykresu Głównego (MAIN)----//
    QChart *Mainchart = new QChart();
    QChartView *MainchartView = new QChartView(Mainchart);
    Mainchart->setAnimationOptions(QChart::NoAnimation);
    MainchartView->setMinimumSize(600, 400);
    Mainchart->addSeries(seriaZad);
    Mainchart->addSeries(seriaRegulowana);
    Mainchart->setTitle("Zadana i Regulowana");

    mainX = new QValueAxis();
    mainY = new QValueAxis();
    mainX->setTitleText("Czas [s]");
    mainY->setTitleText("Wartość");
    mainX->setTickCount(11);
    Mainchart->addAxis(mainX, Qt::AlignBottom);
    Mainchart->addAxis(mainY, Qt::AlignLeft);
    seriaZad->attachAxis(mainX); seriaZad->attachAxis(mainY);
    seriaRegulowana->attachAxis(mainX); seriaRegulowana->attachAxis(mainY);

    MainchartView->setRenderHint(QPainter::Antialiasing);
    ui->horizontalLayout_5->addWidget(MainchartView);

    //----Konfiguracja Wykresu PID----//
    QChart *PIDchart = new QChart();
    QChartView *PIDchartView = new QChartView(PIDchart);
    PIDchart->setAnimationOptions(QChart::NoAnimation);
    PIDchart->addSeries(seriaP);
    PIDchart->addSeries(seriaI);
    PIDchart->addSeries(seriaD);
    PIDchart->setTitle("Składowe sterowania PID");

    pidX = new QValueAxis();
    pidY = new QValueAxis();
    pidX->setTitleText("Czas [s]");
    pidY->setTitleText("Wartość");
    PIDchart->addAxis(pidX, Qt::AlignBottom);
    PIDchart->addAxis(pidY, Qt::AlignLeft);
    seriaP->attachAxis(pidX); seriaP->attachAxis(pidY);
    seriaI->attachAxis(pidX); seriaI->attachAxis(pidY);
    seriaD->attachAxis(pidX); seriaD->attachAxis(pidY);

    PIDchartView->setRenderHint(QPainter::Antialiasing);
    ui->horizontalLayout_4->addWidget(PIDchartView, 1);

    //----Konfiguracja Wykresu Uchybu----//
    QChart *Uchybchart = new QChart();
    QChartView *UchybchartView = new QChartView(Uchybchart);
    Uchybchart->setAnimationOptions(QChart::NoAnimation);
    Uchybchart->addSeries(seriaUchyb);
    Uchybchart->setTitle("Uchyb");

    uchybX = new QValueAxis();
    uchybY = new QValueAxis();
    uchybX->setTitleText("Czas [s]");
    uchybY->setTitleText("Wartość");
    Uchybchart->addAxis(uchybX, Qt::AlignBottom);
    Uchybchart->addAxis(uchybY, Qt::AlignLeft);
    seriaUchyb->attachAxis(uchybX); seriaUchyb->attachAxis(uchybY);

    UchybchartView->setRenderHint(QPainter::Antialiasing);
    ui->horizontalLayout_4->addWidget(UchybchartView, 1);

    //----Konfiguracja Wykresu Regulatora----//
    QChart *Regulatorchart = new QChart();
    QChartView *RegulatorchartView = new QChartView(Regulatorchart);
    Regulatorchart->setAnimationOptions(QChart::NoAnimation);
    Regulatorchart->addSeries(seriaRegulator);
    Regulatorchart->setTitle("Regulator");

    regX = new QValueAxis();
    regY = new QValueAxis();
    regX->setTitleText("Czas [s]");
    regY->setTitleText("Wartość");
    Regulatorchart->addAxis(regX, Qt::AlignBottom);
    Regulatorchart->addAxis(regY, Qt::AlignLeft);
    seriaRegulator->attachAxis(regX); seriaRegulator->attachAxis(regY);

    RegulatorchartView->setRenderHint(QPainter::Antialiasing);
    ui->horizontalLayout_4->addWidget(RegulatorchartView, 1);
    RegulatorchartView->setMinimumSize(0, 300);
}


void MainWindow::onKrokWykonany(double w, double y, double e, double u, int k, double P, double I, double D)
{
    // 1. Parametry okna i czasu
    // Pobieramy wartosc z GUI (5-50s)
    double oknoSekundy = ui->spinBoxOknoczasowe->value();
    double interwalMs = symulator.getInterwalMs();
    double t = k * interwalMs / 1000.0;

    // Obliczamy ile probek powinno byc na ekranie (np. dla 10s i 200ms n=50)
    int n_docelowe = static_cast<int>((oknoSekundy * 1000.0) / interwalMs);

    // 2. Dodanie nowych danych
    seriaZad->append(t, w);
    seriaRegulowana->append(t, y);
    seriaP->append(t, P);
    seriaI->append(t, I);
    seriaD->append(t, D);
    seriaUchyb->append(t, e);
    seriaRegulator->append(t, u);

    // 3. Logika zapominania probek
    int n_aktualne = seriaZad->count();

    if (n_aktualne > n_docelowe)
    {
        // Jesli okno zmalalo, usuwamy 2 stare probki na 1 nowa (rozszerzanie wykresu)
        int doUsuniecia = (n_aktualne - n_docelowe >= 2) ? 2 : 1;

        for (int i = 0; i < doUsuniecia; ++i)
        {
            seriaZad->remove(0);
            seriaRegulowana->remove(0);
            seriaP->remove(0);
            seriaI->remove(0);
            seriaD->remove(0);
            seriaUchyb->remove(0);
            seriaRegulator->remove(0);
        }
    }
    else if (n_aktualne == n_docelowe)
    {
        // Standardowe przesuwanie: usuwamy 1 najstarsza probke
        seriaZad->remove(0);
        seriaRegulowana->remove(0);
        seriaP->remove(0);
        seriaI->remove(0);
        seriaD->remove(0);
        seriaUchyb->remove(0);
        seriaRegulator->remove(0);
    }
    // Jesli n_aktualne < n_docelowe, nic nie usuwamy (sciskanie wykresu)

    // 4. Aktualizacja osi X
    // Zakres od czasu pierwszej dostepnej probki do aktualnego t
    double startTime = seriaZad->at(0).x();
    mainX->setRange(startTime, t);
    pidX->setRange(startTime, t);
    uchybX->setRange(startTime, t);
    regX->setRange(startTime, t);

    // 5. Skalowanie pionowe Y
    dopasujSkalePionowa(mainY, seriaZad, seriaRegulowana);
    dopasujSkalePionowa(pidY, seriaP, seriaI, seriaD);
    dopasujSkalePionowa(uchybY, seriaUchyb);
    dopasujSkalePionowa(regY, seriaRegulator);
}

void MainWindow::dopasujSkalePionowa(QValueAxis *osY, QLineSeries *pierwszaSeria, QLineSeries *drugaSeria, QLineSeries *trzeciaSeria)
{
    double minWartosc = 999999.0;
    double maxWartosc = -999999.0;
    bool czyZnalezionoJakikolwiekPunkt = false;

    // Pobieramy aktualny czas rozpoczecia wykresu z osi X
    double czasStartu = mainX->min();

    // Tworzymy liste pomocnicza, zeby latwiej przejsc przez wszystkie serie w petli
    QList<QLineSeries*> listaSerii;

    if (pierwszaSeria != nullptr)
    {
        listaSerii.append(pierwszaSeria);
    }

    if (drugaSeria != nullptr)
    {
        listaSerii.append(drugaSeria);
    }

    if (trzeciaSeria != nullptr)
    {
        listaSerii.append(trzeciaSeria);
    }

    // Przeszukujemy kazda serie, aby znalezc wartosci min i max widoczne na ekranie
    for (QLineSeries *seria : listaSerii)
    {
        QList<QPointF> punkty = seria->points();

        // Idziemy od konca listy (najnowsze punkty)
        for (int i = punkty.size() - 1; i >= 0; i--)
        {
            // Sprawdzamy, czy punkt miesci sie w widocznym oknie czasowym
            if (punkty[i].x() < czasStartu)
            {
                // Jesli punkt jest starszy niz lewa krawedz wykresu,
                // przerywamy szukanie w tej konkretnej serii.
                break;
            }

            czyZnalezionoJakikolwiekPunkt = true;

            // Sprawdzamy czy wartosc Y punktu jest mniejsza niz dotychczasowe minimum
            if (punkty[i].y() < minWartosc)
            {
                minWartosc = punkty[i].y();
            }

            // Sprawdzamy czy wartosc Y punktu jest wieksza niz dotychczasowe maksimum
            if (punkty[i].y() > maxWartosc)
            {
                maxWartosc = punkty[i].y();
            }
        }
    }

    // Jesli w oknie czasowym znalezlismy jakiekolwiek punkty, ustawiamy skale osi Y
    if (czyZnalezionoJakikolwiekPunkt == true)
    {
        double rozpietosc = maxWartosc - minWartosc;

        // Sprawdzamy, czy sygnal nie jest plaska linia
        if (rozpietosc < 0.001)
        {
            // Dla plaskiego sygnalu dajemy staly margines +/- 1.0
            osY->setRange(minWartosc - 1.0, maxWartosc + 1.0);
        }
        else
        {
            // Obliczamy margines (10% wysokosci wykresu), zeby linia nie dotykala brzegow
            double margines = rozpietosc * 0.1;
            osY->setRange(minWartosc - margines, maxWartosc + margines);
        }
    }
}
void MainWindow::wyczyscWykresy()
{
    // MAIN
    seriaZad->clear();
    seriaRegulowana->clear();

    // PID
    seriaP->clear();
    seriaI->clear();
    seriaD->clear();

    // UCHYB
    seriaUchyb->clear();

    // REGULATOR
    seriaRegulator->clear();

    // Reset osi
    mainX->setRange(0, 10);
    mainY->setRange(-1, 1);

    pidX->setRange(0, 10);
    pidY->setRange(-1, 1);

    uchybX->setRange(0, 10);
    uchybY->setRange(-1, 1);

    regX->setRange(0, 10);
    regY->setRange(-1, 1);
}


MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::on_Sin_Button_clicked()
{
    symulator.setGeneratorTryb(GeneratorSygnalu::SINUS);
}

void MainWindow::on_Square_Button_clicked()
{
    symulator.setGeneratorTryb(GeneratorSygnalu::PROSTOKAT);
}

void MainWindow::on_spinBOX_WzmocK_editingFinished()
{
    symulator.setPID_Kp(ui->spinBOX_WzmocK->value());
    ui->spinBOX_WzmocK->setMinimum(0.0);
    ui->spinBOX_WzmocK->setMaximum(1000.0);
    ui->spinBOX_WzmocK->setSingleStep(0.1);
    ui->spinBOX_WzmocK->setDecimals(5);
}

void MainWindow::on_spinBOX_Amplituda_editingFinished()
{
    symulator.setGeneratorA(ui->spinBOX_Amplituda->value());
    ui->spinBOX_Amplituda->setMinimum(0.0);
    ui->spinBOX_Amplituda->setMaximum(1000.0);
    ui->spinBOX_Amplituda->setSingleStep(0.1);
    ui->spinBOX_Amplituda->setDecimals(5);
}

void MainWindow::on_spinBOX_Czstotliwosc_editingFinished()
{
    symulator.setGeneratorTRZ(ui->spinBOX_Czstotliwosc->value());
    ui->spinBOX_Czstotliwosc->setMinimum(0.01);
    ui->spinBOX_Czstotliwosc->setMaximum(1000.0);
    ui->spinBOX_Czstotliwosc->setSingleStep(0.01);
    ui->spinBOX_Czstotliwosc->setDecimals(5);
}
void MainWindow::on_spinBox_Wypelnienie_editingFinished()
{
    symulator.setGeneratorP(ui->spinBox_Wypelnienie->value());
}


void MainWindow::on_SpinBox_Stala_editingFinished()
{
    symulator.setGeneratorS(ui->SpinBox_Stala->value());
}
void MainWindow::on_spinBOX_Td_editingFinished()
{
    symulator.setPID_Td(ui->spinBOX_Td->value());
    ui->spinBOX_Td->setMinimum(0.0);
    ui->spinBOX_Td->setMaximum(1000.0);
    ui->spinBOX_Td->setSingleStep(0.1);
    ui->spinBOX_Td->setDecimals(5);
}

void MainWindow::on_spinBOX_Ti_editingFinished()
{
    symulator.setPID_Ti(ui->spinBOX_Ti->value());
    ui->spinBOX_Ti->setMinimum(0.0);
    ui->spinBOX_Ti->setMaximum(1000.0);
    ui->spinBOX_Ti->setSingleStep(0.1);
    ui->spinBOX_Ti->setDecimals(5);
}

void MainWindow::on_spinBOX_Interwal_editingFinished()
{
    symulator.setGeneratorTT(ui->spinBOX_Interwal->value());
    ui->spinBOX_Interwal->setMinimum(10);
    ui->spinBOX_Interwal->setMaximum(1000);
    ui->spinBOX_Interwal->setSingleStep(1);
    ui->spinBOX_Interwal->setDecimals(5);


}

void MainWindow::on_spinBoxOknoczasowe_editingFinished()
{
    // Odswiezenie osi X
    if (seriaZad->count() > 0)
    {
        double okno = ui->spinBoxOknoczasowe->value();
        double t_koncowe = seriaZad->points().last().x();
        double t_startowe = t_koncowe - okno;

        if (t_startowe < 0) t_startowe = 0;

        mainX->setRange(t_startowe, t_koncowe);
        pidX->setRange(t_startowe, t_koncowe);
        uchybX->setRange(t_startowe, t_koncowe);
        regX->setRange(t_startowe, t_koncowe);
    }
}


void MainWindow::on_radio_przed_toggled(bool checked)
{
    if (checked)
        symulator.setPID_TypCalki(RegulatorPID::Zew);
}

void MainWindow::on_radio_pod_toggled(bool checked)
{
    if (checked)
        symulator.setPID_TypCalki(RegulatorPID::Wew);
}

void MainWindow::on_Reset_d_clicked()
{
    symulator.setPID_Td(0);
}

void MainWindow::on_Reset_i_clicked()
{
    symulator.setPID_Ti(0);
}

void MainWindow::on_START_Button_clicked()
{
    symulator.start();
}

void MainWindow::on_STOP_Bttun_clicked()
{
    symulator.stop();
}

void MainWindow::on_RESET_Button_clicked()
{
    symulator.reset();

    symulator.setPID_Kp(0);
    ui->spinBOX_WzmocK->setValue(0);
    symulator.setPID_Ti(0);
    ui->spinBOX_Ti->setValue(0);
    symulator.setPID_Td(0);
    ui->spinBOX_Td->setValue(0);
    symulator.setGeneratorA(0);
    ui->spinBOX_Amplituda->setValue(0);
    symulator.setGeneratorTRZ(0);
    ui->spinBOX_Czstotliwosc->setValue(0);
    symulator.setGeneratorTT(200);
    ui->spinBOX_Interwal->setValue(200);
    symulator.setGeneratorS(0);
    ui->SpinBox_Stala->setValue(0);
    symulator.setGeneratorP(0);
    ui->spinBox_Wypelnienie->setValue(0);


    symulator.setPID_TypCalki(RegulatorPID::Zew);
    ui->radio_przed->setChecked(true);
    ui->radio_pod->setChecked(false);
    aktualnyWektorA = {0};
    aktualnyWektorB = {0};
    aktualneOpoznienie = 1;
    aktualnySzum = 0.0;


    symulator.setARX(aktualnyWektorA, aktualnyWektorB, aktualneOpoznienie, aktualnySzum);
    wyczyscWykresy();
}

void MainWindow::ustawARXDane(const std::vector<double> &a,
                              const std::vector<double> &b,
                              int opoznienie,
                              double szum,
                              double uMin, double uMax,
                              double yMin, double yMax,
                              bool aktywne)
{
    aktualnyWektorA = a;
    aktualnyWektorB = b;
    aktualneOpoznienie = opoznienie;
    aktualnySzum = szum;
    arx_uMin = uMin;
    arx_uMax = uMax;
    arx_yMin = yMin;
    arx_yMax = yMax;
    arx_ograniczenia = aktywne;

    symulator.setARX(a, b, opoznienie, szum);

    symulator.setARX_Umin(uMin);
    symulator.setARX_Umax(uMax);
    symulator.setARX_Ymin(yMin);
    symulator.setARX_Ymax(yMax);
    symulator.setARX_Ograniczenia(aktywne);
}

void MainWindow::on_Konf_ARX_Button_clicked()
{
    if (!arxwindow) {
        arxwindow = new ARXwindow(this);

        connect(arxwindow, &ARXwindow::zatwierdzonoARX, this, &MainWindow::ustawARXDane);
    }


    arxwindow->ustawDane(aktualnyWektorA, aktualnyWektorB, aktualneOpoznienie, aktualnySzum,
                         arx_uMin, arx_uMax, arx_yMin, arx_yMax, arx_ograniczenia);

    arxwindow->show();
    arxwindow->raise();
    arxwindow->activateWindow();
}

void MainWindow::on_Zapisz_Button_clicked()
{
    QString sciezka = QFileDialog::getSaveFileName(this, "Zapisz konfigurację",
                                                   "", "JSON (*.json)");
    if (sciezka.isEmpty()) return;
    bool sukces = menedzerKonfig.zapiszKonfiguracje(
        sciezka,
        aktualnyWektorA, aktualnyWektorB, aktualneOpoznienie, aktualnySzum,
        ui->spinBOX_WzmocK->value(), ui->spinBOX_Ti->value(),
        ui->spinBOX_Td->value(), ui->radio_przed->isChecked() ? 0 : 1,
        ui->Sin_Button->isChecked() ? 0 : 1,
        ui->spinBOX_Amplituda->value(),
        ui->spinBOX_Czstotliwosc->value(),
        ui->SpinBox_Stala->value(),
        ui->spinBox_Wypelnienie->value(),
        ui->spinBOX_Interwal->value()
        );

    if (sukces) {
        QMessageBox::information(this ,"sukces", "konfiguracja zapisana");
    } else {
        QMessageBox::warning(this, "blad", "nie udało sie zapisac");
    }
}

void MainWindow::on_Wczytaj_Button_clicked()
{
    QString sciezka = QFileDialog::getOpenFileName(this, "Wczytaj konfigurację",
                                                   "", "JSON (*.json)");
    if (sciezka.isEmpty()) return;

    std::vector<double> a, b;
    int opoznienie;
    double odchylenie, Kp, Ti, Td;
    int typCalki, trybGeneratora;
    double amplituda, StalaSkladniowa, Wypelnienie, czestotliwosc;
    int interwalMs;

    bool sukces = menedzerKonfig.wczytajKonfiguracje(
        sciezka,
        a, b, opoznienie, odchylenie,
        Kp, Ti, Td, typCalki,
        trybGeneratora, amplituda, czestotliwosc, StalaSkladniowa, Wypelnienie,
        interwalMs
        );

    if (sukces) {
        aktualnyWektorA = a;
        aktualnyWektorB = b;
        aktualneOpoznienie = opoznienie;
        aktualnySzum = odchylenie;

        symulator.setARX(a, b, opoznienie, odchylenie);

        ui->spinBOX_WzmocK->setValue(Kp);
        ui->spinBOX_Ti->setValue(Ti);
        ui->spinBOX_Td->setValue(Td);

        if (typCalki == 0) {
            ui->radio_przed->setChecked(true);
            ui->radio_pod->setChecked(false);
        } else {
            ui->radio_przed->setChecked(false);
            ui->radio_pod->setChecked(true);
        }

        if (trybGeneratora == 0) {
            ui->Sin_Button->setChecked(true);
            ui->Square_Button->setChecked(false);
        } else {
            ui->Sin_Button->setChecked(false);
            ui->Square_Button->setChecked(true);
        }

        ui->spinBOX_Amplituda->setValue(amplituda);
        ui->spinBOX_Czstotliwosc->setValue(czestotliwosc);
        ui->SpinBox_Stala->setValue(StalaSkladniowa);
        ui->spinBox_Wypelnienie->setValue(Wypelnienie);
        ui->spinBOX_Interwal->setValue(interwalMs);

        QMessageBox::information(this, "sukces", "konfiguracja wczytana");
    } else {
        QMessageBox::warning(this, "blad", "nie udało sie wczytac");
    }
}




