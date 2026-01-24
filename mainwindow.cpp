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



    //----MAIN----//
    QChart *Mainchart = new QChart();
    QChartView *MainchartView = new QChartView(Mainchart);
    
    // Połączenie z symulatorem – dodawanie punktów w czasie rzeczywistym
    connect(&symulator, &SymulatorUAR::krokWykonany, this,
            [=](double w, double y, double e, double u, int k, double P, double I, double D)
            {
                double t = k * symulator.getInterwalMs() / 1000.0;  // przeliczamy czas na sekundy
                double windowTime = 10.0; // zakres osi X do wyświetlenia
                
                // --- Aktualizacja MAIN CHART ---
                seriaZad->append(t, w);
                seriaRegulowana->append(t, y);
                
                // Automatyczne ustawienie zakresu osi Y z lekkim marginesem
                double minY = w, maxY = w;
                for (const QPointF &p : seriaZad->points())        { minY = std::min(minY, p.y()); maxY = std::max(maxY, p.y()); }
                for (const QPointF &p : seriaRegulowana->points()) { minY = std::min(minY, p.y()); maxY = std::max(maxY, p.y()); }
                
                double margin = (maxY - minY) * 0.1; 
                if (margin < 1e-6) margin = 1.0;
                mainY->setRange(minY - margin, maxY + margin);
                mainX->setRange(t > windowTime ? t - windowTime : 0, t);
                
                // --- Aktualizacja PID ---
                seriaP->append(t, P);
                seriaI->append(t, I);
                seriaD->append(t, D);
                
                auto updateRange = [](QLineSeries* s, double &minVal, double &maxVal){
                    for (const QPointF &pt : s->points()) {
                        minVal = std::min(minVal, pt.y());
                        maxVal = std::max(maxVal, pt.y());
                    }
                };
                
                double minPID = P, maxPID = P;
                updateRange(seriaP, minPID, maxPID);
                updateRange(seriaI, minPID, maxPID);
                updateRange(seriaD, minPID, maxPID);
                
                double pidMargin = (maxPID - minPID) * 0.1; 
                if (pidMargin < 1e-6) pidMargin = 1.0;
                pidY->setRange(minPID - pidMargin, maxPID + pidMargin);
                pidX->setRange(t > windowTime ? t - windowTime : 0, t);
                
                // --- Aktualizacja uchybu ---
                seriaUchyb->append(t, e);
                double minE = e, maxE = e;
                updateRange(seriaUchyb, minE, maxE);
                double eMargin = (maxE - minE) * 0.1; if (eMargin < 1e-6) eMargin = 1.0;
                uchybY->setRange(minE - eMargin, maxE + eMargin);
                uchybX->setRange(t > windowTime ? t - windowTime : 0, t);
                
                // --- Aktualizacja sterowania ---
                seriaRegulator->append(t, u);
                double minU = u, maxU = u;
                updateRange(seriaRegulator, minU, maxU);
                double uMargin = (maxU - minU) * 0.1; if (uMargin < 1e-6) uMargin = 1.0;
                regY->setRange(minU - uMargin, maxU + uMargin);
                regX->setRange(t > windowTime ? t - windowTime : 0, t);
            });
    
    MainchartView->setMinimumSize(600, 400);
    Mainchart->addSeries(seriaZad);
    Mainchart->addSeries(seriaRegulowana);
    Mainchart->setTitle("Zadana i Regulowana");
    
    // Oś X i Y dla głównego wykresu
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
    Mainchart->setAnimationOptions(QChart::SeriesAnimations);
    ui->horizontalLayout_5->addWidget(MainchartView);

    //----PID----//
    QChart *PIDchart = new QChart();
    QChartView *PIDchartView = new QChartView(PIDchart);
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
    PIDchart->setAnimationOptions(QChart::SeriesAnimations);
    ui->horizontalLayout_4->addWidget(PIDchartView, 1);

    //----Uchyb----//    
    QChart *Uchybchart = new QChart();
    QChartView *UchybchartView = new QChartView(Uchybchart);
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
    Uchybchart->setAnimationOptions(QChart::SeriesAnimations);
    
    ui->horizontalLayout_4->addWidget(UchybchartView, 1);

    //----REGULATOR----//
    QChart *Regulatorchart = new QChart();
    QChartView *RegulatorchartView = new QChartView(Regulatorchart);
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
    Regulatorchart->setAnimationOptions(QChart::SeriesAnimations);
    
    ui->horizontalLayout_4->addWidget(RegulatorchartView, 1);
    RegulatorchartView->setMinimumSize(0, 300);
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
    ui->spinBOX_Interwal->setMinimum(1);
    ui->spinBOX_Interwal->setMaximum(9999);
    ui->spinBOX_Interwal->setSingleStep(1);
    ui->spinBOX_Interwal->setDecimals(5);


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

    symulator.setPID_TypCalki(RegulatorPID::ZERO);
    ui->radio_przed->setChecked(false);
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
                              double szum)
{
    aktualnyWektorA = a;
    aktualnyWektorB = b;
    aktualneOpoznienie = opoznienie;
    aktualnySzum = szum;
    symulator.setARX(a, b, opoznienie, szum);
}

void MainWindow::on_Konf_ARX_Button_clicked()
{
    if (!arxwindow) {
        arxwindow = new ARXwindow(this);
    }
    connect(arxwindow, &ARXwindow::zatwierdzonoARX, this, &MainWindow::ustawARXDane);
    arxwindow->ustawDane(aktualnyWektorA, aktualnyWektorB, aktualneOpoznienie, aktualnySzum);
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
        ui->spinBox_StalaSkladniowa->value(),
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
        ui->spinBOX_StalaSkladniowa->setValue(StalaSkladniowa);
        ui->spinBOX_Wypelnienie->setValue(Wypelnienie);
        ui->spinBOX_Interwal->setValue(interwalMs);

        QMessageBox::information(this, "sukces", "konfiguracja wczytana");
    } else {
        QMessageBox::warning(this, "blad", "nie udało sie wczytac");
    }
}


