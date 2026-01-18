#include "mainwindow.h"
#include "Symulator.h"
#include "arxwindow.h"
#include "qvalueaxis.h"
#include "ui_mainwindow.h"

#include <QtCharts/QChart>
#include <QtCharts/QChartView>
#include <QtCharts/QLineSeries>

#include <QDebug>

MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent)
    , ui(new Ui::MainWindow)

    , symulator(GeneratorSygnalu(), RegulatorPID(), ModelARX({0.5}, {0.5}))

{
    ui->setupUi(this);
    this->showMaximized();


    seriaP = new QLineSeries();
    seriaI = new QLineSeries();
    seriaD = new QLineSeries();

    seriaUchyb = new QLineSeries();
    seriaRegulator = new QLineSeries();
    seriaZad = new QLineSeries();
    seriaRegulowana = new QLineSeries();

    //----MAIN----//
    QChart *Mainchart = new QChart();
    QChartView *MainchartView = new QChartView(Mainchart);

    connect(&symulator, &SymulatorUAR::krokWykonany,
            this, [=](double w, double y, double e, double u,
                int k){
                double t = k * symulator.getInterwalMs() / 1000.0;

        // ---- Mainchart ----
        seriaZad->append(t, w);
        seriaRegulowana->append(t, y);

        double minY = std::min(seriaZad->points().first().y(), seriaRegulowana->points().first().y());
        double maxY = std::max(seriaZad->points().first().y(), seriaRegulowana->points().first().y());

        for (auto &p : seriaZad->points()) {
            if (p.y() < minY) minY = p.y();
            if (p.y() > maxY) maxY = p.y();
        }
        for (auto &p : seriaRegulowana->points()) {
            if (p.y() < minY) minY = p.y();
            if (p.y() > maxY) maxY = p.y();
        }

        double margin = (maxY - minY) * 0.1;
        if (margin < 1e-6) margin = 1.0;
        mainY->setRange(minY - margin, maxY + margin);

        double window = 10.0;
        if (t > window)
            mainX->setRange(t - window, t);
        else
            mainX->setRange(0, window);

        // ---- Uchybchart ----
        seriaUchyb->append(t, e);

        double minUchyb = e;
        double maxUchyb = e;

        for (auto &p : seriaUchyb->points()) {
            if (p.y() < minUchyb) minUchyb = p.y();
            if (p.y() > maxUchyb) maxUchyb = p.y();
        }

        double marginU = (maxUchyb - minUchyb) * 0.1;
        if (marginU < 1e-6) marginU = 1.0;

        // załóżmy, że masz w klasie MainWindow uchybX i uchybY
        uchybY->setRange(minUchyb - marginU, maxUchyb + marginU);
        uchybX->setRange(t > window ? t - window : 0, t);

        // ---- Regulatorchart ----
        seriaRegulator->append(t, u);

        double minReg = u;
        double maxReg = u;

        for (auto &p : seriaRegulator->points()) {
            if (p.y() < minReg) minReg = p.y();
            if (p.y() > maxReg) maxReg = p.y();
        }

        double marginR = (maxReg - minReg) * 0.1;
        if (marginR < 1e-6) marginR = 1.0;

        // zakładam, że masz regX i regY jako osie dla wykresu Regulator
        regY->setRange(minReg - marginR, maxReg + marginR);
        regX->setRange(t > window ? t - window : 0, t);
    });

    MainchartView->setMinimumSize(600, 400);

    Mainchart->addSeries(seriaZad);
    Mainchart->addSeries(seriaRegulowana);
    Mainchart->setTitle("Zadana i Regulowana");

    mainX = new QValueAxis();
    mainY = new QValueAxis();

    mainX->setTitleText("Czas [s]");
    mainY->setTitleText("Wartość");

    Mainchart->addAxis(mainX, Qt::AlignBottom);
    Mainchart->addAxis(mainY, Qt::AlignLeft);

    seriaZad->attachAxis(mainX);
    seriaZad->attachAxis(mainY);
    seriaRegulowana->attachAxis(mainX);
    seriaRegulowana->attachAxis(mainY);

    ui->horizontalLayout_5->addWidget(MainchartView);

    //----PID----//
    QChart *PIDchart = new QChart();
    QChartView *PIDchartView = new QChartView(PIDchart);




    PIDchart->addSeries(seriaP);
    PIDchart->addSeries(seriaI);
    PIDchart->addSeries(seriaD);
    PIDchart->setTitle("Sładowe Sterowania PID");

    QValueAxis *axisX = new QValueAxis();
    QValueAxis *axisY = new QValueAxis();

    axisX->setTitleText("Czas [s]");
    axisY->setTitleText("Wartość");

    PIDchart->addAxis(axisX, Qt::AlignBottom);
    PIDchart->addAxis(axisY, Qt::AlignLeft);

    seriaP->attachAxis(axisX);
    seriaP->attachAxis(axisY);
    seriaI->attachAxis(axisX);
    seriaI->attachAxis(axisY);
    seriaD->attachAxis(axisX);
    seriaD->attachAxis(axisY);
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

    seriaUchyb->attachAxis(uchybX);
    seriaUchyb->attachAxis(uchybY);

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

    seriaRegulator->attachAxis(regX);
    seriaRegulator->attachAxis(regY);

    ui->horizontalLayout_4->addWidget(RegulatorchartView, 1);
    RegulatorchartView->setMinimumSize(0, 300);
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
}

void MainWindow::on_spinBOX_Amplituda_editingFinished()
{
    symulator.setGeneratorA(ui->spinBOX_Amplituda->value());
}

void MainWindow::on_spinBOX_Czstotliwosc_editingFinished()
{
    symulator.setGeneratorTRZ(ui->spinBOX_Czstotliwosc->value());
}

void MainWindow::on_spinBOX_Td_editingFinished()
{
    symulator.setPID_Td(ui->spinBOX_Td->value());
}

void MainWindow::on_spinBOX_Ti_editingFinished()
{
    symulator.setPID_Ti(ui->spinBOX_Ti->value());
}

void MainWindow::on_spinBOX_Interwal_editingFinished()
{
    symulator.setGeneratorTT(ui->spinBOX_Interwal->value());
}

void MainWindow::on_radio_przed_toggled(bool checked)
{
    if (checked)
        symulator.setPID_TypCalki(RegulatorPID::Wew);
}

void MainWindow::on_radio_pod_toggled(bool checked)
{
    if (checked)
        symulator.setPID_TypCalki(RegulatorPID::Zew);
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
}

void MainWindow::ustawARXDane(const std::vector<double>& a,
                              const std::vector<double>& b,
                              int opoznienie,
                              double szum)
{

        symulator.setARX(a, b, opoznienie, szum);

}

void MainWindow::on_Konf_ARX_Button_clicked()
{
    if (!arxwindow) {
        arxwindow = new ARXwindow(this);
    }
    connect(arxwindow, &ARXwindow::zatwierdzonoARX,
            this, &MainWindow::ustawARXDane);

    arxwindow->show();
    arxwindow->raise();
    arxwindow->activateWindow();
}

