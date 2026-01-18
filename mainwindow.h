#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include "Symulator.h"
#include "arxwindow.h"
#include "qvalueaxis.h"

#include <QTimer>

#include <QtCharts/QChart>
#include <QtCharts/QChartView>
#include <QtCharts/QLineSeries>

QT_BEGIN_NAMESPACE
namespace Ui {
class MainWindow;
}
QT_END_NAMESPACE

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    MainWindow(QWidget *parent = nullptr);
    ~MainWindow();

private slots:
    void on_Sin_Button_clicked();

    void on_Square_Button_clicked();

    void on_spinBOX_WzmocK_editingFinished();

    void on_spinBOX_Amplituda_editingFinished();

    void on_spinBOX_Czstotliwosc_editingFinished();

    void on_spinBOX_Td_editingFinished();

    void on_spinBOX_Ti_editingFinished();

    void on_spinBOX_Interwal_editingFinished();

    void on_radio_przed_toggled(bool checked);

    void on_radio_pod_toggled(bool checked);

    void on_Reset_d_clicked();

    void on_Reset_i_clicked();

    void on_START_Button_clicked();

    void on_Konf_ARX_Button_clicked();

    void ustawARXDane(const std::vector<double>& a,
                      const std::vector<double>& b,
                      int opoznienie,
                      double szum);


    void on_STOP_Bttun_clicked();

    void on_RESET_Button_clicked();

    void on_TEST_BUTTON_clicked();

private:
    Ui::MainWindow *ui;

    ARXwindow *arxwindow = nullptr;

    SymulatorUAR symulator;

    QLineSeries *seriaP;
    QLineSeries *seriaI;
    QLineSeries *seriaD;
    QLineSeries *seriaUchyb;
    QLineSeries *seriaRegulator;
    QLineSeries *seriaRegulowana;
    QLineSeries *seriaZad;

    QValueAxis *mainX;
    QValueAxis *mainY;

    QValueAxis *uchybX;
    QValueAxis *uchybY;

    QValueAxis *regX;
    QValueAxis *regY;


    QTimer *timer;
    double T; // okres próbkowania [s]
};
#endif // MAINWINDOW_H
