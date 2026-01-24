#include "arxwindow.h"
#include "ui_arxwindow.h"

ARXwindow::ARXwindow(QWidget *parent)
    : QMainWindow(parent)
    , ui(new Ui::ARXwindow)
{
    ui->setupUi(this);

    // Ustawienia minimalnych i maksymalnych wartości SpinBoxów
    ui->SpinBox_WektorA->setMinimum(-1000.0);
    ui->SpinBox_WektorA->setMaximum(1000.0);
    ui->SpinBox_WektorA->setSingleStep(0.01);
    ui->SpinBox_WektorA->setDecimals(4);

    ui->SpinBox_WektorB->setMinimum(-1000.0);
    ui->SpinBox_WektorB->setMaximum(1000.0);
    ui->SpinBox_WektorB->setSingleStep(0.01);
    ui->SpinBox_WektorB->setDecimals(4);

    ui->spinBoxOpoznienie->setMinimum(1);
    ui->spinBoxOpoznienie->setMaximum(100);
    ui->spinBoxOpoznienie->setSingleStep(1);



    ui->doubleSpinBox_sigma->setMinimum(0.0);
    ui->doubleSpinBox_sigma->setMaximum(100.0);
    ui->doubleSpinBox_sigma->setSingleStep(0.01);
    ui->doubleSpinBox_sigma->setDecimals(2);
}
void ARXwindow::ustawDane(const std::vector<double>& wektorA,
                          const std::vector<double>& wektorB,
                          int opoznienie,
                          double sigma)
{
    ui->listWidget_wektorA->clear();
    ui->listWidget_wektorB->clear();

    for (double val : wektorA) {
        QListWidgetItem *item = new QListWidgetItem(QString::number(val, 'f', 4));
        ui->listWidget_wektorA->addItem(item);
    }

    for (double val : wektorB) {
        QListWidgetItem *item = new QListWidgetItem(QString::number(val, 'f', 4));
        ui->listWidget_wektorB->addItem(item);
    }

    ui->spinBoxOpoznienie->setValue(opoznienie);

    ui->doubleSpinBox_sigma->setValue(sigma);
}

ARXwindow::~ARXwindow()
{
    delete ui;
}

void ARXwindow::on_Add_Button_wektorA_clicked()
{
    double val = ui->SpinBox_WektorA->value();

    QListWidgetItem *item = new QListWidgetItem(QString::number(val), ui->listWidget_wektorA);

    ui->listWidget_wektorA->addItem(item);
    ui->SpinBox_WektorA->clear();
}

void ARXwindow::on_Remove_Button_wektorA_clicked()
{
    QListWidgetItem *item = ui->listWidget_wektorA->takeItem(ui->listWidget_wektorA->currentRow());
    delete item;
}

void ARXwindow::on_Add_Button_wektorB_clicked()
{
    double val = ui->SpinBox_WektorB->value();

    QListWidgetItem *item = new QListWidgetItem(QString::number(val), ui->listWidget_wektorB);

    ui->listWidget_wektorB->addItem(item);
    ui->SpinBox_WektorB->clear();
}

void ARXwindow::on_Remove_Button_wektorB_clicked()
{
    QListWidgetItem *item = ui->listWidget_wektorB->takeItem(ui->listWidget_wektorB->currentRow());
    delete item;
}

void ARXwindow::on_Zatwierdz_Button_clicked()
{
    // Pobranie wektorów
    std::vector<double> wektorA, wektorB;
    for (int i = 0; i < ui->listWidget_wektorA->count(); ++i)
        wektorA.push_back(ui->listWidget_wektorA->item(i)->text().toDouble());
    for (int i = 0; i < ui->listWidget_wektorB->count(); ++i)
        wektorB.push_back(ui->listWidget_wektorB->item(i)->text().toDouble());

    // Pobranie spinboxów
    int uchyb = static_cast<int>(ui->spinBoxOpoznienie->value());
    double sygma = ui->doubleSpinBox_sigma->value();

    // Emitowanie sygnału
    emit zatwierdzonoARX(wektorA, wektorB, uchyb, sygma);

    // Zamknięcie okna (opcjonalnie)
    this->close();
}
