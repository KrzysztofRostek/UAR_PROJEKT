#include "arxwindow.h"
#include "ui_arxwindow.h"

ARXwindow::ARXwindow(QWidget *parent)
    : QMainWindow(parent)
    , ui(new Ui::ARXwindow)
{
    ui->setupUi(this);
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
