#ifndef ARXWINDOW_H
#define ARXWINDOW_H

#include <QMainWindow>
//DZIALA

namespace Ui {
class ARXwindow;
}

class ARXwindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit ARXwindow(QWidget *parent = nullptr);
    ~ARXwindow();

signals:
    void zatwierdzonoARX(const std::vector<double>& a,
                         const std::vector<double>& b,
                         int opoznienie,
                         double szum);
private slots:

    void on_Add_Button_wektorA_clicked();

    void on_Remove_Button_wektorA_clicked();

    void on_Add_Button_wektorB_clicked();

    void on_Remove_Button_wektorB_clicked();

    void on_Zatwierdz_Button_clicked();

private:
    Ui::ARXwindow *ui;
};

#endif // ARXWINDOW_H
