#ifndef ALYAWIDGET_H
#define ALYAWIDGET_H

#include <QStackedWidget>

namespace Ui {
class alyaWidget;
}

class alyaWidget : public QStackedWidget
{
    Q_OBJECT
    
public:
    explicit alyaWidget(QWidget *parent = 0);
    ~alyaWidget();
    void saveGUI();
    void loadGUI();
    void newGUI();
    
private slots:

    void on_continueButton_clicked();

    void on_modulesServicesButton_clicked();

    void on_runSimulationButton_clicked();

    void on_pushButton_2_clicked();

    void on_openMesh_clicked();

    void on_saveButton_clicked();

    void on_addModule_clicked();

    void on_delModule_clicked();

    void on_addService_clicked();

    void on_delService_clicked();

private:
    Ui::alyaWidget *ui;
    QString path;
    bool checkBasicData();
};

#endif // ALYAWIDGET_H
