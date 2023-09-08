#ifndef ALYAWIDGET_H
#define ALYAWIDGET_H

#include <QStackedWidget>
#include <QProcess>
#include <QFileDialog>
namespace Ui {
class alyaWidget;
}

class alyaWidget : public QStackedWidget
{
    Q_OBJECT
    
public:
    explicit alyaWidget(QWidget *parent = 0);
    ~alyaWidget();
    void saveGUI(bool saveAs);
    bool loadGUI();
    void newGUI();
    void exportGUI();
    
private slots:

    void on_continueButton_clicked();

    void on_modulesServicesButton_clicked();

    void on_runSimulationButton_clicked();

    void on_pushButton_2_clicked();

    void on_addModule_clicked();

    void on_delModule_clicked();

    void on_addService_clicked();

    void on_delService_clicked();

    void on_importDomainButton_clicked();

    void on_runButton_clicked();
    void printOutput();
    void updateProgress();
    void printError();
    void updateExit();
    void updateExitAlya2pos();

    void on_openResultButton_clicked();

    void on_updateResultButton_clicked();

    void on_modulesServicesButton_2_clicked();

signals:
   void setViewer();

private:
    Ui::alyaWidget *ui;
    QString path;
    QString domainFilePath;
    QString domainFileName;
    QString saveFilePath;
    QString output;
    QProcess *myProcess;
    bool isDomainImported;
    bool checkBasicData();
    bool importDomain(QString path, QString domPath);
    void updateGeneralParamSubGroups();
    void copyAlya2posConf(QString path);
    bool getCodes(QString domPath);
    void exportToAlya(QString dest);
    void emitSetViewer();

public:
   static QDir tmpDir;
};

#endif // ALYAWIDGET_H
