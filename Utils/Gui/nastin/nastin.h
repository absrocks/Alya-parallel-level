#ifndef NASTIN_H
#define NASTIN_H

#include <QToolBox>
#include <QSettings>

namespace Ui {
class Nastin;
}

class Nastin : public QToolBox
{
    Q_OBJECT
    
public:
    explicit Nastin(QWidget *parent = 0);
    std::string saveFile();
    void saveGUIValues(QString filePath);
    void loadGUIValues(QString filePath);
    ~Nastin();
    
private slots:
    void on_nastinRegCombo_currentIndexChanged(const QString &arg1);

    void on_nastinTurbModelCombo1_currentIndexChanged(const QString &arg1);

    void on_addPostProcessButton_clicked();

    void on_delPostProcessButton_clicked();

    void on_addElementSetButton_clicked();

    void on_delElementSetButton_clicked();

    void on_addNodeSetButton_clicked();

    void on_delNodeSetButton_clicked();

    void on_addBoundarySetButton_clicked();

    void on_delBoundarySetButton_clicked();

    void on_addWitnessPointsButton_clicked();

    void on_delWitnessPointsButton_clicked();

    void on_nastinTempCoupCombo_currentIndexChanged(const QString &arg1);

    void on_nastinLeveCoupActiveCheck_clicked(bool checked);

    void on_nastinSurfTensActiveCheck_clicked(bool checked);

    void on_nastinStabilizationCombo1_currentIndexChanged(const QString &arg1);

    void on_nastinTrSubScTimeCheck_clicked(bool checked);

    void on_nastinTrSubScConvCheck_clicked(bool checked);

    void on_nastinLinearitzationCombo_currentIndexChanged(const QString &arg1);

private:
    Ui::Nastin *ui;
};

#endif // NASTIN_H
