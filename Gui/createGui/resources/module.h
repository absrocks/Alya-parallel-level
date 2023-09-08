#ifndef <<defname>>_H
#define <<defname>>_H

#include <QToolBox>
#include <QSettings>
#include "module.h"

namespace Ui {
class <<className>>;
}

class <<className>> : public QToolBox, public Module
{
    Q_OBJECT
    
public:
    explicit <<className>>(QWidget *parent = 0);
    std::string saveFile();
    std::string fileName();
    void saveGUIValues(QString filePath);
    void loadGUIValues(QString filePath);
    void disableModulesServicesParams();
    void enableModule(QString moduleName);
    void enableService(QString serviceName);
    void loadCodes(QMap<QString, QString> boundaryCodes, QMap<QString, QString> nodeCodes);
    ~<<className>>();
    
private slots:
    <<privateSlots>>
//    void on_gravityCheckBox_clicked(bool checked);

//    void on_helpButton_clicked();

private:
    Ui::<<className>> *ui;
    <<privateMethods>>
    void  setValidators();
    QValidator *getValidator(QString validatorName);
};

#endif // <<defname>>_H
