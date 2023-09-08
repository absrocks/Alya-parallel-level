#ifndef MODULE_H
#define MODULE_H
#include <QMap>
#include <QDoubleValidator>
#include <QIntValidator>

namespace Ui {
class Module;
}

class Module
{   
public:
    virtual std::string saveFile() = 0;
    virtual std::string fileName() = 0;
    virtual void saveGUIValues(QString filePath) = 0;
    virtual void loadGUIValues(QString filePath) = 0;
    virtual void loadCodes(QMap<QString, QString> boundaryCodes, QMap<QString, QString> nodeCodes) = 0;
    QDoubleValidator *realValidator;
    QIntValidator *intValidator;
    QRegExpValidator *realseqValidator;
};

#endif // MODULE_H
