#ifndef MODULE_H
#define MODULE_H

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
};

#endif // MODULE_H
