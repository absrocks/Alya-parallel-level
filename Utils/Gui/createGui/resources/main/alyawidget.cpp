#include "alyawidget.h"
#include "ui_alyawidget.h"
#include <QFileDialog>
#include <QTextStream>
#include <sstream>
#include <string>
#include <iostream>
#include <QMessageBox>
#include "module.h"
#include "pqLoadDataReaction.h"
#include <pqSaveStateReaction.h>
#include <pqLoadStateReaction.h>
#include <pqApplicationCore.h>
#include <pqDeleteReaction.h>

QDir alyaWidget::tmpDir = QDir(QDir::temp().homePath() + "/.ParaViewAlya");

alyaWidget::alyaWidget(QWidget *parent) :
    QStackedWidget(parent),
    ui(new Ui::alyaWidget)
{
    ui->setupUi(this);

    this->setCurrentIndex(1);
    this->isDomainImported = false;

    //ui->alyaProblemTabs->setCurrentIndex(2);
}

alyaWidget::~alyaWidget()
{
    delete ui;
}

bool removeFolder(const QString & dirName)
{
    bool result;
    QDir dir(dirName);

    if (dir.exists(dirName)) {
        Q_FOREACH(QFileInfo info, dir.entryInfoList(QDir::NoDotAndDotDot | QDir::System | QDir::Hidden  | QDir::AllDirs | QDir::Files, QDir::DirsFirst)) {
            if (info.isDir()) {
                result = removeFolder(info.absoluteFilePath());
            }
            else {
                result = QFile::remove(info.absoluteFilePath());
            }

            if (!result) {
                return result;
            }
        }
        result = dir.rmdir(dirName);
    }
    return result;
}

bool alyaWidget::checkBasicData() {
    bool valid = true;
    QString message;

    if (valid && ui->problemNameEdit->text() == "") {
        valid = false;
        message = "You need to insert a problem name";
    }
    else if(valid && ui->problemModules->count() == 0) {
        valid = false;
        message = "You need to assign at least one module";
    }
    else if(valid && !this->isDomainImported) {
        valid = false;
        message = "You need to import a problem domain";
    }

    if (!valid) {
        QMessageBox::warning(this,"Warning",message);
    }

    return valid;

}

void alyaWidget::updateGeneralParamSubGroups()
{
     <<disableModulesServicesParams>>
    //disable and hide all module and services edits
    for (int i = 0; i < ui->problemModules->count(); i++) {
        QString text = ui->problemModules->item(i)->text();
        <<enableModule>>
     }
    for (int i = 0; i < ui->problemServices->count(); i++) {
        QString text = ui->problemServices->item(i)->text();
        <<enableService>>
     }
}

void alyaWidget::on_continueButton_clicked()
{
    if (checkBasicData()) {
        //remove all the tabs
        for (int i = ui->alyaProblemTabs->count() - 1; i >= 0; i--) {
            ui->alyaProblemTabs->removeTab(i);
        }
	//activate kernel modules
	<<activeKernelModule>>
	this->updateGeneralParamSubGroups();
        //activate each module tab selected
        for (int i = 0; i < ui->problemModules->count(); i++) {
            <<activeModule>>
        }
        //get the codes from the mesh and loads to the gui
        getCodes(alyaWidget::tmpDir.absolutePath() + "/domain/domain.dom.dat");
        //disable each module tab
        this->setCurrentIndex(0);
    }

}

void alyaWidget::on_modulesServicesButton_clicked()
{
    this->setCurrentIndex(1);
}

void alyaWidget::on_runSimulationButton_clicked()
{
    //Create the alya problem files in the temporal folder
    removeFolder(alyaWidget::tmpDir.absolutePath() + "/problem");
    QDir().mkdir(alyaWidget::tmpDir.absolutePath() + "/problem");
    exportToAlya(alyaWidget::tmpDir.absolutePath() + "/problem");
    this->setCurrentIndex(2);
}

void alyaWidget::on_pushButton_2_clicked()
{
    this->setCurrentIndex(0);
}

bool alyaWidget::importDomain(QString path, QString domPath) {
    bool imported = true;

    if (!path.isNull() && !domPath.isNull()) {
        //copy the file
        removeFolder(path + "/domain");
        QDir().mkdir(path + "/domain");
        QFile::copy(domPath, path + "/domain/" + "domain.dom.dat");

        //open file and load the includes
        QString search = "/";
        QString auxPath = domPath.left(domPath.lastIndexOf(search));
        QFile file(domPath);
        if (!file.open(QIODevice::ReadOnly | QIODevice::Text))
        return imported;
        QTextStream in(&file);
        while (!in.atEnd()) {
            QString line = in.readLine();
            if (line.contains("INCLUDE")) {
                QString fileName;
                if (line.indexOf("/") != -1) {
                    fileName = line.right(line.length() - line.indexOf("/") - 1);
                }
                else {
                    fileName = line.right(line.length() - line.indexOf("INCLUDE") - 7);
                }
                fileName = fileName.trimmed();
                search = ".";
                QString orig = auxPath + "/" + fileName;
                //QString dest = path + "/" + ui->problemNameEdit->text() + extension;
                QString dest = path + "/domain/" + fileName;
                orig = orig.trimmed();
                dest = dest.trimmed();
                QFile copyFile(orig);
                QFile origFile(dest);
                if (origFile.exists()) {
                    origFile.remove();
                }
                bool result = copyFile.copy(dest);
                if (!result) {
                    imported = false;
                    std::cout << "error:" << copyFile.error() << " " << copyFile.errorString().toStdString() << std::endl;
                }
                copyFile.flush();
                copyFile.close();
            }
        }
    }

    return imported;
}

bool alyaWidget::getCodes(QString domPath) {
    bool imported = true;
    if (!domPath.isNull()) {        
        //open file and load the includes
        QString search = "/";
        QString auxPath = domPath.left(domPath.lastIndexOf(search));
        QFile file(domPath);
        if (!file.open(QIODevice::ReadOnly | QIODevice::Text))
            return imported;
        QTextStream in(&file);
        bool startBoundary = false;
        bool endBoundary = false;
        QString content;
        //get the boundary_conditions content
        QRegExp begin("^\\s*\\t*BOUNDARY_CONDITIONS");
        QRegExp finish("^\\s*\\t*END_BOUNDARY_CONDITIONS");
        while (!in.atEnd() && !endBoundary) {
            QString line = in.readLine();            
            if (line.contains(begin)) {
                startBoundary = true;
                line = in.readLine();
            }
            else if (line.contains(finish)) {
                endBoundary = true;
            }
            if (startBoundary && !endBoundary) {
                content = content + line + "#";
            }
        }

       //-If its an include gets the content from the external file
      QRegExp include("^\\s*\\t*INCLUDE");
       if (content.contains(include)) {
           QString fileName;
           if (content.indexOf("/") != -1) {
               fileName = content.right(content.length() - content.indexOf("/") - 1);
           }
           else {
               fileName = content.right(content.length() - content.indexOf("INCLUDE") - 7);
           }
           fileName = fileName.replace("#","");
           fileName = fileName.trimmed();
           search = ".";
           QString orig = auxPath + "/" + fileName;
           QFile file2(orig);
           if (!file2.open(QIODevice::ReadOnly | QIODevice::Text))
               return imported;
           QTextStream inOrig(&file2);
           content = "";
           while (!inOrig.atEnd()) {
               QString line = inOrig.readLine();
               content = content + line + "#";
           }
       }
       //get the on boundaries
       QRegExp onBoundaries("ON_BOUNDARIES(.+)END_ON_BOUNDARIES");
       QRegExp code("\\d+\\s+(\\d+)");
       QMap<QString, QString> onBoundariesList;
       if (onBoundaries.indexIn(content) != -1) {
           QString codes = onBoundaries.cap(1);
           QStringList codesList = codes.split("#");
           for (int i = 0; i < codesList.size(); i++) {
               if (code.indexIn(codesList[i]) != -1) {
                   QString cap = code.cap(1);
                   onBoundariesList[cap] = cap;
               }
           }
       }
       //get the on nodes
       QRegExp onNodes("ON_NODES(.+)END_ON_NODES");
       QMap<QString, QString> onNodesList;
       if (onNodes.indexIn(content) != -1) {
           QString codes = onNodes.cap(1);
           QStringList codesList = codes.split("#");
           for (int i = 0; i < codesList.size(); i++) {
               if (code.indexIn(codesList[i]) != -1) {
                   QString cap = code.cap(1);                  
                   onNodesList[cap] = cap;
               }
           }
       }

       for (int i = 0; i < ui->alyaProblemTabs->count(); i++) {
           QWidget* pWidget= ui->alyaProblemTabs->widget(i);
           dynamic_cast<Module*>(pWidget)->loadCodes(onBoundariesList, onNodesList);
       }
    }
    return imported;
}

void alyaWidget::copyAlya2posConf(QString path) {

    std::string filePath = path.toStdString() + "/" + ui->problemNameEdit->text().toStdString() + ".post.alyadat";
    QFile file(QString::fromStdString(filePath));
    file.open(QIODevice::WriteOnly | QIODevice::Text);
    QTextStream out(&file);
    out << "$-------------------------------------------------------------------\n";
    out << "DATA\n";
    out << "  FORMAT:                   ensight\n";
    out << "  MARK_ELEMENTS:            type\n";
    out << "  ELIMINATE_BOUNDARY_NODES: yes\n";
    out << "  MULTIPLE_FILE:            OFF\n";
    out << "  BOUNDARY:                 ON\n";
    out << "  SUBDOMAINS, ALL\n";
    out << "  END_SUBDOMAINS\n";
    out << "END_DATA\n";
    out << "$-------------------------------------------------------------------\n";

    file.close();
}

void alyaWidget::exportToAlya(QString dest) {

    for (int i = 0; i < ui->alyaProblemTabs->count(); i++) {
        QWidget* pWidget= ui->alyaProblemTabs->widget(i);
        std::string fileData = dynamic_cast<Module*>(pWidget)->saveFile();
        std::string fileName = dynamic_cast<Module*>(pWidget)->fileName();
        std::string filePath = dest.toStdString() + "/" + ui->problemNameEdit->text().toStdString() + fileName;
        QFile file(QString::fromStdString(filePath));
        if (file.open(QFile::WriteOnly | QFile::Truncate)) {
            QTextStream out(&file);
            out << fileData.c_str();
            file.close();
        }
    }
    //export the domain data
    QDir dirDomain(alyaWidget::tmpDir.absolutePath() + "/domain");
    QStringList domainFiles = dirDomain.entryList();
    for (int i = 0; i < domainFiles.size(); ++i) {
        std::cout << "File " << i << domainFiles.at(i).toStdString() << std::endl;
        QString domainFileSour = alyaWidget::tmpDir.absolutePath() + "/domain/" + domainFiles.at(i);
        QString domainFileDest = dest + "/" + domainFiles.at(i);
        if (domainFiles.at(i) == "domain.dom.dat") {
            domainFileDest = dest + "/" + ui->problemNameEdit->text() + ".dom.dat";
        }
        QFile::copy(domainFileSour, domainFileDest);
    }
    //copy the alya2pos conf file, in order to transform alya output to ensight
    copyAlya2posConf(dest);
}

void alyaWidget::newGUI() {

    //remove all the tabs
    for (int i = ui->alyaProblemTabs->count() - 1; i >= 0; i--) {
        ui->alyaProblemTabs->removeTab(i);
    }
    //init modules
    <<newModule>>
    ui->alyaProblemTabs->update();

    //new problem main settings
    ui->problemNameEdit->setText("");
    for (int i = ui->problemModules->count() - 1; i >= 0; i--) {
        int rowIdx = ui->alyaModules->count();
        ui->alyaModules->insertItem(rowIdx,ui->problemModules->item(i)->text());
        //remove the item from the list
        ui->problemModules->takeItem(i);
    }

    for (int i = ui->problemServices->count() - 1; i >= 0; i--) {
        int rowIdx = ui->alyaServices->count();
        ui->alyaServices->insertItem(rowIdx,ui->problemServices->item(i)->text());
        //remove the item from the list
        ui->problemServices->takeItem(i);
    }
    this->isDomainImported = false;
    ui->domainDefinedCheck->setChecked(false);

    this->setCurrentIndex(1);

    removeFolder(alyaWidget::tmpDir.absolutePath());
    alyaWidget::tmpDir.mkpath(".");
}

void alyaWidget::saveGUI(bool saveAs) {
    QDir directory;
    QString filePath;
    if (saveAs) {
        filePath = QFileDialog::getSaveFileName(this, tr("Save as Alya problem"), directory.homePath(),
                                                        "Alya problem (*.alp)");
    }
    else if (saveFilePath.isNull()) {
        filePath = QFileDialog::getSaveFileName(this, tr("Save Alya problem"), directory.homePath(),
                                                        "Alya problem (*.alp)");
        saveFilePath = filePath;
    }
    else {
        filePath = saveFilePath;
    }

    if (!filePath.isNull()) {
        if (!filePath.endsWith(".alp")) {
            filePath = filePath + ".alp";
        }
        //save modules
        for (int i = 0; i < ui->alyaProblemTabs->count(); i++) {
            QWidget* pWidget= ui->alyaProblemTabs->widget(i);
            dynamic_cast<Module*>(pWidget)->saveGUIValues(filePath);
        }

        //save modulesandservices
        QSettings settings(filePath, QSettings::IniFormat);
        settings.beginGroup("modulesServices");
        settings.setValue("problemNameEdit", ui->problemNameEdit->text());
        QString values;
        for (int i = 0; i < ui->problemModules->count(); i++) {
             values = values + ui->problemModules->item(i)->text();
             if (i < (ui->problemModules->count() - 1)) {
                 values = values + "#";
             }
        }
        if (values != "") {
            settings.setValue("problemModules", values);
        }
        values = "";
        for (int i = 0; i < ui->alyaModules->count(); i++) {
             values = values + ui->alyaModules->item(i)->text();
             if (i < (ui->alyaModules->count() - 1)) {
                 values = values + "#";
             }
        }
        if (values != "") {
            settings.setValue("alyaModules", values);
        }
        values = "";
        for (int i = 0; i < ui->problemServices->count(); i++) {
             values = values + ui->problemServices->item(i)->text();
             if (i < (ui->problemServices->count() - 1)) {
                 values = values + "#";
             }
        }
        if (values != "") {
            settings.setValue("problemServices", values);
        }
        values = "";
        for (int i = 0; i < ui->alyaServices->count(); i++) {
             values = values + ui->alyaServices->item(i)->text();
             if (i < (ui->alyaServices->count() - 1)) {
                 values = values + "#";
             }
        }
        if (values != "") {
            settings.setValue("alyaServices", values);
        }
        settings.setValue("isDomainImported",this->isDomainImported);
        settings.setValue("domainDefinedCheck",ui->domainDefinedCheck->isChecked());
        settings.setValue("domainFilePath",this->domainFilePath);
        settings.setValue("path",this->path);
        settings.endGroup();

        //copy the domain data
        QString savePath = filePath.left(filePath.lastIndexOf("/"));
        removeFolder(savePath + "/domain");
        QDir().mkdir(savePath + "/domain");
        QDir dirDomain(alyaWidget::tmpDir.absolutePath() + "/domain");
        QStringList domainFiles = dirDomain.entryList();
        for (int i = 0; i < domainFiles.size(); ++i) {
            QString domainFileSour = alyaWidget::tmpDir.absolutePath() + "/domain/" + domainFiles.at(i);
            QString domainFileDest = savePath + "/domain/" + domainFiles.at(i);
            QFile::copy(domainFileSour, domainFileDest);
        }
    }
}

bool alyaWidget::loadGUI() {
    QDir directory;
    QString filePath = QFileDialog::getOpenFileName(this, tr("Open Alya problem"), directory.homePath(),
                                                    "Alya problem (*.alp)");
    if (!filePath.isNull()) {
        saveFilePath = filePath;
        //load modulesservices
        QSettings settings(filePath, QSettings::IniFormat);
        settings.beginGroup("modulesServices");
        QString vals = settings.value("problemNameEdit","").toString();
        ui->problemNameEdit->setText(vals);
        vals = settings.value("problemModules","").toString();
        ui->problemModules->clear();
        if (vals != "") {
            QStringList valuesList = vals.split("#");
            for (int j = 0; j < valuesList.size(); j++) {
                int rowIdx = ui->problemModules->count();
                ui->problemModules->insertItem(rowIdx,valuesList[j]);
            }
        }

        vals = settings.value("alyaModules","").toString();
        ui->alyaModules->clear();
        if (vals != "") {
            QStringList valuesList = vals.split("#");
            for (int j = 0; j < valuesList.size(); j++) {
                int rowIdx = ui->alyaModules->count();
                ui->alyaModules->insertItem(rowIdx,valuesList[j]);
            }
        }

        vals = settings.value("problemServices","").toString();
        ui->problemServices->clear();
        if (vals != "") {
            QStringList valuesList = vals.split("#");
            for (int j = 0; j < valuesList.size(); j++) {
                int rowIdx = ui->problemServices->count();
                ui->problemServices->insertItem(rowIdx,valuesList[j]);
            }
        }

        vals = settings.value("alyaServices","").toString();
        ui->alyaServices->clear();
        if (vals != "") {
            QStringList valuesList = vals.split("#");
            for (int j = 0; j < valuesList.size(); j++) {
                int rowIdx = ui->alyaServices->count();
                ui->alyaServices->insertItem(rowIdx,valuesList[j]);
            }
        }

        this->updateGeneralParamSubGroups();

        this->isDomainImported = settings.value("isDomainImported",false).toBool();
        ui->domainDefinedCheck->setChecked(settings.value("domainDefinedCheck",false).toBool());
        this->domainFilePath = settings.value("domainFilePath", this->domainFilePath).toString();
        this->path = settings.value("path", "").toString();

        this->setCurrentIndex(1);
        on_continueButton_clicked();

        //load module values
        for (int i = 0; i < ui->alyaProblemTabs->count(); i++) {
            QWidget* pWidget= ui->alyaProblemTabs->widget(i);
            dynamic_cast<Module*>(pWidget)->loadGUIValues(filePath);
        }

        //copy the domain data
        QString savePath = alyaWidget::tmpDir.absolutePath();
        removeFolder(savePath + "/domain");
        QDir().mkdir(savePath + "/domain");
        QDir dirDomain(filePath.left(filePath.lastIndexOf("/")) + "/domain");
        QStringList domainFiles = dirDomain.entryList();
        for (int i = 0; i < domainFiles.size(); ++i) {
            QString domainFileSour = dirDomain.absolutePath() + "/" + domainFiles.at(i);
            QString domainFileDest = savePath + "/domain/" + domainFiles.at(i);
            QFile::copy(domainFileSour, domainFileDest);
        }
    }

    return !filePath.isNull();
}

void alyaWidget::exportGUI() {
    QDir directory;
    if (path.isNull()) {
        path = directory.path();
    }
    path = QFileDialog::getExistingDirectory (this, tr("Directory"), path);
    if ( path.isNull() == false )
    {
        exportToAlya(path);
    }
}

void alyaWidget::on_addModule_clicked()
{
    if (ui->alyaModules->currentItem()!=0 && ui->alyaModules->currentItem()->isSelected()) {
           //add the item to the table
           int rowIdx = ui->problemModules->count();
           ui->problemModules->insertItem(rowIdx,ui->alyaModules->currentItem()->text());

           //remove the item from the list
           ui->alyaModules->takeItem(ui->alyaModules->row(ui->alyaModules->currentItem()));
      }
}

void alyaWidget::on_delModule_clicked()
{
    if (ui->problemModules->currentItem()!=0 && ui->problemModules->currentItem()->isSelected()) {
            //add the item to the table
            int rowIdx = ui->alyaModules->count();
            ui->alyaModules->insertItem(rowIdx,ui->problemModules->currentItem()->text());

            //remove the item from the list
            ui->problemModules->takeItem(ui->problemModules->row(ui->problemModules->currentItem()));
    }

}

void alyaWidget::on_addService_clicked()
{
    if (ui->alyaServices->currentItem()!=0 && ui->alyaServices->currentItem()->isSelected()) {
           //add the item to the table
           int rowIdx = ui->problemServices->count();
           ui->problemServices->insertItem(rowIdx,ui->alyaServices->currentItem()->text());

           //remove the item from the list
           ui->alyaServices->takeItem(ui->alyaServices->row(ui->alyaServices->currentItem()));
      }
}

void alyaWidget::on_delService_clicked()
{
    if (ui->problemServices->currentItem()!=0 && ui->problemServices->currentItem()->isSelected()) {
            //add the item to the table
            int rowIdx = ui->alyaServices->count();
            ui->alyaServices->insertItem(rowIdx,ui->problemServices->currentItem()->text());

            //remove the item from the list
            ui->problemServices->takeItem(ui->problemServices->row(ui->problemServices->currentItem()));
    }
}

void alyaWidget::on_importDomainButton_clicked()
{
    QDir directory;
    QString filePath = QFileDialog::getOpenFileName(this, tr("Open Alya domain"), directory.homePath(),
                                                    "Alya domain (*.dom.dat)");
    this->domainFilePath = filePath;
    if (this->domainFilePath.isNull()) {
        this->isDomainImported = false;
        ui->domainDefinedCheck->setChecked(false);
    }
    else {        
        //Save the domain in the temporary folder
        //import the domain data
        bool imported = imported = importDomain(alyaWidget::tmpDir.absolutePath(), this->domainFilePath);
        if (imported) {
            QMessageBox::information(this,"Import domain", "Domain imported correctly");
        }
        else {
            QMessageBox::information(this,"Import domain", "Domain not imported to tmp folder");
        }
        this->isDomainImported = imported;
        ui->domainDefinedCheck->setChecked(imported);
    }
}

void alyaWidget::on_runButton_clicked()
{
    ui->runButton->setEnabled(false);
    ui->openResultButton->setEnabled(false);
    ui->updateResultButton->setEnabled(false);
    QString        program;
    QStringList    arguments;    
    myProcess = new QProcess( this );
    connect (myProcess, SIGNAL(readyReadStandardOutput()), this, SLOT(printOutput()));
    connect (myProcess, SIGNAL(readyReadStandardOutput()), this, SLOT(updateProgress()));
    connect (myProcess, SIGNAL(readyReadStandardError()), this, SLOT(printError()));
    connect (myProcess, SIGNAL(finished(int)), this, SLOT(updateExit()));

    if (ui->withMPICheck->isChecked()) {
        program = "mpirun";
        QString problemPath = alyaWidget::tmpDir.absolutePath() + "/problem/" + ui->problemNameEdit->text();
        QString numProc = "1";
        if (!ui->procNumberEdit->text().isEmpty()) {
            numProc = ui->procNumberEdit->text();
        }
        arguments << "-np" << numProc << "alyaBin/Alya.x" << problemPath;
    }
    else {
        program = "alyaBin/Alya.x";
        QString problemPath = alyaWidget::tmpDir.absolutePath() + "/problem/" + ui->problemNameEdit->text();
        arguments << problemPath;
    }

    ui->outputRunEdit->setOverwriteMode(false);
    ui->outputRunEdit->setText("");
    this->output = "";
    ui->progressBar->setValue(0);

    myProcess->start(program, arguments);
    myProcess->write("\n");
}

void alyaWidget::printOutput()
{
    //std::cout << "print out" << std::endl;
    QByteArray byteArray = myProcess->readAllStandardOutput();
    QString outputText = QString(byteArray);
    ui->outputRunEdit->append(outputText);
    this->output = this->output + outputText;
}

void alyaWidget::updateProgress()
{
    //std::cout << "print out" << std::endl;
    QByteArray byteArray = myProcess->readAllStandardOutput();
    QString outputText = QString(byteArray);
    ui->progressBar->setValue(ui->progressBar->value() + 1);
}

void alyaWidget::printError()
{
    //ui->e_Log->append("Got to printError()");
    std::cout << "print err" << std::endl;
    QByteArray byteArray = myProcess->readAllStandardError();
    QStringList strLines = QString(byteArray).split("\n");

    foreach (QString line, strLines){
        ui->outputRunEdit->append(line);
    }
}

//When the simulation is finished it launches the alya2pos to obtain ensight result
void alyaWidget::updateExit()
{    
    if (this->output.contains("--| ALYA  FINISHED NORMALLY")) {
        //exec the alya2pos
        QString        program = "alyaBin/alya2pos.x";
        QStringList    arguments;
        arguments << ui->problemNameEdit->text();
        QFileInfo programPath(program);

        myProcess = new QProcess( this );
        QString problemPath = alyaWidget::tmpDir.absolutePath() + "/problem";
        myProcess->setWorkingDirectory(problemPath);
        connect (myProcess, SIGNAL(readyReadStandardOutput()), this, SLOT(printOutput()));
        connect (myProcess, SIGNAL(readyReadStandardError()), this, SLOT(printError()));
        connect (myProcess, SIGNAL(finished(int)), this, SLOT(updateExitAlya2pos()));
        myProcess->start(programPath.absoluteFilePath(), arguments);
        myProcess->write("\n");
    }
    else {
        ui->runButton->setEnabled(true);
        ui->openResultButton->setEnabled(true);
        ui->updateResultButton->setEnabled(true);
        ui->outputRunEdit->append("Finished with errors");
    }
}

void alyaWidget::updateExitAlya2pos()
{
    ui->runButton->setEnabled(true);
    ui->openResultButton->setEnabled(true);
    ui->updateResultButton->setEnabled(true);
    QString dataFile = alyaWidget::tmpDir.absolutePath() + "/problem/" + ui->problemNameEdit->text() + ".ensi.case";
    QFile ensiFile(dataFile);
    if (ensiFile.exists()) {    
        ui->outputRunEdit->append("Finished successfully");
    }
    else {
        ui->outputRunEdit->append("Finished but data not created");
    }
}

void alyaWidget::on_openResultButton_clicked()
{
    QString dataFile = alyaWidget::tmpDir.absolutePath() + "/problem/" + ui->problemNameEdit->text() + ".ensi.case";
    QFile ensiFile(dataFile);
    if (ensiFile.exists()) {
        QStringList    dataToLoad;
        dataToLoad << dataFile;
        pqLoadDataReaction::loadData(dataToLoad);
        emitSetViewer();
    }
    else {
        QMessageBox::warning(this,"Warning","No simulation results found");
    }

}

void alyaWidget::emitSetViewer() {
    emit setViewer();
}

void alyaWidget::on_updateResultButton_clicked()
{
    //save state
    QString dataFile = alyaWidget::tmpDir.absolutePath() + "/problem/" + ui->problemNameEdit->text() + ".pvsm";
    pqSaveStateReaction::saveState(dataFile);
    //delete all
    //pqApplicationCore.instance()->
    pqDeleteReaction::deleteAll();
    //open state
    pqLoadStateReaction::loadState(dataFile);
    //delete state
    emitSetViewer();
}

void alyaWidget::on_modulesServicesButton_2_clicked()
{
    this->setCurrentIndex(1);
}
