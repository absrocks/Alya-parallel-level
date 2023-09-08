#include "alyawidget.h"
#include "ui_alyawidget.h"
#include <QFileDialog>
#include <QTextStream>
#include <sstream>
#include <string>
#include <iostream>
#include <QMessageBox>
#include "module.h"

alyaWidget::alyaWidget(QWidget *parent) :
    QStackedWidget(parent),
    ui(new Ui::alyaWidget)
{
    ui->setupUi(this);

    //modules list
    QStringList modules;
    modules << "alefor" << "chemic" << "codire" << "exmedi" << "fakemo" << "gotita" << "helmoz" << "immbou" << "latbol" << "levels" << "magnet" << "nasedg" << "nastal" << "nastin" << "partis" << "quanty" << "radiat" << "solidz" << "temper" << "turbul" << "wavequ";
    ui->alyaModules->addItems(modules);

    //services list
    QStringList services;
    services << "blas" << "cgns24" << "dodeme" << "gidpos" << "handfp" << "hdfpos" << "parall" << "solmum" << "solpls";
    ui->alyaServices->addItems(services);

    this->setCurrentIndex(0);

    ui->alyaProblemTabs->setCurrentIndex(2);
}

alyaWidget::~alyaWidget()
{
    delete ui;
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

    if (!valid) {
        QMessageBox::warning(this,"Warning",message);
    }

    return valid;

}

void alyaWidget::on_continueButton_clicked()
{
    if (checkBasicData()) {
        //remove all the tabs
        for (int i = ui->alyaProblemTabs->count() - 1; i >= 0; i--) {
            ui->alyaProblemTabs->removeTab(i);
        }
        //activate each module tab selected
        for (int i = 0; i < ui->problemModules->count(); i++) {
            if (ui->problemModules->item(i)->text() == "nastin") {
                ui->alyaProblemTabs->addTab(ui->nastin, "nastin");
            }
        }
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
    this->setCurrentIndex(2);
}

void alyaWidget::on_pushButton_2_clicked()
{
    this->setCurrentIndex(0);
}

void alyaWidget::on_openMesh_clicked()
{
    QString fileName = QFileDialog::getOpenFileName(this, tr("Open File"),
                                                     "",
                                                     tr("Files (*.*)"));
}

void alyaWidget::on_saveButton_clicked()
{
    QDir directory;
    path = QFileDialog::getExistingDirectory (this, tr("Directory"), directory.path());
    if ( path.isNull() == false )
    {
        for (int i = 0; i < ui->alyaProblemTabs->count(); i++) {
            QWidget* pWidget= ui->alyaProblemTabs->widget(i);
            std::string fileData = dynamic_cast<Module*>(pWidget)->saveFile();
            std::string fileName = dynamic_cast<Module*>(pWidget)->fileName();
            std::string filePath = path.toStdString() + "/" + ui->problemNameEdit->text().toStdString() + fileName;
            QFile file(QString::fromStdString(filePath));
            if (file.open(QFile::WriteOnly | QFile::Truncate)) {
                QTextStream out(&file);
                out << fileData.c_str();
                file.close();
            }
        }
        QMessageBox::information(this,"Save problem", "The document has been saved");
    }
}

void alyaWidget::newGUI() {
    //remove all the tabs
    for (int i = ui->alyaProblemTabs->count() - 1; i >= 0; i--) {
        ui->alyaProblemTabs->removeTab(i);
    }
    //init modules
    ui->nastin->deleteLater();
    ui->nastin = new Nastin();

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

    this->setCurrentIndex(1);
}

void alyaWidget::saveGUI() {
    QDir directory;
    QString filePath = QFileDialog::getSaveFileName(this, tr("Save Alya problem"), directory.homePath(),
                                                    "Alya problem (*.alp)");

    //save nastin
    for (int i = 0; i < ui->alyaProblemTabs->count(); i++) {
        QWidget* pWidget= ui->alyaProblemTabs->widget(i);
        dynamic_cast<Module*>(pWidget)->saveGUIValues(filePath);
    }

    //save modulesandservices
    QSettings settings(filePath, QSettings::NativeFormat);
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
    settings.endGroup();
}

void alyaWidget::loadGUI() {
    QDir directory;
    QString filePath = QFileDialog::getOpenFileName(this, tr("Open Alya problem"), directory.homePath(),
                                                    "Alya problem (*.alp)");

    //load modulesservices
    QSettings settings(filePath, QSettings::NativeFormat);
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
    //
    this->setCurrentIndex(1);

    on_continueButton_clicked();

    for (int i = 0; i < ui->alyaProblemTabs->count(); i++) {
        QWidget* pWidget= ui->alyaProblemTabs->widget(i);
        dynamic_cast<Module*>(pWidget)->loadGUIValues(filePath);
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
