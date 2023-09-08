#include "<<moduleName>>.h"
#include "ui_<<moduleName>>.h"
#include <sstream>
#include <string>
#include <iostream>
#include <QMessageBox>
#include <QListWidget>
#include <QTableWidget>
#include <QComboBox>
#include <QCheckBox>
#include <QLineEdit>
#include <QXmlStreamWriter>
#include <QXmlStreamReader>
#include <QFile>

<<className>>::<<className>>(QWidget *parent) :
    QToolBox(parent),
    ui(new Ui::<<className>>)
{
    ui->setupUi(this);
    realValidator = new QDoubleValidator();
    intValidator = new QIntValidator();
    QRegExp rx("^[-+]?[0-9]*\\.?[0-9]+([eE][-+]?[0-9]+)?(\\s*,\\s*[-+]?[0-9]*\\.?[0-9]+([eE][-+]?[0-9]+)?)*$");
    realseqValidator = new QRegExpValidator(rx, this);
    setValidators();
    <<constructorCalls>>
}

<<className>>::~<<className>>()
{
    delete ui;
}

void <<className>>::setValidators() {
    QDoubleValidator *dv = new QDoubleValidator();
    QIntValidator *iv = new QIntValidator();
    <<setValidators>>
}

<<loadSaveGui>>

void <<className>>::saveGUIValues(QString filePath) {
    <<saveGui>>
    QLineEdit* tmpLineEdit;
    QSettings settings(filePath, QSettings::IniFormat);
    settings.beginGroup("<<moduleName>>");
    QList<QComboBox *> combos = this->findChildren<QComboBox *>();
    for (int i=0; i < combos.count(); i++) {
        settings.setValue(combos[i]->objectName(), combos[i]->currentIndex());
    }
    QList<QLineEdit *> lines = this->findChildren<QLineEdit *>();
    for (int i=0; i < lines.count(); i++) {
        settings.setValue(lines[i]->objectName(), lines[i]->text());
    }
    QList<QGroupBox *> groups = this->findChildren<QGroupBox *>();
    for (int i=0; i < groups.count(); i++) {
        if (groups[i]->isCheckable()) {
             settings.setValue(groups[i]->objectName(), groups[i]->isChecked());
        }
    }
    QList<QCheckBox *> checks = this->findChildren<QCheckBox *>();
    for (int i=0; i < checks.count(); i++) {
        settings.setValue(checks[i]->objectName(), checks[i]->isChecked());
    }
    QList<QListWidget *> lists = this->findChildren<QListWidget *>();
    for (int i=0; i < lists.count(); i++) {
        QString values = "";
        for (int j=0; j<lists[i]->count(); j++) {
            values = values + lists[i]->item(j)->text();
            if (j < (lists[i]->count() - 1)) {
                values = values + "#";
            }
        }
        if (values != "") {
             settings.setValue(lists[i]->objectName(), values);
        }
    }
    QList<QTableWidget *> tables = this->findChildren<QTableWidget *>();
    for (int i=0; i < tables.count(); i++) {
        QString values = "";
        for (int j=0; j<tables[i]->rowCount(); j++) {
            for (int k=0; k<tables[i]->columnCount(); k++) {
                 tmpLineEdit = qobject_cast<QLineEdit*>(tables[i]->cellWidget(j,k));
                 if (!tmpLineEdit->text().isEmpty()) {
                     values = values + tmpLineEdit->text();
                 }
                 else {
                     values = values;
                 }
                 values = values + "$" + tmpLineEdit->toolTip();
                 if (k < (tables[i]->columnCount() - 1)) {
                     values = values + "@";
                 }
            }

            if (j < (tables[i]->rowCount() - 1)) {
                values = values + "#";
            }
        }
        if (values != "") {
             settings.setValue(tables[i]->objectName(), values);
        }
    }
    settings.endGroup();

}

QValidator *<<className>>::getValidator(QString validatorName) {
    if (validatorName == "REAL") {
        return realValidator;
    }
    else if (validatorName == "INT") {
        return intValidator;
    }
    else if (validatorName == "REALSEQ") {
        return realseqValidator;
    }
    else {
        return 0;
    }
}

void <<className>>::loadGUIValues(QString filePath) {

    QSettings settings(filePath, QSettings::IniFormat);
    settings.beginGroup("<<moduleName>>");
    QList<QComboBox *> combos = this->findChildren<QComboBox *>();
    for (int i=0; i < combos.count(); i++) {
        int idx = settings.value(combos[i]->objectName(),"0").toInt();
        combos[i]->setCurrentIndex(idx);
    }
    QList<QLineEdit *> lines = this->findChildren<QLineEdit *>();
    for (int i=0; i < lines.count(); i++) {
        QString lineValue = settings.value(lines[i]->objectName(),"").toString();
        lines[i]->setText(lineValue);
    }
    QList<QGroupBox *> groups = this->findChildren<QGroupBox *>();
    for (int i=0; i < groups.count(); i++) {
        if (groups[i]->isCheckable()) {
            bool checked = settings.value(groups[i]->objectName(),false).toBool();
            groups[i]->setChecked(checked);
        }
    }
    QList<QCheckBox *> checks = this->findChildren<QCheckBox *>();
    for (int i=0; i < checks.count(); i++) {
        bool checked = settings.value(checks[i]->objectName(),false).toBool();
        checks[i]->setChecked(checked);
    }
    QList<QListWidget *> lists = this->findChildren<QListWidget *>();
    for (int i=0; i < lists.count(); i++) {
        lists[i]->clear();
        QString vals = settings.value(lists[i]->objectName(),"").toString();
        if (vals != "") {
            QStringList valuesList = vals.split("#");
            for (int j = 0; j < valuesList.size(); j++) {
                int rowIdx = lists[i]->count();
                lists[i]->insertItem(rowIdx,valuesList[j]);
            }
        }
    }
    QList<QTableWidget *> tables = this->findChildren<QTableWidget *>();
    for (int i=0; i < tables.count(); i++) {
        tables[i]->clearContents();
        for (int j=tables[i]->rowCount()-1; j >= 0; --j) {
            tables[i]->removeRow(j);
        }
        QString vals = settings.value(tables[i]->objectName(),"").toString();
        if (vals != "") {
            QStringList valuesList = vals.split("#");
            for (int j = 0; j < valuesList.size(); j++) {
                QStringList valueRowsList = valuesList[j].split("@");
                int rowIdx = tables[i]->rowCount();
                tables[i]->insertRow(rowIdx);
                for (int k = 0; k < valueRowsList.size(); k++) {
                    QStringList valueTooltipList = valueRowsList[k].split("$");
                    QLineEdit *newEdit = new QLineEdit();
                    if (valueTooltipList.size() > 1) {
                        newEdit->setToolTip(valueTooltipList[1]);
                        QValidator *validator = getValidator(valueTooltipList[1]);
                        if (validator != 0) {
                            newEdit->setValidator(validator);
                        }
                    }
                    newEdit->setText(valueTooltipList[0]);
                    tables[i]->setCellWidget(rowIdx, k, newEdit);
                    QTableWidgetItem *newItem = new QTableWidgetItem();
                    tables[i]->setItem(rowIdx, k, newItem);
                }
            }
        }
    }
    settings.endGroup();
    <<loadGui>>
}

void <<className>>::loadCodes(QMap<QString, QString> boundaryCodes, QMap<QString, QString> nodeCodes) {
    <<loadCodes>>
}

std::string <<className>>::fileName() {
    return "<<fileExtension>>";
}

std::string <<className>>::saveFile() {
    std::stringstream file;
    QLineEdit* tmpLineEdit;
    <<uisaves>>
    return file.str();
}



<<privateSlots>>


void <<className>>::disableModulesServicesParams() {

    <<modulesToDisable>>
}

void <<className>>::enableModule(QString moduleName) {
    QList<QGroupBox *> groups = this->findChildren<QGroupBox *>();
    for (int i=0; i < groups.count(); i++) {
        QString name = groups[i]->objectName();
        if (name.contains("sg" + moduleName + "_module") && !groups[i]->isEnabled()) {
            groups[i]->setVisible(true);
            groups[i]->setEnabled(true);
        }
    }
}

void <<className>>::enableService(QString moduleService) {
    QList<QGroupBox *> groups = this->findChildren<QGroupBox *>();
    for (int i=0; i < groups.count(); i++) {
        QString name = groups[i]->objectName();
        if (name.contains("sg" + moduleService + "_service") && !groups[i]->isEnabled()) {
            groups[i]->setVisible(true);
            groups[i]->setEnabled(true);
        }
    }
}
