#include "nastin.h"
#include "ui_nastin.h"
#include <sstream>
#include <string>
#include <iostream>

Nastin::Nastin(QWidget *parent) :
    QToolBox(parent),
    ui(new Ui::Nastin)
{
    ui->setupUi(this);
    QDoubleValidator *validator = new QDoubleValidator(this);
    validator->setNotation(validator->ScientificNotation);
    ui->nastinTempEdit->setValidator(validator);
    //ui->dataWidget->layout()->setSizeConstraint(QLayout::SetFixedSize);
    QGridLayout* dataLayout= new QGridLayout();
    //dataLayout->setGeometry();
    dataLayout->setSizeConstraint(QLayout::SetMinimumSize); //important!
    //ui->page_2->setLayout(dataLayout);
    //ui->page_2->layout()->setSizeConstraint(QLayout::SetFixedSize);
    //ui->widget->setSizePolicy(QSizePolicy::Minimum,QSizePolicy::Minimum);
    //ui->widget->setMinimumSize(200,200);

}

Nastin::~Nastin()
{
    delete ui;
}

void Nastin::saveGUIValues(QString filePath) {

    QSettings settings(filePath, QSettings::NativeFormat);
    settings.beginGroup("nastin");
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
                 if (tables[i]->item(j,k) != 0 && !tables[i]->item(j,k)->text().isEmpty()) {
                     values = values + tables[i]->item(j,k)->text();
                 }
                 else {
                     values = values;
                 }
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

void Nastin::loadGUIValues(QString filePath) {

    QSettings settings(filePath, QSettings::NativeFormat);
    settings.beginGroup("nastin");
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
        tables[i]->clear();
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
                    QTableWidgetItem *newItem = new QTableWidgetItem(valueRowsList[k]);
                    tables[i]->setItem(rowIdx, k, newItem);
                }
            }
        }
    }
    settings.endGroup();

}

std::string Nastin::saveFile() {
    //----------------------------------------------------------------------
    //-------------------------Physical problem part------------------------
    //----------------------------------------------------------------------
    std::stringstream file;
    file << "PHYSICAL_PROBLEM" << std::endl;
    //--------PROBLEM DEFINITION
    file << "PROBLEM_DEFINITION" << std::endl;
    file << "  TEMPORAL_DERIVATIVES:    "<< (ui->nastinTempDeriCheck->isChecked() ? "On" : "Off")   << std::endl;
    file << "  CONVECTIVE_TERM:         "<< (ui->nastinConveTermCheck->isChecked() ? "On" : "Off")   << std::endl;
    file << "  VISCOUS_TERM:            "<< ui->nastinViscoTermcombo->currentText().toStdString()   << std::endl;
    if (ui->nastinRegCombo->currentText() != "Default") {
        file << "  REGIME:                  "<< ui->nastinRegCombo->currentText().toStdString();
        //nastinPressEdit
        if (!ui->nastinTempEdit->text().isEmpty()) {
            file << ", TEMPERATURE= " << ui->nastinTempEdit->text().toStdString();
        }
        if (!ui->nastinPressEdit->text().isEmpty()) {
            file << ", PRESSURE= " << ui->nastinPressEdit->text().toStdString();
        }
        file << std::endl;
    }
    if (!ui->nastinGravNormEdit->text().isEmpty()) {
        file << "  GRAVITY:                  NORM= " << ui->nastinGravNormEdit->text().toStdString();
        file << ", GX= " << ui->nastinGravGXEdit->text().toStdString();
        file << ", GY= " << ui->nastinGravGYEdit->text().toStdString();
        file << ", GZ= " << ui->nastinGravGZEdit->text().toStdString();
        file << std::endl;
    }
    if (!ui->nastinaAxeRotNormEdit->text().isEmpty()) {
        file << "  AXES_ROTATION:                  NORM= " << ui->nastinaAxeRotNormEdit->text().toStdString();
        file << ", OX= " << ui->nastinaAxeRotOXEdit->text().toStdString();
        file << ", OY= " << ui->nastinaAxeRotOYEdit->text().toStdString();
        file << ", OZ= " << ui->nastinaAxeRotOZEdit->text().toStdString();
        file << std::endl;
    }
    if (!ui->nastinaAxeVelNormEdit->text().isEmpty()) {
        file << "  AXES_VELOCITY:                  NORM= " << ui->nastinaAxeVelNormEdit->text().toStdString();
        file << ", VX= " << ui->nastinaAxeVelVXEdit->text().toStdString();
        file << ", VY= " << ui->nastinaAxeVelVYEdit->text().toStdString();
        file << ", VZ= " << ui->nastinaAxeVelVZEdit->text().toStdString();
        file << std::endl;
    }
    if (ui->nastinRotFuncCombo->currentText() != "Default") {
        file << "  ROTATION_FUNCTION:                  "<< ui->nastinRotFuncCombo->currentText().toStdString();
        file << std::endl;
    }
    if (ui->nastinVelFuncCombo->currentText() != "Default") {
        file << "  VELOCITY_FUNCTION:                  "<< ui->nastinVelFuncCombo->currentText().toStdString();
        file << std::endl;
    }
    if (!ui->nastinCentRotEdit1->text().isEmpty()) {
        file << "  CENTER_ROTATION:                  " << ui->nastinCentRotEdit1->text().toStdString();
        file << ", " << ui->nastinCentRotEdit2->text().toStdString();
        file << ", " << ui->nastinCentRotEdit3->text().toStdString();
        file << std::endl;
    }
    if (ui->nastinTempCoupCombo->currentText() != "Default") {
        file << "  TEMPER_COUPLING:                  "<< ui->nastinTempCoupCombo->currentText().toStdString();
        if (!ui->nastinTempCoupBetaEdit->text().isEmpty()) {
            file << ", BETA= " << ui->nastinTempCoupBetaEdit->text().toStdString();
        }
        if (!ui->nastinTempCoupTREdit->text().isEmpty()) {
            file << ", TR= " << ui->nastinTempCoupTREdit->text().toStdString();
        }
        if (!ui->nastinTempCoupGEdit->text().isEmpty()) {
            file << ", G= " << ui->nastinTempCoupGEdit->text().toStdString();
        }
        file << std::endl;
    }
    if (ui->nastinTurbModelCombo1->currentText() != "Default") {
        file << "  TURBULENCE_MODEL:                  "<< ui->nastinTurbModelCombo1->currentText().toStdString();
        if (ui->nastinTurbModelCombo2->currentText() != "Default") {
            file << ", " << ui->nastinTurbModelCombo2->currentText().toStdString();
        }
        if (!ui->nastinTurbModelParamEdit->text().isEmpty()) {
            file << ", PARAMETER= " << ui->nastinTurbModelParamEdit->text().toStdString();
        }
        file << std::endl;
    }
    if (ui->nastinLeveCoupActiveCheck->isChecked()) {
        file << "  LEVELS_COUPLING:    "<< (ui->nastinLeveCoupCheck->isChecked() ? "On" : "Off")   << std::endl;
    }
    if (ui->nastinSurfTensActiveCheck->isChecked()) {
        file << "  SURFACE_TENSION:    "<< (ui->nastinSurfTensCheck->isChecked() ? "On" : "Off");
        if (!ui->nastinSurfTensCoefEdit->text().isEmpty()) {
            file << ", COEFFICIENT= " << ui->nastinSurfTensCoefEdit->text().toStdString();
        }
        file << std::endl;
    }
    file << "END_PROBLEM_DEFINITION" << std::endl;
    //------------------------PROPERTIES
    file << "PROPERTIES" << std::endl;
    if (!ui->nastinTherPressEdit->text().isEmpty()) {
        file << "  THERMODYNAMIC_PRESSURE:                  "<< ui->nastinTherPressEdit->text().toStdString();
        if (ui->nastinTherPressCombo->currentText() != "Default") {
            file << ", " << ui->nastinTherPressCombo->currentText().toStdString();
        }
        file << std::endl;
    }
    if (!ui->nastinAirDensEdit->text().isEmpty()) {
        file << "  AIR_DENSITY=                  "<< ui->nastinAirDensEdit->text().toStdString();
        file << std::endl;
    }
    if (!ui->nastinAirViscEdit->text().isEmpty()) {
        file << "  AIR_VISCOSITY=                  "<< ui->nastinAirViscEdit->text().toStdString();
        file << std::endl;
    }
    file << "END_PROPERTIES" << std::endl;
    file << "END_PHYSICAL_PROBLEM" << std::endl;

    //----------------------------------------------------------------------
    //-------------------------Numerical treatment part---------------------
    //----------------------------------------------------------------------
    file << "NUMERICAL_TREATMENT" << std::endl;
    if (ui->nastinTAUStratCombo->currentText() != "Default") {
        file << "  TAU_STRATEGY:                  "<< ui->nastinTAUStratCombo->currentText().toStdString() << std::endl;
    }
    if (ui->nastinStabilizationCombo1->currentText() != "Default") {
        file << "  STABILIZATION:                  "<< ui->nastinStabilizationCombo1->currentText().toStdString();
        if (ui->nastinStabilizationCombo1_2->currentText() != "Default") {
            file << ", " << ui->nastinStabilizationCombo1_2->currentText().toStdString();
        }
        file << std::endl;
    }
    if(ui->groupBox->isChecked()) {
        file << "  TRACKING_SUBGRID_SCALE:                  ";
        file << (ui->nastinTrSubScTimeCheck->isChecked() ? "Time, " : "");
        file << (ui->nastinTrSubScConvCheck->isChecked() ? "CONVECTION, " : "");
        if (ui->nastinTrSubScConvCheck->isChecked()) {
            file << ui->nastinTrSubScConvCombo->currentText().toStdString() << ", ";
            file << "ITERA= " << ui->nastinTrSubScIteraEdit->text().toStdString() << ", ";
            file << "TOLER= " << ui->nastinTrSubScTolerEdit->text().toStdString() << ", ";
            file << "RELAX= " << ui->nastinTrSubScRelaxEdit->text().toStdString();
                    if (ui->nastinTrSubScTimeCheck->isChecked()) {
                        file << ", ";
                    }
        }
        if (ui->nastinTrSubScTimeCheck->isChecked()) {
            file << "ORDER= " << ui->nastinTrSubScOrderCombo->currentText().toStdString();
        }
        file << std::endl;
    }
    if (ui->nastinDiriCondCombo->currentText() != "Default") {
        file << "  DIRICHLET_CONDITION:                  "<< ui->nastinDiriCondCombo->currentText().toStdString() << std::endl;
    }
    if (ui->nastinIntForCombo->currentText() != "Default") {
        file << "  INTERNAL_FORCES:                  "<< ui->nastinIntForCombo->currentText().toStdString() << std::endl;
    }
    if (ui->nastinElemLenghtCombo->currentText() != "Default") {
        file << "  ELEMENT_LENGTH:                  "<< ui->nastinElemLenghtCombo->currentText().toStdString() << std::endl;
    }
    if (!ui->nastinPenalEdit->text().isEmpty()) {
        file << "  PRE_PENALIZATION:                  VALUE="<< ui->nastinPenalEdit->text().toStdString();
        file << std::endl;
    }
    if (ui->nastinShockCaptCombo->currentText() != "Default") {
        file << "  SHOCK_CAPTURING:                  "<< ui->nastinShockCaptCombo->currentText().toStdString();
        if (!ui->nastinShockCaptValueEdit->text().isEmpty()) {
            file << ", VALUE= " << ui->nastinTurbModelParamEdit->text().toStdString();
        }
        file << std::endl;
    }
    if (ui->nastinUpdPressActiveCheck->isChecked()) {
        file << "  UPDATE_PRESSURE:    "<< (ui->nastinUpdPressCheck->isChecked() ? "On" : "Off")   << std::endl;
    }
    if (ui->nastinTimeIntegCombo->currentText() != "Default") {
        file << "  TIME_INTEGRATION:                  "<< ui->nastinTimeIntegCombo->currentText().toStdString();
        if (!ui->nastinOrderEdit->text().isEmpty()) {
            file << ", ORDER= " << ui->nastinOrderEdit->text().toStdString();
        }
        if (!ui->nastinEulerEdit->text().isEmpty()) {
            file << ", EULER= " << ui->nastinEulerEdit->text().toStdString();
        }
        file << std::endl;
    }
    if (!ui->nastinSafFactorEdit->text().isEmpty()) {
        file << "  SAFETY_FACTOR=                  "<< ui->nastinSafFactorEdit->text().toStdString();
        file << std::endl;
    }
    if (!ui->nastinBemolEdit->text().isEmpty()) {
        file << "  BEMOL=                  "<< ui->nastinBemolEdit->text().toStdString();
        file << std::endl;
    }
    if (!ui->nastinSteStaTolEdit->text().isEmpty()) {
        file << "  STEADY_STATE_TOLERANCE=                  "<< ui->nastinSteStaTolEdit->text().toStdString();
        file << std::endl;
    }
    if (!ui->nastinMaxIntEdit->text().isEmpty()) {
        file << "  MAXIMUM_ITERATIONS=                  "<< ui->nastinMaxIntEdit->text().toStdString();
        if (!ui->nastinStartAtEdit->text().isEmpty()) {
            file << ", START_AT= " << ui->nastinStartAtEdit->text().toStdString();
        }
        file << std::endl;
    }
    if (!ui->nastinConvTolerEdit->text().isEmpty()) {
        file << "  CONVERGENCE_TOLERANCE=                  "<< ui->nastinConvTolerEdit->text().toStdString();
        file << std::endl;
    }
    if (ui->nastinLinearitzationCombo->currentText() != "Default") {
        file << "  LINEARIZATION:                  "<< ui->nastinLinearitzationCombo->currentText().toStdString() << std::endl;
        if (!ui->nastinPicardEdit->text().isEmpty()) {
            file << ", PICARD= " << ui->nastinPicardEdit->text().toStdString();
        }
        file << std::endl;
    }

    if (ui->momentuSolverBox->isChecked()) {
        file << "  MOMENTUM" << std::endl;
        file << "    ALGEBRAIC_SOLVER" << std::endl;
        file << "      SOLVER:                  " << ui->momSolverCombo->currentText().toStdString();
        if (!ui->momSolverKrylovEdit->text().isEmpty()) {
            file << " , KRYLOV=" << ui->momSolverKrylovEdit->text().toStdString();
        }
        if (ui->momSolverCoarseCombo->currentText() != "Default") {
            file << " , COARSE=" << ui->momSolverCoarseCombo->currentText().toStdString();
        }
        file << std::endl;
        if (!ui->momSolverIteraEdit->text().isEmpty()) {
            file << "      CONVERGENCE:                  ITERA=" << ui->momSolverIteraEdit->text().toStdString();
            if (!ui->momSolverTolerEdit->text().isEmpty()) {
                file << " , TOLER=" << ui->momSolverTolerEdit->text().toStdString();
            }
            if (!ui->momSolverRatioEdit->text().isEmpty()) {
                file << " , RATIO=" << ui->momSolverRatioEdit->text().toStdString();
            }
            if (!ui->momSolverAdapCheck->isChecked()) {
                file << " , ADAPTATIVE";
            }
            file << std::endl;
        }
        if (ui->momSolverOutputCombo->currentText() != "Default") {
            file << "  OUTPUT:                  "<< ui->momSolverOutputCombo->currentText().toStdString();
            file << std::endl;
        }
        if (ui->momSolverPrecondCombo->currentText() != "Default") {
            file << "      PRECONDITIONER:                  "<< ui->momSolverPrecondCombo->currentText().toStdString();
            file << std::endl;
        }

        file << "    END_ALGEBRAIC_SOLVER" << std::endl;
        file << "  END_MOMENTUM" << std::endl;
    }

    if (ui->continuitySolverBox->isChecked()) {
        file << "  CONTINUITY" << std::endl;
        file << "    ALGEBRAIC_SOLVER" << std::endl;
        file << "      SOLVER:                  " << ui->conSolverCombo->currentText().toStdString();
        if (!ui->conSolverKrylovEdit->text().isEmpty()) {
            file << " , KRYLOV=" << ui->conSolverKrylovEdit->text().toStdString();
        }
        if (ui->conSolverCoarseCombo->currentText() != "Default") {
            file << " , COARSE=" << ui->conSolverCoarseCombo->currentText().toStdString();
        }
        file << std::endl;
        if (!ui->conSolverIteraEdit->text().isEmpty()) {
            file << "      CONVERGENCE:                  ITERA=" << ui->conSolverIteraEdit->text().toStdString();
            if (!ui->conSolverTolerEdit->text().isEmpty()) {
                file << " , TOLER=" << ui->conSolverTolerEdit->text().toStdString();
            }
            if (!ui->conSolverRatioEdit->text().isEmpty()) {
                file << " , RATIO=" << ui->conSolverRatioEdit->text().toStdString();
            }
            if (!ui->conSolverAdapCheck->isChecked()) {
                file << " , ADAPTATIVE";
            }
        }
        if (ui->conSolverOutputCombo->currentText() != "Default") {
            file << "  OUTPUT:                  "<< ui->conSolverOutputCombo->currentText().toStdString();
            file << std::endl;
        }
        if (ui->conSolverPrecondCombo->currentText() != "Default") {
            file << "  PRECONDITIONER:                  "<< ui->conSolverPrecondCombo->currentText().toStdString();
            file << std::endl;
        }

        file << "    END_ALGEBRAIC_SOLVER" << std::endl;
        file << "  END_CONTINUITY" << std::endl;
    }

    if (ui->hydrostaticSolverBox->isChecked()) {
        file << "  HYDROSTATIC_STATE" << std::endl;
        file << "    ALGEBRAIC_SOLVER" << std::endl;
        file << "      SOLVER:                  " << ui->hydSolverCombo->currentText().toStdString();
        if (!ui->hydSolverKrylovEdit->text().isEmpty()) {
            file << " , KRYLOV=" << ui->hydSolverKrylovEdit->text().toStdString();
        }
        if (ui->hydSolverCoarseCombo->currentText() != "Default") {
            file << " , COARSE=" << ui->hydSolverCoarseCombo->currentText().toStdString();
        }
        file << std::endl;
        if (!ui->hydSolverIteraEdit->text().isEmpty()) {
            file << "      CONVERGENCE:                  ITERA=" << ui->hydSolverIteraEdit->text().toStdString();
            if (!ui->hydSolverTolerEdit->text().isEmpty()) {
                file << " , TOLER=" << ui->hydSolverTolerEdit->text().toStdString();
            }
            if (!ui->hydSolverRatioEdit->text().isEmpty()) {
                file << " , RATIO=" << ui->hydSolverRatioEdit->text().toStdString();
            }
            if (!ui->hydSolverAdapCheck->isChecked()) {
                file << " , ADAPTATIVE";
            }
        }
        if (ui->hydSolverOutputCombo->currentText() != "Default") {
            file << "  OUTPUT:                  "<< ui->hydSolverOutputCombo->currentText().toStdString();
            file << std::endl;
        }
        if (ui->hydSolverPrecondCombo->currentText() != "Default") {
            file << "  PRECONDITIONER:                  "<< ui->hydSolverPrecondCombo->currentText().toStdString();
            file << std::endl;
        }

        file << "    END_ALGEBRAIC_SOLVER" << std::endl;
        file << "  END_HYDROSTATIC_STATE" << std::endl;
    }

    if (ui->schurActiveCheck->isChecked()) {
        file << "  ALGORITHM: SCHUR_COMPLEMENT" << std::endl;
        if (ui->nastinSolverCombo1->currentText() != "Default") {
            file << "    SOLVER:                  "<< ui->nastinSolverCombo1->currentText().toStdString();
            if (ui->nastinSolverCombo2->currentText() != "Default") {
                file << " ," << ui->nastinSolverCombo2->currentText().toStdString();
            }
            file << std::endl;
        }
        if (ui->nastinPrecondCombo->currentText() != "Default") {
            file << "    PRECONDITIONER:                  "<< ui->nastinPrecondCombo->currentText().toStdString();
            file << std::endl;
        }
        if (ui->nastinTauStratCombo->currentText() != "Default") {
            file << "    TAU_STRATEGY:                  "<< ui->nastinTauStratCombo->currentText().toStdString();
            file << std::endl;
        }
        if (ui->nastinElemLengCombo->currentText() != "Default") {
            file << "    ELEMENT_LENGTH:                  "<< ui->nastinElemLengCombo->currentText().toStdString();
            file << std::endl;
        }
        if (ui->nastinCorrectionCombo->currentText() != "Default") {
            file << "    CORRECTION:                  "<< ui->nastinCorrectionCombo->currentText().toStdString();
            file << std::endl;
        }
        file << "  END_ALGORITHM" << std::endl;
    }
    file << "END_NUMERICAL_TREATMENT" << std::endl;


    //----------------------------------------------------------------------
    //-------------------------Output & postprocess-------------------------
    //----------------------------------------------------------------------
    file << "OUTPUT_&_POST_PROCESS" << std::endl;
    ui->postProcessTable->rowCount();
    for (int i = 0; i < ui->postProcessTable->rowCount(); i++) {
        file << "  POSTPROCESS " << ui->postProcessTable->item(i,0)->text().toStdString();
        if (ui->postProcessTable->item(i,1) != 0 && !ui->postProcessTable->item(i,1)->text().isEmpty()) {
            file << ", STEPS=" << ui->postProcessTable->item(i,1)->text().toStdString();
        }
        if (ui->postProcessTable->item(i,2) != 0 && !ui->postProcessTable->item(i,2)->text().isEmpty()) {
            file << " , FILTER=" << ui->postProcessTable->item(i,2)->text().toStdString();
        }
        file << std::endl;
    }

    file << "  ELEMENT_SET" << std::endl;
    for (int i = 0; i < ui->elementSetList->count(); i++) {
        file << "    " << ui->elementSetList->item(i)->text().toStdString() << std::endl;
    }
    file << "  END_ELEMENT_SET" << std::endl;
    file << "  BOUNDARY_SET" << std::endl;
    for (int i = 0; i < ui->boundarySetList->count(); i++) {
        file << "    " << ui->boundarySetList->item(i)->text().toStdString() << std::endl;
    }
    file << "  END_BOUNDARY_SET" << std::endl;
    file << "  NODE_SET" << std::endl;
    for (int i = 0; i < ui->nodeSetList->count(); i++) {
        file << "    " << ui->nodeSetList->item(i)->text().toStdString() << std::endl;
    }
    file << "  END_NODE_SET" << std::endl;
    file << "  WITNESS_POINTS" << std::endl;
    for (int i = 0; i < ui->witnessPointsList->count(); i++) {
        file << "    " << ui->witnessPointsList->item(i)->text().toStdString() << std::endl;
    }
    file << "  END_WITNESS_POINTS" << std::endl;
    if (ui->outputErrSolCombo->currentText() != "Default") {
        file << "  OUTPUT ERROR, SOLUTION ="<< ui->outputErrSolCombo->currentText().toStdString();
        file << std::endl;
    }
    file << "END_OUTPUT_&_POST_PROCESS" << std::endl;

    return file.str();
}

void Nastin::on_nastinRegCombo_currentIndexChanged(const QString &arg1)
{
    if (arg1 == "LOW_MACH") {
        ui->nastinTempEdit->setEnabled(true);
        ui->nastinPressEdit->setEnabled(true);
    }
    else {
        ui->nastinTempEdit->setEnabled(false);
        ui->nastinPressEdit->setEnabled(false);
        ui->nastinTempEdit->setText("");
        ui->nastinPressEdit->setText("");
    }
}

void Nastin::on_nastinTurbModelCombo1_currentIndexChanged(const QString &arg1)
{
    if (arg1 == "LES_MODEL") {
        ui->nastinTurbModelCombo2->setEnabled(true);
        ui->nastinTurbModelParamEdit->setEnabled(true);
    }
    else {
        ui->nastinTurbModelCombo2->setEnabled(false);
        ui->nastinTurbModelParamEdit->setEnabled(false);
        ui->nastinTurbModelCombo2->setCurrentIndex(0);
        ui->nastinTurbModelParamEdit->setText("");
    }

}

void Nastin::on_addPostProcessButton_clicked()
{
    if (ui->postProcessVarslist->currentItem()!=0 && ui->postProcessVarslist->currentItem()->isSelected()) {
        //add the item to the table
        int rowIdx = ui->postProcessTable->rowCount();
        ui->postProcessTable->insertRow(rowIdx);
        QTableWidgetItem *newItem = new QTableWidgetItem(ui->postProcessVarslist->currentItem()->text());
        ui->postProcessTable->setItem(rowIdx, 0, newItem);

        //remove the item from the list
        ui->postProcessVarslist->takeItem(ui->postProcessVarslist->row(ui->postProcessVarslist->currentItem()));

    }
}

void Nastin::on_delPostProcessButton_clicked()
{
    if (ui->postProcessTable->currentItem() != 0 && ui->postProcessTable->currentItem()->isSelected()) {
        //add the item to the list
        int rowIdx = ui->postProcessVarslist->count();
        ui->postProcessVarslist->insertItem(rowIdx,ui->postProcessTable->item(ui->postProcessTable->currentRow(),0)->text());

        //remove the item from the list
        ui->postProcessTable->removeRow(ui->postProcessTable->currentRow());
        //ui->postProcessVarslist->takeItem(ui->postProcessVarslist->row(ui->postProcessVarslist->currentItem()));
    }

}

void Nastin::on_addElementSetButton_clicked()
{
    if (ui->elementSetVarsList->currentItem()!=0 && ui->elementSetVarsList->currentItem()->isSelected()) {
        //add the item to the table
        int rowIdx = ui->elementSetList->count();
        ui->elementSetList->insertItem(rowIdx,ui->elementSetVarsList->currentItem()->text());

        //remove the item from the list
        ui->elementSetVarsList->takeItem(ui->elementSetVarsList->row(ui->elementSetVarsList->currentItem()));

    }
}

void Nastin::on_delElementSetButton_clicked()
{
    if (ui->elementSetList->currentItem()!=0 && ui->elementSetList->currentItem()->isSelected()) {
        //add the item to the table
        int rowIdx = ui->elementSetVarsList->count();
        ui->elementSetVarsList->insertItem(rowIdx,ui->elementSetList->currentItem()->text());

        //remove the item from the list
        ui->elementSetList->takeItem(ui->elementSetList->row(ui->elementSetList->currentItem()));
    }
}

void Nastin::on_addNodeSetButton_clicked()
{
    if (ui->nodeSetVarsList->currentItem()!=0 && ui->nodeSetVarsList->currentItem()->isSelected()) {
        //add the item to the table
        int rowIdx = ui->nodeSetList->count();
        ui->nodeSetList->insertItem(rowIdx,ui->nodeSetVarsList->currentItem()->text());

        //remove the item from the list
        ui->nodeSetVarsList->takeItem(ui->nodeSetVarsList->row(ui->nodeSetVarsList->currentItem()));

    }

}

void Nastin::on_delNodeSetButton_clicked()
{
    if (ui->nodeSetList->currentItem()!=0 && ui->nodeSetList->currentItem()->isSelected()) {
        //add the item to the table
        int rowIdx = ui->nodeSetVarsList->count();
        ui->nodeSetVarsList->insertItem(rowIdx,ui->nodeSetList->currentItem()->text());

        //remove the item from the list
        ui->nodeSetList->takeItem(ui->nodeSetList->row(ui->nodeSetList->currentItem()));
    }
}

void Nastin::on_addBoundarySetButton_clicked()
{
    if (ui->boundarySetVarsList->currentItem()!=0 && ui->boundarySetVarsList->currentItem()->isSelected()) {
        //add the item to the table
        int rowIdx = ui->boundarySetList->count();
        ui->boundarySetList->insertItem(rowIdx,ui->boundarySetVarsList->currentItem()->text());

        //remove the item from the list
        ui->boundarySetVarsList->takeItem(ui->boundarySetVarsList->row(ui->boundarySetVarsList->currentItem()));

    }

}

void Nastin::on_delBoundarySetButton_clicked()
{
    if (ui->boundarySetList->currentItem()!=0 && ui->boundarySetList->currentItem()->isSelected()) {
        //add the item to the table
        int rowIdx = ui->boundarySetVarsList->count();
        ui->boundarySetVarsList->insertItem(rowIdx,ui->boundarySetList->currentItem()->text());

        //remove the item from the list
        ui->boundarySetList->takeItem(ui->boundarySetList->row(ui->boundarySetList->currentItem()));
    }
}

void Nastin::on_addWitnessPointsButton_clicked()
{
    if (ui->witnessPointsVarslist->currentItem()!=0 && ui->witnessPointsVarslist->currentItem()->isSelected()) {
        //add the item to the table
        int rowIdx = ui->witnessPointsList->count();
        ui->witnessPointsList->insertItem(rowIdx,ui->witnessPointsVarslist->currentItem()->text());
        //remove the item from the list
        ui->witnessPointsVarslist->takeItem(ui->witnessPointsVarslist->row(ui->witnessPointsVarslist->currentItem()));
    }
}

void Nastin::on_delWitnessPointsButton_clicked()
{
    if (ui->witnessPointsList->currentItem()!=0 && ui->witnessPointsList->currentItem()->isSelected()) {
        //add the item to the table
        int rowIdx = ui->witnessPointsVarslist->count();
        ui->witnessPointsVarslist->insertItem(rowIdx,ui->witnessPointsList->currentItem()->text());

        //remove the item from the list
        ui->witnessPointsList->takeItem(ui->witnessPointsList->row(ui->witnessPointsList->currentItem()));
    }
}

void Nastin::on_nastinTempCoupCombo_currentIndexChanged(const QString &arg1)
{
    if (arg1 == "BOUSSINESQ") {
        ui->nastinTempCoupBetaEdit->setEnabled(true);
        ui->nastinTempCoupGEdit->setEnabled(true);
        ui->nastinTempCoupTREdit->setEnabled(true);
    }
    else {
        ui->nastinTempCoupBetaEdit->setEnabled(false);
        ui->nastinTempCoupGEdit->setEnabled(false);
        ui->nastinTempCoupTREdit->setEnabled(false);
        ui->nastinTempCoupBetaEdit->setText("");
        ui->nastinTempCoupGEdit->setText("");
        ui->nastinTempCoupTREdit->setText("");
    }
}

void Nastin::on_nastinLeveCoupActiveCheck_clicked(bool checked)
{
    if (checked) {
        ui->nastinLeveCoupCheck->setEnabled(true);
    }
    else {
        ui->nastinLeveCoupCheck->setEnabled(false);
    }
}

void Nastin::on_nastinSurfTensActiveCheck_clicked(bool checked)
{
    if (checked) {
        ui->nastinSurfTensCheck->setEnabled(true);
        ui->nastinSurfTensCoefEdit->setEnabled(true);
    }
    else {
        ui->nastinSurfTensCheck->setEnabled(false);
        ui->nastinSurfTensCoefEdit->setEnabled(false);
        ui->nastinSurfTensCheck->setChecked(false);
        ui->nastinSurfTensCoefEdit->setText("");
    }
}

void Nastin::on_nastinStabilizationCombo1_currentIndexChanged(const QString &arg1)
{
    if (arg1 == "SPLIT_OSS") {
        ui->nastinStabilizationCombo1_2->setEnabled(true);
    }
    else {
        ui->nastinStabilizationCombo1_2->setEnabled(false);
        ui->nastinStabilizationCombo1_2->setCurrentIndex(0);
    }
}

void Nastin::on_nastinTrSubScTimeCheck_clicked(bool checked)
{
    if (checked) {
        ui->nastinTrSubScOrderCombo->setEnabled(true);
    }
    else {
        ui->nastinTrSubScOrderCombo->setEnabled(false);
    }
}

void Nastin::on_nastinTrSubScConvCheck_clicked(bool checked)
{
    if (checked) {
        ui->nastinTrSubScConvCombo->setEnabled(true);
        ui->nastinTrSubScIteraEdit->setEnabled(true);
        ui->nastinTrSubScRelaxEdit->setEnabled(true);
        ui->nastinTrSubScTolerEdit->setEnabled(true);
    }
    else {
        ui->nastinTrSubScConvCombo->setEnabled(false);
        ui->nastinTrSubScIteraEdit->setEnabled(false);
        ui->nastinTrSubScRelaxEdit->setEnabled(false);
        ui->nastinTrSubScTolerEdit->setEnabled(false);
    }
}

void Nastin::on_nastinLinearitzationCombo_currentIndexChanged(const QString &arg1)
{
    if (arg1 == "NEWTON") {
        ui->nastinPicardEdit->setEnabled(true);
    }
    else {
        ui->nastinPicardEdit->setEnabled(false);
        ui->nastinPicardEdit->setText("");
    }
}
