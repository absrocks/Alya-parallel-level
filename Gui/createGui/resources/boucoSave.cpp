//----------------------------------------------------------------------
//-------------------------Boundary_conditions--------------------------
//----------------------------------------------------------------------
file << "BOUNDARY_CONDITIONS" << std::endl;
file << "PARAMETERS" << std::endl;
if (ui->INITIAL_CONDITION300CheckBox->isChecked()) {
    file << "  INITIAL_CONDITION:    ";
              if (ui->ie3000Combo->currentText() != "EMPTY") {
          file << "  " << ui->ie3000Combo->currentText().toStdString();
      }

        if (!ui->ieVALUES3001Edit->text().isEmpty()) {
          file << " ," << "VALUES= " << ui->ieVALUES3001Edit->text().toStdString();
        }

        if (!ui->ie3002Edit->text().isEmpty()) {
          file << " ," << ui->ie3002Edit->text().toStdString();
        }

        if (!ui->ie3003Edit->text().isEmpty()) {
          file << " ," << ui->ie3003Edit->text().toStdString();
        }


    file << std::endl;
}
if (ui->FIX_PRESSURE301CheckBox->isChecked()) {
    file << "  FIX_PRESSURE:    ";
              if (ui->ie3010Combo->currentText() != "EMPTY") {
          file << "  " << ui->ie3010Combo->currentText().toStdString();
      }

        if (!ui->ieON_NODE3011Edit->text().isEmpty()) {
          file << " ," << "ON_NODE= " << ui->ieON_NODE3011Edit->text().toStdString();
        }

        if (!ui->ieVALUE3012Edit->text().isEmpty()) {
          file << " ," << "VALUE= " << ui->ieVALUE3012Edit->text().toStdString();
        }


    file << std::endl;
}
if (ui->VARIATION302CheckBox->isChecked()) {
    file << "  VARIATION:    ";
              if (ui->ie3020Combo->currentText() != "EMPTY") {
          file << "  " << ui->ie3020Combo->currentText().toStdString();
      }


    file << std::endl;
}
    file << "END_PARAMETERS" << std::endl;
    if (ui->bcCNPTable->rowCount() > 0) {
        file << "CODES, NODES, PRESSURE" << std::endl;
        for (int i = 0; i < ui->bcCNPTable->rowCount(); i++) {
            file << "    " << ui->bcCNPTable->item(i,0)->text().toStdString() << " ";
            file << ui->bcCNPTable->item(i,1)->text().toStdString() << std::endl;
        }
        file << "END_CODES" << std::endl;
    }
    if (ui->bcCNTable->rowCount() > 0) {
        file << "CODES, NODES" << std::endl;
        for (int i = 0; i < ui->bcCNTable->rowCount(); i++) {
            QString code = ui->bcCNTable->item(i,0)->text();
            //f1f2f3
            QWidget *f1f2f3 = ui->bcCNTable->cellWidget(i,1);
            QWidget *velocity = ui->bcCNTable->cellWidget(i,2);
            QString timeFunc = ui->bcCNTable->item(i,3)->text();
            QString velBasis = ui->bcCNTable->item(i,4)->text();
            file << "    " << code.toStdString() << " ";

            //f1f2f3
            QComboBox *fCombo = (QComboBox*)f1f2f3->layout()->itemAt(0)->widget();
            if (fCombo->currentText() != "-") {
                file << fCombo->currentText().toStdString();
            }
            fCombo = (QComboBox*)f1f2f3->layout()->itemAt(1)->widget();
            if (fCombo->currentText() != "-") {
                file << fCombo->currentText().toStdString();
            }
            fCombo = (QComboBox*)f1f2f3->layout()->itemAt(2)->widget();
            if (fCombo->currentText() != "-") {
                file << fCombo->currentText().toStdString();
            }
            //velocity
            QLineEdit *veloc = (QLineEdit*)velocity->layout()->itemAt(0)->widget();
            if (veloc->text() != "") {
                file << " " << veloc->text().toStdString();
            }
            veloc = (QLineEdit*)velocity->layout()->itemAt(1)->widget();
            if (veloc->text() != "") {
                file << " " << veloc->text().toStdString();
            }
            veloc = (QLineEdit*)velocity->layout()->itemAt(2)->widget();
            if (veloc->text() != "") {
                file << " " << veloc->text().toStdString();
            }
            //timeFunc only if function is defined
            fCombo = (QComboBox*)ui->bcFTable->cellWidget(i,1)->layout()->itemAt(0)->widget();
            if (fCombo->currentText() != "-") {
                file << " " << timeFunc.toStdString();
            }

            //velBasis
            if (velBasis != "") {
                file << " " << velBasis.toStdString();
            }
            file << std::endl;
        }
        file << "END_CODES" << std::endl;
    }
    if (ui->bcCBTable->rowCount() > 0) {
        file << "CODES, BOUNDARIES" << std::endl;
        for (int i = 0; i < ui->bcCBTable->rowCount(); i++) {
            QString code = ui->bcCBTable->item(i,0)->text();
            QComboBox *typeOfCondition = (QComboBox*)ui->bcCBTable->cellWidget(i,1)->layout()->itemAt(0)->widget();
            QString value = ui->bcCBTable->item(i,2)->text();

            file << "    " << code.toStdString() << " ";
            file << typeOfCondition->currentText().toStdString();
            if (value != "") {
                file << " " << value.toStdString();
            }
            file << std::endl;
        }
        file << "END_CODES" << std::endl;
    }
    if (ui->bcFTable->rowCount() > 0) {
        std::stringstream fileAux;
        bool withFunctions = false;
        fileAux << "FUNCTIONS" << std::endl;
        for (int i = 0; i < ui->bcFTable->rowCount(); i++) {
            QString code = ui->bcFTable->item(i,0)->text();
            QComboBox *type = (QComboBox*)ui->bcFTable->cellWidget(i,1)->layout()->itemAt(0)->widget();
            if (type->currentText() != "-") {
                withFunctions = true;
                fileAux << "    FUNCTION= " << code.toStdString() << ", ";
                fileAux << type->currentText().toStdString();

                //parameters
                QWidget *parameters = ui->bcFTable->cellWidget(i,2);
                QString paramList = "";
                QLineEdit *param = (QLineEdit*)parameters->layout()->itemAt(0)->widget();
                if (param->text() != "") {
                    paramList = paramList + param->text();
                }
                param = (QLineEdit*)parameters->layout()->itemAt(1)->widget();
                if (param->text() != "") {
                    paramList = paramList + " " + param->text();
                }
                param = (QLineEdit*)parameters->layout()->itemAt(2)->widget();
                if (param->text() != "") {
                    paramList = paramList + " " + param->text();
                }
                param = (QLineEdit*)parameters->layout()->itemAt(3)->widget();
                if (param->text() != "") {
                    paramList = paramList + " " + param->text();
                }
                param = (QLineEdit*)parameters->layout()->itemAt(4)->widget();
                if (param->text() != "") {
                    paramList = paramList + " " + param->text();
                }
                param = (QLineEdit*)parameters->layout()->itemAt(5)->widget();
                if (param->text() != "") {
                    paramList = paramList + " " + param->text();
                }
                if (paramList != "") {
                    fileAux << ", PARAMETERS= " << paramList.toStdString();
                }

                fileAux << std::endl << "END_FUNCTION" << std::endl;
            }
        }
        fileAux << "END_FUNCTIONS" << std::endl;
        if (withFunctions) {
            file << fileAux.str();
        }
    }

file << "END_BOUNDARY_CONDITIONS" << std::endl;
