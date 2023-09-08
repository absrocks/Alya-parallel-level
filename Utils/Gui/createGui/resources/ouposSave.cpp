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
    if (ui->outputErrorBox->isChecked()) {
        file << "  OUTPUT ERROR, SOLUTION ="<< ui->solutionCombo->currentText().toStdString();
        file << std::endl;
    }
    file << "END_OUTPUT_&_POST_PROCESS" << std::endl;

