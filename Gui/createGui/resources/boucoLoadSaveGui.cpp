void <<moduleClass>>::loadBoundCondGui(QString filePath) {

    ui->bcCNPTable->clearContents();
    ui->bcCNPTable->setRowCount(0);
    ui->bcCNTable->clearContents();
    ui->bcCNTable->setRowCount(0);
    ui->bcCBTable->clearContents();
    ui->bcCBTable->setRowCount(0);
    ui->bcFTable->clearContents();
    ui->bcFTable->setRowCount(0);

    filePath.replace(".alp", "_<<moduleClass>>_boundCond.xml");
    QFile* file = new QFile(filePath);
    if (file->open(QIODevice::ReadOnly | QIODevice::Text)) {
        QXmlStreamReader xml(file);
        while(!xml.atEnd() && !xml.hasError()) {
            /* Read next element.*/
            QXmlStreamReader::TokenType token = xml.readNext();
            /* If token is just StartDocument, we'll go to next.*/
            if(token == QXmlStreamReader::StartDocument) {
                continue;
            }
            /* If token is StartElement, we'll see if we can read it.*/
            if(token == QXmlStreamReader::StartElement) {
                /* If it's named persons, we'll go to the next.*/
                if(xml.name() == "codesNodesPressureTable" ||
                   xml.name() == "codesNodesTable" ||
                   xml.name() == "codesBoundariesTable" ||
                   xml.name() == "functionsTable") {
                    continue;
                }
                //----------------CodesNodePressure part----------------------------
                if(xml.name() == "codesNodesPressure") {
                    while(!(xml.tokenType() == QXmlStreamReader::EndElement &&
                                xml.name() == "codesNodesPressure")) {
                        if(xml.tokenType() == QXmlStreamReader::StartElement) {
                            //----------------code data----------------
                            if(xml.name() == "code") {                                
                                xml.readNext();
                                if(xml.tokenType() == QXmlStreamReader::Characters) {
                                    int rowIdx = ui->bcCNPTable->rowCount();
                                    ui->bcCNPTable->insertRow(rowIdx);
                                    QTableWidgetItem *newItem = new QTableWidgetItem(xml.text().toString());
                                    newItem->setFlags(Qt::ItemIsSelectable);
                                    ui->bcCNPTable->setItem(rowIdx, 0, newItem);
                                }
                            }
                            //---------------value data----------------
                            if(xml.name() == "value") {
                                xml.readNext();
                                if(xml.tokenType() == QXmlStreamReader::Characters) {
                                    int rowIdx = ui->bcCNPTable->rowCount() - 1;
                                    QTableWidgetItem *newItem = new QTableWidgetItem(xml.text().toString());
                                    ui->bcCNPTable->setItem(rowIdx, 1, newItem);
                                }
                            }
                        }
                        /* ...and next... */
                        xml.readNext();
                    }
                }
                //----------------CodesNode part----------------------------
                if(xml.name() == "codesNodes") {
                    while(!(xml.tokenType() == QXmlStreamReader::EndElement &&
                                xml.name() == "codesNodes")) {
                        if(xml.tokenType() == QXmlStreamReader::StartElement) {
                            //----------------code data----------------
                            if(xml.name() == "code") {
                                xml.readNext();
                                if(xml.tokenType() == QXmlStreamReader::Characters) {
                                    addbcCNRow(xml.text().toString());
                                }
                            }
                            //---------------f1f2f3 data----------------
                            if(xml.name() == "f1" || xml.name() == "f2" || xml.name() == "f3") {
                                int rowIdx = ui->bcCNTable->rowCount() - 1;
                                int widIdx;
                                if(xml.name() == "f1") {
                                    widIdx = 0;
                                }
                                else if(xml.name() == "f2") {
                                    widIdx = 1;
                                }
                                else if(xml.name() == "f3") {
                                    widIdx = 2;
                                }
                                xml.readNext();
                                if(xml.tokenType() == QXmlStreamReader::Characters) {
                                    //f1f2f3
                                    QWidget *f1f2f3 = ui->bcCNTable->cellWidget(rowIdx,1);
                                    //QString velBasis = ui->bcCNTable->item(i,4)->text();
                                    QComboBox *fCombo = (QComboBox*)f1f2f3->layout()->itemAt(widIdx)->widget();
                                    for (int i = 0; i < fCombo->count(); i++) {
                                        if (fCombo->itemText(i) == xml.text().toString()) {
                                            fCombo->setCurrentIndex(i);
                                        }
                                    }

                                }
                            }
                            //---------------v1v2v3 data----------------
                            if(xml.name() == "v1" || xml.name() == "v2" || xml.name() == "v3") {
                                int rowIdx = ui->bcCNTable->rowCount() - 1;
                                int widIdx;
                                if(xml.name() == "v1") {
                                    widIdx = 0;
                                }
                                else if(xml.name() == "v2") {
                                    widIdx = 1;
                                }
                                else if(xml.name() == "v3") {
                                    widIdx = 2;
                                }
                                xml.readNext();
                                if(xml.tokenType() == QXmlStreamReader::Characters && !xml.text().isEmpty()) {
                                    QWidget *velocity = ui->bcCNTable->cellWidget(rowIdx,2);
                                    QLineEdit *fedit = (QLineEdit*)velocity->layout()->itemAt(widIdx)->widget();
                                    fedit->setText(xml.text().toString());
                                }
                            }
                            if(xml.name() == "timeFunc") {
                                int rowIdx = ui->bcCNTable->rowCount() - 1;
                                xml.readNext();
                                if(xml.tokenType() == QXmlStreamReader::Characters && !xml.text().isEmpty()) {
                                    ui->bcCNTable->item(rowIdx,3)->setText(xml.text().toString());
                                }
                            }
                            if(xml.name() == "velBasis") {
                                int rowIdx = ui->bcCNTable->rowCount() - 1;
                                xml.readNext();
                                if(xml.tokenType() == QXmlStreamReader::Characters && !xml.text().isEmpty()) {
                                    ui->bcCNTable->item(rowIdx,4)->setText(xml.text().toString());
                                }
                            }
                        }
                        /* ...and next... */
                        xml.readNext();
                    }
                }
                //----------------Codes Boundaries part----------------------------
                if(xml.name() == "codesBoundaries") {
                    while(!(xml.tokenType() == QXmlStreamReader::EndElement &&
                                xml.name() == "codesBoundaries")) {
                        if(xml.tokenType() == QXmlStreamReader::StartElement) {
                            //----------------code data----------------
                            if(xml.name() == "code") {
                                xml.readNext();
                                if(xml.tokenType() == QXmlStreamReader::Characters) {
                                    addbcCBRow(xml.text().toString());
                                }
                            }
                            //---------------typeOfCondition data----------------
                            if(xml.name() == "typeOfCondition") {
                                int rowIdx = ui->bcCBTable->rowCount() - 1;
                                xml.readNext();
                                if(xml.tokenType() == QXmlStreamReader::Characters && !xml.text().isEmpty()) {
                                    QWidget *fComboWidget = ui->bcCBTable->cellWidget(rowIdx,1);
                                    QComboBox *fCombo = (QComboBox*)fComboWidget->layout()->itemAt(0)->widget();
                                    for (int i = 0; i < fCombo->count(); i++) {
                                        if (fCombo->itemText(i) == xml.text().toString()) {
                                            fCombo->setCurrentIndex(i);
                                        }
                                    }
                                }
                            }
                            //---------------value data----------------
                            if(xml.name() == "value") {
                                int rowIdx = ui->bcCBTable->rowCount() - 1;
                                xml.readNext();
                                if(xml.tokenType() == QXmlStreamReader::Characters && !xml.text().isEmpty()) {
                                    ui->bcCBTable->item(rowIdx,2)->setText(xml.text().toString());
                                }
                            }
                        }
                        /* ...and next... */
                        xml.readNext();
                    }
                }
                //----------------Functions part----------------------------
                QString fCode;
                if(xml.name() == "functions") {
                    while(!(xml.tokenType() == QXmlStreamReader::EndElement &&
                                xml.name() == "functions")) {
                        if(xml.tokenType() == QXmlStreamReader::StartElement) {
                            //----------------code data----------------
                            if(xml.name() == "code") {
                                xml.readNext();
                                if(xml.tokenType() == QXmlStreamReader::Characters) {
                                    fCode = xml.text().toString();
                                }
                            }
                            //---------------type data----------------
                            if(xml.name() == "type" && !fCode.isEmpty()) {
                                int rowIdx = fCode.toInt();
                                xml.readNext();
                                if(xml.tokenType() == QXmlStreamReader::Characters && !xml.text().isEmpty()) {
                                    QWidget *fComboWidget = ui->bcFTable->cellWidget(rowIdx,1);
                                    QComboBox *fCombo = (QComboBox*)fComboWidget->layout()->itemAt(0)->widget();
                                    for (int i = 0; i < fCombo->count(); i++) {
                                        if (fCombo->itemText(i) == xml.text().toString()) {
                                            fCombo->setCurrentIndex(i);
                                        }
                                    }
                                    ui->bcFTable->resizeColumnToContents(1);
                                    ui->bcFTable->resizeRowToContents(rowIdx);
                                }
                            }
                            //---------------parameter data----------------
                            if((xml.name() == "p1" || xml.name() == "p2" || xml.name() == "p3" ||
                               xml.name() == "p4" || xml.name() == "p5" || xml.name() == "p6") &&
                               !fCode.isEmpty()) {
                                int rowIdx = fCode.toInt();
                                int widIdx;
                                if(xml.name() == "p1") {
                                    widIdx = 0;
                                }
                                else if(xml.name() == "p2") {
                                    widIdx = 1;
                                }
                                else if(xml.name() == "p3") {
                                    widIdx = 2;
                                }
                                else if(xml.name() == "p4") {
                                    widIdx = 3;
                                }
                                else if(xml.name() == "p5") {
                                    widIdx = 4;
                                }
                                else if(xml.name() == "p6") {
                                    widIdx = 5;
                                }
                                xml.readNext();
                                if(xml.tokenType() == QXmlStreamReader::Characters && !xml.text().isEmpty()) {
                                    QWidget *parameter = ui->bcFTable->cellWidget(rowIdx,2);
                                    QLineEdit *fedit = (QLineEdit*)parameter->layout()->itemAt(widIdx)->widget();
                                    fedit->setText(xml.text().toString());
                                }
                            }
                        }
                        /* ...and next... */
                        xml.readNext();
                    }
                }
            }
        }
    }
}

void <<moduleClass>>::saveBoundCondGui(QString filePath) {

    filePath.replace(".alp", "_<<moduleClass>>_boundCond.xml");
    QFile output(filePath);
    if (output.open(QFile::WriteOnly | QFile::Truncate)) {
            QXmlStreamWriter stream(&output);
            stream.setAutoFormatting(true);
            stream.writeStartDocument();
            stream.writeStartElement("boundaryConditions");
            stream.writeStartElement("codesNodesPressureTable");
            for (int i = 0; i < ui->bcCNPTable->rowCount(); i++) {
                stream.writeStartElement("codesNodesPressure");
                stream.writeTextElement("code", ui->bcCNPTable->item(i,0)->text());
                stream.writeTextElement("value", (ui->bcCNPTable->item(i,1) != 0) ? ui->bcCNPTable->item(i,1)->text() : "");
                //end codesNodesPressure
                stream.writeEndElement();
            }
            //end codesNodesPressureTable
            stream.writeEndElement();

            stream.writeStartElement("codesNodesTable");
            for (int i = 0; i < ui->bcCNTable->rowCount(); i++) {
                stream.writeStartElement("codesNodes");
                stream.writeTextElement("code", ui->bcCNTable->item(i,0)->text());

                QWidget *f1f2f3 = ui->bcCNTable->cellWidget(i,1);
                QWidget *velocity = ui->bcCNTable->cellWidget(i,2);
                QString timeFunc = (ui->bcCNTable->item(i,3) != 0) ? ui->bcCNTable->item(i,3)->text() : "";
                QString velBasis = (ui->bcCNTable->item(i,4) != 0) ? ui->bcCNTable->item(i,4)->text() : "";

                //f1f2f3
                QComboBox *fCombo = (QComboBox*)f1f2f3->layout()->itemAt(0)->widget();
                stream.writeTextElement("f1", fCombo->currentText());
                fCombo = (QComboBox*)f1f2f3->layout()->itemAt(1)->widget();
                stream.writeTextElement("f2", fCombo->currentText());
                fCombo = (QComboBox*)f1f2f3->layout()->itemAt(2)->widget();
                stream.writeTextElement("f3", fCombo->currentText());

                //velocity
                QLineEdit *veloc = (QLineEdit*)velocity->layout()->itemAt(0)->widget();
                stream.writeTextElement("v1", veloc->text());
                veloc = (QLineEdit*)velocity->layout()->itemAt(1)->widget();
                stream.writeTextElement("v2", veloc->text());
                veloc = (QLineEdit*)velocity->layout()->itemAt(2)->widget();
                stream.writeTextElement("v3", veloc->text());
                //timeFunc
                stream.writeTextElement("timeFunc", timeFunc);
                //velBasis
                stream.writeTextElement("velBasis", velBasis);

                //end codesNodes
                stream.writeEndElement();
            }
            //end codesNodesTable
            stream.writeEndElement();

            stream.writeStartElement("codesBoundariesTable");
            for (int i = 0; i < ui->bcCBTable->rowCount(); i++) {
                stream.writeStartElement("codesBoundaries");
                QString code = ui->bcCBTable->item(i,0)->text();
                stream.writeTextElement("code", code);

                QComboBox *typeOfCondition = (QComboBox*)ui->bcCBTable->cellWidget(i,1)->layout()->itemAt(0)->widget();
                stream.writeTextElement("typeOfCondition", typeOfCondition->currentText());

                QString value = (ui->bcCBTable->item(i,2) != 0) ? ui->bcCBTable->item(i,2)->text() : "";
                stream.writeTextElement("value", value);

                //end codesBoundaries
                stream.writeEndElement();
            }
            //end codesBoundariesTable
            stream.writeEndElement();

            stream.writeStartElement("functionsTable");
            for (int i = 0; i < ui->bcFTable->rowCount(); i++) {
                stream.writeStartElement("functions");

                QString code = ui->bcFTable->item(i,0)->text();
                stream.writeTextElement("code", code);

                QComboBox *type = (QComboBox*)ui->bcFTable->cellWidget(i,1)->layout()->itemAt(0)->widget();
                stream.writeTextElement("type", type->currentText());

                //parameters
                QWidget *parameters = ui->bcFTable->cellWidget(i,2);
                QLineEdit *param = (QLineEdit*)parameters->layout()->itemAt(0)->widget();
                stream.writeTextElement("p1", param->text());
                param = (QLineEdit*)parameters->layout()->itemAt(1)->widget();
                stream.writeTextElement("p2", param->text());
                param = (QLineEdit*)parameters->layout()->itemAt(2)->widget();
                stream.writeTextElement("p3", param->text());
                param = (QLineEdit*)parameters->layout()->itemAt(3)->widget();
                stream.writeTextElement("p4", param->text());
                param = (QLineEdit*)parameters->layout()->itemAt(4)->widget();
                stream.writeTextElement("p5", param->text());
                param = (QLineEdit*)parameters->layout()->itemAt(5)->widget();
                stream.writeTextElement("p6", param->text());
                //end functions
                stream.writeEndElement();
            }
            //end functionsTable
            stream.writeEndElement();
            //end boundaryConditions
            stream.writeEndElement();
            stream.writeEndDocument();
            output.close();
    }
}
