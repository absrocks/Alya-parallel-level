void <<moduleClass>>::insertTimeFunction(int code) {
    //new row
    QWidget *functionType = new QWidget(this);
    QComboBox *f1Combo = new QComboBox();
    f1Combo->addItem("-");
    f1Combo->addItem("PARABOLIC");
    f1Combo->addItem("PERIODIC");
    f1Combo->addItem("DISCRETE");
    QHBoxLayout *layout = new QHBoxLayout;
    layout->addWidget(f1Combo);
    functionType->setLayout(layout);

    QWidget *parameters = new QWidget(this);
    QLineEdit *lineEdit1 = new QLineEdit();
    lineEdit1->setMaximumWidth(60);
    QLineEdit *lineEdit2 = new QLineEdit();
    lineEdit2->setMaximumWidth(60);
    QLineEdit *lineEdit3 = new QLineEdit();
    lineEdit3->setMaximumWidth(60);
    QLineEdit *lineEdit4 = new QLineEdit();
    lineEdit4->setMaximumWidth(60);
    QLineEdit *lineEdit5 = new QLineEdit();
    lineEdit5->setMaximumWidth(60);
    QLineEdit *lineEdit6 = new QLineEdit();
    lineEdit6->setMaximumWidth(60);
    QHBoxLayout *layout2 = new QHBoxLayout;
    layout2->addWidget(lineEdit1);
    layout2->addWidget(lineEdit2);
    layout2->addWidget(lineEdit3);
    layout2->addWidget(lineEdit4);
    layout2->addWidget(lineEdit5);
    layout2->addWidget(lineEdit6);
    parameters->setLayout(layout2);


    ui->bcFTable->insertRow(code);
    QString codeString = QString::number(code);
    QTableWidgetItem *newItem = new QTableWidgetItem(codeString);
    newItem->setFlags(Qt::ItemIsSelectable);
    ui->bcFTable->setItem(code, 0, newItem);
    ui->bcFTable->setColumnWidth(0,40);
    ui->bcFTable->setCellWidget(code,1,functionType);
    ui->bcFTable->setColumnWidth(1,50);
    ui->bcFTable->resizeColumnToContents(1);
    ui->bcFTable->setCellWidget(code,2,parameters);
    ui->bcFTable->setColumnWidth(2,50);
    ui->bcFTable->resizeColumnToContents(2);
    ui->bcFTable->resizeRowToContents(code);

}

void <<moduleClass>>::on_addbcCNPButton_clicked()
{
    if (ui->bcCNPCodesList->currentItem()!=0 && ui->bcCNPCodesList->currentItem()->isSelected()) {
        //add the item to the table
        int rowIdx = ui->bcCNPTable->rowCount();
        ui->bcCNPTable->insertRow(rowIdx);
        QTableWidgetItem *newItem = new QTableWidgetItem(ui->bcCNPCodesList->currentItem()->text());
        newItem->setFlags(Qt::ItemIsSelectable);
        ui->bcCNPTable->setItem(rowIdx, 0, newItem);

        //remove the item from the list
        ui->bcCNPCodesList->takeItem(ui->bcCNPCodesList->row(ui->bcCNPCodesList->currentItem()));

    }
}

void <<moduleClass>>::on_delbcCNPButton_clicked()
{
    if (ui->bcCNPTable->currentItem() != 0 && ui->bcCNPTable->currentItem()->isSelected()) {
        //add the item to the list
        int rowIdx = ui->bcCNPCodesList->count();
        ui->bcCNPCodesList->insertItem(rowIdx,ui->bcCNPTable->item(ui->bcCNPTable->currentRow(),0)->text());

        //remove the item from the list
        ui->bcCNPTable->removeRow(ui->bcCNPTable->currentRow());
        //ui->postProcessVarsList->takeItem(ui->postProcessVarsList->row(ui->postProcessVarsList->currentItem()));
    }
}

void <<moduleClass>>::addbcCNRow(QString codes) {
    //f1f2f3 column widget
    QWidget *f1f2f3 = new QWidget(this);
    QComboBox *f1Combo = new QComboBox();
    f1Combo->addItem("-");
    f1Combo->addItem("0");
    f1Combo->addItem("1");
    QComboBox *f2Combo = new QComboBox();
    f2Combo->addItem("-");
    f2Combo->addItem("0");
    f2Combo->addItem("1");
    QComboBox *f3Combo = new QComboBox();
    f3Combo->addItem("-");
    f3Combo->addItem("0");
    f3Combo->addItem("1");
    QHBoxLayout *layout = new QHBoxLayout;
    layout->addWidget(f1Combo);
    layout->addWidget(f2Combo);
    layout->addWidget(f3Combo);
    f1f2f3->setLayout(layout);

    //velocity column widget
    //f1f2f3 section
    QWidget *velocity = new QWidget(this);
    QLineEdit *lineEdit1 = new QLineEdit();
    lineEdit1->setMaximumWidth(60);
    QLineEdit *lineEdit2 = new QLineEdit();
    lineEdit2->setMaximumWidth(60);
    QLineEdit *lineEdit3 = new QLineEdit();
    lineEdit3->setMaximumWidth(60);
    QHBoxLayout *layout2 = new QHBoxLayout;
    layout2->addWidget(lineEdit1);
    layout2->addWidget(lineEdit2);
    layout2->addWidget(lineEdit3);
    velocity->setLayout(layout2);

    //add the item to the table
    int rowIdx = ui->bcCNTable->rowCount();
    ui->bcCNTable->insertRow(rowIdx);
    QTableWidgetItem *newItem = new QTableWidgetItem(codes);
    newItem->setFlags(Qt::ItemIsSelectable);
    //Code column
    ui->bcCNTable->setItem(rowIdx, 0, newItem);
    ui->bcCNTable->setColumnWidth(0,40);
    //f1f2f3 columnt
    ui->bcCNTable->setCellWidget(rowIdx,1,f1f2f3);
    ui->bcCNTable->resizeColumnToContents(1);
    //velocity column
    ui->bcCNTable->setCellWidget(rowIdx,2,velocity);
    ui->bcCNTable->resizeColumnToContents(2);
    //time function column
    QTableWidgetItem *newItem2 = new QTableWidgetItem("",1001);
    QString rowString = QString::number(rowIdx);
    newItem2->setText(rowString);
    newItem2->setFlags(Qt::ItemIsSelectable);
    insertTimeFunction(rowIdx);
    ui->bcCNTable->setItem(rowIdx, 3, newItem2);
    ui->bcCNTable->setColumnWidth(3,50);
    //velocity basis column
    QTableWidgetItem *newItem3 = new QTableWidgetItem("");
    ui->bcCNTable->setItem(rowIdx, 4, newItem3);
    ui->bcCNTable->setColumnWidth(4,50);
    ui->bcCNTable->resizeRowToContents(rowIdx);
}

void <<moduleClass>>::on_addbcCNButton_clicked()
{
    if (ui->bcCNCodesList->currentItem()!=0 && ui->bcCNCodesList->currentItem()->isSelected()) {
        //get the selected items list
        QList<QListWidgetItem*> itemList = ui->bcCNCodesList->selectedItems();
        QString codes;
        for (int i = 0; i < itemList.size(); ++i) {
           codes = codes + itemList.at(i)->text();
           if (i + 1< itemList.size()) {
               codes = codes + "&";
           }
         }
        addbcCNRow(codes);

        //remove the item from the list
        for (int i = itemList.size() - 1; i >= 0; i--) {
            ui->bcCNCodesList->takeItem(ui->bcCNCodesList->row(itemList.at(i)));
        }

    }

}

void <<moduleClass>>::on_delbcCNButton_clicked()
{
    if (ui->bcCNTable->currentItem() != 0 && ui->bcCNTable->currentItem()->isSelected()) {
        //add the item to the list
        int rowIdx = ui->bcCNCodesList->count();
        QString codes = ui->bcCNTable->item(ui->bcCNTable->currentRow(),0)->text();
        QStringList codesList = codes.split("&");
        for (int i = 0; i < codesList.size(); ++i) {
                   ui->bcCNCodesList->insertItem(rowIdx,codesList.at(i));
         }

        //remove the item from the list
        int currentRow = ui->bcCNTable->currentRow();
        ui->bcCNTable->removeRow(currentRow);
        //remove the item from functions list
        ui->bcFTable->removeRow(currentRow);
        //ui->postProcessVarsList->takeItem(ui->postProcessVarsList->row(ui->postProcessVarsList->currentItem()));
        //recalculate time function id
        int rowCount = ui->bcCNTable->rowCount();
        for (int i = 0; i < rowCount; ++i) {
           QString iString = QString::number(i);
           ui->bcCNTable->item(i,3)->setText(iString);
           ui->bcFTable->item(i,0)->setText(iString);
         }
    }
}

void <<moduleClass>>::addbcCBRow(QString codes) {

    //condition type widget
    QWidget *conditionType = new QWidget(this);
    QComboBox *f1Combo = new QComboBox();
    f1Combo->addItem("2");
    f1Combo->addItem("3");
    QHBoxLayout *layout = new QHBoxLayout;
    layout->addWidget(f1Combo);
    conditionType->setLayout(layout);

    //add the item to the table
    int rowIdx = ui->bcCBTable->rowCount();
    ui->bcCBTable->insertRow(rowIdx);
    QTableWidgetItem *newItem = new QTableWidgetItem(codes);
    newItem->setFlags(Qt::ItemIsSelectable);
    //Code column
    ui->bcCBTable->setItem(rowIdx, 0, newItem);
    ui->bcCBTable->setColumnWidth(0,40);
    //conditionType column
    ui->bcCBTable->setCellWidget(rowIdx,1,conditionType);
    ui->bcCBTable->setColumnWidth(1,50);
    //value column
    QTableWidgetItem *newItem2 = new QTableWidgetItem("");
    ui->bcCBTable->setItem(rowIdx, 2, newItem2);
    ui->bcCBTable->setColumnWidth(2,50);
    ui->bcCBTable->resizeColumnToContents(1);
    ui->bcCBTable->resizeRowToContents(rowIdx);
}

void <<moduleClass>>::on_addbcCBButton_clicked()
{
    if (ui->bcCBCodesList->currentItem()!=0 && ui->bcCBCodesList->currentItem()->isSelected()) {

        addbcCBRow(ui->bcCBCodesList->currentItem()->text());
        //remove the item from the list
        ui->bcCBCodesList->takeItem(ui->bcCBCodesList->row(ui->bcCBCodesList->currentItem()));
    }

}

void <<moduleClass>>::on_delbcCBButton_clicked()
{
    if (ui->bcCBTable->currentItem() != 0 && ui->bcCBTable->currentItem()->isSelected()) {
        //add the item to the list
        int rowIdx = ui->bcCBCodesList->count();
        ui->bcCBCodesList->insertItem(rowIdx,ui->bcCBTable->item(ui->bcCBTable->currentRow(),0)->text());

        //remove the item from the list
        ui->bcCBTable->removeRow(ui->bcCBTable->currentRow());
        //ui->postProcessVarsList->takeItem(ui->postProcessVarsList->row(ui->postProcessVarsList->currentItem()));
    }
}
