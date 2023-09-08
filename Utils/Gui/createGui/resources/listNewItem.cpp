        //COL#colNum
        newEdit = new QLineEdit();
        newEdit->setValidator(#validator);
        newEdit->setToolTip("#tooltip");
        ui->ie00151List->setCellWidget(rowIdx, #colNum, newEdit);
