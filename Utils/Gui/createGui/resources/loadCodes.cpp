    ui->bcCBCodesList->clear();
    foreach (QString code, boundaryCodes.keys()) {
        ui->bcCBCodesList->addItem(code);
    }
    ui->bcCNCodesList->clear();
    ui->bcCNPCodesList->clear();
    foreach (QString code, nodeCodes.keys()) {
         ui->bcCNCodesList->addItem(code);
         ui->bcCNPCodesList->addItem(code);
    }
