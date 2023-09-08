void <<moduleClass>>::on_addPostProcessButton_clicked()
{
    if (ui->postProcessVarsList->currentItem()!=0 && ui->postProcessVarsList->currentItem()->isSelected()) {
        //add the item to the table
        int rowIdx = ui->postProcessTable->rowCount();
        ui->postProcessTable->insertRow(rowIdx);
        QTableWidgetItem *newItem = new QTableWidgetItem(ui->postProcessVarsList->currentItem()->text());
        ui->postProcessTable->setItem(rowIdx, 0, newItem);

        //remove the item from the list
        ui->postProcessVarsList->takeItem(ui->postProcessVarsList->row(ui->postProcessVarsList->currentItem()));

    }
}

void <<moduleClass>>::on_delPostProcessButton_clicked()
{
    if (ui->postProcessTable->currentItem() != 0 && ui->postProcessTable->currentItem()->isSelected()) {
        //add the item to the list
        int rowIdx = ui->postProcessVarsList->count();
        ui->postProcessVarsList->insertItem(rowIdx,ui->postProcessTable->item(ui->postProcessTable->currentRow(),0)->text());

        //remove the item from the list
        ui->postProcessTable->removeRow(ui->postProcessTable->currentRow());
        //ui->postProcessVarsList->takeItem(ui->postProcessVarsList->row(ui->postProcessVarsList->currentItem()));
    }

}

void <<moduleClass>>::on_addElementSetButton_clicked()
{
    if (ui->elementSetVarsList->currentItem()!=0 && ui->elementSetVarsList->currentItem()->isSelected()) {
        //add the item to the table
        int rowIdx = ui->elementSetList->count();
        ui->elementSetList->insertItem(rowIdx,ui->elementSetVarsList->currentItem()->text());

        //remove the item from the list
        ui->elementSetVarsList->takeItem(ui->elementSetVarsList->row(ui->elementSetVarsList->currentItem()));

    }
}

void <<moduleClass>>::on_delElementSetButton_clicked()
{
    if (ui->elementSetList->currentItem()!=0 && ui->elementSetList->currentItem()->isSelected()) {
        //add the item to the table
        int rowIdx = ui->elementSetVarsList->count();
        ui->elementSetVarsList->insertItem(rowIdx,ui->elementSetList->currentItem()->text());

        //remove the item from the list
        ui->elementSetList->takeItem(ui->elementSetList->row(ui->elementSetList->currentItem()));
    }
}

void <<moduleClass>>::on_addNodeSetButton_clicked()
{
    if (ui->nodeSetVarsList->currentItem()!=0 && ui->nodeSetVarsList->currentItem()->isSelected()) {
        //add the item to the table
        int rowIdx = ui->nodeSetList->count();
        ui->nodeSetList->insertItem(rowIdx,ui->nodeSetVarsList->currentItem()->text());

        //remove the item from the list
        ui->nodeSetVarsList->takeItem(ui->nodeSetVarsList->row(ui->nodeSetVarsList->currentItem()));

    }

}

void <<moduleClass>>::on_delNodeSetButton_clicked()
{
    if (ui->nodeSetList->currentItem()!=0 && ui->nodeSetList->currentItem()->isSelected()) {
        //add the item to the table
        int rowIdx = ui->nodeSetVarsList->count();
        ui->nodeSetVarsList->insertItem(rowIdx,ui->nodeSetList->currentItem()->text());

        //remove the item from the list
        ui->nodeSetList->takeItem(ui->nodeSetList->row(ui->nodeSetList->currentItem()));
    }
}

void <<moduleClass>>::on_addBoundarySetButton_clicked()
{
    if (ui->boundarySetVarsList->currentItem()!=0 && ui->boundarySetVarsList->currentItem()->isSelected()) {
        //add the item to the table
        int rowIdx = ui->boundarySetList->count();
        ui->boundarySetList->insertItem(rowIdx,ui->boundarySetVarsList->currentItem()->text());

        //remove the item from the list
        ui->boundarySetVarsList->takeItem(ui->boundarySetVarsList->row(ui->boundarySetVarsList->currentItem()));

    }

}

void <<moduleClass>>::on_delBoundarySetButton_clicked()
{
    if (ui->boundarySetList->currentItem()!=0 && ui->boundarySetList->currentItem()->isSelected()) {
        //add the item to the table
        int rowIdx = ui->boundarySetVarsList->count();
        ui->boundarySetVarsList->insertItem(rowIdx,ui->boundarySetList->currentItem()->text());

        //remove the item from the list
        ui->boundarySetList->takeItem(ui->boundarySetList->row(ui->boundarySetList->currentItem()));
    }
}

void <<moduleClass>>::on_addWitnessPointsButton_clicked()
{
    if (ui->witnessPointsVarsList->currentItem()!=0 && ui->witnessPointsVarsList->currentItem()->isSelected()) {
        //add the item to the table
        int rowIdx = ui->witnessPointsList->count();
        ui->witnessPointsList->insertItem(rowIdx,ui->witnessPointsVarsList->currentItem()->text());
        //remove the item from the list
        ui->witnessPointsVarsList->takeItem(ui->witnessPointsVarsList->row(ui->witnessPointsVarsList->currentItem()));
    }
}

void <<moduleClass>>::on_delWitnessPointsButton_clicked()
{
    if (ui->witnessPointsList->currentItem()!=0 && ui->witnessPointsList->currentItem()->isSelected()) {
        //add the item to the table
        int rowIdx = ui->witnessPointsVarsList->count();
        ui->witnessPointsVarsList->insertItem(rowIdx,ui->witnessPointsList->currentItem()->text());

        //remove the item from the list
        ui->witnessPointsList->takeItem(ui->witnessPointsList->row(ui->witnessPointsList->currentItem()));
    }
}
