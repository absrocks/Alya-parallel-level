void <<moduleClass>>::on_#listNamePlus_clicked()
{
        int rowIdx = ui->#listName->rowCount();
        ui->#listName->insertRow(rowIdx);
        QLineEdit *newEdit;
        <<listColumns>>
        QTableWidgetItem *newItem = new QTableWidgetItem();
        ui->#listName->setItem(rowIdx, 0, newItem);
}

void <<moduleClass>>::on_#listNameMinus_clicked()
{
    if (ui->#listName->currentItem() != 0 && ui->#listName->currentItem()->isSelected()) {
        //remove the item from the list
        ui->#listName->removeRow(ui->#listName->currentRow());
    } else {
        QMessageBox::information(this,"Remove rows", "Click on the number row you want to delete");
    }

}
