          tmpLineEdit = qobject_cast<QLineEdit*>(ui->#listName->cellWidget(i,#rownum));
          if (!tmpLineEdit->text().isEmpty()) {
            file << "#comma #rowname=" << tmpLineEdit->text().toStdString();
          }
