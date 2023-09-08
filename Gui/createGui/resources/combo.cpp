      if (ui-><<comboIdx>>Combo->currentText() != "EMPTY" && !ui-><<comboIdx>>Combo->isHidden()) {
          file << " <<comma>>" << ui-><<comboIdx>>Combo->currentText().toStdString();
      }
