        if (!ui-><<editIdx>>Edit->text().isEmpty() && !ui-><<editIdx>>Edit->isHidden()) {
          file << " <<comma>>" << ui-><<editIdx>>Edit->text().toStdString();
        }
