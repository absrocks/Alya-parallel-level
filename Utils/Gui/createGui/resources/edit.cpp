        if (!ui-><<editIdx>>Edit->text().isEmpty() && !ui-><<editIdx>>Edit->isHidden()) {
          file << " <<comma>>" << "<<editName>>= " << ui-><<editIdx>>Edit->text().toStdString();
        }
