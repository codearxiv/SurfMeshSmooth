//     Copyright (C) 2021 Piotr (Peter) Beben <pdbcas2@gmail.com>
//     See LICENSE included with this distribution.

#include "NormalsDialog.h"
#include "get_field.h"
#include "constants.h"

#include <QMainWindow>
#include <QFormLayout>
#include <QLabel>
#include <QLineEdit>
#include <QDialogButtonBox>
#include <QObject>
#include <QWidget>
#include <QIntValidator>


NormalsDialog::NormalsDialog(QWidget *parent) : QDialog(parent)
{
	form = new QFormLayout(this);
	form->addRow(new QLabel("Approximate mesh vertex normals"));

	buttonBox = new QDialogButtonBox(
				QDialogButtonBox::Ok | QDialogButtonBox::Cancel,
				Qt::Horizontal, this);

	form->addRow(buttonBox);
	connect(buttonBox, SIGNAL(accepted()), this, SLOT(accept()));
	connect(buttonBox, SIGNAL(rejected()), this, SLOT(reject()));

}

