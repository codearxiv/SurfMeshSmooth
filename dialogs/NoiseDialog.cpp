//     Copyright (C) 2021 Piotr (Peter) Beben <pdbcas2@gmail.com>
//     See LICENSE included with this distribution.

#include "NoiseDialog.h"
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


NoiseDialog::NoiseDialog(QWidget *parent) : QDialog(parent)
{
	validator = new QIntValidator(1, int_infinity, this);

	form = new QFormLayout(this);
	form->addRow(new QLabel(
					 "Add random noise to mesh vertex positions"));
	nSweepsLineEdit = new QLineEdit(this);
	nSweepsLineEdit->setValidator(validator);
	nSweepsLineEdit->setText("1");
	form->addRow(QString("Number of noising sweeps:"), nSweepsLineEdit);

	buttonBox = new QDialogButtonBox(
				QDialogButtonBox::Ok | QDialogButtonBox::Cancel,
				Qt::Horizontal, this);

	form->addRow(buttonBox);
	connect(buttonBox, SIGNAL(accepted()), this, SLOT(accept()));
	connect(buttonBox, SIGNAL(rejected()), this, SLOT(reject()));

}


int NoiseDialog::getFields(int& nSweeps) const
{
	bool ok;
	size_t temp;

	ok = get_integer_field(nSweepsLineEdit, validator, temp);
	if(!ok) return -1;
	nSweeps = int(temp);

	return 0;
}

