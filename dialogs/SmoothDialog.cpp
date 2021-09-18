//     Copyright (C) 2021 Piotr (Peter) Beben <pdbcas2@gmail.com>
//     See LICENSE included with this distribution.

#include "SmoothDialog.h"
#include "get_field.h"
#include "constants.h"

#include <QFormLayout>
#include <QLabel>
#include <QLineEdit>
#include <QComboBox>
#include <QDialogButtonBox>
#include <QObject>
#include <QWidget>
#include <QIntValidator>
#include <QDoubleValidator>


SmoothDialog::SmoothDialog(QWidget *parent) : QDialog(parent)
{
    intValidator = new QIntValidator(1, int_infinity, this);
    doubleValidator = new QDoubleValidator(0.0, double_infinity, 2, this);

	form = new QFormLayout(this);
	form->addRow(new QLabel("Smooth the mesh"));

	nSweepsLineEdit = new QLineEdit(this);
    nSweepsLineEdit->setValidator(intValidator);
	nSweepsLineEdit->setText("15");
	nSweepsLineEdit->setToolTip(
				QString("This sets the number of smoothing iterations\n"));
	form->addRow(QString("Number of smoothing iterations:"), nSweepsLineEdit);

	buttonBox = new QDialogButtonBox(
				QDialogButtonBox::Ok | QDialogButtonBox::Cancel,
				Qt::Horizontal, this);

	form->addRow(buttonBox);
	connect(buttonBox, SIGNAL(accepted()), this, SLOT(accept()));
	connect(buttonBox, SIGNAL(rejected()), this, SLOT(reject()));



}

//---------------------------------------------------------

int SmoothDialog::getFields(int& nSweeps) const
{
	bool ok;
	size_t temp;

	ok = get_integer_field(nSweepsLineEdit, intValidator, temp);
	if(!ok) return -1;
	nSweeps = int(temp);

	return 0;
}
