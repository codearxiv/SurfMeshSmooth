//     Copyright (C) 2021 Piotr (Peter) Beben <pdbcas2@gmail.com>
//     See LICENSE included with this distribution.

#include "RandomSurfDialog.h"
#include "get_field.h"
#include "constants.h"

#include <QFormLayout>
#include <QLabel>
#include <QLineEdit>
#include <QDialogButtonBox>
#include <QObject>
#include <QWidget>
#include <QIntValidator>


RandomSurfDialog::RandomSurfDialog(QWidget *parent) : QDialog(parent)
{
	validator = new QIntValidator(1, int_infinity, this);

	form = new QFormLayout(this);
	form->addRow(new QLabel(
					 "Create a new mesh from a random surface"));
	nDimLineEdit = new QLineEdit(this);
	nDimLineEdit->setValidator(validator);
	nDimLineEdit->setText("100");
	form->addRow(QString("Number of vertices along an axis:"), nDimLineEdit);

	buttonBox = new QDialogButtonBox(
				QDialogButtonBox::Ok | QDialogButtonBox::Cancel,
				Qt::Horizontal, this);

	form->addRow(buttonBox);
	connect(buttonBox, SIGNAL(accepted()), this, SLOT(accept()));
	connect(buttonBox, SIGNAL(rejected()), this, SLOT(reject()));

}


bool RandomSurfDialog::getFields(size_t& nDim) const
{

	bool ok = get_integer_field(nDimLineEdit, validator, nDim);
	return ok;
}


