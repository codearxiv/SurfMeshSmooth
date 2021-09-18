//     Copyright (C) 2021 Piotr (Peter) Beben <pdbcas2@gmail.com>
//     See LICENSE included with this distribution.

#ifndef RECONSTRUCTDIALOG_H
#define RECONSTRUCTDIALOG_H

#include "constants.h"

#include <QWidget>
#include <QDialog>
#include <QObject>

QT_BEGIN_NAMESPACE
class QFormLayout;
class QLineEdit;
class QDialogButtonBox;
class QIntValidator;
class QDoubleValidator;
class QComboBox;
QT_END_NAMESPACE

class SmoothDialog : public QDialog
{
	Q_OBJECT

public:
	SmoothDialog(QWidget *parent = nullptr);
	int getFields(int& nSweeps) const;

public slots:

private:
	QFormLayout *form;
	QLineEdit *nSweepsLineEdit;

	QDialogButtonBox *buttonBox;
    QIntValidator *intValidator;
    QDoubleValidator *doubleValidator;

};

#endif // RECONSTRUCTDIALOG_H
