//     Copyright (C) 2021 Piotr (Peter) Beben <pdbcas2@gmail.com>
//     See LICENSE included with this distribution.

#ifndef NOISEDIALOG_H
#define NOISEDIALOG_H

#include <QWidget>
#include <QDialog>
#include <QObject>

QT_BEGIN_NAMESPACE
class QFormLayout;
class QLineEdit;
class QDialogButtonBox;
class QIntValidator;
QT_END_NAMESPACE

class NoiseDialog : public QDialog
{
	Q_OBJECT

public:
	NoiseDialog(QWidget *parent = nullptr);
	int getFields(int& nSweeps) const;

private:
	QFormLayout *form;
	QLineEdit *nSweepsLineEdit;
	QDialogButtonBox *buttonBox;
	QIntValidator *validator;

};

#endif // NOISEDIALOG_H
