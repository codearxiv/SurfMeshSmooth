//     Copyright (C) 2021 Piotr (Peter) Beben <pdbcas2@gmail.com>
//     See LICENSE included with this distribution.

#ifndef SETRANDOMDIALOG_H
#define SETRANDOMDIALOG_H

#include <QWidget>
#include <QDialog>
#include <QObject>

QT_BEGIN_NAMESPACE
class QFormLayout;
class QLineEdit;
class QDialogButtonBox;
class QIntValidator;
QT_END_NAMESPACE

class RandomSurfDialog : public QDialog
{
	Q_OBJECT

public:
    RandomSurfDialog(QWidget *parent = nullptr);
	bool getFields(size_t& nDim) const;

private:
	QFormLayout *form;
	QLineEdit *nDimLineEdit;
	QDialogButtonBox *buttonBox;
	QIntValidator *validator;

};

#endif // SETRANDOMDIALOG_H
