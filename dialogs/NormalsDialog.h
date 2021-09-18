//     Copyright (C) 2021 Piotr (Peter) Beben <pdbcas2@gmail.com>
//     See LICENSE included with this distribution.

#ifndef NORMALSDIALOG_H
#define NORMALSDIALOG_H

#include <QWidget>
#include <QDialog>
#include <QObject>

QT_BEGIN_NAMESPACE
class QFormLayout;
class QLineEdit;
class QDialogButtonBox;
class QIntValidator;
QT_END_NAMESPACE

class NormalsDialog : public QDialog
{
	Q_OBJECT

public:
	NormalsDialog(QWidget *parent = nullptr);

private:
	QFormLayout *form;
	QDialogButtonBox *buttonBox;

};


#endif // NORMALSDIALOG_H
