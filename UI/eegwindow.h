#ifndef EEGWINDOW_H
#define EEGWINDOW_H

#include <QObject>
#include <QMainWindow>
#include "./ui_eegwindow.h"

QT_BEGIN_NAMESPACE
namespace Ui {
class EegWindow;
}
QT_END_NAMESPACE

class eegWindow : public QMainWindow
{
    Q_OBJECT
public:
    explicit eegWindow(QWidget *parent = nullptr);
    ~eegWindow();

signals:

private:
    Ui::EegWindow *ui;
};

#endif // EEGWINDOW_H
