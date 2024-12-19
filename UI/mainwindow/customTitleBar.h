#ifndef CUSTOMTITLEBAR_H
#define CUSTOMTITLEBAR_H

#include <QWidget>
#include <QLabel>
#include <QComboBox>
#include <QToolButton>
#include <QHBoxLayout>
#include <QStyle>
#include <QDockWidget>

class CustomTitleBar : public QWidget
{
    Q_OBJECT

public:
    explicit CustomTitleBar(const QString &title, const QStringList &sources, QDockWidget *dockWidget, QWidget *parent = nullptr);

signals:
    void signalSourceChanged(const QString &source);

private:
    QLabel *titleLabel;
    QComboBox *sourceComboBox;
    QToolButton *closeButton;
    QToolButton *floatButton;

    QDockWidget *m_dockWidget; // Store the dock widget pointer
};

#endif // CUSTOMTITLEBAR_H
