#ifndef CUSTOMTITLEBAR_H
#define CUSTOMTITLEBAR_H

#include <QWidget>
#include <QLabel>
#include <QComboBox>
#include <QToolButton>
#include <QHBoxLayout>
#include <QStyle>
#include <QDockWidget>
#include <QMainWindow>
#include <QEvent>
#include <QCheckBox>

class CustomTitleBar : public QWidget
{
    Q_OBJECT

public:
    explicit CustomTitleBar(const QString &title, const QStringList &sources, QDockWidget *dockWidget, QWidget *parent = nullptr);
    void updateSources(const QStringList &sources);

signals:
    void signalSourceChanged(int sourceIndex);
    void triggerAVisibilityChanged(bool visible);
    void triggerBVisibilityChanged(bool visible);
    void triggerOutVisibilityChanged(bool visible);
    void viewerTypeChanged(int type);

protected:
    bool eventFilter(QObject *obj, QEvent *event) override;

public slots:
    void updateSignalSources(const QStringList& newSources) {
        sourceComboBox->clear();
        sourceComboBox->addItems(newSources);
    }

private:
    QLabel *titleLabel;
    QComboBox *typeComboBox;
    QComboBox *sourceComboBox;
    QToolButton *closeButton;
    QToolButton *floatButton;
    QCheckBox *triggerABox;
    QCheckBox *triggerBBox;
    QCheckBox *triggerOutBox;
    QLabel* triggersLabel;

    QDockWidget *m_dockWidget; // Store the dock widget pointer
};

#endif // CUSTOMTITLEBAR_H
