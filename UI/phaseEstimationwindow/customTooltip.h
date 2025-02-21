#ifndef CUSTOMTOOLTIP_H
#define CUSTOMTOOLTIP_H

#include <QWidget>
#include <QLabel>
#include <QVBoxLayout>

class CustomTooltip : public QWidget {
    Q_OBJECT
public:
    explicit CustomTooltip(QWidget *parent = nullptr) : QWidget(parent, Qt::ToolTip) {
        setWindowFlags(Qt::ToolTip);
        QVBoxLayout *layout = new QVBoxLayout(this);
        label = new QLabel(this);
        layout->addWidget(label);
        setLayout(layout);
    }

    void setText(const QString &text) {
        label->setText(text);
    }

private:
    QLabel *label;
};

#endif // CUSTOMTOOLTIP_H