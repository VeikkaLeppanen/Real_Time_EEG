#include "customTitleBar.h"

CustomTitleBar::CustomTitleBar(const QString &title, const QStringList &sources, QDockWidget *dockWidget, QWidget *parent)
    : QWidget(parent), m_dockWidget(dockWidget)
{
    // Create the title label
    titleLabel = new QLabel(title, this);

    // Create the combo box
    sourceComboBox = new QComboBox(this);
    sourceComboBox->addItems(sources);

    // Connect the combo box signal to the custom signal
    connect(sourceComboBox, &QComboBox::currentTextChanged, this, &CustomTitleBar::signalSourceChanged);

    // Create the close button
    closeButton = new QToolButton(this);
    QIcon closeIcon = style()->standardIcon(QStyle::SP_TitleBarCloseButton);
    closeButton->setIcon(closeIcon);

    // Create the float button
    floatButton = new QToolButton(this);
    QIcon floatIcon = style()->standardIcon(QStyle::SP_TitleBarNormalButton);
    floatButton->setIcon(floatIcon);

    // Connect the buttons to the dock widget's slots using m_dockWidget
    connect(closeButton, &QToolButton::clicked, m_dockWidget, &QWidget::close);
    connect(floatButton, &QToolButton::clicked, [this]() {
        if (m_dockWidget)
        {
            bool isFloating = m_dockWidget->isFloating();
            m_dockWidget->setFloating(!isFloating);
        }
    });

    // Create the layout
    QHBoxLayout *layout = new QHBoxLayout(this);
    layout->setContentsMargins(5, 2, 5, 2); // Adjust margins as needed
    layout->addWidget(titleLabel);
    layout->addStretch();
    layout->addWidget(sourceComboBox);
    layout->addWidget(floatButton);
    layout->addWidget(closeButton);

    setLayout(layout);
}
