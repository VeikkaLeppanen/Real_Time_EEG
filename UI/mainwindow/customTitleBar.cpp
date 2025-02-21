#include "customTitleBar.h"
#include <QCheckBox>

CustomTitleBar::CustomTitleBar(const QString &title, const QStringList &sources, QDockWidget *dockWidget, QWidget *parent)
    : QWidget(parent), m_dockWidget(dockWidget)
{
    // Create the title label
    titleLabel = new QLabel(title, this);

    // Create the type combo box
    typeComboBox = new QComboBox(this);
    typeComboBox->addItem("Signal viewer");
    typeComboBox->addItem("Response monitor");
    typeComboBox->setFixedWidth(120);

    // Create the source combo box
    sourceComboBox = new QComboBox(this);
    sourceComboBox->addItems(sources);
    sourceComboBox->setFixedWidth(150);

    // Set initial source based on viewer number (extracted from title)
    QString numberStr = title.split(" ").last();
    bool ok;
    int viewerNumber = numberStr.toInt(&ok);
    if (ok && viewerNumber > 0 && (viewerNumber - 1) < sources.size()) {
        sourceComboBox->setCurrentIndex(viewerNumber - 1);
    }

    // Connect the combo box signal to the custom signal
    connect(sourceComboBox, QOverload<int>::of(&QComboBox::currentIndexChanged), this, &CustomTitleBar::signalSourceChanged);

    // Create the close button
    closeButton = new QToolButton(this);
    QIcon closeIcon = style()->standardIcon(QStyle::SP_TitleBarCloseButton);
    closeButton->setIcon(closeIcon);

    // Create the float button
    floatButton = new QToolButton(this);
    QIcon floatIcon = style()->standardIcon(QStyle::SP_TitleBarNormalButton);
    floatButton->setIcon(floatIcon);

    // Create the checkboxes
    QLabel* triggersLabel = new QLabel("Triggers:", this);
    triggerABox = new QCheckBox("A", this);
    triggerBBox = new QCheckBox("B", this);
    triggerOutBox = new QCheckBox("Out", this);

    // Style only the checkbox indicators with background colors
    triggerABox->setStyleSheet("QCheckBox { color: blue; }");
    
    triggerBBox->setStyleSheet("QCheckBox { color: green; }");
    
    triggerOutBox->setStyleSheet("QCheckBox { color: skyblue; }");

    // Connect the buttons to the dock widget's slots using m_dockWidget
    connect(closeButton, &QToolButton::clicked, m_dockWidget, &QWidget::close);
    connect(floatButton, &QToolButton::clicked, [this]() {
        if (m_dockWidget)
        {
            bool isFloating = m_dockWidget->isFloating();
            m_dockWidget->setFloating(!isFloating);
        }
    });

    // Connect checkbox signals
    connect(triggerABox, &QCheckBox::toggled, this, &CustomTitleBar::triggerAVisibilityChanged);
    connect(triggerBBox, &QCheckBox::toggled, this, &CustomTitleBar::triggerBVisibilityChanged);
    connect(triggerOutBox, &QCheckBox::toggled, this, &CustomTitleBar::triggerOutVisibilityChanged);

    // Connect the type combo box signal
    connect(typeComboBox, QOverload<int>::of(&QComboBox::currentIndexChanged),
            this, &CustomTitleBar::viewerTypeChanged);

    // Create the layout
    QHBoxLayout *layout = new QHBoxLayout(this);
    layout->setContentsMargins(5, 2, 5, 2);
    layout->addWidget(titleLabel);
    layout->addWidget(typeComboBox);
    layout->addStretch();
    layout->addWidget(triggersLabel);
    layout->addWidget(triggerABox);
    layout->addWidget(triggerBBox);
    layout->addWidget(triggerOutBox);
    layout->addWidget(sourceComboBox);
    layout->addWidget(floatButton);
    layout->addWidget(closeButton);

    setLayout(layout);

    // Add event filter to combo box
    sourceComboBox->installEventFilter(this);
}

bool CustomTitleBar::eventFilter(QObject *obj, QEvent *event)
{
    if (obj == sourceComboBox && event->type() == QEvent::MouseButtonPress) {
        // When combo box is clicked, request updated sources from MainWindow
        if (QMainWindow *mainWindow = qobject_cast<QMainWindow*>(window())) {
            QStringList sources = mainWindow->property("signalSources").toStringList();
            updateSources(sources);
        }
    }
    return QWidget::eventFilter(obj, event);
}

void CustomTitleBar::updateSources(const QStringList &sources)
{
    QString currentText = sourceComboBox->currentText();
    sourceComboBox->clear();
    sourceComboBox->addItems(sources);
    
    // Restore previous selection if it still exists in the new list
    int index = sourceComboBox->findText(currentText);
    if (index >= 0) {
        sourceComboBox->setCurrentIndex(index);
    }
}
