/****************************************************************************
** Meta object code from reading C++ file 'MainWindow.hpp'
**
** Created by: The Qt Meta Object Compiler version 67 (Qt 5.6.0)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "../include/MainWindow.hpp"
#include <QtCore/qbytearray.h>
#include <QtCore/qmetatype.h>
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'MainWindow.hpp' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 67
#error "This file was generated using the moc from 5.6.0. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
struct qt_meta_stringdata_JaalMainWindow_t {
    QByteArrayData data[10];
    char stringdata0[134];
};
#define QT_MOC_LITERAL(idx, ofs, len) \
    Q_STATIC_BYTE_ARRAY_DATA_HEADER_INITIALIZER_WITH_OFFSET(len, \
    qptrdiff(offsetof(qt_meta_stringdata_JaalMainWindow_t, stringdata0) + ofs \
        - idx * sizeof(QByteArrayData)) \
    )
static const qt_meta_stringdata_JaalMainWindow_t qt_meta_stringdata_JaalMainWindow = {
    {
QT_MOC_LITERAL(0, 0, 14), // "JaalMainWindow"
QT_MOC_LITERAL(1, 15, 11), // "resizeEvent"
QT_MOC_LITERAL(2, 27, 0), // ""
QT_MOC_LITERAL(3, 28, 13), // "QResizeEvent*"
QT_MOC_LITERAL(4, 42, 1), // "e"
QT_MOC_LITERAL(5, 44, 17), // "openNewDataDialog"
QT_MOC_LITERAL(6, 62, 19), // "openMeshToolsDialog"
QT_MOC_LITERAL(7, 82, 24), // "openGlobalSettingsDialog"
QT_MOC_LITERAL(8, 107, 21), // "openObjectsListDialog"
QT_MOC_LITERAL(9, 129, 4) // "Quit"

    },
    "JaalMainWindow\0resizeEvent\0\0QResizeEvent*\0"
    "e\0openNewDataDialog\0openMeshToolsDialog\0"
    "openGlobalSettingsDialog\0openObjectsListDialog\0"
    "Quit"
};
#undef QT_MOC_LITERAL

static const uint qt_meta_data_JaalMainWindow[] = {

 // content:
       7,       // revision
       0,       // classname
       0,    0, // classinfo
       6,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       0,       // signalCount

 // slots: name, argc, parameters, tag, flags
       1,    1,   44,    2, 0x08 /* Private */,
       5,    0,   47,    2, 0x08 /* Private */,
       6,    0,   48,    2, 0x08 /* Private */,
       7,    0,   49,    2, 0x08 /* Private */,
       8,    0,   50,    2, 0x08 /* Private */,
       9,    0,   51,    2, 0x08 /* Private */,

 // slots: parameters
    QMetaType::Void, 0x80000000 | 3,    4,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,

       0        // eod
};

void JaalMainWindow::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    if (_c == QMetaObject::InvokeMetaMethod) {
        JaalMainWindow *_t = static_cast<JaalMainWindow *>(_o);
        Q_UNUSED(_t)
        switch (_id) {
        case 0: _t->resizeEvent((*reinterpret_cast< QResizeEvent*(*)>(_a[1]))); break;
        case 1: _t->openNewDataDialog(); break;
        case 2: _t->openMeshToolsDialog(); break;
        case 3: _t->openGlobalSettingsDialog(); break;
        case 4: _t->openObjectsListDialog(); break;
        case 5: _t->Quit(); break;
        default: ;
        }
    }
}

const QMetaObject JaalMainWindow::staticMetaObject = {
    { &QMainWindow::staticMetaObject, qt_meta_stringdata_JaalMainWindow.data,
      qt_meta_data_JaalMainWindow,  qt_static_metacall, Q_NULLPTR, Q_NULLPTR}
};


const QMetaObject *JaalMainWindow::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->dynamicMetaObject() : &staticMetaObject;
}

void *JaalMainWindow::qt_metacast(const char *_clname)
{
    if (!_clname) return Q_NULLPTR;
    if (!strcmp(_clname, qt_meta_stringdata_JaalMainWindow.stringdata0))
        return static_cast<void*>(const_cast< JaalMainWindow*>(this));
    if (!strcmp(_clname, "Ui::JaalMainWindow"))
        return static_cast< Ui::JaalMainWindow*>(const_cast< JaalMainWindow*>(this));
    return QMainWindow::qt_metacast(_clname);
}

int JaalMainWindow::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QMainWindow::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        if (_id < 6)
            qt_static_metacall(this, _c, _id, _a);
        _id -= 6;
    } else if (_c == QMetaObject::RegisterMethodArgumentMetaType) {
        if (_id < 6)
            *reinterpret_cast<int*>(_a[0]) = -1;
        _id -= 6;
    }
    return _id;
}
QT_END_MOC_NAMESPACE
