/****************************************************************************
** Meta object code from reading C++ file 'GlobalSettingsDialog.hpp'
**
** Created by: The Qt Meta Object Compiler version 67 (Qt 5.6.0)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "GlobalSettingsDialog.hpp"
#include <QtCore/qbytearray.h>
#include <QtCore/qmetatype.h>
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'GlobalSettingsDialog.hpp' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 67
#error "This file was generated using the moc from 5.6.0. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
struct qt_meta_stringdata_JGlobalSettingsDialog_t {
    QByteArrayData data[14];
    char stringdata0[195];
};
#define QT_MOC_LITERAL(idx, ofs, len) \
    Q_STATIC_BYTE_ARRAY_DATA_HEADER_INITIALIZER_WITH_OFFSET(len, \
    qptrdiff(offsetof(qt_meta_stringdata_JGlobalSettingsDialog_t, stringdata0) + ofs \
        - idx * sizeof(QByteArrayData)) \
    )
static const qt_meta_stringdata_JGlobalSettingsDialog_t qt_meta_stringdata_JGlobalSettingsDialog = {
    {
QT_MOC_LITERAL(0, 0, 21), // "JGlobalSettingsDialog"
QT_MOC_LITERAL(1, 22, 15), // "openFontsDialog"
QT_MOC_LITERAL(2, 38, 0), // ""
QT_MOC_LITERAL(3, 39, 16), // "openLightsDialog"
QT_MOC_LITERAL(4, 56, 15), // "openFloorDialog"
QT_MOC_LITERAL(5, 72, 16), // "openCameraDialog"
QT_MOC_LITERAL(6, 89, 20), // "openScreenShotDialog"
QT_MOC_LITERAL(7, 110, 11), // "closeDialog"
QT_MOC_LITERAL(8, 122, 13), // "keyPressEvent"
QT_MOC_LITERAL(9, 136, 10), // "QKeyEvent*"
QT_MOC_LITERAL(10, 147, 1), // "e"
QT_MOC_LITERAL(11, 149, 18), // "setBackgroundColor"
QT_MOC_LITERAL(12, 168, 16), // "checkBoundingBox"
QT_MOC_LITERAL(13, 185, 9) // "checkAxis"

    },
    "JGlobalSettingsDialog\0openFontsDialog\0"
    "\0openLightsDialog\0openFloorDialog\0"
    "openCameraDialog\0openScreenShotDialog\0"
    "closeDialog\0keyPressEvent\0QKeyEvent*\0"
    "e\0setBackgroundColor\0checkBoundingBox\0"
    "checkAxis"
};
#undef QT_MOC_LITERAL

static const uint qt_meta_data_JGlobalSettingsDialog[] = {

 // content:
       7,       // revision
       0,       // classname
       0,    0, // classinfo
      10,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       0,       // signalCount

 // slots: name, argc, parameters, tag, flags
       1,    0,   64,    2, 0x08 /* Private */,
       3,    0,   65,    2, 0x08 /* Private */,
       4,    0,   66,    2, 0x08 /* Private */,
       5,    0,   67,    2, 0x08 /* Private */,
       6,    0,   68,    2, 0x08 /* Private */,
       7,    0,   69,    2, 0x08 /* Private */,
       8,    1,   70,    2, 0x08 /* Private */,
      11,    0,   73,    2, 0x08 /* Private */,
      12,    0,   74,    2, 0x08 /* Private */,
      13,    0,   75,    2, 0x08 /* Private */,

 // slots: parameters
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void, 0x80000000 | 9,   10,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,

       0        // eod
};

void JGlobalSettingsDialog::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    if (_c == QMetaObject::InvokeMetaMethod) {
        JGlobalSettingsDialog *_t = static_cast<JGlobalSettingsDialog *>(_o);
        Q_UNUSED(_t)
        switch (_id) {
        case 0: _t->openFontsDialog(); break;
        case 1: _t->openLightsDialog(); break;
        case 2: _t->openFloorDialog(); break;
        case 3: _t->openCameraDialog(); break;
        case 4: _t->openScreenShotDialog(); break;
        case 5: _t->closeDialog(); break;
        case 6: _t->keyPressEvent((*reinterpret_cast< QKeyEvent*(*)>(_a[1]))); break;
        case 7: _t->setBackgroundColor(); break;
        case 8: _t->checkBoundingBox(); break;
        case 9: _t->checkAxis(); break;
        default: ;
        }
    }
}

const QMetaObject JGlobalSettingsDialog::staticMetaObject = {
    { &QDialog::staticMetaObject, qt_meta_stringdata_JGlobalSettingsDialog.data,
      qt_meta_data_JGlobalSettingsDialog,  qt_static_metacall, Q_NULLPTR, Q_NULLPTR}
};


const QMetaObject *JGlobalSettingsDialog::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->dynamicMetaObject() : &staticMetaObject;
}

void *JGlobalSettingsDialog::qt_metacast(const char *_clname)
{
    if (!_clname) return Q_NULLPTR;
    if (!strcmp(_clname, qt_meta_stringdata_JGlobalSettingsDialog.stringdata0))
        return static_cast<void*>(const_cast< JGlobalSettingsDialog*>(this));
    if (!strcmp(_clname, "Ui::GlobalSettingsDialog"))
        return static_cast< Ui::GlobalSettingsDialog*>(const_cast< JGlobalSettingsDialog*>(this));
    return QDialog::qt_metacast(_clname);
}

int JGlobalSettingsDialog::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QDialog::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        if (_id < 10)
            qt_static_metacall(this, _c, _id, _a);
        _id -= 10;
    } else if (_c == QMetaObject::RegisterMethodArgumentMetaType) {
        if (_id < 10)
            *reinterpret_cast<int*>(_a[0]) = -1;
        _id -= 10;
    }
    return _id;
}
QT_END_MOC_NAMESPACE
