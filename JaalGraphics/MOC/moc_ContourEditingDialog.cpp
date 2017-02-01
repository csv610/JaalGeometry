/****************************************************************************
** Meta object code from reading C++ file 'ContourEditingDialog.hpp'
**
** Created by: The Qt Meta Object Compiler version 67 (Qt 5.6.0)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "ContourEditingDialog.hpp"
#include <QtCore/qbytearray.h>
#include <QtCore/qmetatype.h>
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'ContourEditingDialog.hpp' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 67
#error "This file was generated using the moc from 5.6.0. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
struct qt_meta_stringdata_JContourEditingDialog_t {
    QByteArrayData data[7];
    char stringdata0[99];
};
#define QT_MOC_LITERAL(idx, ofs, len) \
    Q_STATIC_BYTE_ARRAY_DATA_HEADER_INITIALIZER_WITH_OFFSET(len, \
    qptrdiff(offsetof(qt_meta_stringdata_JContourEditingDialog_t, stringdata0) + ofs \
        - idx * sizeof(QByteArrayData)) \
    )
static const qt_meta_stringdata_JContourEditingDialog_t qt_meta_stringdata_JContourEditingDialog = {
    {
QT_MOC_LITERAL(0, 0, 21), // "JContourEditingDialog"
QT_MOC_LITERAL(1, 22, 14), // "reparameterize"
QT_MOC_LITERAL(2, 37, 0), // ""
QT_MOC_LITERAL(3, 38, 15), // "curveShortening"
QT_MOC_LITERAL(4, 54, 9), // "smoothing"
QT_MOC_LITERAL(5, 64, 22), // "openPolySimplifyDialog"
QT_MOC_LITERAL(6, 87, 11) // "closeDialog"

    },
    "JContourEditingDialog\0reparameterize\0"
    "\0curveShortening\0smoothing\0"
    "openPolySimplifyDialog\0closeDialog"
};
#undef QT_MOC_LITERAL

static const uint qt_meta_data_JContourEditingDialog[] = {

 // content:
       7,       // revision
       0,       // classname
       0,    0, // classinfo
       5,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       0,       // signalCount

 // slots: name, argc, parameters, tag, flags
       1,    0,   39,    2, 0x08 /* Private */,
       3,    0,   40,    2, 0x08 /* Private */,
       4,    0,   41,    2, 0x08 /* Private */,
       5,    0,   42,    2, 0x08 /* Private */,
       6,    0,   43,    2, 0x08 /* Private */,

 // slots: parameters
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,

       0        // eod
};

void JContourEditingDialog::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    if (_c == QMetaObject::InvokeMetaMethod) {
        JContourEditingDialog *_t = static_cast<JContourEditingDialog *>(_o);
        Q_UNUSED(_t)
        switch (_id) {
        case 0: _t->reparameterize(); break;
        case 1: _t->curveShortening(); break;
        case 2: _t->smoothing(); break;
        case 3: _t->openPolySimplifyDialog(); break;
        case 4: _t->closeDialog(); break;
        default: ;
        }
    }
    Q_UNUSED(_a);
}

const QMetaObject JContourEditingDialog::staticMetaObject = {
    { &QDialog::staticMetaObject, qt_meta_stringdata_JContourEditingDialog.data,
      qt_meta_data_JContourEditingDialog,  qt_static_metacall, Q_NULLPTR, Q_NULLPTR}
};


const QMetaObject *JContourEditingDialog::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->dynamicMetaObject() : &staticMetaObject;
}

void *JContourEditingDialog::qt_metacast(const char *_clname)
{
    if (!_clname) return Q_NULLPTR;
    if (!strcmp(_clname, qt_meta_stringdata_JContourEditingDialog.stringdata0))
        return static_cast<void*>(const_cast< JContourEditingDialog*>(this));
    if (!strcmp(_clname, "Ui::ContourEditingDialog"))
        return static_cast< Ui::ContourEditingDialog*>(const_cast< JContourEditingDialog*>(this));
    return QDialog::qt_metacast(_clname);
}

int JContourEditingDialog::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QDialog::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        if (_id < 5)
            qt_static_metacall(this, _c, _id, _a);
        _id -= 5;
    } else if (_c == QMetaObject::RegisterMethodArgumentMetaType) {
        if (_id < 5)
            *reinterpret_cast<int*>(_a[0]) = -1;
        _id -= 5;
    }
    return _id;
}
QT_END_MOC_NAMESPACE
