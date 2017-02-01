/****************************************************************************
** Meta object code from reading C++ file 'LSystemDialog.hpp'
**
** Created by: The Qt Meta Object Compiler version 67 (Qt 5.5.1)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "LSystemDialog.hpp"
#include <QtCore/qbytearray.h>
#include <QtCore/qmetatype.h>
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'LSystemDialog.hpp' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 67
#error "This file was generated using the moc from 5.5.1. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
struct qt_meta_stringdata_LSystemDialog_t {
    QByteArrayData data[8];
    char stringdata0[76];
};
#define QT_MOC_LITERAL(idx, ofs, len) \
    Q_STATIC_BYTE_ARRAY_DATA_HEADER_INITIALIZER_WITH_OFFSET(len, \
    qptrdiff(offsetof(qt_meta_stringdata_LSystemDialog_t, stringdata0) + ofs \
        - idx * sizeof(QByteArrayData)) \
    )
static const qt_meta_stringdata_LSystemDialog_t qt_meta_stringdata_LSystemDialog = {
    {
        QT_MOC_LITERAL(0, 0, 13), // "LSystemDialog"
        QT_MOC_LITERAL(1, 14, 13), // "keyPressEvent"
        QT_MOC_LITERAL(2, 28, 0), // ""
        QT_MOC_LITERAL(3, 29, 10), // "QKeyEvent*"
        QT_MOC_LITERAL(4, 40, 1), // "e"
        QT_MOC_LITERAL(5, 42, 8), // "genShape"
        QT_MOC_LITERAL(6, 51, 11), // "applyDialog"
        QT_MOC_LITERAL(7, 63, 12) // "cancelDialog"

    },
    "LSystemDialog\0keyPressEvent\0\0QKeyEvent*\0"
    "e\0genShape\0applyDialog\0cancelDialog"
};
#undef QT_MOC_LITERAL

static const uint qt_meta_data_LSystemDialog[] = {

// content:
    7,       // revision
    0,       // classname
    0,    0, // classinfo
    4,   14, // methods
    0,    0, // properties
    0,    0, // enums/sets
    0,    0, // constructors
    0,       // flags
    0,       // signalCount

// slots: name, argc, parameters, tag, flags
    1,    1,   34,    2, 0x08 /* Private */,
    5,    0,   37,    2, 0x08 /* Private */,
    6,    0,   38,    2, 0x08 /* Private */,
    7,    0,   39,    2, 0x08 /* Private */,

// slots: parameters
    QMetaType::Void, 0x80000000 | 3,    4,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,

    0        // eod
};

void LSystemDialog::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    if (_c == QMetaObject::InvokeMetaMethod) {
        LSystemDialog *_t = static_cast<LSystemDialog *>(_o);
        Q_UNUSED(_t)
        switch (_id) {
        case 0:
            _t->keyPressEvent((*reinterpret_cast< QKeyEvent*(*)>(_a[1])));
            break;
        case 1:
            _t->genShape();
            break;
        case 2:
            _t->applyDialog();
            break;
        case 3:
            _t->cancelDialog();
            break;
        default:
            ;
        }
    }
}

const QMetaObject LSystemDialog::staticMetaObject = {
    {   &QDialog::staticMetaObject, qt_meta_stringdata_LSystemDialog.data,
        qt_meta_data_LSystemDialog,  qt_static_metacall, Q_NULLPTR, Q_NULLPTR
    }
};


const QMetaObject *LSystemDialog::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->dynamicMetaObject() : &staticMetaObject;
}

void *LSystemDialog::qt_metacast(const char *_clname)
{
    if (!_clname) return Q_NULLPTR;
    if (!strcmp(_clname, qt_meta_stringdata_LSystemDialog.stringdata0))
        return static_cast<void*>(const_cast< LSystemDialog*>(this));
    if (!strcmp(_clname, "Ui::LSystemDialog"))
        return static_cast< Ui::LSystemDialog*>(const_cast< LSystemDialog*>(this));
    return QDialog::qt_metacast(_clname);
}

int LSystemDialog::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QDialog::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        if (_id < 4)
            qt_static_metacall(this, _c, _id, _a);
        _id -= 4;
    } else if (_c == QMetaObject::RegisterMethodArgumentMetaType) {
        if (_id < 4)
            *reinterpret_cast<int*>(_a[0]) = -1;
        _id -= 4;
    }
    return _id;
}
QT_END_MOC_NAMESPACE
