/****************************************************************************
** Meta object code from reading C++ file 'SingletDialog.hpp'
**
** Created by: The Qt Meta Object Compiler version 67 (Qt 5.5.1)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "SingletDialog.hpp"
#include <QtCore/qbytearray.h>
#include <QtCore/qmetatype.h>
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'SingletDialog.hpp' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 67
#error "This file was generated using the moc from 5.5.1. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
struct qt_meta_stringdata_JSingletDialog_t {
    QByteArrayData data[5];
    char stringdata0[50];
};
#define QT_MOC_LITERAL(idx, ofs, len) \
    Q_STATIC_BYTE_ARRAY_DATA_HEADER_INITIALIZER_WITH_OFFSET(len, \
    qptrdiff(offsetof(qt_meta_stringdata_JSingletDialog_t, stringdata0) + ofs \
        - idx * sizeof(QByteArrayData)) \
    )
static const qt_meta_stringdata_JSingletDialog_t qt_meta_stringdata_JSingletDialog = {
    {
        QT_MOC_LITERAL(0, 0, 14), // "JSingletDialog"
        QT_MOC_LITERAL(1, 15, 14), // "searchSinglets"
        QT_MOC_LITERAL(2, 30, 0), // ""
        QT_MOC_LITERAL(3, 31, 8), // "setColor"
        QT_MOC_LITERAL(4, 40, 9) // "removeAll"

    },
    "JSingletDialog\0searchSinglets\0\0setColor\0"
    "removeAll"
};
#undef QT_MOC_LITERAL

static const uint qt_meta_data_JSingletDialog[] = {

// content:
    7,       // revision
    0,       // classname
    0,    0, // classinfo
    3,   14, // methods
    0,    0, // properties
    0,    0, // enums/sets
    0,    0, // constructors
    0,       // flags
    0,       // signalCount

// slots: name, argc, parameters, tag, flags
    1,    0,   29,    2, 0x08 /* Private */,
    3,    0,   30,    2, 0x08 /* Private */,
    4,    0,   31,    2, 0x08 /* Private */,

// slots: parameters
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,

    0        // eod
};

void JSingletDialog::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    if (_c == QMetaObject::InvokeMetaMethod) {
        JSingletDialog *_t = static_cast<JSingletDialog *>(_o);
        Q_UNUSED(_t)
        switch (_id) {
        case 0:
            _t->searchSinglets();
            break;
        case 1:
            _t->setColor();
            break;
        case 2:
            _t->removeAll();
            break;
        default:
            ;
        }
    }
    Q_UNUSED(_a);
}

const QMetaObject JSingletDialog::staticMetaObject = {
    {   &QDialog::staticMetaObject, qt_meta_stringdata_JSingletDialog.data,
        qt_meta_data_JSingletDialog,  qt_static_metacall, Q_NULLPTR, Q_NULLPTR
    }
};


const QMetaObject *JSingletDialog::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->dynamicMetaObject() : &staticMetaObject;
}

void *JSingletDialog::qt_metacast(const char *_clname)
{
    if (!_clname) return Q_NULLPTR;
    if (!strcmp(_clname, qt_meta_stringdata_JSingletDialog.stringdata0))
        return static_cast<void*>(const_cast< JSingletDialog*>(this));
    if (!strcmp(_clname, "Ui::SingletDialog"))
        return static_cast< Ui::SingletDialog*>(const_cast< JSingletDialog*>(this));
    return QDialog::qt_metacast(_clname);
}

int JSingletDialog::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QDialog::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        if (_id < 3)
            qt_static_metacall(this, _c, _id, _a);
        _id -= 3;
    } else if (_c == QMetaObject::RegisterMethodArgumentMetaType) {
        if (_id < 3)
            *reinterpret_cast<int*>(_a[0]) = -1;
        _id -= 3;
    }
    return _id;
}
QT_END_MOC_NAMESPACE
