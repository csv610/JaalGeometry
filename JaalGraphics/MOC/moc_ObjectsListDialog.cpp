/****************************************************************************
** Meta object code from reading C++ file 'ObjectsListDialog.hpp'
**
** Created by: The Qt Meta Object Compiler version 67 (Qt 5.5.1)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "ObjectsListDialog.hpp"
#include <QtCore/qbytearray.h>
#include <QtCore/qmetatype.h>
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'ObjectsListDialog.hpp' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 67
#error "This file was generated using the moc from 5.5.1. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
struct qt_meta_stringdata_JObjectsListDialog_t {
    QByteArrayData data[9];
    char stringdata0[81];
};
#define QT_MOC_LITERAL(idx, ofs, len) \
    Q_STATIC_BYTE_ARRAY_DATA_HEADER_INITIALIZER_WITH_OFFSET(len, \
    qptrdiff(offsetof(qt_meta_stringdata_JObjectsListDialog_t, stringdata0) + ofs \
        - idx * sizeof(QByteArrayData)) \
    )
static const qt_meta_stringdata_JObjectsListDialog_t qt_meta_stringdata_JObjectsListDialog = {
    {
        QT_MOC_LITERAL(0, 0, 18), // "JObjectsListDialog"
        QT_MOC_LITERAL(1, 19, 7), // "getCell"
        QT_MOC_LITERAL(2, 27, 0), // ""
        QT_MOC_LITERAL(3, 28, 1), // "i"
        QT_MOC_LITERAL(4, 30, 1), // "j"
        QT_MOC_LITERAL(5, 32, 12), // "deleteObject"
        QT_MOC_LITERAL(6, 45, 10), // "saveObject"
        QT_MOC_LITERAL(7, 56, 11), // "closeDialog"
        QT_MOC_LITERAL(8, 68, 12) // "renameObject"

    },
    "JObjectsListDialog\0getCell\0\0i\0j\0"
    "deleteObject\0saveObject\0closeDialog\0"
    "renameObject"
};
#undef QT_MOC_LITERAL

static const uint qt_meta_data_JObjectsListDialog[] = {

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
    1,    2,   39,    2, 0x08 /* Private */,
    5,    0,   44,    2, 0x08 /* Private */,
    6,    0,   45,    2, 0x08 /* Private */,
    7,    0,   46,    2, 0x08 /* Private */,
    8,    0,   47,    2, 0x08 /* Private */,

// slots: parameters
    QMetaType::Void, QMetaType::Int, QMetaType::Int,    3,    4,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,

    0        // eod
};

void JObjectsListDialog::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    if (_c == QMetaObject::InvokeMetaMethod) {
        JObjectsListDialog *_t = static_cast<JObjectsListDialog *>(_o);
        Q_UNUSED(_t)
        switch (_id) {
        case 0:
            _t->getCell((*reinterpret_cast< int(*)>(_a[1])),(*reinterpret_cast< int(*)>(_a[2])));
            break;
        case 1:
            _t->deleteObject();
            break;
        case 2:
            _t->saveObject();
            break;
        case 3:
            _t->closeDialog();
            break;
        case 4:
            _t->renameObject();
            break;
        default:
            ;
        }
    }
}

const QMetaObject JObjectsListDialog::staticMetaObject = {
    {   &QDialog::staticMetaObject, qt_meta_stringdata_JObjectsListDialog.data,
        qt_meta_data_JObjectsListDialog,  qt_static_metacall, Q_NULLPTR, Q_NULLPTR
    }
};


const QMetaObject *JObjectsListDialog::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->dynamicMetaObject() : &staticMetaObject;
}

void *JObjectsListDialog::qt_metacast(const char *_clname)
{
    if (!_clname) return Q_NULLPTR;
    if (!strcmp(_clname, qt_meta_stringdata_JObjectsListDialog.stringdata0))
        return static_cast<void*>(const_cast< JObjectsListDialog*>(this));
    if (!strcmp(_clname, "Ui::ObjectsListDialog"))
        return static_cast< Ui::ObjectsListDialog*>(const_cast< JObjectsListDialog*>(this));
    return QDialog::qt_metacast(_clname);
}

int JObjectsListDialog::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
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
