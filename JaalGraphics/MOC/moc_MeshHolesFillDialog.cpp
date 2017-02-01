/****************************************************************************
** Meta object code from reading C++ file 'MeshHolesFillDialog.hpp'
**
** Created by: The Qt Meta Object Compiler version 67 (Qt 5.6.0)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "../include/MeshHolesFillDialog.hpp"
#include <QtCore/qbytearray.h>
#include <QtCore/qmetatype.h>
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'MeshHolesFillDialog.hpp' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 67
#error "This file was generated using the moc from 5.6.0. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
struct qt_meta_stringdata_JMeshHolesFillDialog_t {
    QByteArrayData data[8];
    char stringdata0[84];
};
#define QT_MOC_LITERAL(idx, ofs, len) \
    Q_STATIC_BYTE_ARRAY_DATA_HEADER_INITIALIZER_WITH_OFFSET(len, \
    qptrdiff(offsetof(qt_meta_stringdata_JMeshHolesFillDialog_t, stringdata0) + ofs \
        - idx * sizeof(QByteArrayData)) \
    )
static const qt_meta_stringdata_JMeshHolesFillDialog_t qt_meta_stringdata_JMeshHolesFillDialog = {
    {
QT_MOC_LITERAL(0, 0, 20), // "JMeshHolesFillDialog"
QT_MOC_LITERAL(1, 21, 11), // "searchHoles"
QT_MOC_LITERAL(2, 33, 0), // ""
QT_MOC_LITERAL(3, 34, 7), // "fillOne"
QT_MOC_LITERAL(4, 42, 7), // "fillAll"
QT_MOC_LITERAL(5, 50, 10), // "refineHole"
QT_MOC_LITERAL(6, 61, 10), // "smoothHole"
QT_MOC_LITERAL(7, 72, 11) // "closeDialog"

    },
    "JMeshHolesFillDialog\0searchHoles\0\0"
    "fillOne\0fillAll\0refineHole\0smoothHole\0"
    "closeDialog"
};
#undef QT_MOC_LITERAL

static const uint qt_meta_data_JMeshHolesFillDialog[] = {

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
       1,    0,   44,    2, 0x08 /* Private */,
       3,    0,   45,    2, 0x08 /* Private */,
       4,    0,   46,    2, 0x08 /* Private */,
       5,    0,   47,    2, 0x08 /* Private */,
       6,    0,   48,    2, 0x08 /* Private */,
       7,    0,   49,    2, 0x08 /* Private */,

 // slots: parameters
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,

       0        // eod
};

void JMeshHolesFillDialog::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    if (_c == QMetaObject::InvokeMetaMethod) {
        JMeshHolesFillDialog *_t = static_cast<JMeshHolesFillDialog *>(_o);
        Q_UNUSED(_t)
        switch (_id) {
        case 0: _t->searchHoles(); break;
        case 1: _t->fillOne(); break;
        case 2: _t->fillAll(); break;
        case 3: _t->refineHole(); break;
        case 4: _t->smoothHole(); break;
        case 5: _t->closeDialog(); break;
        default: ;
        }
    }
    Q_UNUSED(_a);
}

const QMetaObject JMeshHolesFillDialog::staticMetaObject = {
    { &QDialog::staticMetaObject, qt_meta_stringdata_JMeshHolesFillDialog.data,
      qt_meta_data_JMeshHolesFillDialog,  qt_static_metacall, Q_NULLPTR, Q_NULLPTR}
};


const QMetaObject *JMeshHolesFillDialog::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->dynamicMetaObject() : &staticMetaObject;
}

void *JMeshHolesFillDialog::qt_metacast(const char *_clname)
{
    if (!_clname) return Q_NULLPTR;
    if (!strcmp(_clname, qt_meta_stringdata_JMeshHolesFillDialog.stringdata0))
        return static_cast<void*>(const_cast< JMeshHolesFillDialog*>(this));
    if (!strcmp(_clname, "Ui::MeshHolesFillDialog"))
        return static_cast< Ui::MeshHolesFillDialog*>(const_cast< JMeshHolesFillDialog*>(this));
    return QDialog::qt_metacast(_clname);
}

int JMeshHolesFillDialog::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QDialog::qt_metacall(_c, _id, _a);
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
