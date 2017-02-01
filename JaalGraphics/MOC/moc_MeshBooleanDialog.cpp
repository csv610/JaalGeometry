/****************************************************************************
** Meta object code from reading C++ file 'MeshBooleanDialog.hpp'
**
** Created by: The Qt Meta Object Compiler version 67 (Qt 5.6.0)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "../include/MeshBooleanDialog.hpp"
#include <QtCore/qbytearray.h>
#include <QtCore/qmetatype.h>
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'MeshBooleanDialog.hpp' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 67
#error "This file was generated using the moc from 5.6.0. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
struct qt_meta_stringdata_JMeshBooleanDialog_t {
    QByteArrayData data[9];
    char stringdata0[84];
};
#define QT_MOC_LITERAL(idx, ofs, len) \
    Q_STATIC_BYTE_ARRAY_DATA_HEADER_INITIALIZER_WITH_OFFSET(len, \
    qptrdiff(offsetof(qt_meta_stringdata_JMeshBooleanDialog_t, stringdata0) + ofs \
        - idx * sizeof(QByteArrayData)) \
    )
static const qt_meta_stringdata_JMeshBooleanDialog_t qt_meta_stringdata_JMeshBooleanDialog = {
    {
QT_MOC_LITERAL(0, 0, 18), // "JMeshBooleanDialog"
QT_MOC_LITERAL(1, 19, 9), // "showEvent"
QT_MOC_LITERAL(2, 29, 0), // ""
QT_MOC_LITERAL(3, 30, 11), // "QShowEvent*"
QT_MOC_LITERAL(4, 42, 1), // "e"
QT_MOC_LITERAL(5, 44, 7), // "applyOp"
QT_MOC_LITERAL(6, 52, 9), // "loadMeshA"
QT_MOC_LITERAL(7, 62, 9), // "loadMeshB"
QT_MOC_LITERAL(8, 72, 11) // "closeDialog"

    },
    "JMeshBooleanDialog\0showEvent\0\0QShowEvent*\0"
    "e\0applyOp\0loadMeshA\0loadMeshB\0closeDialog"
};
#undef QT_MOC_LITERAL

static const uint qt_meta_data_JMeshBooleanDialog[] = {

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
       1,    1,   39,    2, 0x08 /* Private */,
       5,    0,   42,    2, 0x08 /* Private */,
       6,    0,   43,    2, 0x08 /* Private */,
       7,    0,   44,    2, 0x08 /* Private */,
       8,    0,   45,    2, 0x08 /* Private */,

 // slots: parameters
    QMetaType::Void, 0x80000000 | 3,    4,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,

       0        // eod
};

void JMeshBooleanDialog::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    if (_c == QMetaObject::InvokeMetaMethod) {
        JMeshBooleanDialog *_t = static_cast<JMeshBooleanDialog *>(_o);
        Q_UNUSED(_t)
        switch (_id) {
        case 0: _t->showEvent((*reinterpret_cast< QShowEvent*(*)>(_a[1]))); break;
        case 1: _t->applyOp(); break;
        case 2: _t->loadMeshA(); break;
        case 3: _t->loadMeshB(); break;
        case 4: _t->closeDialog(); break;
        default: ;
        }
    }
}

const QMetaObject JMeshBooleanDialog::staticMetaObject = {
    { &QDialog::staticMetaObject, qt_meta_stringdata_JMeshBooleanDialog.data,
      qt_meta_data_JMeshBooleanDialog,  qt_static_metacall, Q_NULLPTR, Q_NULLPTR}
};


const QMetaObject *JMeshBooleanDialog::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->dynamicMetaObject() : &staticMetaObject;
}

void *JMeshBooleanDialog::qt_metacast(const char *_clname)
{
    if (!_clname) return Q_NULLPTR;
    if (!strcmp(_clname, qt_meta_stringdata_JMeshBooleanDialog.stringdata0))
        return static_cast<void*>(const_cast< JMeshBooleanDialog*>(this));
    if (!strcmp(_clname, "Ui::MeshBooleanDialog"))
        return static_cast< Ui::MeshBooleanDialog*>(const_cast< JMeshBooleanDialog*>(this));
    return QDialog::qt_metacast(_clname);
}

int JMeshBooleanDialog::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
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
