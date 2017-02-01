/****************************************************************************
** Meta object code from reading C++ file 'MeshComponentsDialog.hpp'
**
** Created by: The Qt Meta Object Compiler version 67 (Qt 5.6.0)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "../include/MeshComponentsDialog.hpp"
#include <QtCore/qbytearray.h>
#include <QtCore/qmetatype.h>
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'MeshComponentsDialog.hpp' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 67
#error "This file was generated using the moc from 5.6.0. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
struct qt_meta_stringdata_JMeshComponentsDialog_t {
    QByteArrayData data[9];
    char stringdata0[128];
};
#define QT_MOC_LITERAL(idx, ofs, len) \
    Q_STATIC_BYTE_ARRAY_DATA_HEADER_INITIALIZER_WITH_OFFSET(len, \
    qptrdiff(offsetof(qt_meta_stringdata_JMeshComponentsDialog_t, stringdata0) + ofs \
        - idx * sizeof(QByteArrayData)) \
    )
static const qt_meta_stringdata_JMeshComponentsDialog_t qt_meta_stringdata_JMeshComponentsDialog = {
    {
QT_MOC_LITERAL(0, 0, 21), // "JMeshComponentsDialog"
QT_MOC_LITERAL(1, 22, 16), // "searchComponents"
QT_MOC_LITERAL(2, 39, 0), // ""
QT_MOC_LITERAL(3, 40, 15), // "removeComponent"
QT_MOC_LITERAL(4, 56, 15), // "mergeComponents"
QT_MOC_LITERAL(5, 72, 17), // "select2Components"
QT_MOC_LITERAL(6, 90, 8), // "mergeAll"
QT_MOC_LITERAL(7, 99, 16), // "displayComponent"
QT_MOC_LITERAL(8, 116, 11) // "closeDialog"

    },
    "JMeshComponentsDialog\0searchComponents\0"
    "\0removeComponent\0mergeComponents\0"
    "select2Components\0mergeAll\0displayComponent\0"
    "closeDialog"
};
#undef QT_MOC_LITERAL

static const uint qt_meta_data_JMeshComponentsDialog[] = {

 // content:
       7,       // revision
       0,       // classname
       0,    0, // classinfo
       7,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       0,       // signalCount

 // slots: name, argc, parameters, tag, flags
       1,    0,   49,    2, 0x08 /* Private */,
       3,    0,   50,    2, 0x08 /* Private */,
       4,    0,   51,    2, 0x08 /* Private */,
       5,    0,   52,    2, 0x08 /* Private */,
       6,    0,   53,    2, 0x08 /* Private */,
       7,    0,   54,    2, 0x08 /* Private */,
       8,    0,   55,    2, 0x08 /* Private */,

 // slots: parameters
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,

       0        // eod
};

void JMeshComponentsDialog::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    if (_c == QMetaObject::InvokeMetaMethod) {
        JMeshComponentsDialog *_t = static_cast<JMeshComponentsDialog *>(_o);
        Q_UNUSED(_t)
        switch (_id) {
        case 0: _t->searchComponents(); break;
        case 1: _t->removeComponent(); break;
        case 2: _t->mergeComponents(); break;
        case 3: _t->select2Components(); break;
        case 4: _t->mergeAll(); break;
        case 5: _t->displayComponent(); break;
        case 6: _t->closeDialog(); break;
        default: ;
        }
    }
    Q_UNUSED(_a);
}

const QMetaObject JMeshComponentsDialog::staticMetaObject = {
    { &QDialog::staticMetaObject, qt_meta_stringdata_JMeshComponentsDialog.data,
      qt_meta_data_JMeshComponentsDialog,  qt_static_metacall, Q_NULLPTR, Q_NULLPTR}
};


const QMetaObject *JMeshComponentsDialog::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->dynamicMetaObject() : &staticMetaObject;
}

void *JMeshComponentsDialog::qt_metacast(const char *_clname)
{
    if (!_clname) return Q_NULLPTR;
    if (!strcmp(_clname, qt_meta_stringdata_JMeshComponentsDialog.stringdata0))
        return static_cast<void*>(const_cast< JMeshComponentsDialog*>(this));
    if (!strcmp(_clname, "Ui::MeshComponentsDialog"))
        return static_cast< Ui::MeshComponentsDialog*>(const_cast< JMeshComponentsDialog*>(this));
    return QDialog::qt_metacast(_clname);
}

int JMeshComponentsDialog::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QDialog::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        if (_id < 7)
            qt_static_metacall(this, _c, _id, _a);
        _id -= 7;
    } else if (_c == QMetaObject::RegisterMethodArgumentMetaType) {
        if (_id < 7)
            *reinterpret_cast<int*>(_a[0]) = -1;
        _id -= 7;
    }
    return _id;
}
QT_END_MOC_NAMESPACE
