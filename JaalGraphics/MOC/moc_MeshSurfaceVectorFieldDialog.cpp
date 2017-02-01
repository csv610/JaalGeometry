/****************************************************************************
** Meta object code from reading C++ file 'MeshSurfaceVectorFieldDialog.hpp'
**
** Created by: The Qt Meta Object Compiler version 67 (Qt 5.6.0)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "../include/MeshSurfaceVectorFieldDialog.hpp"
#include <QtCore/qbytearray.h>
#include <QtCore/qmetatype.h>
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'MeshSurfaceVectorFieldDialog.hpp' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 67
#error "This file was generated using the moc from 5.6.0. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
struct qt_meta_stringdata_JMeshSurfaceVectorFieldDialog_t {
    QByteArrayData data[12];
    char stringdata0[187];
};
#define QT_MOC_LITERAL(idx, ofs, len) \
    Q_STATIC_BYTE_ARRAY_DATA_HEADER_INITIALIZER_WITH_OFFSET(len, \
    qptrdiff(offsetof(qt_meta_stringdata_JMeshSurfaceVectorFieldDialog_t, stringdata0) + ofs \
        - idx * sizeof(QByteArrayData)) \
    )
static const qt_meta_stringdata_JMeshSurfaceVectorFieldDialog_t qt_meta_stringdata_JMeshSurfaceVectorFieldDialog = {
    {
QT_MOC_LITERAL(0, 0, 29), // "JMeshSurfaceVectorFieldDialog"
QT_MOC_LITERAL(1, 30, 11), // "getVecField"
QT_MOC_LITERAL(2, 42, 0), // ""
QT_MOC_LITERAL(3, 43, 10), // "clearField"
QT_MOC_LITERAL(4, 54, 12), // "setVecLength"
QT_MOC_LITERAL(5, 67, 21), // "openEdgeAttribsDialog"
QT_MOC_LITERAL(6, 89, 16), // "readFixedFacesID"
QT_MOC_LITERAL(7, 106, 21), // "readFixedFacesVectors"
QT_MOC_LITERAL(8, 128, 14), // "setConstraints"
QT_MOC_LITERAL(9, 143, 11), // "closeDialog"
QT_MOC_LITERAL(10, 155, 14), // "showAllVectors"
QT_MOC_LITERAL(11, 170, 16) // "showVecComponent"

    },
    "JMeshSurfaceVectorFieldDialog\0getVecField\0"
    "\0clearField\0setVecLength\0openEdgeAttribsDialog\0"
    "readFixedFacesID\0readFixedFacesVectors\0"
    "setConstraints\0closeDialog\0showAllVectors\0"
    "showVecComponent"
};
#undef QT_MOC_LITERAL

static const uint qt_meta_data_JMeshSurfaceVectorFieldDialog[] = {

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
       8,    0,   70,    2, 0x08 /* Private */,
       9,    0,   71,    2, 0x08 /* Private */,
      10,    0,   72,    2, 0x08 /* Private */,
      11,    0,   73,    2, 0x08 /* Private */,

 // slots: parameters
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,

       0        // eod
};

void JMeshSurfaceVectorFieldDialog::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    if (_c == QMetaObject::InvokeMetaMethod) {
        JMeshSurfaceVectorFieldDialog *_t = static_cast<JMeshSurfaceVectorFieldDialog *>(_o);
        Q_UNUSED(_t)
        switch (_id) {
        case 0: _t->getVecField(); break;
        case 1: _t->clearField(); break;
        case 2: _t->setVecLength(); break;
        case 3: _t->openEdgeAttribsDialog(); break;
        case 4: _t->readFixedFacesID(); break;
        case 5: _t->readFixedFacesVectors(); break;
        case 6: _t->setConstraints(); break;
        case 7: _t->closeDialog(); break;
        case 8: _t->showAllVectors(); break;
        case 9: _t->showVecComponent(); break;
        default: ;
        }
    }
    Q_UNUSED(_a);
}

const QMetaObject JMeshSurfaceVectorFieldDialog::staticMetaObject = {
    { &QDialog::staticMetaObject, qt_meta_stringdata_JMeshSurfaceVectorFieldDialog.data,
      qt_meta_data_JMeshSurfaceVectorFieldDialog,  qt_static_metacall, Q_NULLPTR, Q_NULLPTR}
};


const QMetaObject *JMeshSurfaceVectorFieldDialog::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->dynamicMetaObject() : &staticMetaObject;
}

void *JMeshSurfaceVectorFieldDialog::qt_metacast(const char *_clname)
{
    if (!_clname) return Q_NULLPTR;
    if (!strcmp(_clname, qt_meta_stringdata_JMeshSurfaceVectorFieldDialog.stringdata0))
        return static_cast<void*>(const_cast< JMeshSurfaceVectorFieldDialog*>(this));
    if (!strcmp(_clname, "Ui::MeshSurfaceVectorFieldDialog"))
        return static_cast< Ui::MeshSurfaceVectorFieldDialog*>(const_cast< JMeshSurfaceVectorFieldDialog*>(this));
    return QDialog::qt_metacast(_clname);
}

int JMeshSurfaceVectorFieldDialog::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
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
