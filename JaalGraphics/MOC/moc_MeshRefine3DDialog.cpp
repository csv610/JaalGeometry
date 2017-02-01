/****************************************************************************
** Meta object code from reading C++ file 'MeshRefine3DDialog.hpp'
**
** Created by: The Qt Meta Object Compiler version 67 (Qt 5.6.0)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "../include/MeshRefine3DDialog.hpp"
#include <QtCore/qbytearray.h>
#include <QtCore/qmetatype.h>
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'MeshRefine3DDialog.hpp' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 67
#error "This file was generated using the moc from 5.6.0. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
struct qt_meta_stringdata_JMeshRefine3DDialog_t {
    QByteArrayData data[6];
    char stringdata0[72];
};
#define QT_MOC_LITERAL(idx, ofs, len) \
    Q_STATIC_BYTE_ARRAY_DATA_HEADER_INITIALIZER_WITH_OFFSET(len, \
    qptrdiff(offsetof(qt_meta_stringdata_JMeshRefine3DDialog_t, stringdata0) + ofs \
        - idx * sizeof(QByteArrayData)) \
    )
static const qt_meta_stringdata_JMeshRefine3DDialog_t qt_meta_stringdata_JMeshRefine3DDialog = {
    {
QT_MOC_LITERAL(0, 0, 19), // "JMeshRefine3DDialog"
QT_MOC_LITERAL(1, 20, 13), // "insertPillows"
QT_MOC_LITERAL(2, 34, 0), // ""
QT_MOC_LITERAL(3, 35, 11), // "refineHex17"
QT_MOC_LITERAL(4, 47, 11), // "refineHex18"
QT_MOC_LITERAL(5, 59, 12) // "genHexBlocks"

    },
    "JMeshRefine3DDialog\0insertPillows\0\0"
    "refineHex17\0refineHex18\0genHexBlocks"
};
#undef QT_MOC_LITERAL

static const uint qt_meta_data_JMeshRefine3DDialog[] = {

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
       1,    0,   34,    2, 0x08 /* Private */,
       3,    0,   35,    2, 0x08 /* Private */,
       4,    0,   36,    2, 0x08 /* Private */,
       5,    0,   37,    2, 0x08 /* Private */,

 // slots: parameters
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,

       0        // eod
};

void JMeshRefine3DDialog::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    if (_c == QMetaObject::InvokeMetaMethod) {
        JMeshRefine3DDialog *_t = static_cast<JMeshRefine3DDialog *>(_o);
        Q_UNUSED(_t)
        switch (_id) {
        case 0: _t->insertPillows(); break;
        case 1: _t->refineHex17(); break;
        case 2: _t->refineHex18(); break;
        case 3: _t->genHexBlocks(); break;
        default: ;
        }
    }
    Q_UNUSED(_a);
}

const QMetaObject JMeshRefine3DDialog::staticMetaObject = {
    { &QDialog::staticMetaObject, qt_meta_stringdata_JMeshRefine3DDialog.data,
      qt_meta_data_JMeshRefine3DDialog,  qt_static_metacall, Q_NULLPTR, Q_NULLPTR}
};


const QMetaObject *JMeshRefine3DDialog::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->dynamicMetaObject() : &staticMetaObject;
}

void *JMeshRefine3DDialog::qt_metacast(const char *_clname)
{
    if (!_clname) return Q_NULLPTR;
    if (!strcmp(_clname, qt_meta_stringdata_JMeshRefine3DDialog.stringdata0))
        return static_cast<void*>(const_cast< JMeshRefine3DDialog*>(this));
    if (!strcmp(_clname, "Ui::MeshRefine3DDialog"))
        return static_cast< Ui::MeshRefine3DDialog*>(const_cast< JMeshRefine3DDialog*>(this));
    return QDialog::qt_metacast(_clname);
}

int JMeshRefine3DDialog::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
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
