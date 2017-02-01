/****************************************************************************
** Meta object code from reading C++ file 'MeshEdgeEditDialog.hpp'
**
** Created by: The Qt Meta Object Compiler version 67 (Qt 5.5.1)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "MeshEdgeEditDialog.hpp"
#include <QtCore/qbytearray.h>
#include <QtCore/qmetatype.h>
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'MeshEdgeEditDialog.hpp' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 67
#error "This file was generated using the moc from 5.5.1. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
struct qt_meta_stringdata_JMeshEdgeEditDialog_t {
    QByteArrayData data[9];
    char stringdata0[98];
};
#define QT_MOC_LITERAL(idx, ofs, len) \
    Q_STATIC_BYTE_ARRAY_DATA_HEADER_INITIALIZER_WITH_OFFSET(len, \
    qptrdiff(offsetof(qt_meta_stringdata_JMeshEdgeEditDialog_t, stringdata0) + ofs \
        - idx * sizeof(QByteArrayData)) \
    )
static const qt_meta_stringdata_JMeshEdgeEditDialog_t qt_meta_stringdata_JMeshEdgeEditDialog = {
    {
        QT_MOC_LITERAL(0, 0, 19), // "JMeshEdgeEditDialog"
        QT_MOC_LITERAL(1, 20, 11), // "getMinEdges"
        QT_MOC_LITERAL(2, 32, 0), // ""
        QT_MOC_LITERAL(3, 33, 11), // "getMaxEdges"
        QT_MOC_LITERAL(4, 45, 6), // "refine"
        QT_MOC_LITERAL(5, 52, 8), // "collapse"
        QT_MOC_LITERAL(6, 61, 9), // "deleteAll"
        QT_MOC_LITERAL(7, 71, 14), // "deleteInternal"
        QT_MOC_LITERAL(8, 86, 11) // "closeDialog"

    },
    "JMeshEdgeEditDialog\0getMinEdges\0\0"
    "getMaxEdges\0refine\0collapse\0deleteAll\0"
    "deleteInternal\0closeDialog"
};
#undef QT_MOC_LITERAL

static const uint qt_meta_data_JMeshEdgeEditDialog[] = {

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

void JMeshEdgeEditDialog::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    if (_c == QMetaObject::InvokeMetaMethod) {
        JMeshEdgeEditDialog *_t = static_cast<JMeshEdgeEditDialog *>(_o);
        Q_UNUSED(_t)
        switch (_id) {
        case 0:
            _t->getMinEdges();
            break;
        case 1:
            _t->getMaxEdges();
            break;
        case 2:
            _t->refine();
            break;
        case 3:
            _t->collapse();
            break;
        case 4:
            _t->deleteAll();
            break;
        case 5:
            _t->deleteInternal();
            break;
        case 6:
            _t->closeDialog();
            break;
        default:
            ;
        }
    }
    Q_UNUSED(_a);
}

const QMetaObject JMeshEdgeEditDialog::staticMetaObject = {
    {   &QDialog::staticMetaObject, qt_meta_stringdata_JMeshEdgeEditDialog.data,
        qt_meta_data_JMeshEdgeEditDialog,  qt_static_metacall, Q_NULLPTR, Q_NULLPTR
    }
};


const QMetaObject *JMeshEdgeEditDialog::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->dynamicMetaObject() : &staticMetaObject;
}

void *JMeshEdgeEditDialog::qt_metacast(const char *_clname)
{
    if (!_clname) return Q_NULLPTR;
    if (!strcmp(_clname, qt_meta_stringdata_JMeshEdgeEditDialog.stringdata0))
        return static_cast<void*>(const_cast< JMeshEdgeEditDialog*>(this));
    if (!strcmp(_clname, "Ui::MeshEdgeEditDialog"))
        return static_cast< Ui::MeshEdgeEditDialog*>(const_cast< JMeshEdgeEditDialog*>(this));
    return QDialog::qt_metacast(_clname);
}

int JMeshEdgeEditDialog::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
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
