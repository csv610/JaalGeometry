/****************************************************************************
** Meta object code from reading C++ file 'MeshDualEdgesDialog.hpp'
**
** Created by: The Qt Meta Object Compiler version 67 (Qt 5.5.1)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "MeshDualEdgesDialog.hpp"
#include <QtCore/qbytearray.h>
#include <QtCore/qmetatype.h>
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'MeshDualEdgesDialog.hpp' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 67
#error "This file was generated using the moc from 5.5.1. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
struct qt_meta_stringdata_JMeshDualEdgesDialog_t {
    QByteArrayData data[8];
    char stringdata0[89];
};
#define QT_MOC_LITERAL(idx, ofs, len) \
    Q_STATIC_BYTE_ARRAY_DATA_HEADER_INITIALIZER_WITH_OFFSET(len, \
    qptrdiff(offsetof(qt_meta_stringdata_JMeshDualEdgesDialog_t, stringdata0) + ofs \
        - idx * sizeof(QByteArrayData)) \
    )
static const qt_meta_stringdata_JMeshDualEdgesDialog_t qt_meta_stringdata_JMeshDualEdgesDialog = {
    {
        QT_MOC_LITERAL(0, 0, 20), // "JMeshDualEdgesDialog"
        QT_MOC_LITERAL(1, 21, 13), // "keyPressEvent"
        QT_MOC_LITERAL(2, 35, 0), // ""
        QT_MOC_LITERAL(3, 36, 10), // "QKeyEvent*"
        QT_MOC_LITERAL(4, 47, 1), // "e"
        QT_MOC_LITERAL(5, 49, 16), // "openAttribDialog"
        QT_MOC_LITERAL(6, 66, 10), // "checkEdges"
        QT_MOC_LITERAL(7, 77, 11) // "deleteEdges"

    },
    "JMeshDualEdgesDialog\0keyPressEvent\0\0"
    "QKeyEvent*\0e\0openAttribDialog\0checkEdges\0"
    "deleteEdges"
};
#undef QT_MOC_LITERAL

static const uint qt_meta_data_JMeshDualEdgesDialog[] = {

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

void JMeshDualEdgesDialog::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    if (_c == QMetaObject::InvokeMetaMethod) {
        JMeshDualEdgesDialog *_t = static_cast<JMeshDualEdgesDialog *>(_o);
        Q_UNUSED(_t)
        switch (_id) {
        case 0:
            _t->keyPressEvent((*reinterpret_cast< QKeyEvent*(*)>(_a[1])));
            break;
        case 1:
            _t->openAttribDialog();
            break;
        case 2:
            _t->checkEdges();
            break;
        case 3:
            _t->deleteEdges();
            break;
        default:
            ;
        }
    }
}

const QMetaObject JMeshDualEdgesDialog::staticMetaObject = {
    {   &QDialog::staticMetaObject, qt_meta_stringdata_JMeshDualEdgesDialog.data,
        qt_meta_data_JMeshDualEdgesDialog,  qt_static_metacall, Q_NULLPTR, Q_NULLPTR
    }
};


const QMetaObject *JMeshDualEdgesDialog::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->dynamicMetaObject() : &staticMetaObject;
}

void *JMeshDualEdgesDialog::qt_metacast(const char *_clname)
{
    if (!_clname) return Q_NULLPTR;
    if (!strcmp(_clname, qt_meta_stringdata_JMeshDualEdgesDialog.stringdata0))
        return static_cast<void*>(const_cast< JMeshDualEdgesDialog*>(this));
    if (!strcmp(_clname, "Ui::MeshDualEdgesDialog"))
        return static_cast< Ui::MeshDualEdgesDialog*>(const_cast< JMeshDualEdgesDialog*>(this));
    return QDialog::qt_metacast(_clname);
}

int JMeshDualEdgesDialog::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
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
