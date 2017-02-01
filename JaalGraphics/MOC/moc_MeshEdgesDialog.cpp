/****************************************************************************
** Meta object code from reading C++ file 'MeshEdgesDialog.hpp'
**
** Created by: The Qt Meta Object Compiler version 67 (Qt 5.6.0)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "MeshEdgesDialog.hpp"
#include <QtCore/qbytearray.h>
#include <QtCore/qmetatype.h>
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'MeshEdgesDialog.hpp' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 67
#error "This file was generated using the moc from 5.6.0. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
struct qt_meta_stringdata_JMeshEdgesDialog_t {
    QByteArrayData data[21];
    char stringdata0[270];
};
#define QT_MOC_LITERAL(idx, ofs, len) \
    Q_STATIC_BYTE_ARRAY_DATA_HEADER_INITIALIZER_WITH_OFFSET(len, \
    qptrdiff(offsetof(qt_meta_stringdata_JMeshEdgesDialog_t, stringdata0) + ofs \
        - idx * sizeof(QByteArrayData)) \
    )
static const qt_meta_stringdata_JMeshEdgesDialog_t qt_meta_stringdata_JMeshEdgesDialog = {
    {
QT_MOC_LITERAL(0, 0, 16), // "JMeshEdgesDialog"
QT_MOC_LITERAL(1, 17, 9), // "showEvent"
QT_MOC_LITERAL(2, 27, 0), // ""
QT_MOC_LITERAL(3, 28, 11), // "QShowEvent*"
QT_MOC_LITERAL(4, 40, 1), // "e"
QT_MOC_LITERAL(5, 42, 13), // "keyPressEvent"
QT_MOC_LITERAL(6, 56, 10), // "QKeyEvent*"
QT_MOC_LITERAL(7, 67, 15), // "countBoundEdges"
QT_MOC_LITERAL(8, 83, 14), // "getHiddenlines"
QT_MOC_LITERAL(9, 98, 11), // "setInternal"
QT_MOC_LITERAL(10, 110, 11), // "setBoundary"
QT_MOC_LITERAL(11, 122, 19), // "getNonManifoldEdges"
QT_MOC_LITERAL(12, 142, 11), // "closeDialog"
QT_MOC_LITERAL(13, 154, 13), // "setNumVisible"
QT_MOC_LITERAL(14, 168, 10), // "checkState"
QT_MOC_LITERAL(15, 179, 12), // "checkDisplay"
QT_MOC_LITERAL(16, 192, 6), // "lookAt"
QT_MOC_LITERAL(17, 199, 12), // "saveBoundary"
QT_MOC_LITERAL(18, 212, 21), // "setDefaultColorMethod"
QT_MOC_LITERAL(19, 234, 20), // "openAttribListDialog"
QT_MOC_LITERAL(20, 255, 14) // "openEditDialog"

    },
    "JMeshEdgesDialog\0showEvent\0\0QShowEvent*\0"
    "e\0keyPressEvent\0QKeyEvent*\0countBoundEdges\0"
    "getHiddenlines\0setInternal\0setBoundary\0"
    "getNonManifoldEdges\0closeDialog\0"
    "setNumVisible\0checkState\0checkDisplay\0"
    "lookAt\0saveBoundary\0setDefaultColorMethod\0"
    "openAttribListDialog\0openEditDialog"
};
#undef QT_MOC_LITERAL

static const uint qt_meta_data_JMeshEdgesDialog[] = {

 // content:
       7,       // revision
       0,       // classname
       0,    0, // classinfo
      16,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       0,       // signalCount

 // slots: name, argc, parameters, tag, flags
       1,    1,   94,    2, 0x08 /* Private */,
       5,    1,   97,    2, 0x08 /* Private */,
       7,    0,  100,    2, 0x08 /* Private */,
       8,    0,  101,    2, 0x08 /* Private */,
       9,    0,  102,    2, 0x08 /* Private */,
      10,    0,  103,    2, 0x08 /* Private */,
      11,    0,  104,    2, 0x08 /* Private */,
      12,    0,  105,    2, 0x08 /* Private */,
      13,    0,  106,    2, 0x08 /* Private */,
      14,    0,  107,    2, 0x08 /* Private */,
      15,    0,  108,    2, 0x08 /* Private */,
      16,    0,  109,    2, 0x08 /* Private */,
      17,    0,  110,    2, 0x08 /* Private */,
      18,    0,  111,    2, 0x08 /* Private */,
      19,    0,  112,    2, 0x08 /* Private */,
      20,    0,  113,    2, 0x08 /* Private */,

 // slots: parameters
    QMetaType::Void, 0x80000000 | 3,    4,
    QMetaType::Void, 0x80000000 | 6,    4,
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
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,

       0        // eod
};

void JMeshEdgesDialog::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    if (_c == QMetaObject::InvokeMetaMethod) {
        JMeshEdgesDialog *_t = static_cast<JMeshEdgesDialog *>(_o);
        Q_UNUSED(_t)
        switch (_id) {
        case 0: _t->showEvent((*reinterpret_cast< QShowEvent*(*)>(_a[1]))); break;
        case 1: _t->keyPressEvent((*reinterpret_cast< QKeyEvent*(*)>(_a[1]))); break;
        case 2: _t->countBoundEdges(); break;
        case 3: _t->getHiddenlines(); break;
        case 4: _t->setInternal(); break;
        case 5: _t->setBoundary(); break;
        case 6: _t->getNonManifoldEdges(); break;
        case 7: _t->closeDialog(); break;
        case 8: _t->setNumVisible(); break;
        case 9: _t->checkState(); break;
        case 10: _t->checkDisplay(); break;
        case 11: _t->lookAt(); break;
        case 12: _t->saveBoundary(); break;
        case 13: _t->setDefaultColorMethod(); break;
        case 14: _t->openAttribListDialog(); break;
        case 15: _t->openEditDialog(); break;
        default: ;
        }
    }
}

const QMetaObject JMeshEdgesDialog::staticMetaObject = {
    { &QDialog::staticMetaObject, qt_meta_stringdata_JMeshEdgesDialog.data,
      qt_meta_data_JMeshEdgesDialog,  qt_static_metacall, Q_NULLPTR, Q_NULLPTR}
};


const QMetaObject *JMeshEdgesDialog::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->dynamicMetaObject() : &staticMetaObject;
}

void *JMeshEdgesDialog::qt_metacast(const char *_clname)
{
    if (!_clname) return Q_NULLPTR;
    if (!strcmp(_clname, qt_meta_stringdata_JMeshEdgesDialog.stringdata0))
        return static_cast<void*>(const_cast< JMeshEdgesDialog*>(this));
    if (!strcmp(_clname, "Ui::MeshEdgesDialog"))
        return static_cast< Ui::MeshEdgesDialog*>(const_cast< JMeshEdgesDialog*>(this));
    return QDialog::qt_metacast(_clname);
}

int JMeshEdgesDialog::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QDialog::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        if (_id < 16)
            qt_static_metacall(this, _c, _id, _a);
        _id -= 16;
    } else if (_c == QMetaObject::RegisterMethodArgumentMetaType) {
        if (_id < 16)
            *reinterpret_cast<int*>(_a[0]) = -1;
        _id -= 16;
    }
    return _id;
}
QT_END_MOC_NAMESPACE
