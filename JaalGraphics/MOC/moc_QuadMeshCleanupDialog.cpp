/****************************************************************************
** Meta object code from reading C++ file 'QuadMeshCleanupDialog.hpp'
**
** Created by: The Qt Meta Object Compiler version 67 (Qt 5.6.0)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "../include/QuadMeshCleanupDialog.hpp"
#include <QtCore/qbytearray.h>
#include <QtCore/qmetatype.h>
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'QuadMeshCleanupDialog.hpp' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 67
#error "This file was generated using the moc from 5.6.0. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
struct qt_meta_stringdata_JQuadMeshCleanupDialog_t {
    QByteArrayData data[11];
    char stringdata0[151];
};
#define QT_MOC_LITERAL(idx, ofs, len) \
    Q_STATIC_BYTE_ARRAY_DATA_HEADER_INITIALIZER_WITH_OFFSET(len, \
    qptrdiff(offsetof(qt_meta_stringdata_JQuadMeshCleanupDialog_t, stringdata0) + ofs \
        - idx * sizeof(QByteArrayData)) \
    )
static const qt_meta_stringdata_JQuadMeshCleanupDialog_t qt_meta_stringdata_JQuadMeshCleanupDialog = {
    {
QT_MOC_LITERAL(0, 0, 22), // "JQuadMeshCleanupDialog"
QT_MOC_LITERAL(1, 23, 14), // "searchSinglets"
QT_MOC_LITERAL(2, 38, 0), // ""
QT_MOC_LITERAL(3, 39, 14), // "searchDoublets"
QT_MOC_LITERAL(4, 54, 14), // "searchDiamonds"
QT_MOC_LITERAL(5, 69, 14), // "removeSinglets"
QT_MOC_LITERAL(6, 84, 14), // "removeDoublets"
QT_MOC_LITERAL(7, 99, 14), // "removeDiamonds"
QT_MOC_LITERAL(8, 114, 14), // "openDualDialog"
QT_MOC_LITERAL(9, 129, 9), // "swapEdges"
QT_MOC_LITERAL(10, 139, 11) // "closeDialog"

    },
    "JQuadMeshCleanupDialog\0searchSinglets\0"
    "\0searchDoublets\0searchDiamonds\0"
    "removeSinglets\0removeDoublets\0"
    "removeDiamonds\0openDualDialog\0swapEdges\0"
    "closeDialog"
};
#undef QT_MOC_LITERAL

static const uint qt_meta_data_JQuadMeshCleanupDialog[] = {

 // content:
       7,       // revision
       0,       // classname
       0,    0, // classinfo
       9,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       0,       // signalCount

 // slots: name, argc, parameters, tag, flags
       1,    0,   59,    2, 0x08 /* Private */,
       3,    0,   60,    2, 0x08 /* Private */,
       4,    0,   61,    2, 0x08 /* Private */,
       5,    0,   62,    2, 0x08 /* Private */,
       6,    0,   63,    2, 0x08 /* Private */,
       7,    0,   64,    2, 0x08 /* Private */,
       8,    0,   65,    2, 0x08 /* Private */,
       9,    0,   66,    2, 0x08 /* Private */,
      10,    0,   67,    2, 0x08 /* Private */,

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

       0        // eod
};

void JQuadMeshCleanupDialog::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    if (_c == QMetaObject::InvokeMetaMethod) {
        JQuadMeshCleanupDialog *_t = static_cast<JQuadMeshCleanupDialog *>(_o);
        Q_UNUSED(_t)
        switch (_id) {
        case 0: _t->searchSinglets(); break;
        case 1: _t->searchDoublets(); break;
        case 2: _t->searchDiamonds(); break;
        case 3: _t->removeSinglets(); break;
        case 4: _t->removeDoublets(); break;
        case 5: _t->removeDiamonds(); break;
        case 6: _t->openDualDialog(); break;
        case 7: _t->swapEdges(); break;
        case 8: _t->closeDialog(); break;
        default: ;
        }
    }
    Q_UNUSED(_a);
}

const QMetaObject JQuadMeshCleanupDialog::staticMetaObject = {
    { &QDialog::staticMetaObject, qt_meta_stringdata_JQuadMeshCleanupDialog.data,
      qt_meta_data_JQuadMeshCleanupDialog,  qt_static_metacall, Q_NULLPTR, Q_NULLPTR}
};


const QMetaObject *JQuadMeshCleanupDialog::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->dynamicMetaObject() : &staticMetaObject;
}

void *JQuadMeshCleanupDialog::qt_metacast(const char *_clname)
{
    if (!_clname) return Q_NULLPTR;
    if (!strcmp(_clname, qt_meta_stringdata_JQuadMeshCleanupDialog.stringdata0))
        return static_cast<void*>(const_cast< JQuadMeshCleanupDialog*>(this));
    if (!strcmp(_clname, "Ui::QuadMeshCleanupDialog"))
        return static_cast< Ui::QuadMeshCleanupDialog*>(const_cast< JQuadMeshCleanupDialog*>(this));
    return QDialog::qt_metacast(_clname);
}

int JQuadMeshCleanupDialog::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QDialog::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        if (_id < 9)
            qt_static_metacall(this, _c, _id, _a);
        _id -= 9;
    } else if (_c == QMetaObject::RegisterMethodArgumentMetaType) {
        if (_id < 9)
            *reinterpret_cast<int*>(_a[0]) = -1;
        _id -= 9;
    }
    return _id;
}
QT_END_MOC_NAMESPACE
