/****************************************************************************
** Meta object code from reading C++ file 'MeshDualsDialog.hpp'
**
** Created by: The Qt Meta Object Compiler version 67 (Qt 5.6.0)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "../include/MeshDualsDialog.hpp"
#include <QtCore/qbytearray.h>
#include <QtCore/qmetatype.h>
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'MeshDualsDialog.hpp' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 67
#error "This file was generated using the moc from 5.6.0. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
struct qt_meta_stringdata_JMeshDualsDialog_t {
    QByteArrayData data[18];
    char stringdata0[198];
};
#define QT_MOC_LITERAL(idx, ofs, len) \
    Q_STATIC_BYTE_ARRAY_DATA_HEADER_INITIALIZER_WITH_OFFSET(len, \
    qptrdiff(offsetof(qt_meta_stringdata_JMeshDualsDialog_t, stringdata0) + ofs \
        - idx * sizeof(QByteArrayData)) \
    )
static const qt_meta_stringdata_JMeshDualsDialog_t qt_meta_stringdata_JMeshDualsDialog = {
    {
QT_MOC_LITERAL(0, 0, 16), // "JMeshDualsDialog"
QT_MOC_LITERAL(1, 17, 12), // "displayDuals"
QT_MOC_LITERAL(2, 30, 0), // ""
QT_MOC_LITERAL(3, 31, 9), // "pickSeeds"
QT_MOC_LITERAL(4, 41, 9), // "nextSeeds"
QT_MOC_LITERAL(5, 51, 10), // "selectDual"
QT_MOC_LITERAL(6, 62, 11), // "getSomeDual"
QT_MOC_LITERAL(7, 74, 10), // "modifyDual"
QT_MOC_LITERAL(8, 85, 8), // "getStyle"
QT_MOC_LITERAL(9, 94, 12), // "checkDisplay"
QT_MOC_LITERAL(10, 107, 12), // "getAllChords"
QT_MOC_LITERAL(11, 120, 15), // "getAllHexSheets"
QT_MOC_LITERAL(12, 136, 11), // "clearChords"
QT_MOC_LITERAL(13, 148, 15), // "resetPrimalMesh"
QT_MOC_LITERAL(14, 164, 6), // "reject"
QT_MOC_LITERAL(15, 171, 13), // "keyPressEvent"
QT_MOC_LITERAL(16, 185, 10), // "QKeyEvent*"
QT_MOC_LITERAL(17, 196, 1) // "e"

    },
    "JMeshDualsDialog\0displayDuals\0\0pickSeeds\0"
    "nextSeeds\0selectDual\0getSomeDual\0"
    "modifyDual\0getStyle\0checkDisplay\0"
    "getAllChords\0getAllHexSheets\0clearChords\0"
    "resetPrimalMesh\0reject\0keyPressEvent\0"
    "QKeyEvent*\0e"
};
#undef QT_MOC_LITERAL

static const uint qt_meta_data_JMeshDualsDialog[] = {

 // content:
       7,       // revision
       0,       // classname
       0,    0, // classinfo
      14,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       0,       // signalCount

 // slots: name, argc, parameters, tag, flags
       1,    0,   84,    2, 0x08 /* Private */,
       3,    0,   85,    2, 0x08 /* Private */,
       4,    0,   86,    2, 0x08 /* Private */,
       5,    0,   87,    2, 0x08 /* Private */,
       6,    0,   88,    2, 0x08 /* Private */,
       7,    0,   89,    2, 0x08 /* Private */,
       8,    0,   90,    2, 0x08 /* Private */,
       9,    0,   91,    2, 0x08 /* Private */,
      10,    0,   92,    2, 0x08 /* Private */,
      11,    0,   93,    2, 0x08 /* Private */,
      12,    0,   94,    2, 0x08 /* Private */,
      13,    0,   95,    2, 0x08 /* Private */,
      14,    0,   96,    2, 0x08 /* Private */,
      15,    1,   97,    2, 0x08 /* Private */,

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
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void, 0x80000000 | 16,   17,

       0        // eod
};

void JMeshDualsDialog::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    if (_c == QMetaObject::InvokeMetaMethod) {
        JMeshDualsDialog *_t = static_cast<JMeshDualsDialog *>(_o);
        Q_UNUSED(_t)
        switch (_id) {
        case 0: _t->displayDuals(); break;
        case 1: _t->pickSeeds(); break;
        case 2: _t->nextSeeds(); break;
        case 3: _t->selectDual(); break;
        case 4: _t->getSomeDual(); break;
        case 5: _t->modifyDual(); break;
        case 6: _t->getStyle(); break;
        case 7: _t->checkDisplay(); break;
        case 8: _t->getAllChords(); break;
        case 9: _t->getAllHexSheets(); break;
        case 10: _t->clearChords(); break;
        case 11: _t->resetPrimalMesh(); break;
        case 12: _t->reject(); break;
        case 13: _t->keyPressEvent((*reinterpret_cast< QKeyEvent*(*)>(_a[1]))); break;
        default: ;
        }
    }
}

const QMetaObject JMeshDualsDialog::staticMetaObject = {
    { &QDialog::staticMetaObject, qt_meta_stringdata_JMeshDualsDialog.data,
      qt_meta_data_JMeshDualsDialog,  qt_static_metacall, Q_NULLPTR, Q_NULLPTR}
};


const QMetaObject *JMeshDualsDialog::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->dynamicMetaObject() : &staticMetaObject;
}

void *JMeshDualsDialog::qt_metacast(const char *_clname)
{
    if (!_clname) return Q_NULLPTR;
    if (!strcmp(_clname, qt_meta_stringdata_JMeshDualsDialog.stringdata0))
        return static_cast<void*>(const_cast< JMeshDualsDialog*>(this));
    if (!strcmp(_clname, "Ui::MeshDualsDialog"))
        return static_cast< Ui::MeshDualsDialog*>(const_cast< JMeshDualsDialog*>(this));
    return QDialog::qt_metacast(_clname);
}

int JMeshDualsDialog::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QDialog::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        if (_id < 14)
            qt_static_metacall(this, _c, _id, _a);
        _id -= 14;
    } else if (_c == QMetaObject::RegisterMethodArgumentMetaType) {
        if (_id < 14)
            *reinterpret_cast<int*>(_a[0]) = -1;
        _id -= 14;
    }
    return _id;
}
QT_END_MOC_NAMESPACE
