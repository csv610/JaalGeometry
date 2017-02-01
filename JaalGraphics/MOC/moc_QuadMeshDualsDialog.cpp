/****************************************************************************
** Meta object code from reading C++ file 'QuadMeshDualsDialog.hpp'
**
** Created by: The Qt Meta Object Compiler version 67 (Qt 5.6.0)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "../include/QuadMeshDualsDialog.hpp"
#include <QtCore/qbytearray.h>
#include <QtCore/qmetatype.h>
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'QuadMeshDualsDialog.hpp' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 67
#error "This file was generated using the moc from 5.6.0. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
struct qt_meta_stringdata_JQuadMeshDualsDialog_t {
    QByteArrayData data[19];
    char stringdata0[240];
};
#define QT_MOC_LITERAL(idx, ofs, len) \
    Q_STATIC_BYTE_ARRAY_DATA_HEADER_INITIALIZER_WITH_OFFSET(len, \
    qptrdiff(offsetof(qt_meta_stringdata_JQuadMeshDualsDialog_t, stringdata0) + ofs \
        - idx * sizeof(QByteArrayData)) \
    )
static const qt_meta_stringdata_JQuadMeshDualsDialog_t qt_meta_stringdata_JQuadMeshDualsDialog = {
    {
QT_MOC_LITERAL(0, 0, 20), // "JQuadMeshDualsDialog"
QT_MOC_LITERAL(1, 21, 13), // "enablePicking"
QT_MOC_LITERAL(2, 35, 0), // ""
QT_MOC_LITERAL(3, 36, 9), // "nextSeeds"
QT_MOC_LITERAL(4, 46, 10), // "selectDual"
QT_MOC_LITERAL(5, 57, 10), // "getNewDual"
QT_MOC_LITERAL(6, 68, 10), // "modifyDual"
QT_MOC_LITERAL(7, 79, 13), // "diceAllChords"
QT_MOC_LITERAL(8, 93, 19), // "displayCyclicChords"
QT_MOC_LITERAL(9, 113, 8), // "getStyle"
QT_MOC_LITERAL(10, 122, 12), // "checkDisplay"
QT_MOC_LITERAL(11, 135, 12), // "getAllChords"
QT_MOC_LITERAL(12, 148, 11), // "clearChords"
QT_MOC_LITERAL(13, 160, 11), // "closeDialog"
QT_MOC_LITERAL(14, 172, 13), // "keyPressEvent"
QT_MOC_LITERAL(15, 186, 10), // "QKeyEvent*"
QT_MOC_LITERAL(16, 197, 1), // "e"
QT_MOC_LITERAL(17, 199, 16), // "getMaxEdgeLength"
QT_MOC_LITERAL(18, 216, 23) // "openBoundaryLayerDialog"

    },
    "JQuadMeshDualsDialog\0enablePicking\0\0"
    "nextSeeds\0selectDual\0getNewDual\0"
    "modifyDual\0diceAllChords\0displayCyclicChords\0"
    "getStyle\0checkDisplay\0getAllChords\0"
    "clearChords\0closeDialog\0keyPressEvent\0"
    "QKeyEvent*\0e\0getMaxEdgeLength\0"
    "openBoundaryLayerDialog"
};
#undef QT_MOC_LITERAL

static const uint qt_meta_data_JQuadMeshDualsDialog[] = {

 // content:
       7,       // revision
       0,       // classname
       0,    0, // classinfo
      15,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       0,       // signalCount

 // slots: name, argc, parameters, tag, flags
       1,    0,   89,    2, 0x08 /* Private */,
       3,    0,   90,    2, 0x08 /* Private */,
       4,    0,   91,    2, 0x08 /* Private */,
       5,    0,   92,    2, 0x08 /* Private */,
       6,    0,   93,    2, 0x08 /* Private */,
       7,    0,   94,    2, 0x08 /* Private */,
       8,    0,   95,    2, 0x08 /* Private */,
       9,    0,   96,    2, 0x08 /* Private */,
      10,    0,   97,    2, 0x08 /* Private */,
      11,    0,   98,    2, 0x08 /* Private */,
      12,    0,   99,    2, 0x08 /* Private */,
      13,    0,  100,    2, 0x08 /* Private */,
      14,    1,  101,    2, 0x08 /* Private */,
      17,    0,  104,    2, 0x08 /* Private */,
      18,    0,  105,    2, 0x08 /* Private */,

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
    QMetaType::Void, 0x80000000 | 15,   16,
    QMetaType::Void,
    QMetaType::Void,

       0        // eod
};

void JQuadMeshDualsDialog::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    if (_c == QMetaObject::InvokeMetaMethod) {
        JQuadMeshDualsDialog *_t = static_cast<JQuadMeshDualsDialog *>(_o);
        Q_UNUSED(_t)
        switch (_id) {
        case 0: _t->enablePicking(); break;
        case 1: _t->nextSeeds(); break;
        case 2: _t->selectDual(); break;
        case 3: _t->getNewDual(); break;
        case 4: _t->modifyDual(); break;
        case 5: _t->diceAllChords(); break;
        case 6: _t->displayCyclicChords(); break;
        case 7: _t->getStyle(); break;
        case 8: _t->checkDisplay(); break;
        case 9: _t->getAllChords(); break;
        case 10: _t->clearChords(); break;
        case 11: _t->closeDialog(); break;
        case 12: _t->keyPressEvent((*reinterpret_cast< QKeyEvent*(*)>(_a[1]))); break;
        case 13: _t->getMaxEdgeLength(); break;
        case 14: _t->openBoundaryLayerDialog(); break;
        default: ;
        }
    }
}

const QMetaObject JQuadMeshDualsDialog::staticMetaObject = {
    { &QDialog::staticMetaObject, qt_meta_stringdata_JQuadMeshDualsDialog.data,
      qt_meta_data_JQuadMeshDualsDialog,  qt_static_metacall, Q_NULLPTR, Q_NULLPTR}
};


const QMetaObject *JQuadMeshDualsDialog::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->dynamicMetaObject() : &staticMetaObject;
}

void *JQuadMeshDualsDialog::qt_metacast(const char *_clname)
{
    if (!_clname) return Q_NULLPTR;
    if (!strcmp(_clname, qt_meta_stringdata_JQuadMeshDualsDialog.stringdata0))
        return static_cast<void*>(const_cast< JQuadMeshDualsDialog*>(this));
    if (!strcmp(_clname, "Ui::QuadMeshDualsDialog"))
        return static_cast< Ui::QuadMeshDualsDialog*>(const_cast< JQuadMeshDualsDialog*>(this));
    return QDialog::qt_metacast(_clname);
}

int JQuadMeshDualsDialog::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QDialog::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        if (_id < 15)
            qt_static_metacall(this, _c, _id, _a);
        _id -= 15;
    } else if (_c == QMetaObject::RegisterMethodArgumentMetaType) {
        if (_id < 15)
            *reinterpret_cast<int*>(_a[0]) = -1;
        _id -= 15;
    }
    return _id;
}
QT_END_MOC_NAMESPACE
