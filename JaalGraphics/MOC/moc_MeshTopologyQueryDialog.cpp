/****************************************************************************
** Meta object code from reading C++ file 'MeshTopologyQueryDialog.hpp'
**
** Created by: The Qt Meta Object Compiler version 67 (Qt 5.6.0)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "MeshTopologyQueryDialog.hpp"
#include <QtCore/qbytearray.h>
#include <QtCore/qmetatype.h>
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'MeshTopologyQueryDialog.hpp' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 67
#error "This file was generated using the moc from 5.6.0. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
struct qt_meta_stringdata_JMeshTopologyQueryDialog_t {
    QByteArrayData data[14];
    char stringdata0[205];
};
#define QT_MOC_LITERAL(idx, ofs, len) \
    Q_STATIC_BYTE_ARRAY_DATA_HEADER_INITIALIZER_WITH_OFFSET(len, \
    qptrdiff(offsetof(qt_meta_stringdata_JMeshTopologyQueryDialog_t, stringdata0) + ofs \
        - idx * sizeof(QByteArrayData)) \
    )
static const qt_meta_stringdata_JMeshTopologyQueryDialog_t qt_meta_stringdata_JMeshTopologyQueryDialog = {
    {
QT_MOC_LITERAL(0, 0, 24), // "JMeshTopologyQueryDialog"
QT_MOC_LITERAL(1, 25, 13), // "getConsistent"
QT_MOC_LITERAL(2, 39, 0), // ""
QT_MOC_LITERAL(3, 40, 13), // "displayMatrix"
QT_MOC_LITERAL(4, 54, 19), // "displayPrimalMatrix"
QT_MOC_LITERAL(5, 74, 17), // "displayDualMatrix"
QT_MOC_LITERAL(6, 92, 11), // "closeDialog"
QT_MOC_LITERAL(7, 104, 10), // "setLapGrid"
QT_MOC_LITERAL(8, 115, 11), // "setLapFonts"
QT_MOC_LITERAL(9, 127, 12), // "setStepLabel"
QT_MOC_LITERAL(10, 140, 12), // "setPointSize"
QT_MOC_LITERAL(11, 153, 14), // "getBettiNumber"
QT_MOC_LITERAL(12, 168, 24), // "openMeshComponentsDialog"
QT_MOC_LITERAL(13, 193, 11) // "getOrphaned"

    },
    "JMeshTopologyQueryDialog\0getConsistent\0"
    "\0displayMatrix\0displayPrimalMatrix\0"
    "displayDualMatrix\0closeDialog\0setLapGrid\0"
    "setLapFonts\0setStepLabel\0setPointSize\0"
    "getBettiNumber\0openMeshComponentsDialog\0"
    "getOrphaned"
};
#undef QT_MOC_LITERAL

static const uint qt_meta_data_JMeshTopologyQueryDialog[] = {

 // content:
       7,       // revision
       0,       // classname
       0,    0, // classinfo
      12,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       0,       // signalCount

 // slots: name, argc, parameters, tag, flags
       1,    0,   74,    2, 0x08 /* Private */,
       3,    0,   75,    2, 0x08 /* Private */,
       4,    0,   76,    2, 0x08 /* Private */,
       5,    0,   77,    2, 0x08 /* Private */,
       6,    0,   78,    2, 0x08 /* Private */,
       7,    0,   79,    2, 0x08 /* Private */,
       8,    0,   80,    2, 0x08 /* Private */,
       9,    0,   81,    2, 0x08 /* Private */,
      10,    0,   82,    2, 0x08 /* Private */,
      11,    0,   83,    2, 0x08 /* Private */,
      12,    0,   84,    2, 0x08 /* Private */,
      13,    0,   85,    2, 0x08 /* Private */,

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

       0        // eod
};

void JMeshTopologyQueryDialog::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    if (_c == QMetaObject::InvokeMetaMethod) {
        JMeshTopologyQueryDialog *_t = static_cast<JMeshTopologyQueryDialog *>(_o);
        Q_UNUSED(_t)
        switch (_id) {
        case 0: _t->getConsistent(); break;
        case 1: _t->displayMatrix(); break;
        case 2: _t->displayPrimalMatrix(); break;
        case 3: _t->displayDualMatrix(); break;
        case 4: _t->closeDialog(); break;
        case 5: _t->setLapGrid(); break;
        case 6: _t->setLapFonts(); break;
        case 7: _t->setStepLabel(); break;
        case 8: _t->setPointSize(); break;
        case 9: _t->getBettiNumber(); break;
        case 10: _t->openMeshComponentsDialog(); break;
        case 11: _t->getOrphaned(); break;
        default: ;
        }
    }
    Q_UNUSED(_a);
}

const QMetaObject JMeshTopologyQueryDialog::staticMetaObject = {
    { &QDialog::staticMetaObject, qt_meta_stringdata_JMeshTopologyQueryDialog.data,
      qt_meta_data_JMeshTopologyQueryDialog,  qt_static_metacall, Q_NULLPTR, Q_NULLPTR}
};


const QMetaObject *JMeshTopologyQueryDialog::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->dynamicMetaObject() : &staticMetaObject;
}

void *JMeshTopologyQueryDialog::qt_metacast(const char *_clname)
{
    if (!_clname) return Q_NULLPTR;
    if (!strcmp(_clname, qt_meta_stringdata_JMeshTopologyQueryDialog.stringdata0))
        return static_cast<void*>(const_cast< JMeshTopologyQueryDialog*>(this));
    if (!strcmp(_clname, "Ui::MeshTopologyQueryDialog"))
        return static_cast< Ui::MeshTopologyQueryDialog*>(const_cast< JMeshTopologyQueryDialog*>(this));
    return QDialog::qt_metacast(_clname);
}

int JMeshTopologyQueryDialog::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QDialog::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        if (_id < 12)
            qt_static_metacall(this, _c, _id, _a);
        _id -= 12;
    } else if (_c == QMetaObject::RegisterMethodArgumentMetaType) {
        if (_id < 12)
            *reinterpret_cast<int*>(_a[0]) = -1;
        _id -= 12;
    }
    return _id;
}
QT_END_MOC_NAMESPACE
