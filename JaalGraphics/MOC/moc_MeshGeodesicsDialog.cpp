/****************************************************************************
** Meta object code from reading C++ file 'MeshGeodesicsDialog.hpp'
**
** Created by: The Qt Meta Object Compiler version 67 (Qt 5.6.0)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "MeshGeodesicsDialog.hpp"
#include <QtCore/qbytearray.h>
#include <QtCore/qmetatype.h>
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'MeshGeodesicsDialog.hpp' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 67
#error "This file was generated using the moc from 5.6.0. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
struct qt_meta_stringdata_JMeshGeodesicsDialog_t {
    QByteArrayData data[11];
    char stringdata0[153];
};
#define QT_MOC_LITERAL(idx, ofs, len) \
    Q_STATIC_BYTE_ARRAY_DATA_HEADER_INITIALIZER_WITH_OFFSET(len, \
    qptrdiff(offsetof(qt_meta_stringdata_JMeshGeodesicsDialog_t, stringdata0) + ofs \
        - idx * sizeof(QByteArrayData)) \
    )
static const qt_meta_stringdata_JMeshGeodesicsDialog_t qt_meta_stringdata_JMeshGeodesicsDialog = {
    {
QT_MOC_LITERAL(0, 0, 20), // "JMeshGeodesicsDialog"
QT_MOC_LITERAL(1, 21, 15), // "getDijkstraPath"
QT_MOC_LITERAL(2, 37, 0), // ""
QT_MOC_LITERAL(3, 38, 10), // "checkState"
QT_MOC_LITERAL(4, 49, 11), // "closeDialog"
QT_MOC_LITERAL(5, 61, 20), // "openNodeAttribDialog"
QT_MOC_LITERAL(6, 82, 20), // "openEdgeAttribDialog"
QT_MOC_LITERAL(7, 103, 12), // "startNewPath"
QT_MOC_LITERAL(8, 116, 17), // "deleteLastSegment"
QT_MOC_LITERAL(9, 134, 9), // "deleteAll"
QT_MOC_LITERAL(10, 144, 8) // "heatFlow"

    },
    "JMeshGeodesicsDialog\0getDijkstraPath\0"
    "\0checkState\0closeDialog\0openNodeAttribDialog\0"
    "openEdgeAttribDialog\0startNewPath\0"
    "deleteLastSegment\0deleteAll\0heatFlow"
};
#undef QT_MOC_LITERAL

static const uint qt_meta_data_JMeshGeodesicsDialog[] = {

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

void JMeshGeodesicsDialog::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    if (_c == QMetaObject::InvokeMetaMethod) {
        JMeshGeodesicsDialog *_t = static_cast<JMeshGeodesicsDialog *>(_o);
        Q_UNUSED(_t)
        switch (_id) {
        case 0: _t->getDijkstraPath(); break;
        case 1: _t->checkState(); break;
        case 2: _t->closeDialog(); break;
        case 3: _t->openNodeAttribDialog(); break;
        case 4: _t->openEdgeAttribDialog(); break;
        case 5: _t->startNewPath(); break;
        case 6: _t->deleteLastSegment(); break;
        case 7: _t->deleteAll(); break;
        case 8: _t->heatFlow(); break;
        default: ;
        }
    }
    Q_UNUSED(_a);
}

const QMetaObject JMeshGeodesicsDialog::staticMetaObject = {
    { &QDialog::staticMetaObject, qt_meta_stringdata_JMeshGeodesicsDialog.data,
      qt_meta_data_JMeshGeodesicsDialog,  qt_static_metacall, Q_NULLPTR, Q_NULLPTR}
};


const QMetaObject *JMeshGeodesicsDialog::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->dynamicMetaObject() : &staticMetaObject;
}

void *JMeshGeodesicsDialog::qt_metacast(const char *_clname)
{
    if (!_clname) return Q_NULLPTR;
    if (!strcmp(_clname, qt_meta_stringdata_JMeshGeodesicsDialog.stringdata0))
        return static_cast<void*>(const_cast< JMeshGeodesicsDialog*>(this));
    if (!strcmp(_clname, "Ui::MeshGeodesicDialog"))
        return static_cast< Ui::MeshGeodesicDialog*>(const_cast< JMeshGeodesicsDialog*>(this));
    return QDialog::qt_metacast(_clname);
}

int JMeshGeodesicsDialog::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
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
