/****************************************************************************
** Meta object code from reading C++ file 'MeshSegmentationDialog.hpp'
**
** Created by: The Qt Meta Object Compiler version 67 (Qt 5.6.0)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "../include/MeshSegmentationDialog.hpp"
#include <QtCore/qbytearray.h>
#include <QtCore/qmetatype.h>
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'MeshSegmentationDialog.hpp' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 67
#error "This file was generated using the moc from 5.6.0. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
struct qt_meta_stringdata_JMeshSegmentationDialog_t {
    QByteArrayData data[10];
    char stringdata0[135];
};
#define QT_MOC_LITERAL(idx, ofs, len) \
    Q_STATIC_BYTE_ARRAY_DATA_HEADER_INITIALIZER_WITH_OFFSET(len, \
    qptrdiff(offsetof(qt_meta_stringdata_JMeshSegmentationDialog_t, stringdata0) + ofs \
        - idx * sizeof(QByteArrayData)) \
    )
static const qt_meta_stringdata_JMeshSegmentationDialog_t qt_meta_stringdata_JMeshSegmentationDialog = {
    {
QT_MOC_LITERAL(0, 0, 23), // "JMeshSegmentationDialog"
QT_MOC_LITERAL(1, 24, 23), // "openMeshGeodesicsDialog"
QT_MOC_LITERAL(2, 48, 0), // ""
QT_MOC_LITERAL(3, 49, 23), // "openMeshPartitionDialog"
QT_MOC_LITERAL(4, 73, 7), // "segment"
QT_MOC_LITERAL(5, 81, 13), // "keyPressEvent"
QT_MOC_LITERAL(6, 95, 10), // "QKeyEvent*"
QT_MOC_LITERAL(7, 106, 1), // "e"
QT_MOC_LITERAL(8, 108, 14), // "getSDFSegments"
QT_MOC_LITERAL(9, 123, 11) // "closeDialog"

    },
    "JMeshSegmentationDialog\0openMeshGeodesicsDialog\0"
    "\0openMeshPartitionDialog\0segment\0"
    "keyPressEvent\0QKeyEvent*\0e\0getSDFSegments\0"
    "closeDialog"
};
#undef QT_MOC_LITERAL

static const uint qt_meta_data_JMeshSegmentationDialog[] = {

 // content:
       7,       // revision
       0,       // classname
       0,    0, // classinfo
       6,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       0,       // signalCount

 // slots: name, argc, parameters, tag, flags
       1,    0,   44,    2, 0x08 /* Private */,
       3,    0,   45,    2, 0x08 /* Private */,
       4,    0,   46,    2, 0x08 /* Private */,
       5,    1,   47,    2, 0x08 /* Private */,
       8,    0,   50,    2, 0x08 /* Private */,
       9,    0,   51,    2, 0x08 /* Private */,

 // slots: parameters
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void, 0x80000000 | 6,    7,
    QMetaType::Void,
    QMetaType::Void,

       0        // eod
};

void JMeshSegmentationDialog::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    if (_c == QMetaObject::InvokeMetaMethod) {
        JMeshSegmentationDialog *_t = static_cast<JMeshSegmentationDialog *>(_o);
        Q_UNUSED(_t)
        switch (_id) {
        case 0: _t->openMeshGeodesicsDialog(); break;
        case 1: _t->openMeshPartitionDialog(); break;
        case 2: _t->segment(); break;
        case 3: _t->keyPressEvent((*reinterpret_cast< QKeyEvent*(*)>(_a[1]))); break;
        case 4: _t->getSDFSegments(); break;
        case 5: _t->closeDialog(); break;
        default: ;
        }
    }
}

const QMetaObject JMeshSegmentationDialog::staticMetaObject = {
    { &QDialog::staticMetaObject, qt_meta_stringdata_JMeshSegmentationDialog.data,
      qt_meta_data_JMeshSegmentationDialog,  qt_static_metacall, Q_NULLPTR, Q_NULLPTR}
};


const QMetaObject *JMeshSegmentationDialog::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->dynamicMetaObject() : &staticMetaObject;
}

void *JMeshSegmentationDialog::qt_metacast(const char *_clname)
{
    if (!_clname) return Q_NULLPTR;
    if (!strcmp(_clname, qt_meta_stringdata_JMeshSegmentationDialog.stringdata0))
        return static_cast<void*>(const_cast< JMeshSegmentationDialog*>(this));
    if (!strcmp(_clname, "Ui::MeshSegmentationDialog"))
        return static_cast< Ui::MeshSegmentationDialog*>(const_cast< JMeshSegmentationDialog*>(this));
    return QDialog::qt_metacast(_clname);
}

int JMeshSegmentationDialog::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QDialog::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        if (_id < 6)
            qt_static_metacall(this, _c, _id, _a);
        _id -= 6;
    } else if (_c == QMetaObject::RegisterMethodArgumentMetaType) {
        if (_id < 6)
            *reinterpret_cast<int*>(_a[0]) = -1;
        _id -= 6;
    }
    return _id;
}
QT_END_MOC_NAMESPACE
