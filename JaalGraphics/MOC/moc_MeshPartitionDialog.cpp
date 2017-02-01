/****************************************************************************
** Meta object code from reading C++ file 'MeshPartitionDialog.hpp'
**
** Created by: The Qt Meta Object Compiler version 67 (Qt 5.6.0)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "../include/MeshPartitionDialog.hpp"
#include <QtCore/qbytearray.h>
#include <QtCore/qmetatype.h>
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'MeshPartitionDialog.hpp' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 67
#error "This file was generated using the moc from 5.6.0. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
struct qt_meta_stringdata_JMeshPartitionDialog_t {
    QByteArrayData data[18];
    char stringdata0[310];
};
#define QT_MOC_LITERAL(idx, ofs, len) \
    Q_STATIC_BYTE_ARRAY_DATA_HEADER_INITIALIZER_WITH_OFFSET(len, \
    qptrdiff(offsetof(qt_meta_stringdata_JMeshPartitionDialog_t, stringdata0) + ofs \
        - idx * sizeof(QByteArrayData)) \
    )
static const qt_meta_stringdata_JMeshPartitionDialog_t qt_meta_stringdata_JMeshPartitionDialog = {
    {
QT_MOC_LITERAL(0, 0, 20), // "JMeshPartitionDialog"
QT_MOC_LITERAL(1, 21, 12), // "optPartition"
QT_MOC_LITERAL(2, 34, 0), // ""
QT_MOC_LITERAL(3, 35, 14), // "applyAlgorithm"
QT_MOC_LITERAL(4, 50, 25), // "openCellsPartitionsDialog"
QT_MOC_LITERAL(5, 76, 25), // "openFacesPartitionsDialog"
QT_MOC_LITERAL(6, 102, 25), // "openEdgesPartitionsDialog"
QT_MOC_LITERAL(7, 128, 25), // "openNodesPartitionsDialog"
QT_MOC_LITERAL(8, 154, 24), // "openNormalClustersDialog"
QT_MOC_LITERAL(9, 179, 24), // "getRegionGrowingClusters"
QT_MOC_LITERAL(10, 204, 22), // "removeZigZagInterfaces"
QT_MOC_LITERAL(11, 227, 14), // "savePartitions"
QT_MOC_LITERAL(12, 242, 8), // "clearAll"
QT_MOC_LITERAL(13, 251, 19), // "getTopologicalDisks"
QT_MOC_LITERAL(14, 271, 13), // "keyPressEvent"
QT_MOC_LITERAL(15, 285, 10), // "QKeyEvent*"
QT_MOC_LITERAL(16, 296, 1), // "e"
QT_MOC_LITERAL(17, 298, 11) // "closeDialog"

    },
    "JMeshPartitionDialog\0optPartition\0\0"
    "applyAlgorithm\0openCellsPartitionsDialog\0"
    "openFacesPartitionsDialog\0"
    "openEdgesPartitionsDialog\0"
    "openNodesPartitionsDialog\0"
    "openNormalClustersDialog\0"
    "getRegionGrowingClusters\0"
    "removeZigZagInterfaces\0savePartitions\0"
    "clearAll\0getTopologicalDisks\0keyPressEvent\0"
    "QKeyEvent*\0e\0closeDialog"
};
#undef QT_MOC_LITERAL

static const uint qt_meta_data_JMeshPartitionDialog[] = {

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
      14,    1,   96,    2, 0x08 /* Private */,
      17,    0,   99,    2, 0x08 /* Private */,

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

       0        // eod
};

void JMeshPartitionDialog::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    if (_c == QMetaObject::InvokeMetaMethod) {
        JMeshPartitionDialog *_t = static_cast<JMeshPartitionDialog *>(_o);
        Q_UNUSED(_t)
        switch (_id) {
        case 0: _t->optPartition(); break;
        case 1: _t->applyAlgorithm(); break;
        case 2: _t->openCellsPartitionsDialog(); break;
        case 3: _t->openFacesPartitionsDialog(); break;
        case 4: _t->openEdgesPartitionsDialog(); break;
        case 5: _t->openNodesPartitionsDialog(); break;
        case 6: _t->openNormalClustersDialog(); break;
        case 7: _t->getRegionGrowingClusters(); break;
        case 8: _t->removeZigZagInterfaces(); break;
        case 9: _t->savePartitions(); break;
        case 10: _t->clearAll(); break;
        case 11: _t->getTopologicalDisks(); break;
        case 12: _t->keyPressEvent((*reinterpret_cast< QKeyEvent*(*)>(_a[1]))); break;
        case 13: _t->closeDialog(); break;
        default: ;
        }
    }
}

const QMetaObject JMeshPartitionDialog::staticMetaObject = {
    { &QDialog::staticMetaObject, qt_meta_stringdata_JMeshPartitionDialog.data,
      qt_meta_data_JMeshPartitionDialog,  qt_static_metacall, Q_NULLPTR, Q_NULLPTR}
};


const QMetaObject *JMeshPartitionDialog::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->dynamicMetaObject() : &staticMetaObject;
}

void *JMeshPartitionDialog::qt_metacast(const char *_clname)
{
    if (!_clname) return Q_NULLPTR;
    if (!strcmp(_clname, qt_meta_stringdata_JMeshPartitionDialog.stringdata0))
        return static_cast<void*>(const_cast< JMeshPartitionDialog*>(this));
    if (!strcmp(_clname, "Ui::MeshPartitionDialog"))
        return static_cast< Ui::MeshPartitionDialog*>(const_cast< JMeshPartitionDialog*>(this));
    return QDialog::qt_metacast(_clname);
}

int JMeshPartitionDialog::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
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
