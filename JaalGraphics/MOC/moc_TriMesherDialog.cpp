/****************************************************************************
** Meta object code from reading C++ file 'TriMesherDialog.hpp'
**
** Created by: The Qt Meta Object Compiler version 67 (Qt 5.6.0)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "../include/TriMesherDialog.hpp"
#include <QtCore/qbytearray.h>
#include <QtCore/qmetatype.h>
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'TriMesherDialog.hpp' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 67
#error "This file was generated using the moc from 5.6.0. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
struct qt_meta_stringdata_JTriMesherDialog_t {
    QByteArrayData data[15];
    char stringdata0[266];
};
#define QT_MOC_LITERAL(idx, ofs, len) \
    Q_STATIC_BYTE_ARRAY_DATA_HEADER_INITIALIZER_WITH_OFFSET(len, \
    qptrdiff(offsetof(qt_meta_stringdata_JTriMesherDialog_t, stringdata0) + ofs \
        - idx * sizeof(QByteArrayData)) \
    )
static const qt_meta_stringdata_JTriMesherDialog_t qt_meta_stringdata_JTriMesherDialog = {
    {
QT_MOC_LITERAL(0, 0, 16), // "JTriMesherDialog"
QT_MOC_LITERAL(1, 17, 17), // "openCleanupDialog"
QT_MOC_LITERAL(2, 35, 0), // ""
QT_MOC_LITERAL(3, 36, 19), // "openIsotropicDialog"
QT_MOC_LITERAL(4, 56, 21), // "openInstantMeshDialog"
QT_MOC_LITERAL(5, 78, 25), // "openQualityDelaunayDialog"
QT_MOC_LITERAL(6, 104, 28), // "openSurfReconstructionDialog"
QT_MOC_LITERAL(7, 133, 24), // "openContourEditingDialog"
QT_MOC_LITERAL(8, 158, 19), // "openHolesFillDialog"
QT_MOC_LITERAL(9, 178, 13), // "rejectNewMesh"
QT_MOC_LITERAL(10, 192, 10), // "genNewMesh"
QT_MOC_LITERAL(11, 203, 14), // "advancingFront"
QT_MOC_LITERAL(12, 218, 24), // "getIntrinsicDelaunayMesh"
QT_MOC_LITERAL(13, 243, 10), // "rejectMesh"
QT_MOC_LITERAL(14, 254, 11) // "closeDialog"

    },
    "JTriMesherDialog\0openCleanupDialog\0\0"
    "openIsotropicDialog\0openInstantMeshDialog\0"
    "openQualityDelaunayDialog\0"
    "openSurfReconstructionDialog\0"
    "openContourEditingDialog\0openHolesFillDialog\0"
    "rejectNewMesh\0genNewMesh\0advancingFront\0"
    "getIntrinsicDelaunayMesh\0rejectMesh\0"
    "closeDialog"
};
#undef QT_MOC_LITERAL

static const uint qt_meta_data_JTriMesherDialog[] = {

 // content:
       7,       // revision
       0,       // classname
       0,    0, // classinfo
      13,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       0,       // signalCount

 // slots: name, argc, parameters, tag, flags
       1,    0,   79,    2, 0x08 /* Private */,
       3,    0,   80,    2, 0x08 /* Private */,
       4,    0,   81,    2, 0x08 /* Private */,
       5,    0,   82,    2, 0x08 /* Private */,
       6,    0,   83,    2, 0x08 /* Private */,
       7,    0,   84,    2, 0x08 /* Private */,
       8,    0,   85,    2, 0x08 /* Private */,
       9,    0,   86,    2, 0x08 /* Private */,
      10,    0,   87,    2, 0x08 /* Private */,
      11,    0,   88,    2, 0x08 /* Private */,
      12,    0,   89,    2, 0x08 /* Private */,
      13,    0,   90,    2, 0x08 /* Private */,
      14,    0,   91,    2, 0x08 /* Private */,

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

       0        // eod
};

void JTriMesherDialog::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    if (_c == QMetaObject::InvokeMetaMethod) {
        JTriMesherDialog *_t = static_cast<JTriMesherDialog *>(_o);
        Q_UNUSED(_t)
        switch (_id) {
        case 0: _t->openCleanupDialog(); break;
        case 1: _t->openIsotropicDialog(); break;
        case 2: _t->openInstantMeshDialog(); break;
        case 3: _t->openQualityDelaunayDialog(); break;
        case 4: _t->openSurfReconstructionDialog(); break;
        case 5: _t->openContourEditingDialog(); break;
        case 6: _t->openHolesFillDialog(); break;
        case 7: _t->rejectNewMesh(); break;
        case 8: _t->genNewMesh(); break;
        case 9: _t->advancingFront(); break;
        case 10: _t->getIntrinsicDelaunayMesh(); break;
        case 11: _t->rejectMesh(); break;
        case 12: _t->closeDialog(); break;
        default: ;
        }
    }
    Q_UNUSED(_a);
}

const QMetaObject JTriMesherDialog::staticMetaObject = {
    { &QDialog::staticMetaObject, qt_meta_stringdata_JTriMesherDialog.data,
      qt_meta_data_JTriMesherDialog,  qt_static_metacall, Q_NULLPTR, Q_NULLPTR}
};


const QMetaObject *JTriMesherDialog::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->dynamicMetaObject() : &staticMetaObject;
}

void *JTriMesherDialog::qt_metacast(const char *_clname)
{
    if (!_clname) return Q_NULLPTR;
    if (!strcmp(_clname, qt_meta_stringdata_JTriMesherDialog.stringdata0))
        return static_cast<void*>(const_cast< JTriMesherDialog*>(this));
    if (!strcmp(_clname, "Ui::TriMesherDialog"))
        return static_cast< Ui::TriMesherDialog*>(const_cast< JTriMesherDialog*>(this));
    return QDialog::qt_metacast(_clname);
}

int JTriMesherDialog::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QDialog::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        if (_id < 13)
            qt_static_metacall(this, _c, _id, _a);
        _id -= 13;
    } else if (_c == QMetaObject::RegisterMethodArgumentMetaType) {
        if (_id < 13)
            *reinterpret_cast<int*>(_a[0]) = -1;
        _id -= 13;
    }
    return _id;
}
QT_END_MOC_NAMESPACE
