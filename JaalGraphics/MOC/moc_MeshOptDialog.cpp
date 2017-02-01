/****************************************************************************
** Meta object code from reading C++ file 'MeshOptDialog.hpp'
**
** Created by: The Qt Meta Object Compiler version 67 (Qt 5.6.0)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "../include/MeshOptDialog.hpp"
#include <QtCore/qbytearray.h>
#include <QtCore/qmetatype.h>
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'MeshOptDialog.hpp' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 67
#error "This file was generated using the moc from 5.6.0. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
struct qt_meta_stringdata_JMeshOptDialog_t {
    QByteArrayData data[14];
    char stringdata0[200];
};
#define QT_MOC_LITERAL(idx, ofs, len) \
    Q_STATIC_BYTE_ARRAY_DATA_HEADER_INITIALIZER_WITH_OFFSET(len, \
    qptrdiff(offsetof(qt_meta_stringdata_JMeshOptDialog_t, stringdata0) + ofs \
        - idx * sizeof(QByteArrayData)) \
    )
static const qt_meta_stringdata_JMeshOptDialog_t qt_meta_stringdata_JMeshOptDialog = {
    {
QT_MOC_LITERAL(0, 0, 14), // "JMeshOptDialog"
QT_MOC_LITERAL(1, 15, 11), // "getOriginal"
QT_MOC_LITERAL(2, 27, 0), // ""
QT_MOC_LITERAL(3, 28, 9), // "nonlinear"
QT_MOC_LITERAL(4, 38, 14), // "setConstraints"
QT_MOC_LITERAL(5, 53, 11), // "closeDialog"
QT_MOC_LITERAL(6, 65, 8), // "untangle"
QT_MOC_LITERAL(7, 74, 13), // "reparamCurves"
QT_MOC_LITERAL(8, 88, 12), // "smoothCurves"
QT_MOC_LITERAL(9, 101, 17), // "openLaplaceDialog"
QT_MOC_LITERAL(10, 119, 15), // "openLloydDialog"
QT_MOC_LITERAL(11, 135, 27), // "openMeanCurvatureFlowDialog"
QT_MOC_LITERAL(12, 163, 18), // "openUntangleDialog"
QT_MOC_LITERAL(13, 182, 17) // "openShapeOpDialog"

    },
    "JMeshOptDialog\0getOriginal\0\0nonlinear\0"
    "setConstraints\0closeDialog\0untangle\0"
    "reparamCurves\0smoothCurves\0openLaplaceDialog\0"
    "openLloydDialog\0openMeanCurvatureFlowDialog\0"
    "openUntangleDialog\0openShapeOpDialog"
};
#undef QT_MOC_LITERAL

static const uint qt_meta_data_JMeshOptDialog[] = {

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

void JMeshOptDialog::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    if (_c == QMetaObject::InvokeMetaMethod) {
        JMeshOptDialog *_t = static_cast<JMeshOptDialog *>(_o);
        Q_UNUSED(_t)
        switch (_id) {
        case 0: _t->getOriginal(); break;
        case 1: _t->nonlinear(); break;
        case 2: _t->setConstraints(); break;
        case 3: _t->closeDialog(); break;
        case 4: _t->untangle(); break;
        case 5: _t->reparamCurves(); break;
        case 6: _t->smoothCurves(); break;
        case 7: _t->openLaplaceDialog(); break;
        case 8: _t->openLloydDialog(); break;
        case 9: _t->openMeanCurvatureFlowDialog(); break;
        case 10: _t->openUntangleDialog(); break;
        case 11: _t->openShapeOpDialog(); break;
        default: ;
        }
    }
    Q_UNUSED(_a);
}

const QMetaObject JMeshOptDialog::staticMetaObject = {
    { &QDialog::staticMetaObject, qt_meta_stringdata_JMeshOptDialog.data,
      qt_meta_data_JMeshOptDialog,  qt_static_metacall, Q_NULLPTR, Q_NULLPTR}
};


const QMetaObject *JMeshOptDialog::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->dynamicMetaObject() : &staticMetaObject;
}

void *JMeshOptDialog::qt_metacast(const char *_clname)
{
    if (!_clname) return Q_NULLPTR;
    if (!strcmp(_clname, qt_meta_stringdata_JMeshOptDialog.stringdata0))
        return static_cast<void*>(const_cast< JMeshOptDialog*>(this));
    if (!strcmp(_clname, "Ui::MeshOptimizationDialog"))
        return static_cast< Ui::MeshOptimizationDialog*>(const_cast< JMeshOptDialog*>(this));
    return QDialog::qt_metacast(_clname);
}

int JMeshOptDialog::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
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
